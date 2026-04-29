module BiomeRAPI

using Biome
using Biome.BIOME4
using Rasters
using DimensionalData

# ─────────────────────────────────────────────────────────────────────────────
# Object store: Julia-side registry so JuliaCall never tries to serialise
# complex structs (PFTClassification, ModelSetup, …) into R objects.
# R holds only an integer handle; Julia owns the actual object.
# ─────────────────────────────────────────────────────────────────────────────
const _STORE = Dict{Int, Any}()
const _STORE_COUNTER = Ref{Int}(0)

function _store!(obj)
    _STORE_COUNTER[] += 1
    id = _STORE_COUNTER[]
    _STORE[id] = obj
    return id          # return plain Int → safe for JuliaCall
end

function _retrieve(id::Int)
    haskey(_STORE, id) || error("No Julia object with handle $(id). Was it already freed?")
    return _STORE[id]
end

function free_handle(id::Int)
    delete!(_STORE, id)
    return nothing
end


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

_to_sym_dict(x) = x
function _to_sym_dict(d::AbstractDict)
    out = Dict{Symbol,Any}()
    for (k, v) in d
        ks = k isa Symbol ? k : Symbol(String(k))
        out[ks] = v
    end
    return out
end

_namedtuple_from_sym_dict(d::AbstractDict{Symbol,Any}) = (; (k => d[k] for k in keys(d))...)

function _normalize_constraints(x)
    x === nothing && return nothing
    d = _to_sym_dict(x)
    out = Dict{Symbol,Any}()
    for (k, v) in d
        v === nothing && continue
        v isa AbstractVector || error("Constraint $(k) must be a length-2 vector. Got $(typeof(v)).")
        length(v) == 2      || error("Constraint $(k) must be length-2. Got length=$(length(v)).")
        out[k] = [Float64(v[1]), Float64(v[2])]
    end
    return _namedtuple_from_sym_dict(out)
end

function _normalize_triplet(x, label::String)
    x === nothing && return nothing
    d = _to_sym_dict(x)
    for needed in (:clt, :prec, :temp)
        haskey(d, needed) || error("$(label) must contain keys clt, prec, temp. Missing: $(needed)")
    end
    return (clt = Float64(d[:clt]), prec = Float64(d[:prec]), temp = Float64(d[:temp]))
end

function _get_string(d::AbstractDict{Symbol,Any}, key::Symbol; required::Bool=false)
    if !haskey(d, key) || d[key] === nothing
        required && error("Missing required field: $(key)")
        return nothing
    end
    return String(d[key])
end


# ─────────────────────────────────────────────────────────────────────────────
# PFT construction
# ─────────────────────────────────────────────────────────────────────────────

"""
    transform_pft_to_structure(spec) -> AbstractPFT

Convert one R `pft()` spec into a concrete Julia PFT by calling the
constructor named by `spec[:type]`.

All parameters default to the type's BASE_DEFAULTS entry (via the
`ctor(; kwargs...)` → `base_pft()` chain already in Biome.jl). The only
special handling needed here is constraints: a partial constraints NamedTuple
from R would wholesale replace the defaults entry, so we pre-merge user
constraint fields over the type defaults before passing them in.
"""
function transform_pft_to_structure(spec)
    d0 = _to_sym_dict(spec)

    typestr      = _get_string(d0, :type; required=true)
    pft_ctor_sym = Symbol(typestr)

    isdefined(Biome, pft_ctor_sym) ||
        error("Biome has no constructor named $(typestr). Check spelling / exports.")
    ctor = getproperty(Biome, pft_ctor_sym)

    kwargs = Dict{Symbol,Any}()

    # --- params (dominance_factor, kk, Emax, etc.) ---
    if haskey(d0, :params) && d0[:params] !== nothing
        for (k, v) in _to_sym_dict(d0[:params])
            kwargs[k] = v
        end
    end

    # --- name ---
    if haskey(d0, :name) && d0[:name] !== nothing && !haskey(kwargs, :name)
        kwargs[:name] = String(d0[:name])
    end

    # --- constraints: find this type's defaults then merge user values over them ---
    # base_pft(AbstractType; kwargs...) merges kwargs over BASE_DEFAULTS, but a
    # partial constraints NamedTuple in kwargs replaces the whole :constraints
    # entry. We therefore resolve the full default constraints for this concrete
    # type and merge field-by-field before passing to the constructor.
    default_constraints = (
        tcm      = [-Inf, +Inf],
        tmin     = [-Inf, +Inf],
        gdd5     = [-Inf, +Inf],
        gdd0     = [-Inf, +Inf],
        twm      = [-Inf, +Inf],
        maxdepth = [-Inf, +Inf],
        swb      = [-Inf, +Inf]
    )
    for (abstract_type, type_defaults) in Biome.BASE_DEFAULTS
        if ctor <: abstract_type
            default_constraints = type_defaults[:constraints]
            break
        end
    end

    user_cons = (haskey(d0, :constraints) && d0[:constraints] !== nothing) ?
        _normalize_constraints(d0[:constraints]) : (;)
    kwargs[:constraints] = merge(default_constraints, user_cons)

    # --- mean_val / sd_val ---
    haskey(d0, :mean_val) && d0[:mean_val] !== nothing &&
        (kwargs[:mean_val] = _normalize_triplet(d0[:mean_val], "mean_val"))
    haskey(d0, :sd_val) && d0[:sd_val] !== nothing &&
        (kwargs[:sd_val] = _normalize_triplet(d0[:sd_val], "sd_val"))

    # ctor(; kwargs...) calls base_pft(AbstractType; kwargs...) internally,
    # which merges all remaining kwargs over BASE_DEFAULTS — so every
    # PFTCharacteristics field not supplied by the user gets its type default.
    nt = _namedtuple_from_sym_dict(kwargs)
    return ctor(; nt...)
end

function transform_pftlist_to_structure(r_list)
    specs = r_list isa AbstractVector ? r_list : collect(r_list)
    return [transform_pft_to_structure(s) for s in specs]
end


# ─────────────────────────────────────────────────────────────────────────────
# PFTClassification factories  →  return handles, never raw structs
# ─────────────────────────────────────────────────────────────────────────────

"""
    make_pftclassification_from_rlist(r_list) -> Int (handle)

Convert an R list of `pft()` specs into a `PFTClassification` and store it.
Returns an integer handle safe to pass back through JuliaCall.
"""
function make_pftclassification_from_rlist(r_list)
    pfts = transform_pftlist_to_structure(r_list)
    obj  = PFTClassification(pfts)
    return _store!(obj)
end

"""
    make_biome4_pftclassification(; T, U) -> Int (handle)

Create the default BIOME4 PFTClassification (13 PFTs) and store it.
Returns an integer handle.
"""
function make_biome4_pftclassification(; T::Type=Float64, U::Type=Int)
    obj = BIOME4.PFTClassification{T,U}()
    return _store!(obj)
end


# ─────────────────────────────────────────────────────────────────────────────
# PFT editing
# ─────────────────────────────────────────────────────────────────────────────

"""
    set_pft_characteristic(handle, pftname, field, value) -> nothing

Edit one constraint or characteristic field of a stored PFTClassification.

- handle  : integer returned by make_biome4_pftclassification / make_pftclassification_from_rlist
- pftname : String, e.g. "BorealEvergreen"
- field   : String or Symbol, e.g. "gdd5" or :tcm
- value   : length-2 numeric vector for constraints, scalar for other fields
"""
function set_pft_characteristic(handle::Int, pftname, field, value)
    pftlist = _retrieve(handle)
    fname   = field isa Symbol ? field : Symbol(String(field))
    Biome.set_characteristic!(pftlist, String(pftname), fname, value)
    return nothing       # explicit nothing → R NULL, no serialisation attempt
end

"""
    apply_pft_edits!(handle, edits) -> nothing

Apply a list of edits (each a dict with keys name/field/value) to a stored PFTClassification.
"""
function apply_pft_edits!(handle::Int, edits)
    pftlist   = _retrieve(handle)
    edits_vec = edits isa AbstractVector ? edits : collect(edits)
    for e in edits_vec
        d = _to_sym_dict(e)
        Biome.set_characteristic!(pftlist, String(d[:name]),
                                  Symbol(String(d[:field])), d[:value])
    end
    return nothing
end

"""
    get_pft_names(handle) -> Vector{String}

Return the names of all PFTs in a stored PFTClassification.
Safe to return to R (plain string vector).
"""
function get_pft_names(handle::Int)
    pftlist = _retrieve(handle)
    return [pft.characteristics.name for pft in pftlist.pft_list]
end

"""
    get_pft_constraint(handle, pftname, field) -> Vector{Float64}

Retrieve one constraint range from a stored PFTClassification.
Returns a length-2 Float64 vector, safe for R.
"""
function get_pft_constraint(handle::Int, pftname::String, field)
    pftlist = _retrieve(handle)
    fname   = field isa Symbol ? field : Symbol(String(field))
    idx     = findfirst(p -> p.characteristics.name == pftname, pftlist.pft_list)
    idx === nothing && error("No PFT named \"$(pftname)\"")
    cons = pftlist.pft_list[idx].characteristics.constraints
    hasproperty(cons, fname) || error("No constraint field :$(fname)")
    v = cons[fname]
    return [Float64(v[1]), Float64(v[2])]
end


# ─────────────────────────────────────────────────────────────────────────────
# Raster helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    raster_from_spec(spec) -> Rasters.Raster

Build a Raster from an R array-spec with keys: values, lon, lat, [fill_value].
"""
function raster_from_spec(spec)
    d = _to_sym_dict(spec)
    for k in (:values, :lon, :lat)
        haskey(d, k) || error("Raster spec missing `$(k)`; got keys: $(collect(keys(d)))")
    end

    A   = Float64.(d[:values])
    lon = Float64.(d[:lon])
    lat = Float64.(d[:lat])

    if haskey(d, :fill_value) && d[:fill_value] !== nothing
        fillv = Float64(d[:fill_value])
        A     = coalesce.(A, fillv)
    end

    nd = ndims(A)
    if nd == 2
        # R: (y, x) → Julia Raster: (x, y)
        return Raster(permutedims(A, (2, 1)), dims=(X(lon), Y(lat)))
    elseif nd == 3
        # R: (y, x, t) → Julia Raster: (x, y, t)
        Axyt = PermutedDimsArray(A, (2, 1, 3))
        return Raster(Axyt, dims=(X(lon), Y(lat), Dim{:time}(1:size(Axyt, 3))))
    else
        error("values must be 2D or 3D. Got ndims=$(nd).")
    end
end


# ─────────────────────────────────────────────────────────────────────────────
# Model runner
# ─────────────────────────────────────────────────────────────────────────────

"""
    make_bounds(lon_min, lon_max, lat_min, lat_max) -> Int (handle)

Build the DimensionalData spatial bounds tuple expected by `execute()` and
store it in the object registry. Returns an integer handle safe for R.

Example from R:
  bounds <- julia_call("BiomeRAPI.make_bounds", 5.0, 11.0, 45.0, 48.0)
"""
function make_bounds(lon_min::Real, lon_max::Real, lat_min::Real, lat_max::Real)
    bounds = (X(Float64(lon_min) .. Float64(lon_max)),
              Y(Float64(lat_min) .. Float64(lat_max)))
    return _store!(bounds)
end

"""
    run_from_r(; model, co2, pft_handle, pft_specs, rasters,
                 fill_value, bounds_handle, outfile, biome_assignment)
               -> String  (path to output NetCDF)

Main entry point called from R via JuliaCall.

- pft_handle    : integer handle from make_biome4_pftclassification() or
                  make_pftclassification_from_rlist()  [preferred]
- pft_specs     : raw R list of pft() specs            [alternative; ignored if pft_handle given]
- bounds_handle : integer handle from make_bounds()    [optional; runs globally if omitted]
"""
function run_from_r(; model::AbstractString            = "BIOME4Model",
                     co2::Real                         = 378.0,
                     pft_handle::Union{Int,Nothing}    = nothing,
                     pft_specs                         = nothing,
                     rasters,
                     fill_value::Real                  = -9999.0,
                     bounds_handle::Union{Int,Nothing} = nothing,
                     outfile::AbstractString           = "out.nc",
                     biome_assignment                  = nothing)

    # Resolve PFT list
    pftlist = if pft_handle !== nothing
        _retrieve(pft_handle)
    elseif pft_specs !== nothing
        make_pftclassification_from_rlist(pft_specs) |> _retrieve
    else
        nothing
    end

    # Build raster objects
    rd          = _to_sym_dict(rasters)
    raster_objs = Dict{Symbol, Any}(k => raster_from_spec(v) for (k, v) in rd)

    # Resolve model
    model_sym = Symbol(model)
    isdefined(Biome, model_sym) || error("Unknown model: $(model)")
    model_obj = getproperty(Biome, model_sym)()

    setup = ModelSetup(model_obj;
        co2              = co2,
        pftlist          = pftlist,
        biome_assignment = biome_assignment,
        fill_value       = fill_value,
        (; raster_objs...)...
    )

    bounds = bounds_handle !== nothing ? _retrieve(bounds_handle) : nothing
    execute(setup; bounds = bounds, outfile = String(outfile))
    return outfile   # plain String → safe for R
end


end # module BiomeRAPI