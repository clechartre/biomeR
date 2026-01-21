module BiomeRAPI

using Biome
using Rasters
using DimensionalData

# Convert keys to Symbols, leave values unchanged
_to_sym_dict(x) = x

function _to_sym_dict(d::AbstractDict)
    out = Dict{Symbol,Any}()
    for (k,v) in d
        ks = k isa Symbol ? k : Symbol(String(k))
        out[ks] = v
    end
    return out
end

# Build a NamedTuple from a Dict{Symbol,Any} at runtime
_namedtuple_from_sym_dict(d::AbstractDict{Symbol,Any}) = (; (k => d[k] for k in keys(d))...)

# Normalize "constraints" format: values should be length-2 vectors of reals
function _normalize_constraints(x)
    x === nothing && return nothing
    d = _to_sym_dict(x)
    out = Dict{Symbol,Any}()

    for (k,v) in d
        if v === nothing
            continue
        end
        # R typically provides c(a,b) -> Vector{<:Real}
        if v isa AbstractVector
            if length(v) != 2
                error("Constraint $(k) must be a length-2 vector like c(min,max). Got length=$(length(v)).")
            end
            out[k] = [Float64(v[1]), Float64(v[2])]
        else
            error("Constraint $(k) must be a vector like c(min,max). Got $(typeof(v)).")
        end
    end

    return _namedtuple_from_sym_dict(out)
end

# Normalize mean_val/sd_val lists: expect clt/prec/temp keys
function _normalize_triplet(x, label::String)
    x === nothing && return nothing
    d = _to_sym_dict(x)

    # accept either Symbol or String keys originally
    for needed in (:clt, :prec, :temp)
        haskey(d, needed) || error("$(label) must contain keys clt, prec, temp. Missing: $(needed)")
    end

    return (clt = Float64(d[:clt]),
            prec = Float64(d[:prec]),
            temp = Float64(d[:temp]))
end

# Extract a string field from a dict-ish object
function _get_string(d::AbstractDict{Symbol,Any}, key::Symbol; required::Bool=false)
    if !haskey(d, key) || d[key] === nothing
        required && error("Missing required field: $(key)")
        return nothing
    end
    return String(d[key])
end


"""
    transform_pft_to_structure(spec) -> AbstractPFT

Convert one R `pft()` spec (passed via JuliaCall) into a concrete Julia PFT by
calling the constructor named by `spec[:type]` in the `Biome` module.

Expected fields (from your R `pft()`):
- type (String): e.g. "BroadleafDeciduousPFT"
- name (String or NULL)
- constraints (list or NULL)
- mean_val (list or NULL)
- sd_val (list or NULL)
- params (list): named list of additional kwargs
"""
function transform_pft_to_structure(spec)
    d0 = _to_sym_dict(spec)

    typestr = _get_string(d0, :type; required=true)
    pft_ctor_sym = Symbol(typestr)

    # Grab the constructor function from Biome
    if !isdefined(Biome, pft_ctor_sym)
        error("Biome has no constructor named $(typestr). Check spelling / exports.")
    end
    ctor = getproperty(Biome, pft_ctor_sym)

    # Collect kwargs from params (R ... stored there)
    kwargs = Dict{Symbol,Any}()

    if haskey(d0, :params) && d0[:params] !== nothing
        pd = _to_sym_dict(d0[:params])
        for (k,v) in pd
            # Leave v as-is; Julia will type it
            kwargs[k] = v
        end
    end

    # Top-level `name` (your pft() stores name separately)
    if haskey(d0, :name) && d0[:name] !== nothing && !haskey(kwargs, :name)
        kwargs[:name] = String(d0[:name])
    end

    # constraints / mean_val / sd_val
    if haskey(d0, :constraints) && d0[:constraints] !== nothing
        kwargs[:constraints] = _normalize_constraints(d0[:constraints])
    end
    if haskey(d0, :mean_val) && d0[:mean_val] !== nothing
        kwargs[:mean_val] = _normalize_triplet(d0[:mean_val], "mean_val")
    end
    if haskey(d0, :sd_val) && d0[:sd_val] !== nothing
        kwargs[:sd_val] = _normalize_triplet(d0[:sd_val], "sd_val")
    end

    # Call constructor with keyword args
    # Convert Dict -> kwargs splat: (; kwargs...) needs a NamedTuple
    nt = _namedtuple_from_sym_dict(kwargs)
    return ctor(; nt...)
end

"""
    transform_pftlist_to_structure(r_list) -> Vector{AbstractPFT}

Convert an R list of `pft()` specs into a Julia vector of concrete PFT objects.
"""
function transform_pftlist_to_structure(r_list)
    # JuliaCall typically passes an R list-of-objects as Vector{Any}
    specs = r_list isa AbstractVector ? r_list : collect(r_list)
    return [transform_pft_to_structure(s) for s in specs]
end

"""
    make_pftclassification_from_rlist(r_list) -> PFTClassification

Convert an R list of `pft()` specs into a `PFTClassification`.
"""
function make_pftclassification_from_rlist(r_list)
    pfts = transform_pftlist_to_structure(r_list)
    return PFTClassification(pfts)
end

"""
    make_biome4_pftclassification(; T=Float64, U=Int)

Create a mutable BIOME4 PFTClassification with explicit numeric types.
"""
function make_biome4_pftclassification(; T::Type=Float64, U::Type=Int)
    return BIOME4.PFTClassification{T,U}()
end

"""
    apply_pft_edits!(pftlist, edits)

edits: Vector of dict-like objects from R, each with:
- name  : "BorealEvergreen"
- field : "gdd5" (or :gdd5)
- value : e.g. c(750,1200)
"""
function apply_pft_edits!(pftlist, edits)
    edits_vec = edits isa AbstractVector ? edits : collect(edits)
    for e in edits_vec
        d = _to_sym_dict(e)
        name  = String(d[:name])
        field = d[:field]
        value = d[:value]
        set_pft_characteristic(pftlist, name, field, value)
    end
    return pftlist
end



"""
    set_pft_characteristic(pftlist, pftname, field, value)

Convenience wrapper around internal `set_characteristic!` or direct mutation.

- pftname: String like "BorealEvergreen"
- field: Symbol or String like :gdd5
- value: usually a 2-element vector [min, max] or a scalar
"""
function set_pft_characteristic(pftlist, pftname, field, value)
    fname = field isa Symbol ? field : Symbol(String(field))
    return Biome.set_characteristic!(pftlist, String(pftname), fname, value)
end


"""
    _nc_raster(path, var; name=var)

Load a NetCDF variable as a Rasters.jl Raster.
"""
function _nc_raster(path::AbstractString, var::AbstractString; name::AbstractString=var)
    return Raster(path, name=var)  # Rasters.jl uses `name` to select the variable in NetCDF
end

"""
    make_modelsetup_from_rasterspecs(; model, co2, pft_specs, rasters, fill_value, biome_assignment)

rasters: dict-like with keys temp/prec/clt/whc/ksat... and each value is an array-spec
(values/lon/lat[/fill_value]).
"""
function make_modelsetup_from_rasterspecs(; model::AbstractString="BaseModel",
                                         co2::Real=378.0,
                                         pft_specs=nothing,
                                         pftlist=nothing,
                                         rasters,
                                         fill_value::Real=-9999.0,
                                         biome_assignment=nothing)

    pftlist = pftlist !== nothing ? pftlist :
          (pft_specs === nothing ? nothing : make_pftclassification_from_rlist(pft_specs))


    rd = _to_sym_dict(rasters)

    raster_objs = Dict{Symbol,Any}()
    for (k, spec) in rd
        raster_objs[k] = raster_from_spec(spec)
    end

    model_sym = Symbol(model)
    !isdefined(Biome, model_sym) && error("Unknown model $(model).")
    model_obj = getproperty(Biome, model_sym)()

    return ModelSetup(model_obj;
        co2=co2,
        pftlist=pftlist,
        biome_assignment=biome_assignment,
        fill_value=fill_value,
        (; raster_objs...)...
    )
end


"""
    raster_from_spec(spec) -> Rasters.Raster

Build a Rasters.jl Raster ONLY from an R "array spec":
spec must contain:
- values : matrix (ny × nx)
- lon    : vector length nx
- lat    : vector length ny
optional:
- fill_value
"""
function raster_from_spec(spec)
    d = _to_sym_dict(spec)

    for k in (:values, :lon, :lat)
        haskey(d, k) || error("Raster spec missing `$(k)`; got keys: $(collect(keys(d)))")
    end

    A   = Float64.(d[:values])   # could be 2D or 3D
    lon = Float64.(d[:lon])
    lat = Float64.(d[:lat])

    # Optional fill
    if haskey(d, :fill_value) && d[:fill_value] !== nothing
        fillv = Float64(d[:fill_value])
        A = coalesce.(A, fillv)
    end

    nd = ndims(A)
    if nd == 2
        # R: (y, x) -> Julia Raster expects (x, y)
        Axy = permutedims(A, (2, 1))
        return Raster(Axy, dims=(X(lon), Y(lat)))

    elseif nd == 3
        # A is (y, x, t) from R
        Axyt = PermutedDimsArray(A, (2, 1, 3))  # (x, y, t) VIEW, no copy
        tdim = Dim{:time}(1:size(Axyt, 3))
        return Raster(Axyt, dims=(X(lon), Y(lat), tdim))
    else
        error("values must be 2D or 3D. Got ndims=$(nd).")
    end
end


function run_from_r(; model::AbstractString="BIOME4Model",
                    co2::Real=378.0,
                    pft_specs=nothing,
                    pftlist=nothing,
                    rasters,
                    fill_value::Real=-9999.0,
                    coordstring::AbstractString="alldata",
                    outfile::AbstractString="out.nc",
                    biome_assignment=nothing)

    # build pftlist if not provided
    if pftlist === nothing
        pftlist = pft_specs === nothing ? nothing : make_pftclassification_from_rlist(pft_specs)
    end

    setup = make_modelsetup_from_rasterspecs(
        model=model, co2=co2, pft_specs=nothing,  # ignore specs now
        rasters=rasters, fill_value=fill_value,
        biome_assignment=biome_assignment
    )

    # overwrite setup.pftlist with our editable BIOME4 list
    setup.pftlist = pftlist

    run!(setup; coordstring=String(coordstring), outfile=String(outfile))
    return outfile
end



end # module