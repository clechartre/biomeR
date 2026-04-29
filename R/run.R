# R/run.R

#' Define a spatial bounding box for model execution
#'
#' Creates a spatial bounds object on the Julia side that can be passed to
#' \code{run_biome()} to restrict execution to a geographic subset.
#'
#' @param lon_min Minimum longitude.
#' @param lon_max Maximum longitude.
#' @param lat_min Minimum latitude.
#' @param lat_max Maximum latitude.
#' @return An integer handle referencing the Julia-side bounds tuple.
#' @export
bbox <- function(lon_min, lon_max, lat_min, lat_max) {
  julia_call("BiomeRAPI.make_bounds",
             as.double(lon_min), as.double(lon_max),
             as.double(lat_min), as.double(lat_max))
}

#' Run Biome.jl from R
#'
#' @param rasters Named list of Raster* objects (from the raster package).
#' @param model Julia model constructor name, e.g. "BaseModel" or "BIOME4Model".
#' @param co2 Numeric CO2 concentration.
#' @param pft_specs Optional: R list of PFT specs created with \code{pft()}.
#' @param pftlist Optional: integer handle from \code{make_biome4_pftlist()}.
#'   Takes precedence over \code{pft_specs} if provided.
#' @param outfile Output NetCDF file path.
#' @param bounds Optional: integer handle from \code{bbox()} to spatially subset the run.
#' @param fill_value Fill value used for NA replacement. Default -9999.
#'
#' @return \code{outfile} invisibly.
#' @export
run_biome <- function(rasters,
                      model      = "BaseModel",
                      co2        = 378.0,
                      pft_specs  = NULL,
                      pftlist    = NULL,
                      outfile    = "out.nc",
                      bounds     = NULL,
                      fill_value = -9999.0) {

  # --- basic checks ---
  if (is.null(rasters) || !is.list(rasters) || is.null(names(rasters)))
    stop("`rasters` must be a named list, e.g. list(temp=..., prec=..., clt=..., whc=..., ksat=...).")
  if (!("temp" %in% names(rasters))) stop("`rasters` must include a `temp` raster.")
  if (!("prec" %in% names(rasters))) stop("`rasters` must include a `prec` raster.")

  # --- convert R rasters to array specs ---
  rasters_spec <- lapply(rasters, r_to_spec, fill_value = fill_value)

  # --- build args ---
  args <- list(
    model         = model,
    co2           = co2,
    rasters       = rasters_spec,
    bounds_handle = if (!is.null(bounds)) as.integer(bounds) else NULL,
    outfile       = outfile,
    fill_value    = fill_value
  )

  if (!is.null(pftlist)) {
    args$pft_handle <- as.integer(pftlist)
  } else if (!is.null(pft_specs)) {
    args$pft_specs <- as_biome_pft_list(pft_specs)
  }

  # --- call Julia ---
  invisible(do.call(julia_call, c(list("BiomeRAPI.run_from_r"), args)))
}