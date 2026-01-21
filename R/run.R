# R/run.R

#' Run Biome.jl from R
#'
#' High-level wrapper that:
#'  1) ensures Julia + BiomeRAPI are initialized (via biome_setup)
#'  2) converts R rasters to array "specs"
#'  3) calls BiomeRAPI.run_from_r() in Julia
#'
#' @param rasters Named list of Raster* objects (from the raster package).
#'   Required names: at least temp, prec. Common: clt, whc, ksat.
#'   Examples:
#'     list(temp=temp_brick, prec=prec_brick, clt=clt_brick, whc=whc_layer, ksat=ksat_layer)
#' @param model Julia model constructor name as string. e.g. "BaseModel", "BIOME4Model".
#' @param co2 Numeric CO2 concentration.
#' @param pft_specs Optional: R list of PFT specs created with pft(), e.g. list(pft(...), pft(...)).
#' @param pftlist Optional: a Julia PFTClassification object (e.g. from make_biome4_pftclassification()).
#'   If provided, it takes precedence over pft_specs.
#' @param outfile Output NetCDF path.
#' @param coordstring Coordinate bounding box string or "alldata".
#' @param fill_value Fill value used for NA -> numeric.
#' @param biome Optional: result from biome_setup(). If NULL, will call biome_setup().
#' @param project_dir Path to the Biome.jl project directory (needed if biome is NULL).
#' @param ... Passed through to biome_setup() (e.g. julia_bin, installJulia, pkg_check).
#'
#' @return outfile (invisibly), as a string.
#' @export
run_biome <- function(rasters,
                      model = "BaseModel",
                      co2 = 378.0,
                      pft_specs = NULL,
                      pftlist = NULL,
                      outfile = "out.nc",
                      coordstring = "alldata",
                      fill_value = -9999.0,
                      biome = NULL,
                      project_dir = NULL,
                      ...) {

  # --- basic checks ---
  if (is.null(rasters) || !is.list(rasters) || is.null(names(rasters))) {
    stop("`rasters` must be a *named* list, e.g. list(temp=..., prec=..., clt=..., whc=..., ksat=...).")
  }
  if (!("temp" %in% names(rasters))) stop("`rasters` must include a `temp` raster.")
  if (!("prec" %in% names(rasters))) stop("`rasters` must include a `prec` raster.")

  # --- ensure Julia/Biome is initialized ---
  if (is.null(biome)) {
    if (is.null(project_dir)) {
      stop("Provide either `biome` (result of biome_setup()) or `project_dir` (Biome.jl project path).")
    }
    biome <- biome_setup(project_dir = project_dir, ...)
  }

  # --- convert R rasters -> specs expected by Julia ---
  rasters_spec <- lapply(rasters, r_to_spec, fill_value = fill_value)

  # --- choose PFT input mode ---
  # If pftlist is provided, pass it through as-is.
  # Else pass pft_specs (R list of pft() objects) or NULL.
  args <- list(
    model = model,
    co2 = co2,
    rasters = rasters_spec,
    coordstring = coordstring,
    outfile = outfile,
    fill_value = fill_value
  )

  if (!is.null(pftlist)) {
    args$pftlist <- pftlist
  } else {
    args$pft_specs <- pft_specs
  }

  # --- call Julia ---
  # Julia_function wrapper uses need_return. We want the outfile string back.
  res <- do.call(biome$api$run_from_r, c(args, list(need_return = "R")))

  invisible(res)
}
