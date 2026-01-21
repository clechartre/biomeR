# R/biome_setup.R

# internal helper: generate a callable wrapper for one Julia function
# The two functions julia_function and julia_pkg_import 
# have been taken from diffeqr, the R API for DifferentialEquations.jl
julia_function <- function(func_name, pkg_name = "Main",
                           env = emptyenv()){
  fname <- paste0(pkg_name, ".", func_name)
  force(fname)
  f <- function(...,
                need_return = c("R", "Julia", "None"),
                show_value = FALSE){
    if (!isTRUE(env$initialized)) {
      env$setup()
    }
    JuliaCall::julia_do.call(func_name = fname, list(...),
                             need_return = match.arg(need_return),
                             show_value = show_value)
  }
  force(f)
  env[[func_name]] <- f
}

julia_pkg_import <- function(pkg_name, func_list){
  env <- new.env(parent = emptyenv())
  env$setup <- function(...){
    JuliaCall::julia_setup(...)
    JuliaCall::julia_library(pkg_name)
    env$initialized <- TRUE
  }
  for (fname in func_list) {
    julia_function(func_name = fname,
                   pkg_name = pkg_name,
                   env = env)
  }
  env
}

#' Initialize Julia + Biome backend for biomeR
#'
#' @param project_dir Path to the Julia project to activate (Biome.jl dev checkout).
#' @param julia_bin Optional path to Julia home (JULIA_HOME)
#' @param installJulia Whether to let JuliaCall install Julia
#' @param pkg_check Whether to run Pkg.instantiate()
#' @return A list with $api and (optionally) $biome handles.
#' @export
biome_setup <- function(project_dir = NULL,
                        julia_bin = NULL,
                        installJulia = FALSE,
                        pkg_check = TRUE,
                        ...) {

  # 1) start Julia
  JuliaCall::julia_setup(
    installJulia = installJulia,
    JULIA_HOME = julia_bin,
    ...
  )

  # 2) activate project (dev-mode)
  if (!is.null(project_dir)) {
    JuliaCall::julia_eval(sprintf(
      'using Pkg; Pkg.activate(raw"%s"); %s',
      project_dir,
      if (pkg_check) "Pkg.instantiate()" else "nothing"
    ))
  }

  # 3) load Biome
  JuliaCall::julia_library("Biome")

  # 4) load our Julia API shim shipped with biomeR
  api_file <- system.file("julia", "r_api.jl", package = "biomeR")
  if (api_file == "") {
    stop("Could not find inst/julia/r_api.jl in the installed biomeR package.")
  }
  JuliaCall::julia_eval(sprintf('include(raw"%s"); using .BiomeRAPI', api_file))

  # 5) create R handles to call Julia functions
  env_api <- new.env(parent = emptyenv())
  env_api$initialized <- TRUE
  env_api$setup <- function() invisible(TRUE)  # already setup above

  api_functions <- c(
    "make_biome4_pftclassification",
    "set_pft_characteristic",
    "apply_pft_edits_bang",
    "raster_from_spec",
    "make_modelsetup_from_rasterspecs",
    "run_from_r"
  )
  api <- julia_pkg_import("Main.BiomeRAPI", api_functions)

  # (Optional) expose a small subset of Biome constructors if you want
  # You can also skip this entirely and keep only $api.
  env_biome <- new.env(parent = emptyenv())
  env_biome$initialized <- TRUE
  env_biome$setup <- function() invisible(TRUE)

  biome_functions <- c(
    "BaseModel",
    "BIOME4Model",
    "BIOMEDominanceModel",
    "KoppenModel",
    "WissmannModel", 
    "ThornthwaiteModel", 
    "TrollPfaffenModel"
  )
  biome <- julia_pkg_import("Biome", biome_functions)

  list(api = api, biome = biome)
}
