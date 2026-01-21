# R/pft.R

#' Create a Plant Functional Type (PFT) specification
#'
#' This constructs an R-side "PFT spec" that biomeR sends to Julia to build a
#' concrete Biome.jl PFT object.
#'
#' @param type Character. Julia constructor name, e.g. "BroadleafDeciduousPFT".
#' @param name Optional character name for this PFT.
#' @param constraints Optional named list of length-2 numeric vectors
#'   (e.g. list(gdd5 = c(750, 1200), tmin = c(-Inf, Inf))).
#' @param mean_val Optional list with keys clt, prec, temp (numeric).
#' @param sd_val Optional list with keys clt, prec, temp (numeric).
#' @param ... Additional parameters passed as keyword arguments to the Julia
#'   constructor (e.g. kk = 0.3, dominance_factor = 9).
#'
#' @return An object of class `BiomePFT`.
#' @export
pft <- function(type,
                name = NULL,
                constraints = NULL,
                mean_val = NULL,
                sd_val = NULL,
                ...) {

  if (!is.character(type) || length(type) != 1 || is.na(type) || type == "")
    stop("`type` must be a single, non-empty character string.", call. = FALSE)

  x <- list(
    type = type,
    name = name,
    constraints = constraints,
    mean_val = mean_val,
    sd_val = sd_val,
    params = list(...)
  )
  class(x) <- c("BiomePFT", "list")
  x
}

#' @export
print.BiomePFT <- function(x, ...) {
  cat("<BiomePFT>\n")
  cat("  type: ", x$type, "\n", sep = "")
  if (!is.null(x$name)) cat("  name: ", x$name, "\n", sep = "")
  if (!is.null(x$constraints)) cat("  constraints: ", paste(names(x$constraints), collapse = ", "), "\n", sep = "")
  if (!is.null(x$params) && length(x$params) > 0) cat("  params: ", paste(names(x$params), collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' Test if an object is a BiomePFT
#' @param x Any object
#' @return TRUE/FALSE
#' @export
is_pft <- function(x) inherits(x, "BiomePFT")


#' Create a BIOME4 PFT edit instruction
#'
#' Convenience helper to build edits you later apply to a BIOME4 PFTClassification.
#'
#' @param name Character. PFT name used on the Julia side (e.g. "BorealEvergreen").
#' @param field Character or symbol-like. Field to edit (e.g. "gdd5").
#' @param value Numeric scalar or length-2 numeric vector.
#'
#' @return A list with class `BiomePFTEdit`.
#' @export
pft_edit <- function(name, field, value) {
  if (!is.character(name) || length(name) != 1 || name == "")
    stop("`name` must be a single non-empty character string.", call. = FALSE)
  if (!(is.character(field) || is.symbol(field)) || length(field) != 1)
    stop("`field` must be a single string (recommended) or symbol.", call. = FALSE)

  x <- list(
    name = name,
    field = as.character(field),
    value = value
  )
  class(x) <- c("BiomePFTEdit", "list")
  x
}

#' @export
print.BiomePFTEdit <- function(x, ...) {
  cat("<BiomePFTEdit>\n")
  cat("  name : ", x$name, "\n", sep = "")
  cat("  field: ", x$field, "\n", sep = "")
  cat("  value: ", paste(x$value, collapse = ", "), "\n", sep = "")
  invisible(x)
}
