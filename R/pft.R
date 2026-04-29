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

drop_nulls <- function(x) {
  if (!is.list(x) || inherits(x, "data.frame")) return(x)
  x <- x[!vapply(x, is.null, logical(1))]
  for (nm in names(x)) {
    if (is.list(x[[nm]]) && !inherits(x[[nm]], "data.frame")) {
      x[[nm]] <- drop_nulls(x[[nm]])
    }
  }
  x
}

validate_pft <- function(x) {
  if (!is.list(x)) stop("PFT must be a list.", call. = FALSE)

  if (is.null(x$type) || !is.character(x$type) || length(x$type) != 1 || is.na(x$type) || x$type == "") {
    stop("`type` must be a single, non-empty character string.", call. = FALSE)
  }

  if (!is.null(x$name) && (!is.character(x$name) || length(x$name) != 1 || is.na(x$name) || x$name == "")) {
    stop("`name` must be NULL or a single, non-empty character string.", call. = FALSE)
  }

  if (!is.null(x$constraints)) {
    if (!is.list(x$constraints)) {
      stop("`constraints` must be a list.", call. = FALSE)
    }
    bad <- names(x$constraints)[!vapply(x$constraints, function(v) {
      is.numeric(v) && length(v) == 2
    }, logical(1))]
    if (length(bad)) {
      stop(
        "Each constraint must be a numeric vector of length 2: ",
        paste(bad, collapse = ", "),
        call. = FALSE
      )
    }
  }

  class(x) <- unique(c("BiomePFT", class(x)))
  x
}

#' Modify a PFT specification
#' @export
modify_pft <- function(pft, ...) {
  if (!inherits(pft, "BiomePFT")) {
    stop("`pft` must be created with `pft()`.", call. = FALSE)
  }

  updates <- list(...)
  if (is.null(names(updates)) || any(names(updates) == "")) {
    stop("All updates in `modify_pft()` must be named.", call. = FALSE)
  }

  for (nm in names(updates)) {
    if (nm %in% names(pft)) {
      pft[[nm]] <- updates[[nm]]
    } else {
      if (is.null(pft$params)) pft$params <- list()
      pft$params[[nm]] <- updates[[nm]]
    }
  }

  pft <- drop_nulls(pft)
  validate_pft(pft)
}

#' Build a list of PFT specifications
#' @export
pft_list <- function(...) {
  x <- list(...)
  if (!all(vapply(x, inherits, logical(1), "BiomePFT"))) {
    stop("All elements of `pft_list()` must be created with `pft()`.", call. = FALSE)
  }
  structure(x, class = c("BiomePFTList", "list"))
}

as_biome_pft <- function(x) {
  if (!inherits(x, "BiomePFT")) {
    stop("`x` must inherit from `BiomePFT`.", call. = FALSE)
  }
  drop_nulls(unclass(x))
}

as_biome_pft_list <- function(x) {
  if (inherits(x, "BiomePFT")) return(list(as_biome_pft(x)))
  if (!is.list(x)) stop("`x` must be a PFT or a list of PFTs.", call. = FALSE)

  lapply(x, function(z) {
    if (inherits(z, "BiomePFT")) return(as_biome_pft(z))
    if (is.list(z)) return(drop_nulls(z))
    stop("Invalid PFT specification.", call. = FALSE)
  })
}

#' Create the default BIOME4 PFT list
#'
#' Instantiates the 13 BIOME4 plant functional types on the Julia side and
#' returns an opaque handle (integer) that can be passed to \code{set_pft_constraint()}
#' and \code{run_biome()}.
#'
#' @param biome Result of \code{biome_setup()}.
#' @return An integer handle referencing the Julia-side PFTClassification.
#' @export
make_biome4_pftlist <- function(biome) {
  julia_call("BiomeRAPI.make_biome4_pftclassification")
}

#' Edit one constraint of a BIOME4 PFT
#'
#' @param pftlist Integer handle returned by \code{make_biome4_pftlist()}.
#' @param name Character. PFT name, e.g. \code{"CoolConifer"}.
#' @param field Character. Constraint field, e.g. \code{"gdd5"}, \code{"tcm"}.
#' @param value Numeric vector of length 2: \code{c(min, max)}.
#' @return Invisibly \code{NULL}.
#' @export
set_pft_constraint <- function(pftlist, name, field, value) {
  if (!is.numeric(pftlist) || length(pftlist) != 1L)
    stop("`pftlist` must be an integer handle from `make_biome4_pftlist()`.", call. = FALSE)
  if (!is.numeric(value) || length(value) != 2L)
    stop("`value` must be a numeric vector of length 2, e.g. c(min, max).", call. = FALSE)
  julia_call("BiomeRAPI.set_pft_characteristic",
             as.integer(pftlist), name, field, value)
  invisible(NULL)
}

#' Get the names of all PFTs in a PFT list
#'
#' @param pftlist Integer handle returned by \code{make_biome4_pftlist()}.
#' @return Character vector of PFT names.
#' @export
get_pft_names <- function(pftlist) {
  julia_call("BiomeRAPI.get_pft_names", as.integer(pftlist))
}

# Optional helper for applying BIOME4-style edits to a local PFT list
#' @export
apply_pft_edits <- function(pftlist, edits) {
  if (!is.list(edits)) stop("`edits` must be a list of BiomePFTEdit objects.", call. = FALSE)
  for (ed in edits) {
    if (!inherits(ed, "BiomePFTEdit")) {
      stop("Each edit must be created with `pft_edit()`.", call. = FALSE)
    }
    pftlist <- lapply(pftlist, function(p) {
      if (!is.list(p)) return(p)
      if (!is.null(p$name) && identical(p$name, ed$name)) {
        if (is.null(p$params)) p$params <- list()
        p$params[[ed$field]] <- ed$value
      }
      p
    })
  }
  pftlist
}