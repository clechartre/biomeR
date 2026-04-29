#' Assign a biome from a score matrix
#' @export
assign_biome <- function(scores, ...) {
  UseMethod("assign_biome")
}

#' @export
assign_biome.default <- function(scores, pft_names = NULL, ties.method = "first", ...) {
  m <- as.matrix(scores)
  idx <- max.col(m, ties.method = ties.method)

  if (!is.null(pft_names)) return(pft_names[idx])
  if (!is.null(colnames(m))) return(colnames(m)[idx])
  idx
}

#' @export
assign_biome.function <- function(scores, fun, ...) {
  fun(scores, ...)
}