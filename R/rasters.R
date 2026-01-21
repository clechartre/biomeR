# R/rasters.R

#' Handle Rasters to pass them on the the model setup.
#'
#' This starts from a RasterLayer, RasterBrick, and RasterSack and formats it properly
#'
#' @param r Input raster
#' @param fill_value Value to replace missing data. Default to -9999
#' @return Raster value and metadata as a list
#' @export
r_to_spec <- function(r, fill_value = -9999) {
  nx <- ncol(r)
  ny <- nrow(r)
  nt <- raster::nlayers(r)

  # This preserves the raster's row/col layout (and bands, if any)
  if (nt == 1) {
    A <- raster::as.matrix(r)              # [ny, nx]
    A[is.na(A)] <- fill_value
  } else {
    A <- raster::as.array(r)               # [ny, nx, nt]
    A[is.na(A)] <- fill_value
  }

  lon <- raster::xFromCol(r, 1:nx)
  lat <- raster::yFromRow(r, 1:ny)

  list(values = A, lon = lon, lat = lat, fill_value = fill_value)
}
