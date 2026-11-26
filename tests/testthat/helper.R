rng_state <- function() {
  if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  }
}

habitat_ids <- function(rast) {
  which(!is.na(raster::getValues(rast)))
}