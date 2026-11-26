# Needs from helper.R: rng_state(), habitat_ids()
library(here)
library(data.table)
library(raster)
library(purrr)
library(testthat)

skip_if_not_installed("raster")
skip_if_not_installed("withr")
skip_if_not_installed("checkmate")

source(here::here("Model", "src", "habitat_cells.R"))
source(here::here("tests", "testthat", "helper.R"))
source(here::here("Model", "parameters.R"))
list.files(here::here("Model", "src"), full.names = TRUE) %>%
  walk(source)

# ---- fixtures ----------------------------------------------------------------

mod_par$steps_pre_frag <- 4
mod_par$steps_post_frag <- 6
states <- single_run(
  mod_par = mod_par,
  var_par = var_par[1, ],
  species_par = species_par,
  switch = switch,
  sim_id = "test",
  seeds = derive_seeds(origin = 42L, which = "model"),
  record_steps = c("pre_fragmentation", "post_fragmentation")
)

st_full <- states$pre_fragmentation
st_frag <- states$post_fragmentation

ls_full <- st_full$grid
ls_frag <- st_frag$grid

hab_full <- habitat_ids(ls_full)
hab_frag <- habitat_ids(ls_frag)

# ---- input handling ------------------------------------------------------------

test_that("RasterLayer and model_state give the same result", {
  expect_identical(habitat_cells(st_full), habitat_cells(ls_full))
  expect_identical(habitat_cells(st_full, 10, seed = 1L),
                   habitat_cells(ls_full, 10, seed = 1L))
})

test_that("other inputs are rejected", {
  expect_error(habitat_cells(matrix(1, 3, 3)), "RasterLayer")
  expect_error(habitat_cells(NULL), "RasterLayer")
  expect_error(habitat_cells(list(grid = ls_full)), "RasterLayer")
  expect_error(habitat_cells(states), "RasterLayer")
  expect_error(habitat_cells("grid"), "RasterLayer")
})

# ---- all habitat cells (n_cells = Inf) -----------------------------------------

test_that("default returns all habitat cells, ascending, without a seed", {
  out_full <- habitat_cells(ls_full)
  expect_type(out_full, "integer")
  expect_identical(as.integer(out_full), as.integer(hab_full))
  expect_false(is.unsorted(out_full))
  expect_null(attr(out_full, "seed", exact = TRUE))
  expect_identical(habitat_cells(ls_full, n_cells = Inf), out_full)
})

test_that("a continuous landscape returns all cells", {
  expect_identical(as.integer(habitat_cells(ls_full)), seq_len(raster::ncell(ls_full)))
})

test_that("the correct number of habitat cells is returned", {
  expect_identical(length(habitat_cells(ls_full)), length(hab_full))
  expect_identical(length(habitat_cells(ls_frag)), length(hab_frag))
})

test_that("habitat is derived from the grid, not from a cached cells element", {
  stale <- st_full
  stale$cells$habitat <- integer(0)   # empty, but the grid is unchanged
  expect_identical(as.integer(habitat_cells(stale)), as.integer(habitat_ids(ls_full)))
})

# ---- random draws ----------------------------------------------------------------

test_that("a draw returns n distinct habitat cells", {
  out <- habitat_cells(st_frag, n_cells = 10, seed = 7L)
  expect_type(out, "integer")
  expect_length(out, 10)
  expect_false(anyDuplicated(out) > 0)
  expect_true(all(out %in% hab_frag))
})

test_that("n_cells equal to the habitat count is a draw like any other", {
  out <- habitat_cells(st_frag, n_cells = length(hab_frag), seed = 1L)
  expect_setequal(out, hab_frag)
  expect_true(is.unsorted(out))        # random order, not ascending
  expect_identical(attr(out, "seed", exact = TRUE), 1L)
})

test_that("a single habitat cell is returned as itself", {
  # Regression test: sample(x, n) with a length-one x samples from 1:x
  one <- raster::raster(nrows = 10, ncols = 10)
  one[5] <- 1
  expect_identical(as.integer(habitat_cells(one, n_cells = 1, seed = 1L)), 5L)
  expect_identical(as.integer(habitat_cells(one)), 5L)
})

test_that("draws only ever return habitat cells on a sparse grid", {
  sparse <- raster::raster(nrows = 10, ncols = 10)
  sparse[c(4, 25, 61, 99)] <- 1
  for (s in 1:20) {
    expect_true(all(habitat_cells(sparse, n_cells = 2, seed = s) %in% c(4L, 25L, 61L, 99L)))
  }
})

# ---- seeding -----------------------------------------------------------------------

test_that("same seed gives the same draw, different seeds do not", {
  a <- habitat_cells(st_frag, 10, seed = 5L)
  b <- habitat_cells(st_frag, 10, seed = 5L)
  c <- habitat_cells(st_frag, 10, seed = 6L)
  expect_identical(a, b)
  expect_false(identical(sort(as.integer(a)), sort(as.integer(c))))
})

test_that("the seed attribute is the integer seed that was used", {
  expect_identical(attr(habitat_cells(st_frag, 10, seed = 5L), "seed", exact = TRUE), 5L)
  expect_identical(attr(habitat_cells(st_frag, 10, seed = 5), "seed", exact = TRUE), 5L)   # double input
})

test_that("a seedless draw records a fresh seed that reproduces it", {
  a <- habitat_cells(st_frag, 10)
  s <- attr(a, "seed", exact = TRUE)
  expect_type(s, "integer")
  expect_false(is.na(s))

  b <- habitat_cells(st_frag, 10, seed = s)
  expect_identical(a, b)
})

test_that("a seedless draw follows the global seed", {
  set.seed(1); a <- habitat_cells(st_frag, 10)
  set.seed(1); b <- habitat_cells(st_frag, 10)
  expect_identical(a, b)

  # consecutive seedless draws advance the stream
  set.seed(1)
  x <- habitat_cells(st_frag, 10)
  y <- habitat_cells(st_frag, 10)
  expect_false(identical(x, y))
  expect_false(identical(attr(x, "seed", exact = TRUE), attr(y, "seed", exact = TRUE)))
})

test_that("a seeded draw neither consumes nor depends on the global RNG", {
  set.seed(1)
  before <- rng_state()
  a <- habitat_cells(st_frag, 10, seed = 5L)
  expect_identical(rng_state(), before)

  set.seed(999)
  b <- habitat_cells(st_frag, 10, seed = 5L)
  expect_identical(a, b)
})

test_that("returning all habitat cells does not touch the global RNG", {
  set.seed(1)
  before <- rng_state()
  habitat_cells(st_frag)
  expect_identical(rng_state(), before)
})

# ---- validation ------------------------------------------------------------------------

test_that("n_cells above the number of habitat cells is an error", {
  expect_error(habitat_cells(st_full, n_cells = length(hab_full) + 1), "cannot exceed")
  expect_error(habitat_cells(st_frag, n_cells = length(hab_full)), "cannot exceed")
})

test_that("invalid n_cells values are rejected", {
  expect_error(habitat_cells(st_frag, n_cells = 0))
  expect_error(habitat_cells(st_frag, n_cells = -3))
  expect_error(habitat_cells(st_frag, n_cells = NA))
  expect_error(habitat_cells(st_frag, n_cells = NULL))          # NULL must not silently mean "all"
  expect_error(habitat_cells(st_frag, n_cells = 2.5))           # needs the whole-number check
  expect_error(habitat_cells(st_frag, n_cells = "10"))
  expect_error(habitat_cells(st_frag, n_cells = c(5, 10)))
})

# ---- behaviour across states with different habitat -------------------------------------

test_that("the same seed does not give the same cells in grids with different habitat", {
  a <- habitat_cells(ls_full, 10, seed = 1L)
  b <- habitat_cells(ls_frag, 10, seed = 1L)
  expect_false(identical(sort(as.integer(a)), sort(as.integer(b))))
})

test_that("a design drawn from the later grid is valid in the earlier grid", {
  design <- habitat_cells(ls_frag, 10, seed = 1L)
  expect_true(all(design %in% habitat_ids(ls_frag)))
  expect_true(all(design %in% habitat_cells(ls_full)))
})

# ---- attribute behaviour (documented, not desired) ---------------------------------------

test_that("subsetting a draw drops the seed attribute", {
  design <- habitat_cells(st_frag, 10, seed = 3L)
  expect_identical(attr(design, "seed", exact = TRUE), 3L)
  expect_null(attr(design[1:5], "seed", exact = TRUE))
})
