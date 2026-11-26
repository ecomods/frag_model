# Needs from helper.R: rng_state(), habitat_ids()
skip_if_not_installed("raster")
skip_if_not_installed("data.table")
skip_if_not_installed("withr")
library(here)
library(data.table)
library(dplyr)
library(purrr)
library(raster)
library(testthat)
library(withr)

source(here::here("R", "sample_cells.R"))
source(here::here("R", "toroidal_clump.R"))
source(here::here("tests", "testthat", "helper.R"))
source(here::here("Model", "parameters.R"))
list.files(here::here("Model", "src"), full.names = TRUE) %>%
  walk(source)


mod_par$steps_pre_frag <- 4
mod_par$steps_post_frag <- 6
states <- single_run(
  mod_par = mod_par,
  var_par = var_par[1, ],
  species_par = species_par,
  switch = switch,
  sim_id = "test",
  seeds = derive_seeds(origin = 42L, which = "model"),
  record_steps = c("pre_fragmentation", "post_fragmentation", "final")
)

# ---- argument validation -----------------------------------------------------

test_that("invalid method and format are rejected", {
  st <- states$pre_fragmentation
  expect_error(sample_cells(st, method = "chessboard"))
  expect_error(sample_cells(st, format = "tidy"))
})

test_that("random sampling requires n_samples", {
  st <- states$pre_fragmentation
  expect_error(sample_cells(st, "random", seed = 1L), "n_samples must be provided")
})

test_that("n_samples larger than the habitat gives a clear error", {
  # Habitat shrinks after the modification, so a fixed n_samples can exceed it.
  # The roxygen text assumes an error; adjust both if the function should cap and warn.
  st <- fragment(states$pre_fragmentation, habitat = 0.01, fragmentation = 0.5, seed = 1)
  too_large <- length(habitat_ids(st$grid)) + 1L
  expect_error(sample_cells(st, "random", n_samples = too_large), "n_samples cannot exceed")
})

test_that("n_samples is not used (and warns) for other methods", {
  st <- states$pre_fragmentation
  expect_warning(a <- sample_cells(st, "all", n_samples = 5L), "n_samples")
  expect_warning(b <- sample_cells(st, "checkerboard", n_samples = 5L), "n_samples")
  expect_identical(a, sample_cells(st, "all"))
  expect_identical(b, sample_cells(st, "checkerboard"))
})

test_that("n_samples equal to the habitat count returns all cells as a random design", {
  st <- states$final
  hab <- habitat_ids(st$grid)
  out <- sample_cells(st, "random", n_samples = length(hab), seed = 1L)
  expect_equal(nrow(out), length(hab))
  expect_setequal(out$cell_id, hab)
  expect_true(all(out$samp_method == "random"))
  expect_identical(unique(out$seed_sampling), 1L)
})

# ---- sampling methods --------------------------------------------------------

test_that("method 'all' returns every habitat cell once", {
  st <- states$post_fragmentation
  out <- sample_cells(st, "all")
  expect_setequal(out$cell_id, habitat_ids(st$grid))
  expect_equal(nrow(out), length(habitat_ids(st$grid)))
  expect_equal(out$sample_id, seq_len(nrow(out)))
  expect_true(all(out$samp_method == "all"))
})

test_that("method 'checkerboard' picks habitat cells of even row + column parity", {
  st <- states$post_fragmentation
  out <- sample_cells(st, "checkerboard")
  expect_true(all(out$cell_id %in% habitat_ids(st$grid)))
  expect_true(all((out$x_loc + out$y_loc) %% 2 == 0))

  rc <- raster::rowColFromCell(st$grid, habitat_ids(st$grid))
  expect_equal(nrow(out), sum((rc[, 1] + rc[, 2]) %% 2 == 0))
})

test_that("deterministic methods ignore the seed", {
  st <- states$final
  selection <- head(habitat_ids(st$grid), 10)
  expect_identical(sample_cells(st, "all", seed = 1L), sample_cells(st, "all", seed = 2L))
  expect_identical(sample_cells(st, "checkerboard", seed = 1L),
                   sample_cells(st, "checkerboard", seed = 2L))
  expect_identical(sample_cells(st, "fixed", cell_ids = selection, seed = 1L),
                   sample_cells(st, "fixed", cell_ids = selection, seed = 2L))
})

test_that("method 'random' returns n distinct habitat cells", {
  st <- states$final
  out <- sample_cells(st, "random", n_samples = 15L, seed = 7L)
  expect_equal(nrow(out), 15L)
  expect_false(anyDuplicated(out$cell_id) > 0)
  expect_true(all(out$cell_id %in% habitat_ids(st$grid)))
})

# ---- method 'fixed' -----------------------------------------------------------

test_that("method 'fixed' returns exactly the supplied cells in the supplied order", {
  st <- states$post_fragmentation
  hab <- habitat_ids(st$grid)
  selection <- rev(hab[round(seq(1, length(hab), length.out = 4))])
  out <- sample_cells(st, "fixed", cell_ids = selection)
  expect_equal(out$cell_id, selection)
  expect_equal(out$sample_id, seq_along(selection))
  expect_true(all(out$samp_method == "fixed"))
})

test_that("method 'fixed' requires cell_ids, and non-null cell_ids raises warning when not used with method 'fixed'", {
  st <- states$post_fragmentation
  expect_error(sample_cells(st, "fixed"), "cell_ids")
  expect_warning(sample_cells(st, "all", cell_ids = head(habitat_ids(st$grid), 3)), "cell_ids")
  expect_warning(sample_cells(st, "random", n_samples = 5L, cell_ids = head(habitat_ids(st$grid), 3)), "cell_ids")
})

test_that("method 'fixed' errors if any cell is not habitat", {
  st <- states$post_fragmentation
  non_habitat <- setdiff(seq_len(raster::ncell(st$grid)), habitat_ids(st$grid))
  selection <- c(head(habitat_ids(st$grid), 3), head(non_habitat, 2))
  expect_error(sample_cells(st, "fixed", cell_ids = selection), "not habitat")
  expect_error(sample_cells(st, "fixed", cell_ids = selection), "2 ")   # reports how many
})

test_that("method 'fixed' rejects duplicated and out-of-range cell_ids", {
  st <- states$post_fragmentation
  hab <- habitat_ids(st$grid)
  expect_error(sample_cells(st, "fixed", cell_ids = c(hab[1], hab[1])))
  expect_error(sample_cells(st, "fixed", cell_ids = raster::ncell(st$grid) + 1L))
  expect_error(sample_cells(st, "fixed", cell_ids = 0L))
})

test_that("a design drawn from the later state is valid in earlier states", {
  
  pre  <- states$pre_fragmentation
  post <- states$post_fragmentation

  design <- sample_cells(post, "random", n_samples = 10L, seed = 1L)$cell_id
  out_pre <- sample_cells(pre, "fixed", cell_ids = design)
  expect_equal(out_pre$cell_id, design)

  # the other direction fails if the design uses cells that were lost
  lost_cells <- head(which(is.na(post$grid[])), 5)
  expect_error(sample_cells(post, "fixed", cell_ids = lost_cells), "not habitat")
})

test_that("a design from habitat_cells() is applied unchanged across recorded states", {
  fin <- states$final
  design <- habitat_cells(fin, n_cells = min(10L, length(habitat_ids_of(fin))), seed = 1L)
  outs <- lapply(states[c("pre_fragmentation", "post_fragmentation", "final")],
                 sample_cells, method = "fixed", cell_ids = design)
  for (o in outs) {
    expect_identical(o$cell_id, as.integer(design))
    expect_identical(unique(o$seed_sampling), 1L)
  }
})

# ---- seeding: random method ---------------------------------------------------

test_that("same seed gives identical random samples, different seeds do not", {
  st <- states$final
  a <- sample_cells(st, "random", n_samples = 15L, seed = 5L)
  b <- sample_cells(st, "random", n_samples = 15L, seed = 5L)
  c <- sample_cells(st, "random", n_samples = 15L, seed = 6L)
  expect_identical(a, b)
  expect_false(identical(sort(a$cell_id), sort(c$cell_id)))
})

test_that("random sampling records the seed that was used", {
  st <- states$final
  out <- sample_cells(st, "random", n_samples = 15L, seed = 5L)
  expect_type(out$seed_sampling, "integer")
  expect_identical(unique(out$seed_sampling), 5L)
})

test_that("seedless random sampling records a fresh seed that reproduces the sample", {
  st <- states$final
  a <- sample_cells(st, "random", n_samples = 15L)
  expect_false(anyNA(a$seed_sampling))
  expect_length(unique(a$seed_sampling), 1L)

  b <- sample_cells(st, "random", n_samples = 15L, seed = a$seed_sampling[1])
  expect_identical(a, b)
})

test_that("seedless random sampling follows the global seed", {
  st <- states$final

  set.seed(1); a <- sample_cells(st, "random", n_samples = 15L)
  set.seed(1); b <- sample_cells(st, "random", n_samples = 15L)
  expect_identical(a, b)

  # consecutive seedless calls advance the stream and give different draws
  set.seed(1)
  x <- sample_cells(st, "random", n_samples = 15L)
  y <- sample_cells(st, "random", n_samples = 15L)
  expect_false(identical(x$seed_sampling[1], y$seed_sampling[1]))
})

test_that("a seeded call neither consumes nor depends on the global RNG", {
  st <- states$final

  set.seed(1)
  before <- rng_state()
  a <- sample_cells(st, "random", n_samples = 15L, seed = 5L)
  expect_identical(rng_state(), before)

  set.seed(999)
  b <- sample_cells(st, "random", n_samples = 15L, seed = 5L)
  expect_identical(a, b)
})

test_that("seed_sampling is reported only when a seed was used", {
  st <- states$final
  hab <- habitat_ids(st$grid)    # all habitat cells, no seed attribute
  n <- min(5L, length(hab))
  design <- habitat_cells(st, n_cells = n, seed = 3L)

  expect_false("seed_sampling" %in% names(sample_cells(st, "all", seed = 1L)))
  expect_false("seed_sampling" %in% names(sample_cells(st, "checkerboard", seed = 1L)))
  expect_false("seed_sampling" %in% names(sample_cells(st, "fixed", cell_ids = hab, seed = 1L)))
  expect_identical(unique(sample_cells(st, "fixed", cell_ids = design)$seed_sampling), 3L)
  expect_false("seed_sampling" %in% names(
    sample_cells(st, "fixed", cell_ids = as.integer(design))   # same cells, attribute stripped
  ))
  expect_true("seed_sampling" %in% names(sample_cells(st, "random", n_samples = n, seed = 1L)))
})

test_that("works with a seed from derive_seeds()", {
  st <- states$post_fragmentation
  seed <- derive_seeds(42L, "sampling")$sampling
  a <- sample_cells(st, "random", n_samples = 15L, seed = seed)
  b <- sample_cells(st, "random", n_samples = 15L, seed = derive_seeds(42L)$sampling)
  expect_identical(a, b)    # same seed whichever group was requested
  expect_identical(unique(a$seed_sampling), seed)
})

# ---- output formats -----------------------------------------------------------

test_that("wide format has one row per sampled cell and species columns with (only numeric) abundances", {
  st <- states$post_fragmentation
  present_species <- unique(st$ss_abund$species_id)
  out <- sample_cells(st, "all", format = "wide")
  sp_cols <- grep("^sp_", names(out), value = TRUE)

  expect_equal(nrow(out), length(habitat_ids(st$grid)))
  expect_equal(sp_cols, sprintf("sp_%d", sort(present_species)))      # sorted, no sp_NA
  expect_false(anyNA(out[sp_cols]))                      # missing abundances filled with 0
  expect_true(all(unlist(out[sp_cols]) >= 0))
})

test_that("long format keeps cells without any species (as NA species rows)", {
  st <- states$final
  hab <- habitat_ids(st$grid)
  rc <- raster::rowColFromCell(st$grid, hab[1])
  # remove all species from the first habitat cell, so it should appear in the output with NA species
  st$ss_abund <- st$ss_abund[!(st$ss_abund$x_loc == rc[1, 1] & st$ss_abund$y_loc == rc[1, 2]), ]

  out <- sample_cells(st, "all", format = "long")
  expect_true(all(hab %in% out$cell_id))
  expect_equal(sum(out$cell_id == hab[1]), 1L)
  expect_true(is.na(out$species_id[out$cell_id == hab[1]]))
  expect_true(all(c("species_id", "n") %in% names(out)))
})

test_that("wide and long format describe the same abundances", {
  st <- states$final
  long <- sample_cells(st, "all", format = "long")
  wide <- sample_cells(st, "all", format = "wide")
  expect_equal(sum(wide[grep("^sp_", names(wide))]), sum(long$n, na.rm = TRUE))
})

# ---- metadata and patch information -------------------------------------------

test_that("metadata is attached to every row", {
  st <- states$final
  out <- sample_cells(st, "checkerboard")
  expect_true(all(out$sim_id == st$meta$sim_id))
  expect_true(all(out$step == st$step))
  expect_true(all(out$step_label == st$step_label))
  expect_true(all(out$samp_method == "checkerboard"))
})

test_that("patch_size equals the number of habitat cells sharing a patch_id", {
  st <- states$post_fragmentation
  out <- sample_cells(st, "all")      # all habitat cells, so counts are complete
  counts <- table(out$patch_id)
  expect_equal(out$patch_size, as.integer(counts[as.character(out$patch_id)]))
})

test_that("sampling is not altered by the order of ss_abund rows", {
  st <- states$post_fragmentation
  st_shuffled <- st
  st_shuffled$ss_abund <- withr::with_seed(3L, st$ss_abund[sample.int(nrow(st$ss_abund))])
  a <- sample_cells(st, "all", format = "wide")
  b <- sample_cells(st_shuffled, "all", format = "wide")
  expect_equal(a[order(a$cell_id), ], b[order(b$cell_id), ], ignore_attr = TRUE)
})

test_that("the cell_id column carries no attributes", {
  st <- states$final
  expect_null(attributes(sample_cells(st, "random", n_samples = 5L, seed = 1L)$cell_id))
  expect_null(attributes(sample_cells(st, "fixed", cell_ids = habitat_cells(st, 5, seed = 1L))$cell_id))
})
