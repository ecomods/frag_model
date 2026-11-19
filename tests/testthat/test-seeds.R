library(here)
library(testthat)
source(here("Model/src/seeds.R"))

all_processes <- c(
  "landscape", "community", "distribution",
  "spin_up", "fragmentation", "post_frag",
  "sampling", "analysis"
)

# ---- Structure of the return value -----------------------------------------

test_that("Return value is a named list of integer seeds with an origin attribute", {
  s <- derive_seeds(origin = 42, which = "all")

  expect_type(s, "list")
  expect_named(s, all_processes)
  expect_true(all(vapply(s, is.integer, logical(1))))
  expect_false(anyNA(unlist(s)))
  expect_identical(attr(s, "origin"), 42L)
})

test_that("Seeds lie in the 28-bit range implied by the 7-hex-character truncation", {
  for (o in c(1, 42, 12345, 99999)) {
    v <- unlist(derive_seeds(o, "all"))
    expect_true(all(v >= 0L & v <= 268435455L))
  }
})

test_that("Seeds can be passed to set.seed() and give reproducible streams", {
  s <- derive_seeds(42, "landscape")$landscape
  expect_no_error(set.seed(s))

  set.seed(s); a <- runif(3)
  set.seed(s); b <- runif(3)
  expect_identical(a, b)
})

# ---- Selection of processes via `which` ------------------------------------

test_that("Group keys expand to the documented processes, in order", {
  expect_named(derive_seeds(42, "all"), all_processes)
  expect_named(derive_seeds(42, "init"), c("landscape", "community", "distribution"))
  expect_named(
    derive_seeds(42, "model"),
    c("landscape", "community", "distribution", "spin_up", "fragmentation", "post_frag")
  )
  expect_named(derive_seeds(42, "workflow"), c("sampling", "analysis"))
})

test_that("Single process names return exactly that process", {
  s <- derive_seeds(42, "spin_up")
  expect_named(s, "spin_up")
  expect_length(s, 1L)
})

test_that("Groups and names can be combined; names follow first appearance", {
  expect_named(
    derive_seeds(42, c("workflow", "spin_up")),
    c("sampling", "analysis", "spin_up")
  )
  expect_named(
    derive_seeds(42, c("sampling", "init")),
    c("sampling", "landscape", "community", "distribution")
  )
})

test_that("Overlapping selections do not produce duplicates", {
  expect_named(
    derive_seeds(42, c("init", "landscape")),
    c("landscape", "community", "distribution")
  )
  expect_named(
    derive_seeds(42, c("init", "model")),
    c("landscape", "community", "distribution", "spin_up", "fragmentation", "post_frag")
  )
  expect_named(derive_seeds(42, c("spin_up", "spin_up")), "spin_up")
  # "all" combined with anything else is still just "all"
  expect_identical(derive_seeds(42, c("all", "sampling")), derive_seeds(42, "all"))
})

# ---- Reproducibility: seed depends only on origin and process --------------

test_that("A process gets the same seed regardless of selection, order or specification", {
  full <- derive_seeds(42, "all")

  selections <- list(
    "spin_up",
    c("analysis", "spin_up", "landscape"),
    c("model"),
    c("workflow", "post_frag"),
    c("spin_up", "init")
  )
  for (w in selections) {
    sel <- derive_seeds(42, w)
    expect_identical(unlist(sel), unlist(full[names(sel)]))
  }
})

test_that("Reordering `which` only reorders the output, not the seed values", {
  a <- derive_seeds(42, c("landscape", "sampling", "spin_up"))
  b <- derive_seeds(42, c("spin_up", "landscape", "sampling"))
  expect_identical(unlist(a)[sort(names(a))], unlist(b)[sort(names(b))])
})

test_that("Shared processes get the same seed in 'init' and 'model'", {
  init  <- derive_seeds(42, "init")
  model <- derive_seeds(42, "model")
  for (p in names(init)) expect_identical(init[[p]], model[[p]])
})

test_that("Seeds are reproducible with the same origin", {
  seed1 <- derive_seeds(origin = 12345, which = "all")
  seed2 <- derive_seeds(origin = 12345, which = "model")
  seed3 <- derive_seeds(origin = 12345, which = "workflow")
  seed4 <- derive_seeds(origin = 12345, which = "fragmentation")
  expect_equal(seed1$landscape, seed2$landscape)
  expect_equal(seed1$sampling, seed3$sampling)
  expect_equal(seed1$fragmentation, seed4$fragmentation)
  expect_equal(seed2$fragmentation, seed4$fragmentation)
  expect_identical(derive_seeds(12345, "all"), derive_seeds(12345, "all"))
})

test_that("Integer and double origins are equivalent (incl. values that paste() would format as 1e+05)", {
  expect_identical(derive_seeds(42), derive_seeds(42L))
  expect_identical(derive_seeds(1e5), derive_seeds(100000L))
})

# ---- Distinctness ----------------------------------------------------------

test_that("Seeds are distinct for each process", {
  seeds <- list(
    derive_seeds(origin = 12345, which = "all"),
    derive_seeds(origin = 54321, which = "model"),
    derive_seeds(origin = 67890, which = "workflow")
  )
  for (s in seeds) expect_length(unique(unlist(s)), length(s))
  expect_length(unique(unlist(seeds)), sum(lengths(seeds)))
})

test_that("Seeds are distinct across processes for many origins", {
  # 28-bit seeds: collisions are possible in principle, but not for these origins
  for (o in 1:200) {
    v <- unlist(derive_seeds(o, "all"))
    expect_length(unique(v), length(all_processes))
  }
})

test_that("The same process gets different seeds for different origins", {
  v <- vapply(1:50, function(o) derive_seeds(o, "landscape")$landscape, integer(1))
  expect_length(unique(v), 50L)

  a <- derive_seeds(1, "all")
  b <- derive_seeds(2, "all")
  expect_false(any(unlist(a) == unlist(b)))
})

# ---- Hashing method --------------------------------------------------------

test_that("Seeds follow the documented method: first 7 hex chars of md5('<origin>_<process>')", {
  # NOTE: digest() is called with its default serialize = TRUE, as in derive_seeds().
  # If you switch to serialize = FALSE, update this reference calculation too.
  expected <- strtoi(
    substr(digest::digest("42_spin_up", algo = "md5"), 1, 7),
    base = 16L
  )
  expect_identical(derive_seeds(42, "spin_up")$spin_up, expected)
})

# ---- Input validation ------------------------------------------------------

test_that("Invalid `origin` is rejected", {
  expect_error(derive_seeds(0))
  expect_error(derive_seeds(-5))
  expect_error(derive_seeds(1.5))
  expect_error(derive_seeds(NA_integer_))
  expect_error(derive_seeds("42"))
  expect_error(derive_seeds(c(1, 2)))
  expect_error(derive_seeds(NULL))
})

test_that("Boundary origins are accepted", {
  expect_no_error(derive_seeds(1))
  expect_no_error(derive_seeds(.Machine$integer.max))
})

test_that("Invalid `which` is rejected", {
  expect_error(derive_seeds(42, "typo"))
  expect_error(derive_seeds(42, c("init", "typo")))   # unknown entry mixed with a valid group
  expect_error(derive_seeds(42, character(0)))
  expect_error(derive_seeds(42, NA_character_))
  expect_error(derive_seeds(42, 1))
  expect_error(derive_seeds(42, NULL))
})

test_that("Matching is case-insensitive", {
  expect_identical(derive_seeds(42, "All"), derive_seeds(42, "all"))
  expect_identical(derive_seeds(42, "SPIN_UP"), derive_seeds(42, "spin_up"))
  expect_identical(derive_seeds(42, c("Init", "SAMPLING")), derive_seeds(42, c("init", "sampling")))
})

test_that("Unique partial matches are accepted and return full names", {
  expect_named(derive_seeds(42, "land"), "landscape")
  expect_named(derive_seeds(42, c("frag", "post")), c("fragmentation", "post_frag"))
  expect_named(derive_seeds(42, "al"), all_processes)
  expect_identical(derive_seeds(42, "dist")$distribution, derive_seeds(42, "all")$distribution)
})

test_that("Exact matches take precedence over partial ones", {
  expect_named(derive_seeds(42, "model"), names(derive_seeds(42, "model")))
  expect_length(derive_seeds(42, "model"), 6L)
})

test_that("Ambiguous or unknown abbreviations are rejected", {
  expect_error(derive_seeds(42, "s"), "ambiguous")
  expect_error(derive_seeds(42, "a"), "ambiguous")
  expect_error(derive_seeds(42, "typo"), "unknown")
  expect_error(derive_seeds(42, c("init", "typo")), "typo")
  expect_error(derive_seeds(42, ""))
})
