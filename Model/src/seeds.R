#' Derive seeds for each process from an origin seed
#'
#' Deterministically derives one integer seed per process from a single origin
#' seed, so that every process of the framework can be seeded independently
#' while the whole set stays reproducible from the origin alone.
#'
#' @details
#' **Method.** For each process, the string `"<origin>_<process>"` (e.g.
#' `"42_spin_up"`) is hashed with MD5 via [digest::digest()]. The first 7
#' hexadecimal characters of the hash are converted to an integer, giving a
#' seed in the range 0 to 268,435,455 (28 bits), which is safely within R's
#' integer range and valid for [set.seed()].
#'
#' **Reproducibility.** A process's seed is a pure function of the origin and
#' the process name only. It does not depend on which other processes are
#' requested, on the order in which they are given, or on how they are
#' specified (by name, via a group key, or both). For a given origin, a
#' process therefore always receives the same seed, e.g. `derive_seeds(42,
#' "spin_up")$spin_up` equals `derive_seeds(42, "all")$spin_up`. Different
#' origins or different process names yield different hash inputs and hence
#' (practically always) different seeds. Because the seed is truncated to 28
#' bits, distinctness is statistically very likely but not mathematically
#' guaranteed: a collision between two processes for a given origin is
#' possible in principle, with a probability on the order of 1e-7.
#'
#' **Groups.** Group keys are only a convenience for selecting several
#' processes at once and have no influence on the seeds. Processes that
#' belong to several groups (`"landscape"`, `"community"`, `"distribution"`
#' appear in both `"init"` and `"model"`) receive the same seed in either
#' context. Duplicates arising from overlapping selections are removed.
#'
#' @param origin A single positive integer (>= 1) from which all seeds are
#'   derived.
#' @param which Character vector selecting the processes for which to derive
#'   seeds. Matching is case-insensitive and accepts unique abbreviations (e.g. `"land"` 
#'   for `"landscape"`); exact matches take precedence, and ambiguous
#'   abbreviations (e.g. `"s"`) are an error. Returned names are always the
#'   full process names. 
#'   Can contain process names, group keys, or a combination of both:
#'   \itemize{
#'     \item process names: `"landscape"`, `"community"`, `"distribution"`,
#'       `"spin_up"`, `"fragmentation"`, `"post_frag"`, `"sampling"`,
#'       `"analysis"`
#'     \item group keys: `"init"` (landscape, community, distribution),
#'       `"model"` (landscape, community, distribution, spin_up,
#'       fragmentation, post_frag), `"workflow"` (sampling, analysis), and
#'       `"all"` (every process)
#'   }
#'   Defaults to `"all"`.
#'
#' @importFrom checkmate assert_int assert_character
#' @importFrom digest digest
#'
#' @return A named list of integer seeds, one element per selected process
#'   (in order of first appearance, without duplicates), with the origin
#'   stored in the attribute `"origin"`.
#'
#' @examples
#' derive_seeds(42)
#' derive_seeds(42, which = "init")
#' derive_seeds(42, which = c("workflow", "spin_up"))
#'
#' # the seed of a process does not depend on the selection
#' identical(
#'   derive_seeds(42, "spin_up")$spin_up,
#'   derive_seeds(42, "all")$spin_up
#' )
derive_seeds <- function(
    origin, 
    which = "all"
  ) {

  checkmate::assert_int(origin, lower = 1L)
  checkmate::assert_character(which, any.missing = FALSE, min.len = 1L)

  origin <- as.integer(origin)

  groups <- list(
    init = c("landscape", "community", "distribution"),
    model = c("landscape", "community", "distribution", "spin_up", "fragmentation", "post_frag"),
    workflow = c("sampling", "analysis")
  )
  groups$all <- unique(unlist(groups, use.names = FALSE))

  lookup <- c(groups, stats::setNames(as.list(groups$all), groups$all))
  keys <- names(lookup)

  # case-insensitive exact or unique-partial matching against groups and process names
  idx <- pmatch(tolower(which), keys, duplicates.ok = TRUE)

  if (anyNA(idx)) {
    bad <- which[is.na(idx)]
    ambiguous <- vapply(
      tolower(bad), function(b) sum(startsWith(keys, b)) > 1L, logical(1)
    )
    stop(
      "Cannot match `which` entries: ",
      paste0("'", bad, "' (", ifelse(ambiguous, "ambiguous", "unknown"), ")", collapse = ", "),
      ".\nValid keys: ", paste(keys, collapse = ", "), call. = FALSE
    )
  }

  processes <- unique(unlist(lookup[idx], use.names = FALSE))

  hash_strings <- paste(origin, processes, sep = "_")
  hash_hex <- vapply(
    hash_strings, 
    digest::digest, character(1), algo = "md5", USE.NAMES = FALSE
  )
  seed_values <- strtoi(substr(hash_hex, 1, 7), base = 16)

  seeds <- as.list(stats::setNames(seed_values, processes))

  attr(seeds, "origin") <- origin
  return(seeds)

}
