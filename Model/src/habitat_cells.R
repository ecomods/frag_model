#' Get the habitat cells of a state or landscape
#'
#' Returns the cell IDs of all habitat cells (non-\code{NA} cells of the grid),
#' or a random subset of them. The result can be used directly as
#' \code{cell_ids} for \code{\link{sample_cells}()} with
#' \code{method = "fixed"}, and is used internally for
#' \code{method = "random"}.
#'
#' @param x A \code{model_state} object (its \code{grid} is used) or a
#'   \code{RasterLayer}. Habitat cells are the cells that are not \code{NA}.
#'
#' @param n_cells Number of habitat cells to draw. Either a whole number
#'   between 1 and the number of habitat cells, or \code{Inf} (default), which
#'   returns all habitat cells without a random draw. Values above the number of
#'   habitat cells give an error. A value equal to the number of habitat cells
#'   is a draw like any other: all cells are returned in random order and the
#'   seed is recorded.
#'
#' @param seed Optional integer seed for the draw; ignored if
#'   \code{n_cells = Inf}. If \code{NULL}, a seed is drawn from the current
#'   random number stream (so results follow \code{set.seed()}) and used
#'   instead.
#'
#' @details
#' Habitat is always derived from the grid itself, not from the cached
#' \code{cells$habitat} of a state, so the result cannot drift from the
#' landscape. Cell IDs are raster cell numbers, counted row by row from the top
#' left, and are returned in ascending order when \code{n_cells = Inf}.
#'
#' The draw runs under \code{withr::with_seed()}: a call with a \code{seed}
#' leaves the global random number state untouched, and a seedless call
#' consumes one integer from the global stream. In both cases passing the
#' recorded seed back reproduces the draw for the same grid. The same seed does
#' \emph{not} give the same cells for grids with different habitat, because the
#' draw is made from the cells that are habitat in the given grid.
#'
#' @return An integer vector of cell IDs. For a random draw (finite
#'   \code{n_cells}) it carries the attribute \code{seed}, the integer seed that
#'   was actually used. \code{sample_cells()} reports this seed in the column
#'   \code{seed_sampling} when the vector is passed as \code{cell_ids}. The
#'   attribute is lost when the vector is subsetted, sorted or otherwise
#'   modified, so pass it on unchanged. For \code{n_cells = Inf} there is no
#'   attribute.
#'
#' @seealso \code{\link{sample_cells}}, \code{\link{derive_seeds}}
#'
#' @examples
#' r <- raster::raster(nrows = 10, ncols = 10)
#' raster::values(r) <- ifelse(seq_len(100) %% 3 == 0, NA, 1)
#'
#' habitat_cells(r)                            # all habitat cells
#' design <- habitat_cells(r, n_cells = 5, seed = 1L)
#' attr(design, "seed")                        # 1
#'
#' @importFrom raster getValues
#' @importFrom withr with_seed
#' @importFrom checkmate assert_number assert_count
#'
#' @export
habitat_cells <- function(
  x,
  n_cells = +Inf,
  seed = NULL
) {

  grid <- if (inherits(x, "model_state")) {
    x$grid
  } else if (inherits(x, "RasterLayer")) {
    x
  } else {
    stop("x must be a single model_state or RasterLayer object")
  } 

  habitat_ids <- which(!is.na(raster::getValues(grid)))

  checkmate::assert_number(n_cells, lower = 1)      
  # accepts Inf, rejects NULL, NA, NaN, negative numbers, and non-numeric values
  if (is.infinite(n_cells)) return(habitat_ids)     # all habitat cells
  checkmate::assert_count(n_cells)                  # a whole number
  if (n_cells > length(habitat_ids)) {
    stop("n_cells cannot exceed the number of habitat cells", call. = FALSE)
  }

  used_seed <- if (is.null(seed)) {
    sample.int(.Machine$integer.max, 1L)
  } else {
    as.integer(seed)
  }

  selection <- withr::with_seed(
    used_seed,
    habitat_ids[sample.int(length(habitat_ids), n_cells)]
  )

  attr(selection, "seed") <- used_seed
  return(selection)
}