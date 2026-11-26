#' Sample habitat cells from an ABM state
#'
#' Extracts a subset of habitat cells from a model state object and returns
#' per-cell species abundances along with spatial coordinates, patch
#' information, and simulation metadata. Cells can be sampled exhaustively
#' (\code{"all"}), randomly (\code{"random"}), in a deterministic checkerboard
#' pattern (\code{"checkerboard"}), or from a user-supplied set of cell IDs
#' (\code{"fixed"}).
#'
#' @param full_state A (\code{model_state}) object. A named list representing
#'   the ABM state at one recorded step. It must contain:
#'   \itemize{
#'     \item \code{grid}: a \code{RasterLayer} defining habitat (non-\code{NA})
#'       and non-habitat (\code{NA}) cells.
#'     \item \code{ss_abund}: a data table with columns
#'       \code{x_loc}, \code{y_loc}, \code{species_id}, and \code{n}, giving
#'       species abundances per cell.
#'     \item \code{step}, \code{step_label}: the recorded step and its label.
#'     \item \code{meta}: a named list of simulation metadata with the elements
#'       \code{sim_id}, \code{grid_size}, \code{ac_amount},
#'       \code{habitat}, \code{fragmentation}, \code{niche_breadth},
#'       \code{edge_effect}, \code{dispersal_type}, \code{dispersal_kernel},
#'       \code{dispersal_ratio}, and \code{dispersal_dist}.
#'   }
#'
#' @param method Sampling method used to select habitat cells. One of:
#'   \describe{
#'     \item{\code{"all"}}{All habitat cells are returned.}
#'     \item{\code{"random"}}{\code{n_samples} habitat cells are drawn at random
#'       (without replacement).}
#'     \item{\code{"checkerboard"}}{Habitat cells whose row and column indices
#'       sum to an even number (deterministic checkerboard pattern).}
#'     \item{\code{"fixed"}}{The habitat cells given in \code{cell_ids}, in the
#'       order supplied.}
#'   }
#'
#' @param n_samples Number of habitat cells to draw. Required for 
#'   \code{method = "random"}; ignored, with a warning, for other methods. 
#'   A whole number between 1 and the number of habitat cells; values outside
#'   this range give an error. A value equal to the number of habitat cells
#'   returns all cells in random order, with the seed recorded. 
#'
#' @param cell_ids Integer vector of cell IDs (cell numbers of \code{grid}) to
#'   sample. Required for \code{method = "fixed"}; ignored, with a warning, for
#'   other methods. The IDs must be unique, and every one must be a habitat cell
#'   in \code{full_state}, otherwise an error is raised. Designs are best drawn
#'   with \code{\link{habitat_cells}()}: if the vector still carries the
#'   \code{seed} attribute set by that function, the seed is added to the output
#'   in the column \code{seed_sampling}.
#'
#' @param format Output format. Either \code{"wide"} (one row per sampled cell,
#'   one column per species) or \code{"long"} (one row per cell-species
#'   combination).
#'
#' @param seed Optional integer seed for \code{method = "random"}; ignored for
#'   other methods. If \code{NULL}, a seed is drawn from the current random
#'   number stream (so results follow \code{set.seed()}) and used instead. The
#'   seed that was actually used is returned in the column \code{seed_sampling}.
#'   Note that the argument is ignored for \code{method = "fixed"} as well; in
#'   this case \code{seed_sampling} is reported only if \code{cell_ids} has a
#'   \code{seed} attribute (as set by \code{\link{habitat_cells}()}).
#'
#' @details
#' Cell coordinates (\code{x_loc}, \code{y_loc}) refer to row and column indices
#' of the underlying raster grid, as returned by
#' \code{raster::rowColFromCell()}.
#'
#' Patch identifiers and patch sizes are computed using
#' \code{toroidal_clump()}, which assigns patch IDs while accounting for the
#' toroidal nature of the grid.
#'
#' In wide format, missing species abundances are filled with zeros.
#'
#' \strong{Reproducibility.} Random cell selection runs under
#' \code{withr::with_seed()}, so a call with a \code{seed} leaves the global
#' random number state untouched. If \code{seed} is \code{NULL}, one integer is
#' drawn from the global stream to serve as the seed. In both cases passing the
#' reported \code{seed_sampling} back as \code{seed} reproduces the sample for
#' the same state. The same seed does \emph{not} give the same cells in states
#' with different habitat (for example before and after landscape modification),
#' because the draw is made from the cells that are habitat in the given state.
#' To use one set of cells across several recorded steps, draw it once from the
#' state with the least habitat, with \code{habitat_cells(state, n_cells,
#' seed)}, and pass the result as \code{cell_ids} to \code{method = "fixed"} for
#' the other states.
#'
#' @return A data frame (tibble) containing sampled cells and associated data.
#'   Columns, in order: \code{sim_id}, \code{seed_sampling} (only for
#'   \code{method = "random"} and if \code{cell_ids} has a \code{seed} attribute
#'   for \code{method = "fixed"}): the integer seed used for the cell selection,
#'   then simulation metadata ( \code{grid_size}, \code{ac_amount},
#'   \code{habitat}, \code{fragmentation}, \code{niche_breadth},
#'   \code{edge_effect}, \code{dispersal_type}, \code{dispersal_kernel},
#'   \code{dispersal_ratio}, \code{dispersal_dist}), the columns describing the
#'   recorded step ( \code{step}, \code{step_label}), \code{samp_method}, then
#'   the sample columns (\code{sample_id}, \code{cell_id}, \code{x_loc},
#'   \code{y_loc}, \code{patch_id}, \code{patch_size}), then species abundances.
#'   The last part depends on \code{format}:
#'   \itemize{
#'     \item \code{"wide"}: one row per sampled cell, with species abundances
#'       spread across columns prefixed with \code{"sp_"}.
#'     \item \code{"long"}: one row per sampled cell and species
#'       (\code{species_id}, \code{n}); cells without any species appear once
#'       with \code{NA}.
#'   }
#'
#' @seealso \code{\link{derive_seeds}} for the seed of the sampling process,
#'   \code{\link{habitat_cells}} for retrieving habitat cell IDs, and
#'   \code{\link{toroidal_clump}} for patch identification.
#'
#' @importFrom rlang arg_match
#' @importFrom raster getValues rowColFromCell freq ncell
#' @importFrom dplyr left_join select any_of
#' @importFrom tidyr pivot_wider
#' @importFrom tibble tibble
#' @importFrom withr with_seed
#' @importFrom checkmate assert_integerish
#' @importFrom magrittr %>%
#' @importFrom stats setNames
#' @importFrom utils head
#'
#' @export
sample_cells <- function(
  full_state,
  method = c("all", "checkerboard", "random", "fixed"),
  n_samples = NULL,
  cell_ids = NULL,
  format = c("wide", "long"),
  seed = NULL
) {
  method <- rlang::arg_match(method)
  format <- rlang::arg_match(format)

  if (!is.null(cell_ids) && method != "fixed") {
    warning("cell_ids argument is only used for fixed sampling.\nIf you want to sample specific cell_ids, use method = 'fixed'.")
  }
  if (!is.null(n_samples) && method != "random") {
    warning("n_samples argument is only used for random sampling.\nIf you want to sample a specific number of random cells, use method = 'random'.")
  }

  # Extract grid and species abundances per cell from state
  grid            <- full_state$grid
  ss_abund        <- full_state$ss_abund

  # Extract IDs of all habitat cells
  habitat_all <- habitat_cells(x = full_state)

  used_seed <- NULL # will be set only for random sampling and, 
                    # if cell_ids have a seed attribute, for fixed sampling

  # Select cells to be sampled
  sampled_cells <- switch(
    method,

    all = habitat_all,

    checkerboard = 
      habitat_all[((habitat_all - 1) %/% ncol(grid) + 
        (habitat_all - 1) %% ncol(grid)) %% 2 == 0],
      
    random = {
      if (is.null(n_samples)) {
        stop("n_samples must be provided for random sampling")
      } 
      checkmate::assert_count(n_samples, positive = TRUE)
      if (n_samples > length(habitat_all)) {
        stop("n_samples cannot exceed the number of habitat cells")
      } 
      selection <- habitat_cells(
        x = full_state,
        n_cells = n_samples,
        seed = seed
      )
      used_seed <- attr(selection, "seed", exact = TRUE)
      as.integer(selection)
    },

    fixed = {
      if (is.null(cell_ids)) {
        stop("cell_ids must be provided for fixed sampling")
      }
      checkmate::assert_integerish(
        cell_ids, lower = 1, upper = raster::ncell(grid),
        any.missing = FALSE, unique = TRUE, min.len = 1
      )
      bad <- setdiff(cell_ids, habitat_all)
      if (length(bad) > 0) {
        stop(
          sprintf("%d supplied cell(s) are not habitat in this state (e.g. cell(s) with ID(s) %s)",
                  length(bad), toString(utils::head(bad, 5))), call. = FALSE)
      }
      used_seed <- attr(cell_ids, "seed", exact = TRUE)
      as.integer(cell_ids)
    }
  )

  # Assess coordinates of sample cells
  sample_coords <- raster::rowColFromCell(grid, sampled_cells)

  samples <- tibble::tibble(
    sample_id = seq_len(length(sampled_cells)),
    cell_id = sampled_cells,
    x_loc = sample_coords[, 1],
    y_loc = sample_coords[, 2]
  )

  # Assess patch information
  clumped <- toroidal_clump(grid, directions = 4)
  patch_freq <- raster::freq(clumped, useNA = "no")
  patch_ids <- patch_freq[, 1]
  patch_sizes <- patch_freq[, 2]
  # Create named vector for patches: name = patch ID, value = patch size
  patches <- setNames(patch_sizes, patch_ids)
  # Add patch info to samples df
  samples$patch_id <- clumped[samples$cell_id]
  samples$patch_size <- unname(patches[as.character(samples$patch_id)])

  # Add static metadata columns
  meta <- full_state$meta
  meta_df <- tibble::tibble(
    sim_id           = meta$sim_id,
    seed_sampling    = used_seed,
    grid_size        = meta$grid_size,
    ac_amount        = meta$ac_amount,
    habitat          = meta$habitat,
    fragmentation    = meta$fragmentation,
    niche_breadth    = meta$niche_breadth,
    edge_effect      = meta$edge_effect,
    dispersal_type   = meta$dispersal_type,
    dispersal_kernel = meta$dispersal_kernel,
    dispersal_ratio  = meta$dispersal_ratio,
    dispersal_dist   = meta$dispersal_dist,
    step             = full_state$step,
    step_label       = full_state$step_label,
    samp_method      = method
  )

  # Merge metadata with samples df (metadata first, then sample-specific columns)
  meta_samp <- cbind(meta_df, samples)

  # Merge samples df and species abundances per cell
  out_long <- dplyr::left_join(meta_samp, ss_abund, by = c("x_loc", "y_loc"))

  # Return long format
  if (format == "long") return(out_long)

  # If needed, reformat and return wide format
  out_wide <- out_long %>%
    tidyr::pivot_wider(
      names_from = species_id,
      values_from = n,
      values_fill = 0,
      names_prefix = "sp_",
      names_sort = TRUE
    ) %>%
    dplyr::select(-dplyr::any_of("sp_NA"))  # Remove sp_NA if present

  return(out_wide)
}
