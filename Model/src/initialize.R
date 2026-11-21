################ Simulation initialization Function ############################

# The function combines 3 functions ("generate_grid", "generate_agents", and
# "distribute_agent") together to simplify the simulation initialization process.
# The function then returns the simulation grid, the updated agents list and
# the agents distribution grid in the form of a list

initialize <- function(
  grid_size,
  ac_amount,
  n_species,
  n_pop,
  species_par,
  k_inter,
  k_intra,
  niche_breadth = NULL,
  seeds,
  species_specific = FALSE,
  state_type = c("core_state", "full_state")
) {

  state_type <- match.arg(state_type)
  
  if (is.null(seeds$landscape) || is.null(seeds$community) || is.null(seeds$distribution)) {
    stop("Seeds for landscape, community, and distribution must be provided in the `seeds` list.")
  }

  # For model-detached use
  if (is.null(species_par) && !species_specific) {
    species_par <- data.frame(
      species_id = 1:n_species,
      n_value = seq(from = 0, to = 1, length.out = n_species)
    )
  }

  grid <- fbm_fft(
    gr_size = grid_size,
    ac_amount = ac_amount,
    seed = seeds$landscape
  )

  cells <- list(
    habitat = 1:grid_size^2,
    edges = NULL
  )

  agents_init <- generate_agents(
    n_species = n_species,
    n_pop = n_pop,
    seed = seeds$community
  )

  agents <- distribute_agent(
    agents = agents_init,
    grid = grid,
    species_par = species_par,
    k_inter = k_inter,
    k_intra = k_intra,
    nb = niche_breadth,
    species_specific = species_specific,
    seed = seeds$distribution
  )

  state_out <- record_step(
    grid = grid,
    agents = agents,
    step = 0L,
    cells = cells
  )

  if (state_type == "core_state") {
    return(state_out)
  } else if (state_type == "full_state") {
    meta <- build_meta(
      sim_id = NA,
      run_date = Sys.Date(),
      seeds = seeds,
      grid_size = grid_size,
      ac_amount = ac_amount,
      habitat = 1,
      fragmentation = NA,
      niche_breadth = niche_breadth,
      edge_effect = NA,
      dispersal_type = NA,
      dispersal_kernel = NA,
      dispersal_ratio = NA,
      dispersal_dist = NA
    )
    full_state <- record_state(
      model_state = state_out,
      meta = meta,
      step_label = "start"
    )
  }
}
