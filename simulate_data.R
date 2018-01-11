simulate_stop_dataset <- function(df, post, noise=0) {
  # Simulates a new dataset using the stop model where each race-department threshold with a random error ~ N(0,se=noise). 
  message(sprintf("Simulating stop dataset with noise %2.3f", noise))
  
  # extract police department, race, and number of stops
  nraces <- length(unique(df$driver_race))
  n      <- integer(nrow(df))
  
  # initialize vectors
  searches_noisy    <- n
  hits_noisy        <- n
  thresh_noisy      <- n
  
  # extract posterior means for parameters
  phi    <- colMeans(post$phi)
  delta    <- colMeans(post$delta)
  t        <- df$thresholds
  
  # loop over precincts. For each precinct, randomly sample from
  for(d in unique(df$location_variable)){
    # get parameters for precinct. 
    idxs = df$location_variable == d
    
    precinct_thresholds = t[idxs]
    precinct_phi = phi[idxs]
    precinct_delta = delta[idxs]
    base_rates = df$base_population[idxs]
    stopifnot(abs(sum(base_rates) - 1) < 1e-6)
    total_searches = sum(df$num_searches[idxs])
    
    # draw the races of a whole bunch of randomly sampled pedestrians from the base population. 
    searches = c()
    r = c()
    signal = NULL
    n_draws = 1e4
    while (sum(searches) < total_searches) {
      multinomial_draw = rmultinom(n_draws, 1, base_rates) 
      new_r = integer(n_draws)
      for(r_i in 1:nraces){  # race order is NOT necessarily White-Black-Hispanic here, it is whatever order  
        new_r[multinomial_draw[r_i,] == 1] = r_i
      }
      t_noisy    <- inv.logit(logit(precinct_thresholds[new_r]) + rnorm(n_draws, mean=0, sd=noise))
      signal_new = draw_from_signal_distribution(n_draws, phi = precinct_phi[new_r], delta = precinct_delta[new_r], sigma_g = 1) %>% 
        sample_frac() # shuffle
      p_guilt    <- signal_to_p(signal_new$signal_samples, phi = precinct_phi[new_r], delta = precinct_delta[new_r], sigma_g = 1)
      r = c(r, new_r)
      signal = bind_rows(signal, signal_new)
      searches = c(searches, p_guilt > t_noisy)
      n_draws = n_draws*10
    }
    searches = which(searches)[1:total_searches]

    r = factor(r, levels = 1:nraces)
    
    searches_by_race = table(r[searches])
    hits_by_race = table(r[searches][signal$guilty[searches] == 'guilty'])

    searches_noisy[idxs] = searches_by_race
    hits_noisy[idxs] = hits_by_race
  }
  
  # return simulated dataset along with aggregate search & hit rates
  df <- df %>% 
    mutate(t_baseline   = t,
           num_searches = searches_noisy,
           num_hits     = hits_noisy)
  return(df)
}


simulate_search_dataset <- function(df, post, noise=0) {
  # Simulates a new dataset using the search model where each race-department threshold with a random error ~ N(0,se=noise). 

  # extract police department, race, and number of stops
  nraces <- length(unique(df$driver_race))
  n      <- integer(nrow(df))
  
  # initialize vectors
  searches_noisy    <- n
  hits_noisy        <- n
  thresh_noisy      <- n
  
  # extract posterior means for parameters
  phi_r    <- colMeans(post$phi_r)
  phi_d    <- colMeans(post$phi_d)
  delta_r <- colMeans(post$delta_r)
  delta_d <- colMeans(post$delta_d)
  t        <- df$thresholds
  
  for (i in 1:nrow(df)) {
    r <- as.integer(df$driver_race[i])
    d <- as.integer(df$location_variable[i])
    nr_stops = df$num_stops[i]
    
    phi    = inv.logit(phi_d[d] + phi_r[r])
    delta = exp(delta_r[r] + delta_d[d])
    
    # draw probability of posessing contraband for each stop in the race-department pair
    signal = draw_from_signal_distribution(nr_stops, phi = phi, delta = delta, sigma_g = 1)
    p_guilt    <- signal_to_p(signal$signal_samples, phi = phi, delta = delta, sigma_g = 1)
    
    # perturb threshold:  t'_rd = t_rd + N(mu=0, se=noise)
    t_noisy    <- inv.logit(rep(logit(t[i]), nr_stops) + rnorm(nr_stops, mean=0, sd=noise))
    
    # update results
    searches_noisy[i] <- sum(p_guilt > t_noisy)
    hits_noisy[i]     <- sum((signal$guilty == 'guilty') & (p_guilt > t_noisy))
    thresh_noisy[i]   <- mean(t_noisy)
    
    
  }
  
  # return simulated dataset along with aggregate search & hit rates
  df <- df %>% 
    mutate(t_baseline   = t,
           t_noisy      = thresh_noisy,
           num_searches = searches_noisy,
           num_hits     = hits_noisy,
           search_rate  = num_searches/num_stops,
           hit_rate     = num_hits/num_searches)
  return(df)
}

run_on_simulated_data = function(noise_levels, decision, rerun_models = TRUE){
  model = gsub('.stan', '', ifelse(decision == 'search_decision', new_frisk_model_name, new_stop_model_name))
  for(noise in noise_levels){
    message(sprintf("Running simulation for noise level %2.3f, decision %2.3f", noise, decision))
    filestring = sprintf('simulated_stop_and_frisk_noise_%2.3f_%s', noise, decision)
    output_filename = paste0(base_output_dir, sprintf('%s_%s.RData', filestring, model))
    if(rerun_models){
      if(decision == 'stop_decision'){
        load(paste0(base_output_dir, 'nyc_stop_and_frisk_multinomial_stop_and_frisk.RData'))
        stops = simulate_stop_dataset(obs, post, noise = noise)
      }else{
        load(paste0(base_output_dir, 'nyc_stop_and_frisk_search_decision_model_mixture.RData'))
        stops = simulate_search_dataset(obs, post, noise = noise)
      }
      save(stops, file = paste0(base_input_dir, filestring, '.RData'))
      run_threshold_test(filestring, paste0(model, '.stan'))
    }
  }
}

plot_all_noise_levels = function(decision, noise_levels = c(0, 0.05, 0.1, 0.2, 0.5, 1)) {
  model = gsub('.stan', '', ifelse(decision == 'search_decision', new_frisk_model_name, new_stop_model_name))
  df = list()
  for(noise in noise_levels){
    filestring = sprintf('simulated_stop_and_frisk_noise_%2.3f_%s', noise, decision)
    output_filename = paste0(base_output_dir, sprintf('%s_%s.RData', filestring, model))
    load(output_filename)
    df[[length(df)+1]] = get_single_threshold_from_state(obs) %>%
      mutate(noise = noise)
  }
  df = bind_rows(df)
  if(decision == 'stop_decision'){
    load(paste0(base_output_dir, 'nyc_stop_and_frisk_multinomial_stop_and_frisk.RData'))
  }else{
    load(paste0(base_output_dir, 'nyc_stop_and_frisk_search_decision_model_mixture.RData'))
  }
  true_thresholds = get_single_threshold_from_state(obs)
  ggsave(sprintf('figures/%s_simulation.pdf', decision),
    ggplot(df) + geom_line(aes(x=noise, y=weighted_threshold , color=driver_race)) +
      geom_hline(data = true_thresholds, mapping = aes(yintercept = weighted_threshold, color = driver_race), linetype=2) +
      our_theme() + theme(legend.position = c(0.8,0.65)) +
      scale_x_continuous(expression('\nThreshold heterogeneity (' * sigma * ')')) + 
      scale_y_continuous('Inferred threshold', labels = scales::percent) +
      scale_color_manual('', values = c('White' = 'blue', 'Black' = 'black', 'Hispanic' = 'red')),
    width = 6, height = 3.5)
}

# lines below are useful if you want to run several simulations in parallel using a bash script. 
# noise_level = as.numeric(commandArgs(trailingOnly=TRUE)[1])
# decision = commandArgs(trailingOnly=TRUE)[2]
# run_on_simulated_data(c(noise_level), decision) -
