source('mixture_analysis.R')

get_stan_data = function(model_file, obs){
  if(model_file == 'working_countrywide_mixture.stan'){
    stan_data = with(obs, list(
      N = nrow(obs),
      D = length(unique(location_variable)),
      R = length(unique(driver_race)),
      K = length(unique(state)),
      state = as.integer(state),
      d = as.integer(location_variable),
      r = as.integer(driver_race),
      n = num_stops,
      s = num_searches,
      h = num_hits
    ))
  }else if(model_file %in% c(old_stop_model_name, new_stop_model_name)){
    # add a couple asserts to make sure data is ordered correctly. 
    error_message = 'Data is not formatted correctly for the stop model. Every location should have data for all races, and dataframe should be sorted by location_variable.'
    if(!(as.character(all_equal(obs, obs %>% arrange(location_variable), ignore_row_order = FALSE)) == 'TRUE')){stop(error_message)}
    if(!all(obs %>% group_by(location_variable) %>% summarise(n = n()) %>% .$ n %>% unique() == length(unique(obs$driver_race)))){stop(error_message)}
    if(!all(obs %>% group_by(driver_race) %>% summarise(n = n()) %>% .$ n %>% unique() == length(unique(obs$location_variable)))){stop(error_message)}
    # package data for Stan. 
    stan_data = with(obs, list(
      N = nrow(obs),
      D = length(unique(location_variable)),
      R = length(unique(driver_race)),
      d = as.integer(location_variable),
      r = as.integer(driver_race),
      base_population_proportions = base_population,
      s = num_searches,
      h = num_hits
    ))
  }else if(model_file %in% c(old_frisk_model_name, new_frisk_model_name)){
    stan_data = with(obs, list(
      N = nrow(obs),
      D = length(unique(location_variable)),
      R = length(unique(driver_race)),
      d = as.integer(location_variable),
      r = as.integer(driver_race),
      n = num_stops,
      s = num_searches,
      h = num_hits
    ))
  }else{
    message("Error: not a valid model name.")
    stopifnot(FALSE)
  }
  return(stan_data)
}

get_stan_initialization = function(model_file, stan_data){
  if(model_file %in% c(new_frisk_model_name, new_stop_model_name)){
    init_fn = 'random' # checked. No special initialization needed. 
  }else if(model_file %in% c(old_frisk_model_name, old_stop_model_name)){
    initializer <- function(num_obs, num_races, num_locations) {	
      # force immediate evaluation of arguments
      force(num_obs); force(num_races); force(num_locations);
      function() {
        list(sigma_t        = runif(num_races, 0.05, 0.25),
             t_r            = runif(num_races, -3, -1),  
             t_i_raw        = runif(num_obs, -0.75, 0.75),
             phi_r          = runif(num_races, -3.75, -3.25),
             lambda_r       = runif(num_races, 1.5, 2.5), 
             phi_d_raw      = runif(num_locations-1, -0.25, 0.25),
             lambda_d_raw   = runif(num_locations-1, -0.25, 0.25)
        )
      }
    }
    init_fn <- initializer(stan_data$N, stan_data$R, stan_data$D)
  }else if(model_file == 'working_countrywide_mixture.stan'){ # this gives an example of how to fit the model on the national traffic stops data. 
    initializer <- function(num_obs, num_races, num_depts, num_states) {	
      # force immediate evaluation of arguments
      force(num_obs); force(num_races); force(num_depts); force(num_states);
      function() {
        list(sigma_t        = runif(1, 0.1, 1),
             phi_r          = runif(num_races, -5, -3),
             phi_state          = runif(num_states, -.25, -.25),
             sigma_phi          = runif(num_states, 0.01, .1),
             phi_d      = runif(num_depts, -0.25, 0.25),
             delta_r       = runif(num_races, -.3, .3), 
             delta_state       = runif(num_states, -.3, .3),
             sigma_delta          = runif(num_states, .01, .1),
             delta_d   = runif(num_depts, -0.25, .25)
        )
      }
    }
    init_fn <- initializer(stan_data$N, stan_data$R, stan_data$D, stan_data$K)
  }else{
    message("Error: not a valid model name.")
    stopifnot(FALSE)
  }
  return(init_fn)
}

run_mixture_mcmc <- function(stops, output_filename, iter = 5000, warmup = NULL, chains = 5, adapt_delta = 0.95, max_treedepth = 12, sample_from_prior = FALSE, verbose = FALSE, simulation=FALSE, model_file = 'mixture_model.stan') {
  # checked. 
  if (is.null(warmup)) {
    if (sample_from_prior) {
      warmup = 0
    } else {
      warmup = ceiling(iter/2)
    }
  }
  warmup = ceiling(iter / 2)
  obs = stops

  # check model identifiability
  r <- length(unique(obs$driver_race))
  d <- length(unique(obs$location_variable))
  if (! (r >= 3 & d >= 5)) {
    stop('Not enough departments to constrain estimates')
  }
  
  # Package data for Stan
  stan_data = get_stan_data(model_file, obs)

  # set up parameter initialization
  init_fn = get_stan_initialization(model_file, stan_data)
  
  # fit the model
  my_model = stan_model(paste0('stan_models/', model_file))
  t1 = now()
  message('Starting sampling now!')
  fit <- sampling(my_model, data = stan_data, iter=iter, init = init_fn, chains=chains, cores=chains, refresh = 50, warmup = warmup, control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth, adapt_engaged = !sample_from_prior), verbose = verbose, diagnostic_file = paste0(output_filename, '_diag.txt'))
  seconds_required = as.numeric(now() - t1, units = "secs")
  post = rstan::extract(fit)
  s = summary(fit)
  Rhat = max(s$summary[,'Rhat'], na.rm = TRUE)
  message(sprintf("RHat is %2.3f", Rhat))
  if((model_file == 'model_flat.stan') | (model_file == 'multinomial_model_flat.stan')){
    obs = get_thresholds_from_post_old_model(post, obs);
  }else{
    obs$thresholds = get_thresholds_from_post(post, obs)
  }
  print(get_single_threshold_from_state(obs))
  save(file=paste0(output_filename, '.RData'), obs, post, fit, seconds_required, Rhat)
}

run_threshold_test <- function(file_prefix, model_file){
  #Checked. This actually runs the threshold test. 
  input_file = paste0(base_input_dir, file_prefix, '.RData')
  stopifnot(file.exists(input_file))
  stopifnot(file.exists(paste0('stan_models/', model_file)))
  message(sprintf("Loading %s", input_file))
  load(input_file)
  model_name = gsub('.stan', '', model_file)
  out_name = paste0(base_output_dir, sprintf('%s_%s', file_prefix, model_name))
  
  message(sprintf('Running threshold test on %s: %i locations, %i races, saving results to %s', 
                  file_prefix, 
                  length(unique(stops$location_variable)), 
                  length(unique(stops$driver_race)), 
                  out_name))
  output = run_mixture_mcmc(stops, out_name, iter=5000, chains=5, adapt_delta=.9, max_treedepth=12, model_file=model_file)
}






