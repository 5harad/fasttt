source('constants_and_libraries.R')

signal_to_p = function(x, phi, delta, sigma_g){
  #Checked. Converts x -> p. 
  p = phi * dnorm(x, delta, sigma_g) / (phi * dnorm(x, delta, sigma_g) + (1 - phi) * dnorm(x, 0, 1));
  return(p)
}

draw_from_signal_distribution = function(n, phi, delta, sigma_g){
  #Checked. Draws samples from the x distribution. 
  n_guilty = rbinom(1, size = n, prob = phi)
  n_innocent = n - n_guilty
  signal_samples = c(rnorm(n = n_guilty, mean = delta, sd = sigma_g), 
                     rnorm(n = n_innocent, mean = 0, sd = 1))
  return(data.frame(signal_samples, guilty = c(rep('guilty', n_guilty), rep('innocent', n_innocent))))
}

plot_probability_histogram = function(phi, delta, sigma_g){
  #Checked. 
  n = 1e6
  sampled_signal = draw_from_signal_distribution(n, phi = phi, delta = delta, sigma_g = sigma_g)
  
  p_samples = signal_to_p(sampled_signal$signal_samples,
                          phi = phi, delta = delta, sigma_g = sigma_g)
  p = ggplot() + geom_density(aes(x = p_samples))
  print(p)
}


get_thresholds_from_post = function(post, obs, quantile_to_report = NULL){
  # Returns the threshold for each race-precinct pair. 
  # Checked. If quantile_to_report is NULL, returns the mean; otherwise, returns the quantile_to_report.  
  thresholds = c()
  for(i in 1:dim(post$t_i)[2]){
    all_draws = signal_to_p(post$t_i[,i], phi = post$phi[,i], delta = post$delta[,i], sigma_g = 1)
    if(is.null(quantile_to_report)){
      thresholds = c(thresholds, 
                     mean(all_draws))
    }else{
      thresholds = c(thresholds, 
                     quantile(all_draws, quantile_to_report))
    }
    
  }
  return(thresholds)
}

get_thresholds_from_post_old_model <- function(post, obs){
  # Checked (old model).
  obs$thresholds    <- colMeans(post$t_i)
  ci    <- apply(post$t_i, 2, function(x) quantile(x, probs = c(0.025, 0.975))) 
  obs$low_ci  = ci[1,]
  obs$high_ci = ci[2,]
  return(obs)  
}

search_rate_ppc <- function(obs, post, filename, ylim=0.03) {
  # Checked. 
  obs$pred_search_rate = colMeans(post$search_rate)
  p = ggplot(data=obs, aes(x=pred_search_rate, y=pred_search_rate-search_rate)) +
    geom_point(aes(size=num_stops, color=driver_race), alpha = 0.8) + scale_size_area(max_size=10) +
    scale_x_continuous('\nPredicted frisk rate', labels=percent)+
    scale_y_continuous('Frisk rate prediction error\n', labels=percent, limits=c(-ylim, ylim)) +
    geom_abline(slope=0, intercept=0, linetype='dashed') +
    scale_color_manual(values=c('blue','black','red')) +
    guides(size=FALSE) + 
    our_theme()
  ggsave(plot = p, 
         filename = filename, 
         width = 8, 
         height = 4)
}

hit_rate_ppc <- function(obs, post, filename, ylim=0.3) {
  # Checked. 
  obs$pred_hit_rate = colMeans(post$hit_rate)
  p = ggplot(data=obs, aes(x=pred_hit_rate, y=hit_rate-pred_hit_rate)) +
    geom_point(aes(size=num_stops, color=driver_race), alpha=0.8) + scale_size_area(max_size=10) +
    scale_x_continuous('\nPredicted hit rate', labels=percent) +
    scale_y_continuous('Hit rate prediction error\n', labels=percent, limits = c(-ylim, ylim)) +
    geom_abline(slope=0, intercept=0, linetype='dashed') +
    scale_color_manual(values=c('blue','black','red')) +
    guides(size=FALSE) + 
    our_theme()
  ggsave(plot = p, 
         filename = filename, 
         width = 8, 
         height = 4)
}

make_multinomial_search_rate_ppc = function(obs, post, size_column, filename){
  # Checked. 
  obs = obs %>% as.data.frame()
  obs$size_column = obs[,size_column]
  obs$theta = colMeans(post$theta)
  obs$relative_search_probability = obs$theta
  obs = obs %>% 
    group_by(location_variable) %>% 
    mutate(true_fraction_of_searches = num_searches / sum(num_searches), 
           predicted_fraction_of_searches = relative_search_probability / sum(relative_search_probability)) %>%
    ungroup()
  
  p = ggplot(obs) + 
    geom_point(aes(x = true_fraction_of_searches, 
                   y = predicted_fraction_of_searches - true_fraction_of_searches, 
                   color = driver_race, 
                   group = driver_race, 
                   size = size_column)) + 
    xlab("True fraction of CPW stops") + 
    ylab("Prediction error") + 
    theme_bw(base_size=13) + scale_size_area(max_size=10) +
    scale_x_continuous(limits=c(0, 1), labels = scales::percent, expand=c(0,0)) + 
    scale_y_continuous(limits=c(-.1, .1), labels = scales::percent, expand=c(0,0)) +
    geom_abline(slope = 0, intercept = 0) + 
    scale_color_manual(values=c('blue','black','red')) + 
    guides(size=FALSE) + 
    our_theme()
  ggsave(plot = p, 
         filename = filename, 
         width = 8, 
         height = 4)
}

plot_signal <- function(obs, post, filename, xlim = c(.001,0.5), ymax = 5, for_paper = FALSE){
  # Checked. Very specific to stop-and-frisk data: weights by stops-per-precinct. 
  races = levels(obs$driver_race)  
  stopifnot(sum(races == c('White', 'Black', 'Hispanic')) == 3)
  x = seq(0.0001, 0.9999, 0.0001)
  
  phi = colMeans(post$phi)
  delta = colMeans(post$delta)
  
  threshold = get_single_threshold_from_state(obs)
  stops_per_precinct = get_stops_per_precinct() %>% 
    mutate(total_stops = total_stops / sum(total_stops)) %>%
    rename(weighting_for_precinct = total_stops)
  
  obs = obs %>% left_join(stops_per_precinct, by = 'location_variable')
  n = 5e6
  black_samples = c()
  white_samples = c()
  hispanic_samples = c()
  
  black_signal = c()
  white_signal = c()
  hispanic_signal = c()
  
  black_guilty = c()
  white_guilty = c()
  hispanic_guilty = c()
  
  for(i in 1:nrow(obs)){
    race = obs$driver_race[i]
    r_i = match(race, levels(obs$driver_race))
    n_stops = round(n * obs$weighting_for_precinct[i])
    signal_df = draw_from_signal_distribution(n_stops, phi = phi[i], delta = delta[i], sigma_g = 1)
    guilty_indicator = signal_df$guilty
    signal_samples = signal_df$signal_samples
    p_samples = signal_to_p(signal_samples, phi = phi[i], delta = delta[i], sigma_g = 1)
    
    if(race == 'White'){
      white_samples = c(white_samples, p_samples)
      white_signal = c(white_signal, signal_samples)
      white_guilty = c(white_guilty, guilty_indicator)
    }
    if(race == 'Hispanic'){
      hispanic_samples = c(hispanic_samples, p_samples)
      hispanic_signal = c(hispanic_signal, signal_samples)
      hispanic_guilty = c(hispanic_guilty, guilty_indicator)
    }
    if(race == 'Black'){
      black_samples = c(black_samples, p_samples)
      black_signal = c(black_signal, signal_samples)
      black_guilty = c(black_guilty, guilty_indicator)
    }
  }
  p_plt <- ggplot() + 
    geom_vline(data = threshold, aes(xintercept = weighted_threshold, color=driver_race), linetype="dashed")+
    geom_density(aes(x = white_samples), color = 'blue') +
    geom_density(aes(x = black_samples), color = 'black') +
    geom_density(aes(x = hispanic_samples), color = 'red') +
    scale_x_continuous('Likelihood of carrying a weapon', labels=percent, expand = c(0, 0), breaks = c(0, .2, .4)) +
    scale_y_continuous('\nDensity', breaks = c(), expand = c(0, 0)) + 
    coord_cartesian(ylim=c(0, 10), xlim = xlim) + 
    scale_color_manual(values = c('blue', 'black', 'red')) + 
    our_theme() + 
    theme(legend.position=c(.7,.7))
  p_plt
  df = data.frame(guilty = c(white_guilty, black_guilty, hispanic_guilty), 
                  signal = c(white_signal, black_signal, hispanic_signal), 
                  race = c(rep('White', length(white_guilty)), rep('Black', length(black_guilty)), rep('Hispanic', length(hispanic_guilty)))) 
  signal_plt <- ggplot(df) + 
    geom_density(aes(x = signal, y = ..count.., group = guilty)) + 
    scale_x_continuous('Signal', expand = c(0, 0)) + 
    scale_y_continuous('\nDensity', breaks = c(), expand = c(0, 0)) + 
    facet_grid(.~race) + 
    our_theme() + 
    theme(legend.position=c(.8,.8))
  if(for_paper){
    ggsave(p_plt, filename = filename, width = 5, height = 3)
  }else{
    pdf(filename, width = 5, height = 4)
    grid.arrange(p_plt, signal_plt, nrow = 2)
    dev.off()
  }
}

get_stops_per_precinct = function(){
  # Checked. Just reads in a dataframe of stops per precinct. 
  stops_per_precinct = suppressMessages(read_csv(stops_per_precinct_filename)) %>%
    as.data.frame() %>%
    mutate(location_variable = factor(location_variable))
  return(stops_per_precinct)
}


get_single_threshold_from_state = function(obs){
  # this infers a single threshold from a state. Checked. 
  # First read in stops per location to use as weighting. 
  obs = as.data.frame(obs)
  stops_per_location = get_stops_per_precinct()
  stops_per_location$location_variable = as.character(stops_per_location$location_variable)
  obs$location_variable = as.character(obs$location_variable)
  obs = obs %>% left_join(stops_per_location, by = c('location_variable'))
  
  # Now set the column that you want to compute the weighted mean of. 
  aggregate_thresholds = obs %>% 
    group_by(driver_race) %>% 
    summarise(weighted_threshold = weighted.mean(thresholds, total_stops))
  return(aggregate_thresholds %>% select(driver_race, weighted_threshold))
}

get_threshold_CIs_for_state = function(obs, post){
  # This function computes the CIs for thresholds for each race group (aggregated over precincts), weighting by total stops in each precinct. 
  # Checked. Verified that this returns the same mean threshold as get_single_threshold_from_state.
  # first compute weights. 
  obs = as.data.frame(obs)
  stops_per_location = get_stops_per_precinct()
  stops_per_location$location_variable = as.character(stops_per_location$location_variable)
  obs$location_variable = as.character(obs$location_variable)
  obs = obs %>% left_join(stops_per_location, by = c('location_variable'))
  
  # now compute the thresholds in probability space (not signal space). 
  all_thresholds_in_p_space = NULL
  for(i in 1:dim(post$t_i)[2]){# loop over precincts. 
    all_draws = signal_to_p(post$t_i[,i], phi = post$phi[,i], delta = post$delta[,i], sigma_g = 1)
    all_thresholds_in_p_space = cbind(all_thresholds_in_p_space, all_draws)
  }
  # now loop over race groups, computing weighted CIs for each group. 
  CIs = NULL
  for(race in unique(obs$driver_race)){
    idxs = which(obs$driver_race == race)
    weights = obs$total_stops[idxs]
    weights = weights / sum(weights)
    weighted_samples = 0
    for(i in 1:length(idxs)){ 
      weighted_samples = weighted_samples + weights[i] * all_thresholds_in_p_space[,idxs[i]]
    }
    CIs = rbind(CIs, data.frame(driver_race = race, 
                                low_CI = quantile(weighted_samples, .025), 
                                mean_threshold = mean(weighted_samples), 
                                high_CI = quantile(weighted_samples, .975), 
                                stringsAsFactors = FALSE))
  }
  return(CIs)
}


make_threshold_plot <- function(obs, var, filename, size_column, breaks, limits, log_scale = TRUE) {
  #Checked. 
  obs = data.frame(obs)
  obs$size_column = obs[,size_column] # this column is used to size circles. 
  df <- obs %>%
    filter(driver_race == 'White') %>%
    right_join(obs %>% filter(driver_race != 'White'), by='location_variable') %>%
    rename(white_threshold=thresholds.x, minority_threshold = thresholds.y, minority_n = size_column.y, minority = driver_race.y)
  
  p <- ggplot(df) +
    geom_point(aes(x=white_threshold, y=minority_threshold, size=minority_n), shape=1, alpha=0.6) +
    facet_grid(.~minority) +
    geom_abline(slope=1, intercept=0, linetype='dashed') +
    scale_size_area(max_size=10) +
    guides(size=FALSE, color=FALSE) + 
    our_theme()
  if(log_scale){
    p = p + 
      scale_x_log10(paste('White'   , var),  expand=c(0,0), limits = limits, breaks = breaks) +
      scale_y_log10(paste('Minority', var), expand=c(0,0), limits = limits, breaks = breaks) 
  }else{
    p = p + 
      scale_x_continuous(paste('White'   , var),  expand=c(0,0), limits = limits, breaks = breaks) +
      scale_y_continuous(paste('Minority', var), expand=c(0,0), limits = limits, breaks = breaks) 
  }
  message(sprintf("Made threshold plot and saved to %s", filename))
  ggsave(p, filename = filename, width = 8, height = 4)
  return(p)
}

get_timing_statistics = function(fit){
  # Checked. Returns timing statistics for a fit. 
  chain_times = as.data.frame(get_elapsed_time(fit))
  total_time = sum(chain_times$warmup + chain_times$sample)
  max_time = max(chain_times$warmup + chain_times$sample)
  n_eff_df = data.frame(summary(fit)$summary)
  n_eff = min(n_eff_df$n_eff)
  sampler_params = get_sampler_params(fit)
  n_chains = length(sampler_params)
  total_iter = 0
  total_leapfrog_steps = 0
  for(i in 1:n_chains){
    chain_results = as.data.frame(sampler_params[[i]])
    total_iter = total_iter + nrow(chain_results)
    total_leapfrog_steps = total_leapfrog_steps + sum(chain_results$n_leapfrog__)
  }
  time_per_neff = total_time / n_eff 
  iter_per_neff = total_iter / n_eff
  leapfrog_per_iter = total_leapfrog_steps / total_iter
  time_per_leapfrog = total_time / total_leapfrog_steps
  
  df = data.frame(time_per_neff, 
                  iter_per_neff, 
                  leapfrog_per_iter, 
                  time_per_leapfrog, 
                  total_iter, 
                  total_time, 
                  max_time,
                  total_leapfrog_steps, 
                  n_eff)
  return(df)
}

# Compare runtime and results on old and new datasets. Checked. 
compare_old_and_new_results = function(filenames, old_model_name, new_model_name){ 
  comparison_table = NULL
  message(sprintf("Comparing %s to %s in terms of timing and correlation", old_model_name, new_model_name))
  for(filename in filenames){
    message(sprintf("Analyzing %s", filename))
    load(paste0(base_output_dir, sprintf('%s_%s.RData', filename, gsub('.stan', '', new_model_name))))
    new_Rhat = Rhat
    new_obs = obs
    new_timing_statistics = get_timing_statistics(fit)
    load(paste0(base_output_dir, sprintf('%s_%s.RData', filename, gsub('.stan', '', old_model_name))))
    old_Rhat = Rhat
    old_obs = obs
    old_timing_statistics = get_timing_statistics(fit)
    combined_df = left_join(old_obs, new_obs, by = c('location_variable', 'driver_race'))
    correlation = cor(combined_df$thresholds.x, combined_df$thresholds.y)
    message(sprintf('Old RHat: %2.5f; new RHat: %2.5f; correlation: %2.5f', old_Rhat, new_Rhat, correlation))
    timing_statistics_ratio = old_timing_statistics / new_timing_statistics
    stopifnot((nrow(old_obs) == nrow(new_obs)) & (nrow(new_obs) == nrow(combined_df)))
    comparison_table = rbind(comparison_table, old_timing_statistics %>% mutate(filename = paste('old', filename)))
    comparison_table = rbind(comparison_table, new_timing_statistics %>% mutate(filename = paste('new', filename)))
    comparison_table = rbind(comparison_table, timing_statistics_ratio %>% mutate(filename = paste('ratio', filename)))
    
  }
  comparison_table = comparison_table %>% 
    select(filename, time_per_neff, iter_per_neff, leapfrog_per_iter, time_per_leapfrog, total_time, max_time)
  print(data.frame(comparison_table))
}

#Makes a graph of sensitivity to perturbation. Checked. 
run_perturbation_tests = function(){
  message(sprintf("Running perturbation tests."))
  all_thresholds = NULL
  for(perturbation in c(white_population_perturbations, 1)){
    if(perturbation == 1){
      model_name = sprintf('nyc_stop_and_frisk_multinomial_stop_and_frisk.RData', perturbation)
    }else{
      model_name = sprintf('nyc_stop_and_frisk_white_perturbation_%2.5f_multinomial_stop_and_frisk.RData', perturbation)
    }
    model_fn = paste0(base_output_dir, model_name)
    if(!file.exists(model_fn)){message(sprintf("%s does not exist; continuing", model_fn)); next;}
    load(model_fn)
    obs$num_stops = obs$race_base_pop
    message(sprintf('Perturbation %2.3f, RHat %2.3f, seconds required %2.3f', perturbation, Rhat, seconds_required))
    thresholds = get_single_threshold_from_state(obs)
    print(thresholds)
    all_thresholds = rbind(all_thresholds, thresholds %>% mutate(perturbation = perturbation))
  }
  
  p = ggplot(all_thresholds %>% filter(perturbation >= .5)) + 
    geom_line(aes(x = perturbation, 
                  y = weighted_threshold, 
                  group = driver_race, 
                  color = driver_race)) + 
    scale_color_manual(values = c('blue', 'black', 'red')) + 
    our_theme() + 
    theme(legend.position=c(.8,.7)) + 
    scale_y_continuous(labels = scales::percent) + 
    scale_x_continuous(expand = c(0, 0), limits = c(.5, 2.1)) + 
    xlab(expression(k[white])) + 
    ylab("Aggregate threshold")
  print(p)
  ggsave(plot = p, 
         filename = 'figures/stop_rate_robustness.pdf', 
         width = 6, 
         height = 4)
}


