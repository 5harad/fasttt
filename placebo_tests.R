# This script runs placebo tests, replacing race with a column we would not expect to affect the outcome
# like weekday or month. 

make_placebo_filename = function(column_to_use_as_placebo, decision){
  # small helper function -- make filename of placebo data. 
  return(sprintf('placebo_%s_stop_and_frisk_%s.RData', column_to_use_as_placebo, decision))
}

run_placebo_tests = function(decision){
  # Fit the threshold models. Checked. 
  stopifnot(decision %in% c('stop_decision', 'search_decision'))
  for(column_to_use_as_placebo in c('weekday', 'month')){
    filename = make_placebo_filename(column_to_use_as_placebo, decision)
    make_stop_and_frisk_dataframe(decision, column_to_use_as_placebo = column_to_use_as_placebo, filename = filename)
    if(decision == 'search_decision'){run_threshold_test(gsub('.RData', '', filename), new_frisk_model_name)}
    else{run_threshold_test(gsub('.RData', '', filename), new_stop_model_name)}
  }
}

make_placebo_plots = function(decision){
  # makes the placebo plots for stop or search decision. Checked. 
  stopifnot(decision %in% c('stop_decision', 'search_decision'))
  df_to_plot = NULL
  month_names = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'June', 'July', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
  for(placebo_category in c('weekday', 'month')){
    model = gsub('.stan', '', ifelse(decision == 'search_decision', new_frisk_model_name, new_stop_model_name))
    results_fn = sprintf("%s/%s", base_output_dir, sprintf("placebo_%s_stop_and_frisk_%s_%s.RData", placebo_category, decision, model))
    message(sprintf("Loading placebo data for %s, %s\nPath is %s", placebo_category, decision, results_fn))
    load(results_fn)
    message(sprintf("Rhat is %s", Rhat))

    # aggregate across precincts. 
    mean_and_CIs = get_threshold_CIs_for_state(obs, post)
    mean_and_CIs$driver_race = as.character(mean_and_CIs$driver_race)
    mean_and_CIs$placebo_category = placebo_category
    
    # reformat names nicely. 
    if(placebo_category == 'weekday'){mean_and_CIs$driver_race = str_sub(mean_and_CIs$driver_race, 1, 3)}
    if(placebo_category == 'month'){
      mean_and_CIs$driver_race = plyr::mapvalues(mean_and_CIs$driver_race, 
                                               c(1,2,3,4,5,6,7,8,9,10,11,12), 
                                               month_names)
    }
    print(mean_and_CIs)
    df_to_plot = rbind(df_to_plot, mean_and_CIs)
  }

  df_to_plot$driver_race = factor(df_to_plot$driver_race, 
                                  levels = c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun', 
                                             month_names))
  df_to_plot$placebo_category = Hmisc::capitalize(df_to_plot$placebo_category)
  ymax = ifelse(decision == 'search_decision', .1, .01)
  p = ggplot(df_to_plot) + 
    geom_errorbar(aes(x = driver_race, ymin = low_CI, ymax = high_CI, width = 0)) + 
    geom_point(aes(x = driver_race, y = mean_threshold)) + 
    facet_grid(.~placebo_category, scales = 'free_x', space = 'free_x') + 
    ylab("Inferred threshold") + xlab("") + 
    scale_y_continuous(labels = scales::percent, limits = c(0, ymax), expand = c(0, 0)) + 
    our_theme() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(p, filename = paste0(base_code_dir, sprintf('figures/placebo_test_%s.pdf', decision)), width = 8, height = 3)
}


