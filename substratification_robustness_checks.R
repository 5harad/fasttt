# This function fits models on subsets of the data. 

all_substratification_thresholds_fn = paste0(base_code_dir, 'all_substratification_thresholds.csv')

get_filter_fxns = function(stratification_cat){
  # returns the functions applied to the stop dataframe to filter for a category. Checked. 
  stopifnot(stratification_cat %in% c('suspect.sex', 'suspect.age', 'location.housing', 'time', 'year'))
  if(stratification_cat == 'suspect.sex'){
    filter_fxns = c(function(x){return(x == 'female')}, 
                    function(x){return(x == 'male')})
  }
  if(stratification_cat == 'suspect.age'){
    filter_fxns = c(function(x){return(x < 18)}, 
                    function(x){return((x >= 18) & (x < 25))},
                    function(x){return((x >= 25) & (x < 35))},
                    function(x){return((x >= 35) & (x <= 45))},
                    function(x){return((x > 45))})
  }
  if(stratification_cat == 'location.housing'){
    filter_fxns = c(function(x){return(x == 'housing')}, 
                    function(x){return(x == 'neither')},
                    function(x){return(x == 'transit')})
  }
  if(stratification_cat == 'time'){
    get_hour = function(x){return(hour(strptime(x, '%H:%M')))}
    filter_fxns = c(function(x){get_hour(x) < 6},
                    function(x){(get_hour(x) >= 6) & (get_hour(x) < 12)},
                    function(x){(get_hour(x) >= 12) & (get_hour(x) < 18)},
                    function(x){(get_hour(x) >= 18)})
  }
  if(stratification_cat == 'year'){
    filter_fxns = c(function(x){x == 2008}, 
                    function(x){x == 2009}, 
                    function(x){x == 2010}, 
                    function(x){x == 2011}, 
                    function(x){x == 2012})
    
  }
  return(filter_fxns)
}

get_category_names = function(stratification_cat){
  # returns the names for each subgroup category. Needs to be in same order as filter functions above. Checked.
  if(stratification_cat == 'suspect.sex'){
    filenames = c('sex_female','sex_male')
  }
  if(stratification_cat == 'suspect.age'){
    filenames = c('age_under_18', 'age_18_25', 'age_25_35', 'age_35_45', 'age_over_45')
  }
  if(stratification_cat == 'location.housing'){
    filenames = c('location_public_housing', 'location_neither', 'location_transit')
  }
  if(stratification_cat == 'time'){
    filenames = c('time_midnight_6AM', 'time_6AM_12PM', 'time_12PM_6PM', 'time_6PM_midnight')
  }
  if(stratification_cat == 'year'){
    filenames = c('year_2008', 'year_2009', 'year_2010', 'year_2011', 'year_2012')
  }
  return(filenames)
}

make_substratification_filename = function(subgroup_name, decision){
  return(sprintf('substratification_%s_stop_and_frisk_%s.RData', subgroup_name, decision))
}


run_all_stratification_models = function(){
  # runs all substratification for a set of decisions (ie, stop_decision or search_decision) and categories to substratify by. 
  # Checked. 
  decisions = c('stop_decision', 'search_decision')
  stratification_categories = c('suspect.sex', 'suspect.age', 'time', 'year')
  for(decision in decisions){
    for(stratification_category in stratification_categories){
      # for each category, loop over the subsets of that category we want to analyze. 
      # for example, if category were "suspect.sex", we would loop over subgroups "male" and "female". 
      filter_fxns = get_filter_fxns(stratification_category)
      subgroup_names = get_category_names(stratification_category)
      for(i in 1:length(subgroup_names)){
        filename = make_substratification_filename(subgroup_names[i], decision)
        stops = make_stop_and_frisk_dataframe(decision,
                                              col_to_filter_on = stratification_category, 
                                              filter_fxn = filter_fxns[[i]], 
                                              filename = filename)
        
        if(decision == 'search_decision'){
          run_threshold_test(gsub('.RData', '', filename), new_frisk_model_name)
        }else{
          run_threshold_test(gsub('.RData', '', filename), new_stop_model_name)
        }
      }
    }
  }
}

make_combined_plot_of_substratifications = function(all_substratification_thresholds, decision_name){
  # make a combined plot of all the thresholds for each substratification group. 
  # all_substratification_thresholds is the dataframe that has all thresholds. 
  # Checked. 
  message(sprintf("Making substratification plot for %s", decision_name))
  factor_levels = c('Under 18', '18-25', '25-35', '35-45', 'Over 45', 'Female', 'Male', 
                    '12AM-6AM', '6AM-12PM', '12PM-6PM', '6PM-12AM', '2008', '2009', '2010', '2011', '2012')

  # rename a few things for nice display. 
  all_substratification_thresholds$subgroup_name = 
    plyr::mapvalues(all_substratification_thresholds$subgroup_name, 
                    c('under-18', 'over-45', 'midnight-6AM', '6PM-midnight'),
                    c('under 18', 'over 45', '12AM-6AM', '6PM-12AM'))
                    
  all_substratification_thresholds$subgroup_name = Hmisc::capitalize(all_substratification_thresholds$subgroup_name)
  all_substratification_thresholds = filter(all_substratification_thresholds, 
                                            decision == decision_name, 
                                            substratification_category != 'location') %>% 
    mutate(substratification_category = Hmisc::capitalize(substratification_category), 
           subgroup_name = factor(subgroup_name, levels = factor_levels), 
           driver_race = factor(driver_race, levels = c('White', 'Black', 'Hispanic'))) # set factor levels in correct order. 
  ymax = max(all_substratification_thresholds$weighted_threshold) * 1.1
  p = ggplot(all_substratification_thresholds) + 
    geom_line(aes(x = subgroup_name, y = weighted_threshold, color = driver_race, group = driver_race)) + 
    facet_grid(.~substratification_category, scales = 'free_x', space = 'free_x') + 
    scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(-.001, ymax)) + 
    our_theme() + 
    xlab("") + ylab("Inferred threshold") + 
    scale_color_manual(values=c('blue','black','red')) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(plot.margin = unit(c(.2,.2,0,.2), "cm"))
  
  fn = paste0(base_code_dir, sprintf('figures/all_substratifications_%s.pdf', decision_name))
  ggsave(p, filename = fn, width = 8, height = 3)
}

get_thresholds_for_category = function(category, decision){
  # returns a dataframe of thresholds for a category and decision. Checked. 
  output_files = list.files(base_output_dir)
  output_files = output_files[grepl('.RData', output_files) & 
                                grepl('substratification', output_files) & 
                                grepl(decision, output_files) & 
                                grepl(category, output_files)]
  thresholds = NULL
  for(fn in output_files){
    load(paste0(base_output_dir, fn))
    # extract name of substratification. This is messy :(
    val = str_split(str_split(fn, '_stop_and_frisk')[[1]][1], category)[[1]][2]
    val = str_sub(val, 2, -1)
    val = gsub('_', '-', val)
    message(sprintf('substratification model %s, value %s, Rhat is %2.3f', fn, val, Rhat))
    thresholds_for_subcategory = get_single_threshold_from_state(obs)
    thresholds_for_subcategory = thresholds_for_subcategory %>% 
      mutate(substratification_category = category, 
             decision = decision, 
             subgroup_name = val)
    thresholds = rbind(thresholds, thresholds_for_subcategory)
  }
  
  if(is.null(thresholds)){
    message(sprintf("No models for %s and %s", category, decision))
    return(NULL)
  }
  
  return(thresholds)
  
}

collect_all_substratification_thresholds_into_dataframe = function(){
  # combine all the substratification thresholds into a single dataframe. 
  # Convenience method because loading all the models in takes a while. Checked. 
  all_substratification_thresholds = NULL
  for(decision in c('stop_decision', 'search_decision')){
    for(category in c('sex', 'age', 'year', 'time', 'location')){
      message(sprintf('Loading substratification thresholds for %s %s', decision, category))
      thresholds = get_thresholds_for_category(category = category, decision = decision)
      all_substratification_thresholds = rbind(all_substratification_thresholds, thresholds)
    }
  }
  stopifnot(all_substratification_thresholds %>% group_by(decision, subgroup_name) %>% summarise(n = n()) %>% .$n == 3)
  write_csv(all_substratification_thresholds, all_substratification_thresholds_fn)
}



