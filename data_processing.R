read_in_stop_and_frisk_data = function(){
  # Checked. 
  load(stop_and_frisk_data_path)
  stops = filter(stops, year >= 2008, year <= 2012)
  stops$suspect.race = as.character(stops$suspect.race) 
  stops$suspect.race = plyr::mapvalues(stops$suspect.race, 
                                       c('black', 'black hispanic', 'white', 'white hispanic'), 
                                       c('Black', 'Hispanic', 'White', 'Hispanic'))
  
  stops = filter(stops, suspect.race %in% c('Black', 'Hispanic', 'White')) %>%
    rename(driver_race = suspect.race, 
           location_variable = precinct)
  stops$location_variable = as.numeric(as.character(stops$location_variable))
  return(stops)
}

add_placebo_column = function(stops, column_to_use_as_placebo){
  # Replaces race with a different column which we do not expect to have an effect. Checked. 
  stopifnot(column_to_use_as_placebo %in% c('weekday', 'month', 'quarter'))
  if(column_to_use_as_placebo == 'month'){
    stops$driver_race = month(strptime(stops$date, '%Y-%m-%d'))
  }
  if(column_to_use_as_placebo == 'quarter'){
    stops$driver_race = quarter(strptime(stops$date, '%Y-%m-%d'))
  }
  if(column_to_use_as_placebo == 'weekday'){
    stops$driver_race = weekdays(strptime(stops$date, '%Y-%m-%d'))
  }
  stops = filter(stops, !is.na(driver_race))
  return(stops)
}
make_stop_and_frisk_dataframe = function(analysis_to_conduct, 
                                         white_population_counterfactual_perturbation = 1, 
                                         col_to_filter_on = NULL,
                                         filter_fxn = NULL, 
                                         filename = NULL, 
                                         column_to_use_as_placebo = NULL){
  #Creates the dataframes to analyze. Checked. 
  # analysis_to_conduct: either stop_decision or search_decision. search_decision refers to frisk model. 
  # white_population_counterfactual_perturbation: 1 by default. Robustness check for stop model. 
  # col_to_filter_on: if not NULL, apply filter_fxn to this column and only take stops returning TRUE. 
  # filename: the filename to save the data under. 
  # column_to_use_as_placebo: used to create a placebo dataset. 
  
  # Data source for stop dataset: 
  # https://fusiontables.google.com/DataSource?dsrcid=767562#rows:id=1
  # http://johnkeefe.net/nyc-police-precinct-and-census-data
  # http://www.census.gov/prod/cen2010/doc/pl94-171.pdf
  stopifnot(analysis_to_conduct %in% c('stop_decision', 'search_decision'))
  # make sure we're not running multiple weird tests at once. 
  stopifnot((!is.null(col_to_filter_on)) + 
              (!is.null(column_to_use_as_placebo)) + 
              (white_population_counterfactual_perturbation != 1) <= 1)
  
  # read in the data. 
  stops = read_in_stop_and_frisk_data() %>% as.data.frame()
  
  # if we want to run the substratification analysis, filter dataset first. 
  # we make no other modifications to the pipeline.
  if(!is.null(col_to_filter_on)){ 
    message(sprintf("Filtering on column %s", col_to_filter_on))
    stops = stops[filter_fxn(stops[,col_to_filter_on]),]
  }
  
  # if we want to run placebo analysis, change race column to weekday or month. 
  if(!is.null(column_to_use_as_placebo)){
    message("Creating placebo dataset!")
    stopifnot(!is.null(filename))
    stops = add_placebo_column(stops, column_to_use_as_placebo)
  }
  
  # Now create the dataframe. 
  if(analysis_to_conduct == 'stop_decision'){
    if(white_population_counterfactual_perturbation != 1){
      message(sprintf("Perturbing the white population by %2.3f", white_population_counterfactual_perturbation))
    }
    census_pop_by_race_and_precinct = read_csv(paste0(base_code_dir, 'NYC_Blocks_2010CensusData_Plus_Precincts.csv')) %>%
      group_by(precinct) %>% summarise(Hispanic = sum(P0020002), 
                                       White = sum(P0020005), 
                                       Black = sum(P0020006)) %>%
      mutate(White = White * white_population_counterfactual_perturbation) %>%
      gather(key = 'driver_race', value = 'race_base_pop', Hispanic:Black) %>%
      rename(location_variable = precinct)
    # compute summary stats for stops. 
    stops = stops %>% 
      mutate(stopped_because_suspected_weapon = (suspected.crime == 'cpw') & (!is.na(suspected.crime))) %>%
      group_by(location_variable, driver_race) %>%
      summarise(num_searches = sum(stopped_because_suspected_weapon), 
                num_hits = sum(stopped_because_suspected_weapon & found.weapon)) %>% 
      ungroup()
    
    # if we're not using a placebo column, add the race 
    if(is.null(column_to_use_as_placebo)){
      stops = stops %>% 
      left_join(census_pop_by_race_and_precinct, by = c("location_variable", "driver_race")) %>% 
      ungroup() 
    }else{ # otherwise, we just have a constant race base population. 
      stops$race_base_pop = 1000
    }

    #normalize race base pop
    stops = stops %>% 
      group_by(location_variable) %>% 
      mutate(total_pop_in_precinct = sum(race_base_pop), 
             base_population = race_base_pop / total_pop_in_precinct) %>% 
      ungroup() %>%
      arrange(location_variable, driver_race) %>%
      mutate(searches_per_capita = num_searches / race_base_pop, 
             hit_rate = num_hits / num_searches) %>%
      filter(total_pop_in_precinct >= 1000)
  }else{
    stops = stops %>% 
      group_by(location_variable, driver_race) %>%
      summarise(num_stops = n(), 
                num_searches = sum(frisked.bc.weapons), 
                num_hits = sum(frisked.bc.weapons & found.weapon)) %>% 
      ungroup() %>%
      mutate(search_rate = num_searches / num_stops, 
             hit_rate = num_hits / num_searches)
  }
  
  
  stops = filter(stops, location_variable != 22) # filter out central park b/c no good population data.
  stops = stops %>% mutate(location_variable = as.factor(location_variable))
  
  # change race column to a factor. If it is really race (not a placebo or something else) reorder levels. 
  if(is.null(column_to_use_as_placebo)){
    stops$driver_race = factor(stops$driver_race, levels = c('White', 'Black', 'Hispanic'))
  }else{
    stops$driver_race = factor(stops$driver_race)
  }
  
  if(is.null(filename)){
    if(analysis_to_conduct == 'stop_decision'){
      if(white_population_counterfactual_perturbation == 1){
        filename = 'nyc_stop_and_frisk.RData'
      }
      else{
        filename = sprintf('nyc_stop_and_frisk_white_perturbation_%2.5f.RData', white_population_counterfactual_perturbation)
      }
    }else{
      filename = 'nyc_stop_and_frisk_search_decision.RData'
    }
  }
  filename = paste0(base_input_dir, filename)
  message(sprintf("Saving dataframe with %i rows to %s", nrow(stops), filename))
  save(stops, file = filename)
  message("Head of dataframe")
  print(head(stops %>% as.data.frame()))
  return(stops)
}

create_stops_per_precinct_df = function(){
  load(paste0(base_input_dir, 'nyc_stop_and_frisk_search_decision.RData'))
  stops = stops %>% 
    group_by(location_variable) %>%
    summarise(total_stops = sum(num_stops)) %>% 
    ungroup()
  write_csv(stops, stops_per_precinct_filename)
}

make_all_stop_and_frisk_data = function(){
  message("Making all stop and frisk data.")
  make_stop_and_frisk_dataframe('stop_decision')
  for(perturbation in white_population_perturbations){
    make_stop_and_frisk_dataframe('stop_decision', white_population_counterfactual_perturbation = perturbation)
  }
  make_stop_and_frisk_dataframe('search_decision')
  create_stops_per_precinct_df()
}
