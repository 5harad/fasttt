# Looking at frisk decision. See mixture_analysis.R. 
message("\n\nComparing old and new frisk results!")
compare_old_and_new_results(filenames = c('nyc_stop_and_frisk_search_decision'), 
                            old_model_name = old_frisk_model_name, 
                            new_model_name = new_frisk_model_name)
new_model = 'nyc_stop_and_frisk_search_decision_model_mixture.RData'
load(paste0(base_output_dir, new_model))
message("\n\nRemaking signal, PPC, and threshold plots for frisk data!")
hit_rate_ppc(obs, post, filename = 'figures/search_decision_hit_rate_ppc.pdf')
search_rate_ppc(obs, post, filename = 'figures/search_decision_search_rate_ppc.pdf')
make_threshold_plot(obs, 
                    var = 'frisk threshold', 
                    filename = 'figures/stop_and_frisk_search_thresholds.pdf', 
                    size_column = 'num_stops', 
                    breaks = c(.02, .05, .1), 
                    limits = c(.01, .2))


# Analysis of stop rates. 
message(sprintf("Analyzing stop rates!"))
compare_old_and_new_results(filenames = c('nyc_stop_and_frisk'), 
                            old_model_name = old_stop_model_name, 
                            new_model_name = new_stop_model_name)
load(paste0(base_output_dir, 'nyc_stop_and_frisk_multinomial_stop_and_frisk.RData'))
obs$num_stops = obs$race_base_pop # we size circles (for stop threshold plots and PPCs) by the size of the race population in the precinct
message("\n\nRemaking signal, PPC, and threshold plots for stop data!")
make_threshold_plot(obs, 
                    var = 'stop threshold', 
                    filename = 'figures/stop_and_frisk_stop_thresholds.pdf', 
                    size_column = 'race_base_pop', 
                    breaks = c(.01, .1, .5), 
                    limits = c(.0002, 1))

hit_rate_ppc(obs, post, 'figures/stop_decision_hit_rate_ppc.pdf')
make_multinomial_search_rate_ppc(obs, post, size_column = 'race_base_pop', 'figures/stop_decision_search_rate_ppc.pdf')
run_perturbation_tests() # this function comes from mixture_analysis.R. 

# Make placebo robustness plots. Only do this for search decision (doesn't make sense for stop decision). 
message("Making placebo plots!")
source('placebo_tests.R')
make_placebo_plots('search_decision')

# Make substratification robustness plots. 
message("Making substratification plots!")
source('substratification_robustness_checks.R')
collect_all_substratification_thresholds_into_dataframe()
all_substratification_thresholds = read_csv(all_substratification_thresholds_fn)
make_combined_plot_of_substratifications(all_substratification_thresholds, 'search_decision')
make_combined_plot_of_substratifications(all_substratification_thresholds, 'stop_decision')

message('Making heterogeneous threshold plots!')
source('simulate_data.R')
plot_all_noise_levels('search_decision')
plot_all_noise_levels('stop_decision')
message("Successfully remade all plots and tables.")



