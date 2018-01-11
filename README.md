# fast-threshold-test

This repo contains code to reproduce the analysis in the paper ["Fast Threshold Tests for Detecting Discrimination"](https://arxiv.org/abs/1702.08536) by Emma Pierson, Sam Corbett-Davies, and Sharad Goel (AISTATS 2018).  

All scripts should be run using redo_all_analysis.R; do not attempt to run the other scripts directly. Set the appropriate TRUE/FALSE flags in order to run the analyses you need to run. To reproduce results: 

1. Download the stop-and-frisk data, available at https://5harad.com/data/sqf.RData. Also download NYC_Blocks_2010CensusData_Plus_Precincts.csv, available at http://johnkeefe.net/nyc-police-precinct-and-census-data, and store it the same folder where you store the code.  
2. Set the paths in constants_and_libraries.R in order to suit your system. You will need to create an input and an output directory where you want to store the fitted models -- probably not in the same directory where you store the code, because the models get quite large. 
3. Set the appropriate TRUE/FALSE flags in redo_all_analysis.R to rerun the analyses you want to run. Note that some analyses must precede others. For example, you cannot produce figures in the paper without first fitting the Stan models. 

Please contact emmap1@cs.stanford.edu with any comments or questions about the code. 

## R scripts

**redo_all_analysis.R**: the script for rerunning all the analysis.

**data_processing.R**: processes data into dataframes suited for Stan models. 

**run_mixture_mcmc.R**: runs the Stan models.

**make_plots_for_paper.R**: makes the plots for the paper. 

**mixture_analysis.R**: various helper methods for the analysis.

**constants_and_libraries.R**: paths, constants, and libraries. Edit this file to suit your system.

**placebo_tests.R**: run placebo tests. 

**substratification_robustness_checks.R**: run substratification robustness checks. 

**simulate_data.R**: run simulated data robustness checks. 

## Stan models

A note on naming conventions in the Stan models: throughout the codebase, we refer to the frisk decision as the search_decision and the stop decision as the stop_decision. driver_race refers to the race of the stopped pedestrian. These namings were chosen for consistency with the original threshold models to which we compare. In both the frisk model and the stop model, num_searches refers to the number of cases in which a pedestrian was above the inferred threshold.

**model_mixture.stan**: this is the new model which uses discriminant distributions to model the frisk decision. 

**multinomial_stop_and_frisk.stan**: this is the new model which uses discriminant distributions to model the stop decision. 

**model_flat.stan**: this is the old model which uses beta distributions to model the frisk decision.

**multinomial_model_flat.stan**: this is the old model which uses beta distributions to model the stop decision.

**working_countrywide_mixture.stan**: this is the model we fit to the national traffic stop dataset. Data to fit this model is available at https://openpolicing.stanford.edu/data/; we do not analyze this model in detail in the paper, but include it as an illustration of how one might fit a multistate model. 


