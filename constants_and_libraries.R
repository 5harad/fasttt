library(rstan)
library(dplyr)
library(parallel)
library(readr)
library(tidyr)
library(gridExtra)
library(scales)
library(xtable)
library(lubridate)
library(stringr)
library(boot)

# CHANGE THESE PATHS AS NECESSARY TO SUIT YOUR SYSTEM.
# If the appropriate folders don't exist, the code will throw an error. 
base_code_dir = '~/fast-threshold-test/' # the root directory for code. 
base_input_dir = '/share/data/policing-data/mixture_model/input/' # the path where input data for the Stan models is stored. 
base_output_dir = '/share/data/policing-data/mixture_model/output/' # the path where fitted Stan models are stored.
stop_and_frisk_data_path = '/share/data/algobias-data/sqf.RData' # the path for the original stop-and-frisk data. 
stopifnot(file.exists(base_code_dir), 
          file.exists(base_input_dir), 
          file.exists(base_output_dir), 
          file.exists(stop_and_frisk_data_path))

stops_per_precinct_filename = paste0(base_input_dir, 'stops_per_precinct.csv')
setwd(base_code_dir)
white_population_perturbations = c(.2, .3, .4, .5, .6, .7, .8, .9, 1.1, 1.2, 1.5, 1.8, 2, 3, 4, 5, 10)

new_frisk_model_name = 'model_mixture.stan' # checked. 
old_frisk_model_name = 'model_flat.stan' # checked
new_stop_model_name = 'multinomial_stop_and_frisk.stan' # checked. 
old_stop_model_name = 'multinomial_model_flat.stan' # checked. 

our_theme <- function() {
  # Set the base size
  theme_bw(base_size=15) +
    theme(
      # Remove the title
      plot.title=element_blank(),
      # Make the background white
      panel.background=element_rect(fill='white', colour='white'),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      # Minimize margins
      plot.margin=unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
      panel.margin=unit(0.25, "lines"),
      # Tiny space between axis labels and tick labels
      axis.title.x=element_text(margin=ggplot2::margin(t=6.0)),
      axis.title.y=element_text(margin=ggplot2::margin(r=6.0)),
      # Simplify the legend
      legend.key=element_blank(),
      legend.background=element_rect(fill='transparent'),
      legend.title=element_blank()
    )
}
