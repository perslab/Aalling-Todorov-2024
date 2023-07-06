# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed. # nolint
library(future)
library(future.callr)
plan(callr)

# Set target options:
tar_option_set(
#   packages = c("tibble"), # packages that your targets need to run
  format = "qs", # default storage format,
  error = "null",
  retrieval = "worker",
  storage = "worker"
  # Set other options as needed.
)

source('_targets_define.R')

obj_to_load = make_da_results %>%
tar_select_names(starts_with("da_results_02") | contains('pldf'))

tar_load(obj_to_load)

run_list = list(
  tarchetypes::tar_combine(combined_da_results_obob5v5, 
                            make_da_results %>% 
                            tar_select_targets(starts_with("da_results_02")) %>% 
                            tar_select_targets(ends_with("obob5v5"))),
    tarchetypes::tar_combine(combined_da_results_obob14v14, 
                          make_da_results %>% 
                          tar_select_targets(starts_with("da_results_02")) %>% 
                          tar_select_targets(ends_with("obob14v14"))),
    tarchetypes::tar_combine(combined_da_results_obobBL6d5, 
                          make_da_results %>% 
                          tar_select_targets(starts_with("da_results_02")) %>% 
                          tar_select_targets(ends_with("obobBL6d5"))),
    tarchetypes::tar_combine(combined_da_results_obobFGF1BL6d5, 
                          make_da_results %>% 
                          tar_select_targets(starts_with("da_results_02")) %>% 
                          tar_select_targets(ends_with("obobFGF1BL6d5"))),
    tarchetypes::tar_combine(combined_pldf_obob5v5, 
                             make_da_results %>% 
                             tar_select_targets(contains(tar_meta(fields = NULL) %>% 
                              filter(str_detect(name, 'pldf')) %>% 
                              filter(str_detect(name, 'obob5v5')) %>% 
                              pull(name)))),
    tarchetypes::tar_combine(combined_pldf_obob14v14, 
                             make_da_results %>% 
                             tar_select_targets(contains(tar_meta(fields = NULL) %>% 
                              filter(str_detect(name, 'pldf')) %>% 
                              filter(str_detect(name, 'obob14v14')) %>% 
                              pull(name)))),
    tarchetypes::tar_combine(combined_pldf_obobBL6d5, 
                          make_da_results %>% 
                          tar_select_targets(contains(tar_meta(fields = NULL) %>% 
                          filter(str_detect(name, 'pldf')) %>% 
                          filter(str_detect(name, 'obobBL6d5')) %>% 
                          pull(name)))),
    tarchetypes::tar_combine(combined_pldf_obobFGF1BL6d5, 
                          make_da_results %>% 
                          tar_select_targets(contains(tar_meta(fields = NULL) %>% 
                          filter(str_detect(name, 'pldf')) %>% 
                          filter(str_detect(name, 'obobFGF1BL6d5')) %>% 
                          pull(name))))
)

run_list

