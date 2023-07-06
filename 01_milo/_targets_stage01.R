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


run_list = list(    
  tar_target(exp_labelled_other_path, paste0(Sys.getenv("PROJECT_DIR"), '/00_preprocessing/_targets/objects/exp_labelled_other'), format = "file"),
  tar_target(exp_labelled_neuron_path, paste0(Sys.getenv("PROJECT_DIR"), '/00_preprocessing/_targets/objects/exp_RElabelled_neuron'), format = "file"),
  tar_target(exp_labelled_other, qs::qread(exp_labelled_other_path)),
  tar_target(exp_labelled_neuron, qs::qread(exp_labelled_neuron_path)),
  make_split_objs,
  make_da_results)

run_list

