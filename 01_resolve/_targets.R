# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed. # nolint
library(future)
library(future.callr)
library(tidyverse)
plan(callr)

# Set target options:
tar_option_set(
  packages = c("tidyverse"), # packages that your targets need to run
  format = "qs", # default storage format,
  error = "null",
  retrieval = "worker",
  storage = "worker"
  # Set other options as needed.
)

PROJECT_DIR = Sys.getenv('PROJECT_DIR')
source(paste0(PROJECT_DIR, "/00_preprocessing/preprocessing.R"))
source(paste0(PROJECT_DIR, "/code/resolve2xe/LoadResolveBaysor.R"))

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore",
        clustermq.ssh.timeout=36000,
        clustermq.worker.timeout=36000,
        clustermq.error.timeout=36000
        )

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
# tar_source()
# source("other_functions.R") # Source other scripts as needed. # nolint

ingest_tibble = qs::qread('00_ingest_tibble.qs')
stage_01 = tar_map(
  values = ingest_tibble,
  names = sample_name,
  tar_target(sample_path,
             paste0(PROJECT_DIR, '/data/resolve/xe/32810-1377-slide3_', sample_name, '-1/'), 
             format='file'),
  tar_target(obj, 
             sample_path %>% 
             LoadResolveBaysor %>%
             AddMetaData(sample_name, col.name='sample_name') %>%
             AddMetaData(treatment, col.name='treatment') %>%
             AddMetaData(strain, col.name='strain') %>%
             AddMetaData(strain, col.name='time'),
             packages = c("Seurat", "tidyverse"),
             priority=1)
)

merge_samples = list(
  tar_target(obj_merged, 
             merge(obj_A1, list(obj_A2, obj_B1, obj_B2, obj_C1, obj_C2, obj_D1, obj_D2)),
             packages = c("Seurat", "tidyverse"))
)
qc_filters = list(
  tar_target(obj_merged_sct,
             obj_merged %>% filter_down_cells %>% sc_transform_resolve,
             packages = c("Seurat", "tidyverse"))
)
stage_01 = list(stage_01, merge_samples, qc_filters)


run_list = list(
  stage_01
)

run_list
