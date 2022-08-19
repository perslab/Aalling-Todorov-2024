# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
# library(future)
# library(future.callr)
library(tarchetypes)
library(tidyr)
# plan(callr)

# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
#   packages = c("tibble"), # packages that your targets need to run
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore")

#globalsmax
options(future.globals.maxSize= (200 * 1000)*1024^2) #200GB

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Load the R scripts with your custom functions:
lapply(list.files("R", full.names = TRUE, recursive = TRUE), source)
# source("other_functions.R") # Source other scripts as needed. # nolint

source(paste0("preprocessing.R"))

tar_option_set(packages = c("readr", "dplyr", "ggplot2"))

values <- tibble(
  campbell_subset = rlang::syms(c("campbell_sct_neurons_00", "campbell_sct_other_00")),
  exp_subset = rlang::syms(c("exp_neurons_00", "exp_other_00")),
  names = c("neurons", "other")
)
transfer_campbell_labels_pipeline = tar_map(
  values = values,
  names = "names",
  tar_target(campbell_sct_sub_01, sc_transform_campbell(campbell_subset)),
  tar_target(exp_sub_01, sc_transform_fgf1(exp_subset)),
  tar_target(transfer_anchors_campbell_sub, find_anchors_campbell(exp_sub_01, campbell_sct_sub_01)),
  tar_target(exp_sub_02, transfer_campbell_labels(exp_sub_01, campbell_sct_sub_01, transfer_anchors_campbell_sub)),
  tar_target(exp_sub_03, run_pca(exp_sub_02)),
  tar_target(exp_sub_04, run_umap(exp_sub_03))
)
transfer_campbell_labels_pipeline = list(transfer_campbell_labels_pipeline)


list(    
  tar_target(exp22_path, paste0(PROJECT_DIR, 'data/received/SCOP_2022_0189/Output/data/aggregated-filtered/20220622_mouse/20220622_mouse_rna-seurat.rds'), format = "file"),
  tar_target(exp23_path, paste0(PROJECT_DIR, 'data/received/SCOP_2022_0189/Output/data/aggregated-filtered/20220623_mouse/20220623_mouse_rna-seurat.rds'), format = "file"),
  tar_target(exp24_path, paste0(PROJECT_DIR, 'data/received/SCOP_2022_0189/Output/data/aggregated-filtered/20220624_mouse/20220624_mouse_rna-seurat.rds'), format = "file"),
  tar_target(exp22, readRDS(exp22_path)),
  tar_target(exp23, readRDS(exp23_path)),
  tar_target(exp24, readRDS(exp24_path)),
  tar_target(exp_list_00, make_exp_list(exp22, exp23, exp24)),
  tar_target(features, select_integration_featuers(exp_list_00)),
  tar_target(exp_list_01, prep_sc_transform(exp_list_00, features)),
  tar_target(anchors, find_integration_anchors(exp_list_01, features)),
  tar_target(exp_00, integrate_data(anchors)),
  tar_target(exp_01, run_pca(exp_00)),
  tar_target(exp_02, run_umap(exp_01)), # this is the completly sctransformed data from all 3
  tar_target(meta_path, paste0(PROJECT_DIR, 'data/meta/metadata.csv'), format = "file"),
  tar_target(meta, read_metadata(meta_path)),
  tar_target(exp_03, add_meta_to_exp(exp_02, meta)),
  tar_target(campbell, load_campbell()),
  tar_target(campbell_sct, sc_transform_campbell(campbell)),
  tar_target(transfer_anchors_campbell, find_anchors_campbell(exp_03, campbell_sct)),
  tar_target(exp_04, transfer_campbell_labels(exp_03, campbell, transfer_anchors_campbell)),
  tar_target(exp_neurons_00, subset_neurons(exp_04)),
  tar_target(exp_other_00, subset_other(exp_04)),
  tar_target(campbell_sct_neurons_00, subset_neurons(campbell_sct)),
  tar_target(campbell_sct_other_00, subset_other(campbell_sct)),
  transfer_campbell_labels_pipeline
)
