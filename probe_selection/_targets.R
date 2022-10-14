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

# Set target options:
tar_option_set(
#   packages = c("tibble"), # packages that your targets need to run
  format = "qs" # default storage format,
  # Set other options as needed.
)

source(paste0("../00_preprocessing/preprocessing.R"))
source(paste0("../01_milo/milo.R"))
source(paste0("prep_probe_selection.R"))

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore")
options(future.globals.maxSize= (200 * 1000)*1024^2) #200GB
Sys.setenv(BASILISK_EXTERNAL_DIR = paste0(Sys.getenv("PROJECT_DIR"), '/.basilisk'))

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
# tar_source()
# source("other_functions.R") # Source other scripts as needed. # nolint

# Replace the target list below with your own:
list(
  tar_target(exp_labelled_other_path, paste0(Sys.getenv("PROJECT_DIR"), '/00_preprocessing/_targets/objects/exp_labelled_other'), format = "file"),
  tar_target(exp_labelled_neuron_path, paste0(Sys.getenv("PROJECT_DIR"), '/00_preprocessing/_targets/objects/exp_RElabelled_neuron'), format = "file"),
  tar_target(pldf_obob5v5_path, paste0(Sys.getenv("PROJECT_DIR"), '/01_milo/_targets/objects/combined_pldf_obob5v5'), format = "file"),
  tar_target(exp_labelled_other, qs::qread(exp_labelled_other_path)),
  tar_target(exp_labelled_neuron, qs::qread(exp_labelled_neuron_path)),
  tar_target(pldf_obob5v5, 
             qs::qread(pldf_obob5v5_path)),
  tar_target(exp_labelled_other_ngo,
             exp_labelled_other %>%
               prep_obj(pldf_obob5v5) %>%
               make_new_seurat_ngo %>%
               sc_transform_fgf1),
  tar_target(exp_labelled_neuron_ngo,
             exp_labelled_neuron %>%
               prep_obj(pldf_obob5v5) %>%
               make_new_seurat_ngo %>%
               sc_transform_fgf1),
  tar_target(exp_ngo,
             merge_neuron_and_other(exp_labelled_neuron_ngo, exp_labelled_other_ngo) %>%
             sc_transform_fgf1),
  tar_target(sce_ngo,
             convert_obj_to_sce(exp_ngo)),
  # tar_target(basilisk_env,
  #            make_basilisk_env()),
  # tar_target(basilisk_roundtrip,
  #            run_basilisk_convert(sce_ngo, '/scratch/nmq407/sce_ngo_h5ad', basilisk_env))
  tar_target(sct_var_feats_10k,
             get_sct_variable_features(exp_ngo, 10000)),
  tar_target(exp_labelled_other_ngo_10kf,
             subset_seurat_obj_genes(exp_labelled_other_ngo, sct_var_feats_10k) %>%
             sc_transform_campbell %>% run_sct_chaser),
  tar_target(exp_labelled_neuron_ngo_10kf,
             subset_seurat_obj_genes(exp_labelled_neuron_ngo, sct_var_feats_10k) %>%
             sc_transform_campbell %>% run_sct_chaser)
             #####working here to make objects that wont change in nfeats
)
