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
  format = "qs", # default storage format,
  error = "continue",
  # Set other options as needed.
)

source(paste0("../00_preprocessing/preprocessing.R"))
source(paste0("../01_milo/milo.R"))
source(paste0("genebasisr.R"))

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore")
options(future.globals.maxSize= (200 * 1000)*1024^2) #200GB

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
# tar_source()
# source("other_functions.R") # Source other scripts as needed. # nolint

# Replace the target list below with your own:
list(
  tar_target(sce_ngo_path, paste0(Sys.getenv("PROJECT_DIR"), '/probe_selection/_targets/objects/sce_ngo'), format = "file"),
  tar_target(sce_ngo_00, qs::qread(sce_ngo_path)),
  tar_target(sce_ngo_01, gbr_retain_informative_genes(sce_ngo_00)),
  tar_target(genes_stat_100,
             gbr_genes_stat(sce_ngo_01, 100),
             cue = tar_cue(mode = "never")), ### remember to get rid of this cue when you can let it run
  tar_target(sce_ngo_02, gbr_retain_informative_genes(sce_ngo_00, n=10000)),
  # tar_target(genes_stat_100_10k,
  #            gbr_genes_stat(sce_ngo_02, 100)),
  tar_target(genes_100,
             get_selected_genes(genes_stat_100)),
  tar_target(ctm_00_100,
             gbr_get_celltype_mapping(sce_ngo_02, genes_stat_100$gene, 'polar_label'))
  
  # tar_target(sce_ngo_path, paste0(Sys.getenv("PROJECT_DIR"), '/probe_selection/_targets/objects/exp_ngo'), format = "file"),
  # tar_target(sce_ngo_00, qs::qread(sce_ngo_path))
)
