# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
# library(future)
# library(future.callr)
library(tarchetypes)
# plan(callr)

# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  packages = c("Seurat", "SingleCellExperiment", "scRNAseq", "edgeR", "tidyverse", "ggplot2"), # packages that your targets need to run
  format = "qs",
  error = "null",
  retrieval = "worker",
  storage = "worker"
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

source(paste0("../00_preprocessing/preprocessing.R"))
source("../00_preprocessing/mapCamp.R")
source("../00_preprocessing/build_edger.R")
source("../00_preprocessing/splitwrapper.R")
source("bmp_preprocessing.R")


list(    
  tar_target(exp_bmp_path, paste0(PROJECT_DIR, 'data/received/SCOP_2022_0214/Output/data//aggregated-filtered//221115_brain/221115_brain_rna-seurat.rds'), format = "file"),
  tar_target(bmp_meta_path, paste0(PROJECT_DIR, 'data/meta/0214_scnRNA_metadata_NA041122.xlsx'), format='file'),
  tar_target(exp_bmp, readRDS(exp_bmp_path), cue=tar_cue(file=FALSE)),
  tar_target(bmp_meta_00, readxl::read_xlsx(bmp_meta_path)),
  tar_target(exp_bmp_00, 
             exp_bmp %>%
             bmp_combine_astrocytes),
  tar_target(exp_bmp_01, add_meta_to_bmp_obj(exp_bmp_00, bmp_meta_00)),
  tar_target(exp_bmp_01_all, mapCamp(exp_bmp_01, class="all")),
  tar_target(exp_bmp_01_other, mapCamp(exp_bmp_01, class="other")),
  tar_target(exp_fgf1_other_path, paste0(PROJECT_DIR, '00_preprocessing/_targets/objects/exp_labelled_other'), format = "file"),
  tar_target(exp_fgf1_other, 
             qs::qread(exp_fgf1_other_path) %>%
             add_umap_model),
  tar_target(exp_bmp_fgf1_other,
             merge(exp_bmp_01_other, y=exp_fgf1_other) %>%
             sc_transform_fgf1),
  tar_target(exp_bmp_fgf1_other_forb,
             merge(exp_bmp_01_other, y=exp_fgf1_other) %>%
             sc_transform_fgf1_forb),
  tar_target(exp_bmp_fgf1_sct_batchvar_all,
             merge(exp_bmp_01_other, y=exp_fgf1_other) %>%
             sc_transform_fgf1_bmp1_batchvar_all),
  tar_target(exp_bmp_fgf1_sct_batchvar_forb,
             merge(exp_bmp_01_other, y=exp_fgf1_other) %>%
             sc_transform_fgf1_bmp1_batchvar_forb),
  tar_target(exp_list_00, make_exp_list(exp_bmp_01_other, exp_fgf1_other), cue=tar_cue(file=FALSE)),
  tar_target(features, select_integration_featuers(exp_list_00), cue=tar_cue(file=FALSE)),
  tar_target(exp_bmp_fgf1_other_integrated, 
             exp_list_00 %>%
                prep_sc_transform(., features) %>%
                find_integration_anchors(., features) %>%
                integrate_data %>%
                run_sct_chaser,
             cue=tar_cue(file=FALSE)),
  tar_target(bmp1_fgf1_mapped,
             exp_bmp_fgf1_other %>%
             project_bmp1_into_fgf1_unimod %>%
             set_unimod_embedding)
)
