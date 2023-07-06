# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(tidyverse)
library(targets)
library(tarchetypes) # Load other packages as needed. # nolint
# library(future)
# library(future.callr)
# plan(callr)

# Set target options:
tar_option_set(
  packages = c("Seurat", "miloR", "SingleCellExperiment", "tidyverse", "ggplot2"), # packages that your targets need to run
  format = "qs", # default storage format,
  error = "null",
  retrieval = "worker",
  storage = "worker"
  # Set other options as needed.
)


source("../00_preprocessing/splitwrapper.R")
source(paste0("../00_preprocessing/preprocessing.R"))
source(paste0("../01_milo/milo.R"))
source(paste0("../01_milo/milo_plotting.R"))
source(paste0("../00_bmp_preprocessing/bmp_preprocessing.R"))

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
# tar_source()
# source("other_functions.R") # Source other scripts as needed. # nolint


clusters_tibble = qs::qread('clusters_tibble_bmp1_bmp_only.qs')
make_split_objs_bmp1_only = tar_map(
  values = clusters_tibble,
  names = label,
  tar_target(obj,
             object %>%
             subset_exp_by_time(day) %>%
             single_split(cluster) %>%
             sc_transform_fgf1_nobatch %>%
             prep_obj_for_milo)
)
clusters_tibble = qs::qread('clusters_tibble_bmp1_unimod.qs')
make_split_objs_unimod = tar_map(
  values = clusters_tibble,
  names = label,
  tar_target(obj,
             object %>%
             subset_exp_by_time(day) %>%
             single_split(cluster) %>%
             project_bmp1_into_fgf1_unimod %>%
             set_unimod_embedding %>%
             prep_obj_for_milo)
)
clusters_tibble = qs::qread('clusters_tibble_bmp1_sct.qs')
make_split_objs_sct = tar_map(
  values = clusters_tibble,
  names = label,
  tar_target(obj,
             object %>%
             subset_exp_by_time(day) %>%
             single_split(cluster) %>%
             sc_transform_fgf1 %>%
             prep_obj_for_milo)
)
clusters_tibble = qs::qread('clusters_tibble_bmp1_sct_forb.qs')
make_split_objs_sct_forb = tar_map(
  values = clusters_tibble,
  names = label,
  tar_target(obj,
             object %>%
             subset_exp_by_time(day) %>%
             single_split(cluster) %>%
             sc_transform_fgf1_forb %>%
             prep_obj_for_milo)
)
clusters_tibble = qs::qread('clusters_tibble_bmp1_sct_batchvar_all.qs')
make_split_objs_sct_batchvar_all = tar_map(
  values = clusters_tibble,
  names = label,
  tar_target(obj,
             object %>%
             subset_exp_by_time(day) %>%
             single_split(cluster) %>%
             sc_transform_fgf1_bmp1_batchvar_all %>%
             prep_obj_for_milo)
)
clusters_tibble = qs::qread('clusters_tibble_bmp1_sct_batchvar_forb.qs')
make_split_objs_sct_batchvar_forb = tar_map(
  values = clusters_tibble,
  names = label,
  tar_target(obj,
             object %>%
             subset_exp_by_time(day) %>%
             single_split(cluster) %>%
             sc_transform_fgf1_bmp1_batchvar_forb %>%
             prep_obj_for_milo)
)


clusters_tibble = qs::qread('clusters_tibble_bmp1.qs')
make_split_milos = tar_map(
  values = clusters_tibble,
  names = label,
  tar_target(milo, 
             convert_obj_to_sce(seurat_obj) %>%
             make_milo,
             priority=1),
  tar_target(design_df,
             make_design_df(milo),
             priority=1),
  tar_target(mm,
             make_model_matrix_bmp1(design_df),
             priority=1)
  
)


da_recipe = qs::qread('da_recipe_bmp1.qs')
make_da_results = tar_map(
  values = da_recipe,
  names = label,
  tar_target(da_results_00,
             get_da_results(milo_obj=milo_obj,
                            model_matrix=mm, design_df=design_df, model_contrasts_str=contrast)),
  tar_target(da_results_01,
             add_polarity_to_da_results(da_results_00, spatial_fdr_cutoff=0.1) %>%
             annotate_nhoods(milo_obj, col="labels")),
  tar_target(nhm,
             make_nhm(milo_obj),
             cue=tar_cue(file=FALSE)),
  tar_target(da_results_02,
             annotate_nhood_counts(da_results_01, nhm) %>%
             combine_splits(., nhm)),
  tar_target(nonzero_mean_df,
             make_logFC_df(da_results_02) %>%
             make_logFC_cells_df(nhm, .) %>%
             make_nzdf(.),
             cue=tar_cue(file=FALSE)),
  tar_target(nhood_summary,
             make_nhood_summary(da_results_02, nhm) %>% check_split_validity()),
  tar_target(barcode_df,
             make_barcode_polarity_df(nonzero_mean_df, nhood_summary)),
  tar_target(pldf,
             make_pldf(barcode_df, milo_obj))

)


combination_recipe = qs::qread('combination_recipe_bmp1.qs')
stage_02 = tar_eval(
  values = combination_recipe,
  tar_combine(output_name,
              make_da_results %>% tar_select_targets(all_of(da_names)))
)




stage_01 = list(
  tar_target(bmp1_fgf1_mapped_path, paste0(Sys.getenv("PROJECT_DIR"), '/00_bmp_preprocessing/_targets/objects/bmp1_fgf1_mapped'), format = "file"),
  tar_target(bmp1_fgf1_mapped, qs::qread(bmp1_fgf1_mapped_path)),
  tar_target(bmp1_only_mapped,
             subset(x = bmp1_fgf1_mapped, subset = batch == 'Batch 4') %>%
             sc_transform_fgf1_nobatch),
  make_split_objs_bmp1_only,
  make_split_objs_unimod,
  make_split_objs_sct,
  make_split_objs_sct_forb,
  make_split_objs_sct_batchvar_all,
  make_split_objs_sct_batchvar_forb,
  make_split_milos,
  make_da_results
  )



run_list = list(
  stage_01,
  stage_02
)

run_list
