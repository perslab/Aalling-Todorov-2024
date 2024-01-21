# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed. # nolint
library(tidyverse)
library(future)
library(future.callr)
plan(callr)


# Set target options:
tar_option_set(
  packages = c("tidyverse", "Seurat", "miloR"), # packages that your targets need to run
  format = "qs", # default storage format,
  error = "null",
  retrieval = "worker",
  storage = "worker",
  priority = 0.1,
  memory = "transient",
  garbage_collection = TRUE,
  cue = tar_cue(
      mode = c("thorough", "always", "never"),
      command = TRUE,
      depend = TRUE,
      format = TRUE,
      repository = TRUE,
      iteration = TRUE,
      file = TRUE
)
)

source("../00_preprocessing/splitwrapper.R")
source(paste0("../00_preprocessing/preprocessing.R"))
source(paste0("../01_milo/milo.R"))
source(paste0("../01_milo/milo_plotting.R"))
source(paste0("../00_cellbender/cellbender_preprocess.R"))
source("../00_cellbender/process_seurat.R")
source('../01_milo_cellbender/milo_cellbender.R')

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore",
        clustermq.ssh.timeout=36000,
        clustermq.worker.timeout=36000,
        clustermq.error.timeout=36000,
        clustermq.ssh.log='clustermq_sshlog.log'
        )

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
# tar_source()
# source("other_functions.R") # Source other scripts as needed. # nolint


#lvl1
clusters_tibble_lvl1 = qs::qread('clusters_tibble_lvl1.qs')
make_split_objs_lvl1 = tar_map(
  values = clusters_tibble_lvl1,
  names = label,
  tar_target(obj,
             object %>%
             set_labels_to_lvl1 %>%
             prep_obj_for_milo_cb_v01 %>%
             set_batch_to_lane %>% # do not set batch to lane for cluster splits, will error on design or model matrix. reset later
             prep_obj_for_milo_cb_v01(set_orig.batch = FALSE) %>%
             subset_exp_by_time(day) %>%
             single_split(cluster) %>%
             reconsitute_rna_seurat %>%
#              process_seurat(method = "integrate", 
#                             batch ="batch", # batch to batch
#                             dims = 30, res = 0.8,
#                             k.weight = 40)
             process_seurat(method = "integrate", 
                            batch = batch,
                            dims = dims, res = 0.8,
                            k.weight = k.weight,
                            k.anchor = k.anchor) %>% 
             reset_orig.batch %>% #reset batch for downstream applications
             prep_obj_for_milo_cb_v01,
             priority=0.99
             ),
    tar_target(batchy_score_tibble,
               obj %>%
               calc_batchy_score(obj=.,
                            batch = batch,
                        dims = dims,
                        reduction = 'pca',
                        k.param = 40) %>% 
                `[[`('logp_bxsqp') %>% sum %>% tibble(name=label, batchy_score=.),
               priority=0.999)
)


stage_01 = list(
  tar_target(exp_labelled_other_path, paste0(Sys.getenv("PROJECT_DIR"), '/00_cellbender/_targets/objects/obj_cb_other_01'), format = "file"),
  tar_target(exp_labelled_neuron_path, paste0(Sys.getenv("PROJECT_DIR"), '/00_cellbender/_targets/objects/obj_cb_neuron_00'), format = "file"),
  tar_target(exp_labelled_other, qs::qread(exp_labelled_other_path)),
  tar_target(exp_labelled_neuron, qs::qread(exp_labelled_neuron_path)),
  make_split_objs_lvl1
  )


stage_02a = list(
  tar_target(clusters_tibble_path, 'clusters_tibble_lvl1.qs', format = "file"),
  tar_target(clusters_tibble_obj, qs::qread(clusters_tibble_path))
)

combination_recipe = qs::qread('combination_recipe.qs')
stage_02b = tar_eval(
  values = combination_recipe,
  tar_combine(output_name,
              make_split_objs_lvl1 %>% tar_select_targets(any_of(da_names_idx)))
)
combination_recipe = qs::qread('combination_recipe.qs')
stage_02c = tar_map(
  values = combination_recipe,
  names = old.label,
  tar_target(combined_batchy_score,
             combination_tibble_obj  %>%
                 rename(label = name) %>%
                 left_join(clusters_tibble_obj, by='label')
            )
)



stage_02 = list(stage_02a, stage_02b, stage_02c)


run_list = list(
  stage_01,
  stage_02,
  list()
)

run_list
