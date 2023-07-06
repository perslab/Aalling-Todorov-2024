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

source("wgcna_functions.R")

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multiprocess")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
# tar_source()
# source("other_functions.R") # Source other scripts as needed. # nolint


milo_clusters_tibble = qs::qread('milo_clusters_tibble.qs')
stage_01 = tar_map(
  values = milo_clusters_tibble,
  names = label,
  tar_target(obj_path, paste0(Sys.getenv("PROJECT_DIR"), milo_obj_path), format = "file"),
  tar_target(obj,
             qs::qread(obj_path),
             deployment='main'),
  tar_target(expr_wgcna,
             obj %>%
             build_wgcna(gamma=20),
             deployment='main'),
  tar_target(sft,
             expr_wgcna %>%
             find_threshold(cluster, corfnc="cor"),
             deployment='worker'),
  tar_target(modules_wgcna,
             sft %>%
             generate_tom(expr_wgcna, sftpwr = ., method = "average", cortype=1) %>% 
             find_modules(expr_wgcna, tom = .,  modsize=30, deepsplit = c(4)) %>% 
             merge_modules(expr_wgcna, ., thresh = 0.1),
             deployment='worker'),
  tar_target(hubgenes,
             get_hubgenes(modules_wgcna, expr_wgcna, label),
             deployment='main'),
  tar_target(hubgenes_summary,
             make_hubgenes_summary(hubgenes),
             deployment='main'),
  tar_target(ME,
             calc_eigengenes(obj, modules_wgcna, expr_wgcna, label),
             deployment='worker'),
  tar_target(gost_result,
             make_gost(hubgenes, label),
             deployment='worker'),
  tar_target(top_table,
             fit_limma_me(ME) %>%
             get_sig_mods),
  tar_target(hubgenes_summary_tt,
             make_hubgenes_tt_summary(hubgenes, top_table),
             deployment='main'),
  tar_target(gost_tt_summary,
             make_gost_tt_summary(gost_result, top_table))
)

combination_recipe = qs::qread('combination_recipe.qs')
stage_02 = tar_eval(
  values = combination_recipe,
  tar_combine(output_name,
              stage_01 %>% tar_select_targets(any_of(obj_names)))
)


run_list = list(
  stage_01,
  stage_02
)

run_list
