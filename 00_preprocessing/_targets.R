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
  format = "qs" # default storage format
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
source("mapCamp.R")
source("build_edger.R")
source("splitwrapper.R")
source("gprofiler_requests.R")

tar_option_set(packages = c("readr", "dplyr", "ggplot2"))

tcl_values <- tibble(
  subset_class = c("neuron", "other"),
  names = c("neuron", "other")
)
transfer_campbell_labels_pipeline = tar_map(
  values = tcl_values,
  names = "names",
  tar_target(exp_subset,
             subset_exp(exp_04, class=subset_class),
             cue=tar_cue(file=FALSE)),
  tar_target(exp_subset_sct, sc_transform_fgf1(exp_subset), cue=tar_cue(file=FALSE)),
  tar_target(exp_labelled,
             mapCamp(exp_subset_sct, class=subset_class))
  
)
transfer_campbell_labels_pipeline = list(transfer_campbell_labels_pipeline)


subset_values = tibble(
  obj_to_subset = rlang::syms(c("exp_RElabelled_neuron", "exp_RElabelled_neuron", "exp_labelled_other", "exp_labelled_other")),
  strain = c("obob", "BL6", "obob", "BL6"),
  names = c("neuron_obob", "neuron_BL6", "other_obob", "other_BL6")
)
make_strain_subsets = tar_map(
  values = subset_values,
  names = 'names',
  tar_target(exp, subset_exp_by_strain(obj_to_subset, strain))
)

design_edger_obob = c("treatment", "time")
batch_edger_obob = c("batch")
contrasts_list_obob = c("groupFGF1.Day5-groupVeh_PF.Day5",
                   "groupFGF1.Day14-groupVeh_PF.Day14")
deg_values <- tibble(
  input_objs_deg = rlang::syms(c("exp_other_obob", "exp_other_obob", "exp_neuron_obob", "exp_neuron_obob")),
  split_by_col = c("labels", "predicted.id", "labels", "predicted.id"),
  names = c("obob_other_labels", "obob_other_predicted.id", "obob_neuron_labels", "obob_neuron_predicted.id")
)
find_degs = tar_map(
  values = deg_values,
  names = "names",
  tar_target(edger,
             splitwrapper(input_objs_deg, split.by=split_by_col) %>% 
             map(~build_edger(.x, 10, design = design_edger_obob, batch = batch_edger_obob)) %>%
             add_cluster_names(),
             cue=tar_cue(file=FALSE)),
  tar_target(qlf,
             get_qlf(edger, contrasts_list_obob)),
  tar_target(top_tags,
             get_toptags(qlf))
)

gpdf_values = tibble(
  tt_obj = rlang::syms(c("top_tags_obob_neuron_labels", "top_tags_obob_other_labels")),
  names = c("neuron_obob", "other_obob")
)
get_gpdfs = tar_map(
  values = gpdf_values,
  names = 'names',
  tar_target(gpdf, get_gosts_df(tt_obj))
)

bind_gpdfs = list(
    tar_target(gpdfb_other, bind_gosts_dfs(gpdf_other_obob)),
    tar_target(gpdfb_neuron, bind_gosts_dfs(gpdf_neuron_obob))
)

bind_top_tags = list(
    tar_target(ttb_neuron, bind_top_tags_list(top_tags_obob_neuron_labels)),
    tar_target(ttb_other, bind_top_tags_list(top_tags_obob_other_labels))
)


list(    
  tar_target(exp22_path, paste0(PROJECT_DIR, 'data/received/SCOP_2022_0189/Output/data/aggregated-filtered/20220622_mouse/20220622_mouse_rna-seurat.rds'), format = "file"),
  tar_target(exp23_path, paste0(PROJECT_DIR, 'data/received/SCOP_2022_0189/Output/data/aggregated-filtered/20220623_mouse/20220623_mouse_rna-seurat.rds'), format = "file"),
  tar_target(exp24_path, paste0(PROJECT_DIR, 'data/received/SCOP_2022_0189/Output/data/aggregated-filtered/20220624_mouse/20220624_mouse_rna-seurat.rds'), format = "file"),
  tar_target(exp22, readRDS(exp22_path), cue=tar_cue(file=FALSE)),
  tar_target(exp23, readRDS(exp23_path), cue=tar_cue(file=FALSE)),
  tar_target(exp24, readRDS(exp24_path), cue=tar_cue(file=FALSE)),
  tar_target(exp_list_00, make_exp_list(exp22, exp23, exp24), cue=tar_cue(file=FALSE)),
  tar_target(features, select_integration_featuers(exp_list_00), cue=tar_cue(file=FALSE)),
  tar_target(exp_list_01, prep_sc_transform(exp_list_00, features), cue=tar_cue(file=FALSE)),
  tar_target(anchors, find_integration_anchors(exp_list_01, features), cue=tar_cue(file=FALSE)),
  tar_target(exp_00, integrate_data(anchors), cue=tar_cue(file=FALSE)),
  tar_target(exp_01, run_pca(exp_00), cue=tar_cue(file=FALSE)),
  tar_target(exp_02, run_umap(exp_01), cue=tar_cue(file=FALSE)), # this is the completly sctransformed data from all 3
  tar_target(meta_path, paste0(PROJECT_DIR, 'data/meta/metadata.csv'), format = "file"),
  tar_target(meta, read_metadata(meta_path)),
  tar_target(exp_03, add_meta_to_exp(exp_02, meta), cue=tar_cue(file=FALSE)),
  tar_target(exp_04, mapCamp(exp_03, class="all")),
  transfer_campbell_labels_pipeline,
  tar_target(exp_RElabelled_neuron, do_neuron_cluster_surgery(exp_labelled_neuron)),
  make_strain_subsets,
  find_degs,
  get_gpdfs,
  bind_gpdfs,
  bind_top_tags
)
