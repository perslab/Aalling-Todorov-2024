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
  packages = c("tidyverse", "ggplot2", "patchwork", "Seurat", "miloR"), # packages that your targets need to run
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
source('milo_cellbender.R')

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


#class
clusters_tibble_class = qs::qread('clusters_tibble_class.qs')
make_split_objs_class = tar_map(
  values = clusters_tibble_class,
  names = label,
  tar_target(obj,
             object %>%
             set_labels_to_class %>%
             prep_obj_for_milo_cb_v01 %>%
#              set_batch_to_lane %>%
#              prep_obj_for_milo_cb_v01(set_orig.batch = FALSE) %>%
             subset_exp_by_time(day) %>%
             single_split(cluster) %>%
             reconsitute_rna_seurat %>%
             process_seurat(method = "integrate", 
                            batch ="batch",
                            k.anchor=25,
                            dims = 30, res = 0.8),
             priority=0.99),
  tar_target(milo, 
             convert_obj_to_sce(obj) %>%
             make_milo(k=40, d=30, prop=0.1),
             priority=0.8),
  tar_target(milo_index_tibble, 
             make_milo_index_tibble(milo, cluster),
             priority=0.8),
  tar_target(design_df,
             make_design_df_cb(milo),
             priority=0.8),
  tar_target(mm,
             make_model_matrix(design_df),
             priority=0.8)
  
)
da_recipe_class = qs::qread('da_recipe_class.qs')
make_da_results_class = tar_map(
  values = da_recipe_class,
  names = label,
  tar_target(da_results_00,
             get_da_results(milo_obj=milo_obj,
                            model_matrix=mm, design_df=design_df, model_contrasts_str=contrast)),
  tar_target(da_results_01,
             add_polarity_to_da_results(da_results_00, spatial_fdr_cutoff=0.1) %>%
             annotate_nhoods(milo_obj, col="labels")),
  tar_target(nhm,
             make_nhm(milo_obj)),
  tar_target(da_results_02,
             annotate_nhood_counts(da_results_01, nhm) # %>%
             # combine_splits(., nhm)
            ),
  tar_target(da_results_02_idx,
             milo_index_obj %>%
             select(-any_of('labels')) %>%
             left_join(da_results_02, ., by='Nhood'),
             packages=c("tidyverse", "ggplot2", "patchwork", "Seurat")),
  tar_target(nonzero_mean_df,
             make_logFC_df(da_results_02) %>%
             make_logFC_cells_df(nhm, .) %>%
             make_nzdf(.)),
  tar_target(nhood_summary,
             make_nhood_summary(da_results_02, nhm) %>% check_split_validity()),
  tar_target(barcode_df,
             make_barcode_polarity_df(nonzero_mean_df, nhood_summary)),
  tar_target(pldf,
             make_pldf(barcode_df, milo_obj))

)

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
             run_until_success(function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                              batch = "batch",
                                              dims = 15, res = 0.8,
                                              k.weight = 100, k.anchor=25),
                               function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                              batch = "batch",
                                              dims = 15, res = 0.8,
                                              k.weight = 50, k.anchor=25),
                               function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                              batch = "batch",
                                              dims = 15, res = 0.8,
                                              k.weight = 20, k.anchor=25),
                               function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                              batch = "batch",
                                              dims = 15, res = 0.8,
                                              k.weight = 10, k.anchor=25),
                               function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                              batch = "orig.batch",
                                              dims = 15, res = 0.8,
                                              k.weight = 100, k.anchor=25),
                               function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                              batch = "orig.batch",
                                              dims = 15, res = 0.8,
                                              k.weight = 50, k.anchor=25),
                               function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                              batch = "orig.batch",
                                              dims = 15, res = 0.8,
                                              k.weight = 20, k.anchor=25),
                               function() subset(eval(.)) %>% Seurat::SCTransform(assay='RNA',
                                               method="glmGamPoi",
                                               vars.to.regress= 'batch',
                                               vst.flavor="v2",
                                               verbose=TRUE) %>% run_sct_chaser,
                               function() subset(eval(.)) %>% Seurat::SCTransform(assay='RNA',
                                               method="glmGamPoi",
                                               vars.to.regress= 'orig.batch',
                                               vst.flavor="v2",
                                               verbose=TRUE) %>% run_sct_chaser
                              ) %>% #reset batch for downstream applications
             reset_orig.batch %>% 
             prep_obj_for_milo_cb_v01,
             priority=0.99
             ),
  tar_target(milo, 
             convert_obj_to_sce(obj) %>%
             make_milo(k=40, d=15, prop=0.1),
             priority=0.8),
  tar_target(milo_index_tibble, 
             make_milo_index_tibble(milo, cluster),
             priority=0.8),
  tar_target(design_df,
             make_design_df(milo),
             priority=0.8),
  tar_target(mm,
             make_model_matrix(design_df),
             priority=0.8)
  
)
da_recipe_lvl1 = qs::qread('da_recipe_lvl1.qs')
make_da_results_lvl1 = tar_map(
  values = da_recipe_lvl1,
  names = label,
  tar_target(da_results_00,
             get_da_results(milo_obj=milo_obj,
                            model_matrix=mm, design_df=design_df, model_contrasts_str=contrast)),
  tar_target(da_results_01,
             add_polarity_to_da_results(da_results_00, spatial_fdr_cutoff=0.1) %>%
             annotate_nhoods(milo_obj, col="labels")),
  tar_target(nhm,
             make_nhm(milo_obj)),
  tar_target(da_results_02,
             annotate_nhood_counts(da_results_01, nhm) # %>%
             # combine_splits(., nhm)
            ),
  tar_target(da_results_02_idx,
             milo_index_obj %>%
             select(-any_of('labels')) %>%
             left_join(da_results_02, ., by='Nhood'),
             packages=c("tidyverse", "ggplot2", "patchwork", "Seurat")),
  tar_target(nonzero_mean_df,
             make_logFC_df(da_results_02) %>%
             make_logFC_cells_df(nhm, .) %>%
             make_nzdf(.)),
  tar_target(nhood_summary,
             make_nhood_summary(da_results_02, nhm) %>% check_split_validity()),
  tar_target(barcode_df,
             make_barcode_polarity_df(nonzero_mean_df, nhood_summary)),
  tar_target(pldf,
             make_pldf(barcode_df, milo_obj))

)

#labels_chunk
clusters_tibble_labels_chunk = qs::qread('clusters_tibble_labels_chunk.qs')
make_split_objs_labels_chunk = tar_map(
  values = clusters_tibble_labels_chunk,
  names = label,
  tar_target(obj,
             object %>%
             set_labels_to_labels_chunk %>%
             prep_obj_for_milo_cb_v01 %>%
             set_batch_to_lane %>%
             prep_obj_for_milo_cb_v01(set_orig.batch = FALSE) %>%
             subset_exp_by_time(day) %>%
             single_split(cluster) %>%
             reconsitute_rna_seurat %>%
              run_until_success(function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                              batch = "batch",
                                              dims = 15, res = 0.8,
                                              k.weight = 100, k.anchor=25),
                               function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                              batch = "batch",
                                              dims = 15, res = 0.8,
                                              k.weight = 50, k.anchor=25),
                               function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                              batch = "batch",
                                              dims = 15, res = 0.8,
                                              k.weight = 20, k.anchor=25),
                               function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                              batch = "batch",
                                              dims = 15, res = 0.8,
                                              k.weight = 10, k.anchor=25),
                               function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                              batch = "orig.batch",
                                              dims = 15, res = 0.8,
                                              k.weight = 100, k.anchor=25),
                               function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                              batch = "orig.batch",
                                              dims = 15, res = 0.8,
                                              k.weight = 50, k.anchor=25),
                               function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                              batch = "orig.batch",
                                              dims = 15, res = 0.8,
                                              k.weight = 20, k.anchor=25),
                               function() subset(eval(.)) %>% Seurat::SCTransform(assay='RNA',
                                               method="glmGamPoi",
                                               vars.to.regress= 'batch',
                                               vst.flavor="v2",
                                               verbose=TRUE) %>% run_sct_chaser,
                               function() subset(eval(.)) %>% Seurat::SCTransform(assay='RNA',
                                               method="glmGamPoi",
                                               vars.to.regress= 'orig.batch',
                                               vst.flavor="v2",
                                               verbose=TRUE) %>% run_sct_chaser
                              ) %>%
             reset_orig.batch %>%
             prep_obj_for_milo_cb_v01,
             priority=0.99),
  tar_target(milo, 
             convert_obj_to_sce(obj) %>%
             make_milo(k=40, d=15, prop=0.1),
             priority=0.8),
  tar_target(milo_index_tibble, 
             make_milo_index_tibble(milo, cluster),
             priority=0.8),
  tar_target(design_df,
             make_design_df_cb(milo),
             priority=0.8),
  tar_target(mm,
             make_model_matrix(design_df),
             priority=0.8)
  
)
da_recipe_labels_chunk = qs::qread('da_recipe_labels_chunk.qs')
make_da_results_labels_chunk = tar_map(
  values = da_recipe_labels_chunk,
  names = label,
  tar_target(da_results_00,
             get_da_results(milo_obj=milo_obj,
                            model_matrix=mm, design_df=design_df, model_contrasts_str=contrast)),
  tar_target(da_results_01,
             add_polarity_to_da_results(da_results_00, spatial_fdr_cutoff=0.1) %>%
             annotate_nhoods(milo_obj, col="labels")),
  tar_target(nhm,
             make_nhm(milo_obj)),
  tar_target(da_results_02,
             annotate_nhood_counts(da_results_01, nhm) # %>%
             # combine_splits(., nhm)
            ),
  tar_target(da_results_02_idx,
             milo_index_obj %>%
             select(-any_of('labels')) %>%
             left_join(da_results_02, ., by='Nhood'),
             packages=c("tidyverse", "ggplot2", "patchwork", "Seurat")),
  tar_target(nonzero_mean_df,
             make_logFC_df(da_results_02) %>%
             make_logFC_cells_df(nhm, .) %>%
             make_nzdf(.)),
  tar_target(nhood_summary,
             make_nhood_summary(da_results_02, nhm) %>% check_split_validity()),
  tar_target(barcode_df,
             make_barcode_polarity_df(nonzero_mean_df, nhood_summary)),
  tar_target(pldf,
             make_pldf(barcode_df, milo_obj))

)



make_split_objs = list(
    make_split_objs_class,
    make_split_objs_lvl1,
    make_split_objs_labels_chunk
)

make_da_results = list(
    make_da_results_class,
    make_da_results_lvl1,
    make_da_results_labels_chunk
)


stage_01 = list(
  tar_target(exp_labelled_other_path, paste0(Sys.getenv("PROJECT_DIR"), '/00_cellbender/_targets/objects/obj_cb_other_01'), format = "file"),
  tar_target(exp_labelled_neuron_path, paste0(Sys.getenv("PROJECT_DIR"), '/00_cellbender/_targets/objects/obj_cb_neuron_00'), format = "file"),
  tar_target(exp_labelled_other, qs::qread(exp_labelled_other_path)),
  tar_target(exp_labelled_neuron, qs::qread(exp_labelled_neuron_path)),
  make_split_objs,
  make_da_results
  )


combination_recipe = qs::qread('combination_recipe.qs')
stage_02 = tar_eval(
  values = combination_recipe,
  tar_combine(output_name,
              make_da_results %>% tar_select_targets(any_of(da_names)))
#   tar_combine(output_name_idx,
#               make_da_results %>% tar_select_targets(any_of(da_names_idx))) ### no idea why this doesn't work as a list of targets!
)
combination_recipe = qs::qread('combination_recipe.qs')
stage_02b = tar_eval(
  values = combination_recipe,
  tar_combine(output_name_idx,
              make_da_results %>% tar_select_targets(any_of(da_names_idx)))
)
combination_recipe = qs::qread('combination_recipe_pldf.qs')
stage_02c = tar_eval(
  values = combination_recipe,
  tar_combine(output_name,
              make_da_results %>% tar_select_targets(any_of(da_names)))
)

stage_02 = list(stage_02, stage_02b, stage_02c)


restored_recipe = qs::qread("restored_recipe.qs")
stage_03 = tar_map(
    values = restored_recipe,
    names = restored_suffix,
    tar_target(restored_df,
               make_nh_restored_tibble(output_obj.fgf1,
                                       output_obj.BL6,
                                       comparison=restored_suffix)),
    tar_target(restored_summary, make_nh_restored_summary_tibble(restored_df, comparison=restored_suffix))
)
restored_df_names = restored_recipe %>% pull(restored_suffix) %>% paste0("restored_df_", .)
restored_summary_names = restored_recipe %>% pull(restored_suffix) %>% paste0("restored_summary_", .)
stage_03 = list(
  stage_03,
  tar_combine(all_restored_df,
              stage_03 %>%tar_select_targets(any_of(restored_df_names))),
  tar_combine(all_restored_summary,
              stage_03 %>%tar_select_targets(any_of(restored_summary_names)))
)

annotate_nh_grouping_recipe = qs::qread("annotate_nh_grouping_recipe.qs")
stage_04 = tar_map(
  values = annotate_nh_grouping_recipe,
  names = da_results_nhg_output_suffix,
  tar_target(da_results_nhg, 
             annotate_summary_groupings(da_results_obj, restored_df_obj)),
  tar_target(da_results_nhg_idx,
             milo_index_tibble %>%
             select(-any_of('labels')) %>%
             left_join(da_results_nhg, ., by='Nhood'),
             packages=c("tidyverse", "ggplot2", "patchwork", "Seurat")),
  tar_target(nhgc,
             nhg2cell(nhm_obj, da_results_nhg))
)


deg_restored_recipe = qs::qread("deg_restored_recipe.qs")
stage_05 = tar_map(
  values = deg_restored_recipe,
  names = deg_output_suffix,
  tar_target(deg,
             get_seurat_nhg_markers_nhBulkDESeq(seurat_obj, nhm_obj, da_results_nhg_obj,
                                                nhood_grouping, group_a, group_b=group_b, tag=deg_output_suffix),
             packages=c("tidyverse", "Seurat", "Matrix")),
  tar_target(deg_ensmus,
             add_gsea_cols_to_seurat_marker_results(deg)),
  tar_target(gsea,
             run_all_fgsea(deg_ensmus, padj_cutoff=0.10, tag=deg_output_suffix)),
  tar_target(gost,
             make_gost(deg_ensmus)),
  tar_target(gost_result,
             gost %>% `[[`('result') %>% mutate(tag = deg_output_suffix))
)
summary_plot_recipe = qs::qread("summary_plot_recipe.qs")
stage_05 = list(
  stage_05,
  tar_map(
    values = summary_plot_recipe,
    names = deg_output_suffix,
    # tar_target(deg_summary_plot_gSCT,
    #            make_summary_deg_plot(global_seurat_obj, nhgc_obj, deg_ensmus_obj, nhood_grouping, name),
    #            packages=c("tidyverse", "ggplot2", "patchwork", "Seurat")),
    # tar_target(deg_summary_file_gSCT,
    #            save_summary_plot(deg_summary_plot_gSCT, c("gSCT_", deg_output_suffix)),
    #            packages=c("tidyverse", "ggplot2", "patchwork", "Seurat"),
    #            format="file"),
# 
# 
# 
    tar_target(deg_summary_plot_gSCT,
               make_summary_deg_plot(global_seurat_obj, nhgc_obj, deg_ensmus_obj, nhood_grouping, name) %>%
               save_summary_plot(., c("gSCT_", deg_output_suffix)),
               packages=c("tidyverse", "ggplot2", "patchwork", "Seurat"),
               priority=0.0001),
    tar_target(deg_summary_plot_reSCT,
               make_summary_deg_plot(seurat_obj, nhgc_obj, deg_ensmus_obj, nhood_grouping, name) %>%
               save_summary_plot(c("reSCT_", deg_output_suffix)),
               packages=c("tidyverse", "ggplot2", "patchwork", "Seurat"),
               priority=0.0001)
    # tar_target(deg_summary_plot_gSCT,
    #            make_and_save_summary_plot(global_seurat_obj, nhgc_obj, deg_ensmus_obj, nhood_grouping, name, 'gSCT_', deg_output_suffix),
    #            packages=c("tidyverse", "ggplot2", "patchwork", "Seurat"),
    #            format="file")
  )
)



combination_recipe_degs = qs::qread("combined_recipe_degs.qs")
stage_06 = tar_eval(
  values = combination_recipe_degs,
  tar_combine(output_name,
              stage_05 %>% tar_select_targets(any_of(results_name)))
)
stage_07 = list(
  tar_target(combined_deg_seurat_formatted,
             combined_deg_seurat %>%
             dplyr::filter(p_val_adj < 0.05) %>%
             rehydrate_deg_tag),
  tar_target(combined_gsea_seurat_formatted,
             rehydrate_deg_tag(combined_gsea_seurat)),
  tar_target(combined_gost_seurat_formatted,
             rehydrate_deg_tag(combined_gost_seurat))
)

stage_xx_recipe = qs::qread('stage_xx_recipe.qs')
stage_xx = tar_eval(
  values = stage_xx_recipe,
  tar_combine(output_name,
              stage_04 %>% tar_select_targets(any_of(da_results_nhg_names)))
)



run_list = list(
  stage_01,   
  stage_02,
  stage_03,
  stage_04,
  stage_05,
  stage_06,
  stage_07,
  stage_xx,
  list()

)

run_list
