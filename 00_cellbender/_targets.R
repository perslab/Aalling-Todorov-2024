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
  packages = c("tidyverse", "Seurat"), # packages that your targets need to run
  format = "qs", # default storage format,
  error = "null",
  retrieval = "worker",
  storage = "worker",
  memory = "transient",
  garbage_collection = TRUE,
  workspace_on_error = TRUE,
  priority=0.5
  # Set other options as needed.
)

source(paste0("../00_preprocessing/preprocessing.R"))
source(paste0("../01_milo/milo.R"))
source(paste0("cellbender_preprocess.R"))
source("process_seurat.R")
source("map_ref.R")
source("../01_milo_cellbender/milo_cellbender.R")


# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore",
        clustermq.ssh.timeout=36000,
        clustermq.worker.timeout=36000,
        clustermq.error.timeout=36000,
        clustermq.ssh.log='clustermq_sshlog.log'
        )


stage_01 = list(
  tar_target(path_to_cb, 
             paste0(Sys.getenv("PROJECT_DIR"), 'data/cellbender/cellbender_h5'), 
             format = "file",
             packages=c("tidyverse", "Seurat", "scCustomize")),
  tar_target(path_to_sample_meta, 
             paste0(Sys.getenv("PROJECT_DIR"), 'data/meta/metadata.csv'), 
             format = "file",
             packages=c("tidyverse", "Seurat")),
  tar_target(obj_cb_00,
             path_to_cb %>%
             make_mat %>%
             fix_mat_colnames %>%
             make_obj_from_mat),
  tar_target(cb_metadata,
             make_metadata(obj_cb_00, path_to_cb, path_to_sample_meta)),
  tar_target(obj_cb_01,
             AddMetaData(obj_cb_00, cb_metadata) %>%
             filter_down_cells),
  tar_target(obj_cb_02,
             obj_cb_01 %>%
             reconsitute_rna_seurat %>%
             process_seurat(method = "integrate", batch ="Index.10x_SCOP", dims = 30, res = 0.5, k.anchor=25)),
    tar_target(obj_cb_02_sct,
             obj_cb_01 %>%
             reconsitute_rna_seurat %>%
             sc_transform_fgf1 %>% run_sct_chaser)
  )

stage_02 = list(
    tar_target(path_to_class_markers, 
               paste0('marker_major_class.txt'), 
               format = "file"),
    tar_target(obj_cb_class,
               annotate_by_level(obj_cb_02,
                                 counts_assay='RNA', 
                                 graph_name=NULL, 
                                 classification_data_path=path_to_class_markers, 
                                 annotation_col='class'),
               packages=c("tidyverse", "Seurat", "CellAnnotatoR")),
    tar_target(path_to_other_markers, 
               paste0('non_neuron_markers.txt'), 
               format = "file"),

    tar_target(obj_cb_other_00,
               obj_cb_class %>% 
               subset(subset = class == 'other') %>%
               reconsitute_rna_seurat %>%
               process_seurat(method = "integrate", batch ="Index.10x_SCOP", dims = 30, res = 0.8, k.anchor=25) %>%
               annotate_by_level(counts_assay='RNA', 
                                 graph_name=NULL, 
                                 classification_data_path=path_to_other_markers, 
                                 annotation_col='labels_lvl1') %>%
               subset(subset = labels_lvl1 != 'neuron') %>%
               reconsitute_rna_seurat %>%
               process_seurat(method = "integrate", batch ="Index.10x_SCOP", dims = 30, res = 0.8, k.anchor=25),
               packages=c("tidyverse", "Seurat", "CellAnnotatoR")),
    tar_target(label_chunks_other_path, paste0(Sys.getenv("PROJECT_DIR"), '00_cellbender/cluster_surgery/labels_chunk_other.qs'), format = "file"),
    tar_target(label_chunks_other,
               qs::qread(label_chunks_other_path)),
    tar_target(obj_cb_other_01,
               obj_cb_other_00 %>% 
               AddMetaData(label_chunks_other),
               packages=c("tidyverse", "Seurat", "CellAnnotatoR")),
    tar_target(obj_cb_other_02,
               obj_cb_other_01 %>%
               subset(subset = labels_chunk != 'g_drop') %>%
               reconsitute_rna_seurat %>%
               process_seurat(method = "integrate", batch ="Index.10x_SCOP", dims = 30, res = 0.8, k.anchor=25) %>%
               annotate_by_level(counts_assay='RNA', 
                                 graph_name=NULL, 
                                 classification_data_path=path_to_other_markers, 
                                 annotation_col='labels_lvl1'),
               packages=c("tidyverse", "Seurat", "CellAnnotatoR")),
    tar_target(obj_cb_neuron,
               obj_cb_class %>% 
               subset(subset = class == 'neuron') %>%
               reconsitute_rna_seurat %>%
               process_seurat(method = "integrate", batch ="Index.10x_SCOP", dims = 30, res = 0.8, k.anchor=25) %>%
               annotate_by_level(counts_assay='RNA', 
                                 graph_name=NULL, 
                                 classification_data_path=path_to_other_markers, 
                                 annotation_col='labels_lvl1_mg') %>%
               subset(subset = labels_lvl1_mg == 'neuron') %>%
               reconsitute_rna_seurat %>%
               process_seurat(method = "integrate", batch ="Index.10x_SCOP", dims = 30, res = 0.8, k.anchor=25) %>%
               map_ref_merged_clusters(column_name='labels_lvl1', ref="/projects/amj/my_projects/reference_ARCDVC/20231010_full_integration/1_dat/1_seurat_obj/20231027_ARC_VMH_filtered.qs")  %>%
               AddMetaData(.,
                           .[[]] %>% select(labels_lvl1) %>% rename(labels_lvl2 = labels_lvl1)),
               packages=c("tidyverse", "Seurat", "CellAnnotatoR")),
# neuron cluster surgery
    tar_target(label_chunks_neuron_path, paste0(Sys.getenv("PROJECT_DIR"), '00_cellbender/cluster_surgery/labels_chunk_neuron_v01.qs'), format = "file"),
    tar_target(label_chunks_neuron,
               qs::qread(label_chunks_neuron_path)),
   tar_target(obj_cb_neuron_00,
               obj_cb_neuron %>%
               AddMetaData(label_chunks_neuron),
               packages=c("tidyverse", "Seurat", "CellAnnotatoR")
              ),
    tar_target(obj_cb_neuron_01,
               obj_cb_neuron_00 %>%
               AddMetaData(label_chunks_neuron) %>%
               subset(subset = labels_chunk != 'n_drop') %>%
               reconsitute_rna_seurat %>%
               process_seurat(method = "integrate", batch ="Index.10x_SCOP", dims = 30, res = 0.5) %>%
               map_ref_merged_clusters(column_name='labels_lvl1', ref="/projects/amj/my_projects/reference_ARCDVC/20231010_full_integration/1_dat/1_seurat_obj/20231027_ARC_VMH_filtered.qs"),
               packages=c("tidyverse", "Seurat", "CellAnnotatoR")
              ),
    tar_target(obj_cb_tany,
               obj_cb_other_01 %>% 
               subset(subset = labels_lvl1 ==  'Tanycytes') %>%
               reconsitute_rna_seurat %>%
               process_seurat(method = "integrate",
                              nfeats = 2000,
                              batch ="batch",
                              dims = 20, res = 0.5) %>%
               annotate_by_level(counts_assay='RNA', 
                                 graph_name=NULL, 
                                 classification_data_path=path_to_tany_markers, 
                                 annotation_col='labels_lvl2'),
               packages=c("tidyverse", "Seurat", "CellAnnotatoR")),
    tar_target(path_to_tany_markers, 
               paste0('tany_markers.txt'), 
               format = "file"),
    tar_target(obj_cb_other_03,
               obj_cb_tany %>% 
               `[[` %>% 
               select(labels_lvl2) %>%
               AddMetaData(obj_cb_other_02, .) %>% 
               set_empty_labels_lvl2_to_lvl1,
               packages=c("tidyverse", "Seurat"))
)

stage_03 = list(
    tar_target(obj_cb_neuron_sv4,
               obj_cb_neuron %>% seurat_v5_to_v4),
    tar_target(obj_cb_neuron_00_sv4,
               obj_cb_neuron_00 %>% seurat_v5_to_v4),
    tar_target(obj_cb_neuron_01_sv4,
               obj_cb_neuron_01 %>% seurat_v5_to_v4),
    tar_target(obj_cb_other_00_sv4,
               obj_cb_other_00 %>% seurat_v5_to_v4),
    tar_target(obj_cb_other_02_sv4,
               obj_cb_other_02 %>% seurat_v5_to_v4),
    tar_target(obj_cb_tany_sv4,
               obj_cb_tany %>% seurat_v5_to_v4)

)

stage_04 = list(
    tar_target(exp_labelled_other,
               obj_cb_other_01,
               priority=0.1),
        tar_target(exp_labelled_neuron,
               obj_cb_neuron_00,
               priority=0.1)
)

clusters_tibble_lvl1 = qs::qread('../01_milo_cellbender/clusters_tibble_lvl1.qs') #%>% 
#     filter(cluster %in% c('Astrocytes', "Oligodendrocytes", "Microglia"))
stage_04a = tar_map(
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
#                             dims = 30, res = 0.8, k.anchor=25,
#                             k.weight = 40)
             run_until_success(function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                              batch = "batch",
                                              dims = 30, res = 0.8, k.anchor=25,
                                              k.weight = 100),
                               function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                              batch = "batch",
                                              dims = 30, res = 0.8, k.anchor=25,
                                              k.weight = 50),
                               function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                              batch = "batch",
                                              dims = 30, res = 0.8, k.anchor=25,
                                              k.weight = 20),
                               function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                              batch = "batch",
                                              dims = 30, res = 0.8, k.anchor=25,
                                              k.weight = 10),
                               function() subset(eval(.)) %>% Seurat::SCTransform(assay='RNA',
                                               method="glmGamPoi",
                                               vars.to.regress= 'batch',
                                               vst.flavor="v2",
                                               verbose=TRUE),
                               function() subset(eval(.)) %>% Seurat::SCTransform(assay='RNA',
                                               method="glmGamPoi",
                                               vars.to.regress= 'orig.batch',
                                               vst.flavor="v2",
                                               verbose=TRUE)
                              ) %>% #reset batch for downstream applications
             reset_orig.batch %>% 
             prep_obj_for_milo_cb_v01,
             priority=0.0001
             )
)

stage_04 = list(stage_04, stage_04a)
    
    
# xx_stage_03 = list(
#     tar_target(obj_cb_other_00_all,
#                obj_cb_class %>% 
# #                subset(subset = class == 'other') %>%
#.               reconsitute_rna_seurat %>%
#                process_seurat(method = "integrate", batch ="Index.10x_SCOP", dims = 30, res = 0.75, k.anchor=25, k.weight=25) %>%
#                annotate_by_level(counts_assay='RNA', 
#                                  graph_name=NULL, 
#                                  classification_data_path=path_to_other_markers, 
#                                  annotation_col='labels_lvl1'),
#                packages=c("tidyverse", "Seurat", "CellAnnotatoR")),
#     tar_target(obj_cb_neuron_2,
#                obj_cb_neuron %>%
#                map_ref(ref = "/projects/amj/my_projects/reference_ARCDVC/20231010_full_integration/1_dat/1_seurat_obj/20231027_ARC_VMH_filtered.qs"),
#                packages=c("tidyverse", "Seurat", "CellAnnotatoR")),
#     tar_target(obj_cb_neuron_3,
#                obj_cb_neuron %>%
#                transfer_labels_via_map_query(ref_path = "/projects/amj/my_projects/reference_ARCDVC/20231010_full_integration/1_dat/1_seurat_obj/20231027_ARC_VMH_filtered.qs"),
#                packages=c("tidyverse", "Seurat", "CellAnnotatoR")),
# tar_target(obj_cb_neuron_4,
#                obj_cb_neuron %>%
#                map_camp_seurat_cluster(ref_path = "/projects/amj/my_projects/reference_ARCDVC/20231010_full_integration/1_dat/1_seurat_obj/20231027_ARC_VMH_filtered.qs"),
#                packages=c("tidyverse", "Seurat", "CellAnnotatoR")),
# 
# )


run_list = list(
  stage_01,
  stage_02,
  stage_03,
    stage_04
)

run_list
