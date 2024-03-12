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
  packages = c("tidyverse", "Seurat"), # packages that your targets need to run
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

PROJECT_DIR = Sys.getenv('PROJECT_DIR')

source("../00_preprocessing/splitwrapper.R")
source(paste0("../00_preprocessing/preprocessing.R"))
source(paste0("../01_milo/milo.R"))
source(paste0("../01_milo/milo_plotting.R"))
source(paste0("../00_cellbender/cellbender_preprocess.R"))
source("../00_cellbender/process_seurat.R")
source('../01_milo_cellbender/milo_cellbender.R')
source("../code/resolve2xe/LoadResolveBaysor.R")
source("../01_resolve/resolve.R")
source("resolve_cellbender.R")


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

ingest_tibble = qs::qread('00_ingest_tibble.qs')
stage_01a = tar_map(
  values = ingest_tibble,
  names = sample_name,
  tar_target(sample_path,
             paste0(PROJECT_DIR, '/data/resolve/xe/32810-1377-slide3_', sample_name, '-1/'), 
             format='file'),
  tar_target(obj, 
             sample_path %>% 
             LoadResolveBaysor %>%
             AddMetaData(sample_name, col.name='sample_name') %>%
             AddMetaData(treatment, col.name='treatment') %>%
             AddMetaData(strain, col.name='strain') %>%
             AddMetaData(strain, col.name='time'),
             priority=1)
)





stage_01b = list(
      tar_target(
        nhgc_d5_path, 
        '../01_milo_cellbender/_targets_MAST/objects/combined_nhgc_Day5.obob5v5__v__Day5.obobBL6d5',
        format = "file"),
      tar_target(
        neuron_path, 
        '../00_cellbender/_targets/objects/obj_cb_neuron_00',
        format = "file"),
      tar_target(
        other_path, 
        '../00_cellbender/_targets/objects/obj_cb_other_01',
        format = "file"),
      tar_target(
        tany_path, 
        '../00_cellbender/_targets/objects/obj_cb_tany',
        format = "file"),
      tar_target(obj_d5_00,
                 make_obj_d5(nhgc_d5_path, neuron_path, other_path, tany_path),
                 packages=c("tidyverse", "Seurat")),
      tar_target(obj_d5_01,
                 obj_d5_00 %>%
                 set_labels_to_lvl2 %>%
                 prep_obj_for_milo_cb_v01 %>%
                 set_batch_to_lane %>% # do not set batch to lane for cluster splits, will error on design or model matrix. reset later
                 prep_obj_for_milo_cb_v01(set_orig.batch = FALSE) %>%
                 reconsitute_rna_seurat %>%
                 add_cell_class_and_polar_label %>%
                 process_seurat(method = "integrate", batch ="orig.batch", dims = 30, res = 0.5, k.anchor=25, k.weight=100),
                 packages=c("tidyverse", "Seurat")),
      tar_target(obj_d5_other,
                 obj_d5_00 %>% 
                 subset(subset = class == 'other') %>%
                 set_labels_to_lvl2 %>%
                 prep_obj_for_milo_cb_v01 %>%
                 set_batch_to_lane %>% # do not set batch to lane for cluster splits, will error on design or model matrix. reset later
                 prep_obj_for_milo_cb_v01(set_orig.batch = FALSE) %>%
                 reconsitute_rna_seurat %>%
                 add_cell_class_and_polar_label %>%
                 process_seurat(method = "integrate", batch ="batch", dims = 30, res = 0.8, k.anchor=25),
                 packages=c("tidyverse", "Seurat")),
      tar_target(obj_d5_neuron,
                 obj_d5_00 %>% 
                 subset(subset = class == 'neuron') %>%
                 set_labels_to_lvl2 %>%
                 prep_obj_for_milo_cb_v01 %>%
                 set_batch_to_lane %>% # do not set batch to lane for cluster splits, will error on design or model matrix. reset later
                 prep_obj_for_milo_cb_v01(set_orig.batch = FALSE) %>%
                 reconsitute_rna_seurat %>%
                 add_cell_class_and_polar_label %>%
                 process_seurat(method = "integrate", batch ="orig.batch", dims = 30, res = 0.8, k.anchor=25),
                 packages=c("tidyverse", "Seurat"))

)

stage_01c = list(
    tar_target(obj_merged, 
               merge(obj_A1, list(obj_A2, obj_B1, obj_B2, obj_C1, obj_C2, obj_D1, obj_D2)) %>% JoinLayers
               ),
    tar_target(xe_class_features,
                 c("Agrp", "Aqp4", "Bmp4", "Cfap299", "Cntn4", "Cntn5",
                   "Grm8", "Hs3st4", "Kcnip4", "Pdgfra", "Pdzrn3", "Plp1", 
                   "Pomc", "Rax", "Rbfox1", "Slc1a2", "Slc7a11", "Tenm2",
                   "Trpm3", "Zfp804b")),
    tar_target(resolve_major_classes_markers_path,
               'resolve_major_classes_markers_v2.txt',
               format='file'),
      tar_target(xenium_genes_all,
            get_xe_genes_all_v02(obj_merged)
            ),
      tar_target(xenium_genes,
             get_xe_genes_v02(obj_merged %>% set_default_assay('Xenium'),
                              obj_d5_01 %>% set_default_assay('RNA'))
            ),
    tar_target(xe_obj_00,
                 obj_merged %>% 
                 filter_down_cells_v02 %>% 
                 process_xenium(xenium_genes=xenium_genes) %>%
                 annotate_by_level_v02(counts_assay='Xenium', 
                           graph_name=NULL, 
                           classification_data_path=resolve_major_classes_markers_path, 
                           annotation_col='class'),
                 packages=c("tidyverse", "Seurat", "CellAnnotatoR")),
     tar_target(high_conf_cells_class,
                get_high_conf_cells(xe_obj_00, cutoff_score=0.9)),
     tar_target(xe_obj_01,
                xe_obj_00 %>%
                    subset(cells = high_conf_cells_class) %>%
                    AddMetaData(xe_obj_00 %>% `[[` %>% mutate(cell_class = class)) %>%
                    reclass_by_gene_hilo_v02('Plp1', 20, 20, 'neuron', 'other') %>%
                    reclass_by_gene_hilo_v02('Agrp', 25, 25, 'other', 'neuron') %>%
                    reclass_by_gene_hilo_v02('Pomc', 25, 25, 'other', 'neuron') %>%
                    process_xenium(xenium_genes=xenium_genes)),
      tar_target(obj_d5_other_01,
                 obj_d5_other %>% 
                 reconsitute_rna_seurat %>%
                 process_seurat(method = "integrate", features=xenium_genes, batch ="batch", dims = 30, res = 0.8, k.anchor=25),
                 packages=c("tidyverse", "Seurat")),
      tar_target(obj_d5_neuron_01,
                 obj_d5_neuron %>%
                 reconsitute_rna_seurat %>%
                 process_seurat(method = "integrate", features=xenium_genes, batch ="orig.batch", dims = 30, res = 0.8, k.anchor=25),
                 packages=c("tidyverse", "Seurat"))
)


stage_04 = list(
#     tar_target(xe_obj_cca_td_class,
#                transfer_data_cca_00_v02(obj_merged_sct, obj_fgf1, 'cell_class') %>%
#               #  reclass_by_gene('Agrp', 1, 'neuron') %>%
#               #  reclass_by_gene('Pomc', 1, 'neuron')
#               reclass_by_gene_hilo('Agrp', 1, 10, 'other', 'neuron') %>%
#               reclass_by_gene_hilo('Pomc', 1, 10, 'other', 'neuron')
#             ),
    tar_target(xe_obj_cca_td_neuron,
               xe_obj_01 %>% split_cell_class('neuron') %>% process_xenium(xenium_genes=xenium_genes)),
    tar_target(xe_obj_cca_td_neuron_2s,
            xe_obj_cca_td_neuron %>%
            drop_blah_class %>% 
            process_xenium(xenium_genes=xenium_genes) %>%
            transfer_data_cca_00_v02(obj_d5_neuron_01, 'polar_label') %>%
            labels_from_polar_label %>%
#             relabel_by_gene_hilo_v02('Agrp', 1, 1000, 'Agrp') %>%
            drop_blah_labels %>% 
            process_xenium(xenium_genes=xenium_genes)
        ),
    tar_target(xe_obj_cca_td_neuron_labels,
               transfer_data_cca_00_v02(xe_obj_cca_td_neuron, obj_d5_neuron_01, "labels") %>%
               relabel_by_gene_hilo_v02('Agrp', 1, 10, 'Agrp')
               ),
    tar_target(xe_obj_cca_td_other,
               xe_obj_01 %>% split_cell_class('other') %>% process_xenium(xenium_genes=xenium_genes)
               ),
    tar_target(xe_obj_cca_td_other_2s,
            xe_obj_cca_td_other %>%
            drop_blah_class %>% 
            process_xenium(xenium_genes=xenium_genes) %>%
            transfer_data_cca_00_v02(obj_d5_other_01, 'polar_label') %>%
            labels_from_polar_label
            ),
    tar_target(xe_obj_cca_td_other_labels,
               transfer_data_cca_00_v02(xe_obj_cca_td_other, obj_d5_other_01, 'labels')
               )
)


stage_05_neuron = list(
  tar_target(neuron_labels,
             obj_d5_neuron_01 %>% `[[` %>% distinct(labels) %>% pull(labels)),
  tar_target(xe_obj_cca_td_neuron_labels_polar_label_metadata,
             transfer_data_cca_00_polar_label_v02(xe_obj_cca_td_neuron_labels, obj_d5_neuron_01, xenium_genes, neuron_labels) %>%
             `[[`,
             pattern = map(neuron_labels)
             ),
  tar_target(xe_obj_cca_td_neuron_labels_polar_label_markers,
             transfer_data_cca_00_polar_label_v02(xe_obj_cca_td_neuron_labels, obj_d5_neuron_01, xenium_genes, neuron_labels) %>%
             resolve_find_all_markers(idents='polar_label'),
             pattern = map(neuron_labels)
             ),
  tar_target(xe_obj_cca_td_neuron_2s_labels_polar_label_markers,
            transfer_data_cca_00_polar_label_v02(xe_obj_cca_td_neuron_2s, obj_d5_neuron_01, xenium_genes, neuron_labels) %>%
            resolve_find_all_markers(idents='polar_label'),
            pattern = map(neuron_labels)
            )
)
stage_05_other = list(
  tar_target(other_labels,
             obj_d5_other_01 %>% `[[` %>% distinct(labels) %>% pull(labels)),
  tar_target(xe_obj_cca_td_other_labels_polar_label_metadata,
             transfer_data_cca_00_polar_label_v02(xe_obj_cca_td_other_labels, obj_d5_other_01, xenium_genes, other_labels) %>%
             `[[`,
             pattern = map(other_labels)
             ),
  tar_target(xe_obj_cca_td_other_labels_polar_label_markers,
             transfer_data_cca_00_polar_label_v02(xe_obj_cca_td_other_labels, obj_d5_other_01, xenium_genes, other_labels) %>%
             resolve_find_all_markers(idents='polar_label'),
             pattern = map(other_labels)
             ),
  tar_target(xe_obj_cca_td_other_2s_labels_polar_label_markers,
            transfer_data_cca_00_polar_label_v02(xe_obj_cca_td_other_2s, obj_d5_other_01, xenium_genes, other_labels) %>%
            resolve_find_all_markers(idents='polar_label'),
            pattern = map(other_labels)
            )
)
stage_05_both = list(
  tar_combine(neuron_polar_label_metadata,
              stage_05_neuron %>% tar_select_targets(starts_with("xe_obj_cca_td_neuron_labels_polar_label_metadata"))),
  tar_combine(neuron_polar_label_markers,
              stage_05_neuron %>% tar_select_targets(starts_with("xe_obj_cca_td_neuron_labels_polar_label_markers"))),
  tar_combine(neuron_polar_label_markers_2s,
              stage_05_neuron %>% tar_select_targets(starts_with("xe_obj_cca_td_neuron_2s_labels_polar_label_markers"))),
  tar_combine(other_polar_label_metadata,
              stage_05_other %>% tar_select_targets(starts_with("xe_obj_cca_td_other_labels_polar_label_metadata"))),
  tar_combine(other_polar_label_markers,
              stage_05_other %>% tar_select_targets(starts_with("xe_obj_cca_td_other_labels_polar_label_markers"))),
  tar_combine(other_polar_label_markers_2s,
              stage_05_other %>% tar_select_targets(starts_with("xe_obj_cca_td_other_2s_labels_polar_label_markers"))),              
  tar_target(polar_label_meta,
             bind_rows(neuron_polar_label_metadata,
                       other_polar_label_metadata)
             ),
  tar_target(polar_label_markers,
            bind_rows(neuron_polar_label_markers,
                      other_polar_label_markers)
             ),
  tar_target(polar_label_markers_2s,
            bind_rows(neuron_polar_label_markers_2s,
                      other_polar_label_markers_2s)
             ),              
  tar_target(xe_obj_cca_td_3s,
             merge(xe_obj_cca_td_neuron_labels, xe_obj_cca_td_other_labels) %>%
             AddMetaData(metadata = polar_label_meta[c('polar_label', 'polar_label_prediction.score.max')]) %>%
             process_xenium(xenium_genes=xenium_genes)
             ),
  tar_target(xe_obj_cca_td_2s,
            merge(xe_obj_cca_td_neuron_2s, xe_obj_cca_td_other_2s) %>%
            process_xenium(xenium_genes=xenium_genes)
            )
)


stage_04_uni = list(
#     tar_target(xe_obj_cca_td_class,
#                transfer_data_cca_00_unimodal_v02(obj_merged_sct, obj_fgf1, 'cell_class') %>%
#               #  reclass_by_gene('Agrp', 1, 'neuron') %>%
#               #  reclass_by_gene('Pomc', 1, 'neuron')
#               reclass_by_gene_hilo('Agrp', 1, 10, 'other', 'neuron') %>%
#               reclass_by_gene_hilo('Pomc', 1, 10, 'other', 'neuron')
#             ),
    tar_target(xe_obj_cca_td_neuron_uni,
               xe_obj_01 %>% split_cell_class('neuron') %>% process_xenium(xenium_genes=xenium_genes)),
    tar_target(xe_obj_cca_td_neuron_2s_uni,
            xe_obj_cca_td_neuron_uni %>%
            drop_blah_class %>% 
            process_xenium(xenium_genes=xenium_genes) %>%
            transfer_data_cca_00_unimodal_v02(obj_d5_neuron_01, 'polar_label') %>%
            labels_from_polar_label %>%
#             relabel_by_gene_hilo_v02('Agrp', 1, 1000, 'Agrp') %>%
            drop_blah_labels %>% 
            process_xenium(xenium_genes=xenium_genes)
        ),
    tar_target(xe_obj_cca_td_neuron_labels_uni,
               transfer_data_cca_00_unimodal_v02(xe_obj_cca_td_neuron_uni, obj_d5_neuron_01, "labels") %>%
               relabel_by_gene_hilo_v02('Agrp', 1, 10, 'Agrp')
               ),
    tar_target(xe_obj_cca_td_other_uni,
               xe_obj_01 %>% split_cell_class('other') %>% process_xenium(xenium_genes=xenium_genes)
               ),
    tar_target(xe_obj_cca_td_other_2s_uni,
            xe_obj_cca_td_other_uni %>%
            drop_blah_class %>% 
            process_xenium(xenium_genes=xenium_genes) %>%
            transfer_data_cca_00_unimodal_v02(obj_d5_other_01, 'polar_label') %>%
            labels_from_polar_label
            ),
    tar_target(xe_obj_cca_td_other_labels_uni,
               transfer_data_cca_00_unimodal_v02(xe_obj_cca_td_other_uni, obj_d5_other_01, 'labels')
               )
)
stage_05_neuron_uni = list(
  tar_target(neuron_labels_uni,
             obj_d5_neuron_01 %>% `[[` %>% distinct(labels) %>% pull(labels)),
  tar_target(xe_obj_cca_td_neuron_labels_polar_label_metadata_uni,
             transfer_data_cca_00_polar_label_unimodal_v02(xe_obj_cca_td_neuron_labels, obj_d5_neuron_01, xenium_genes, neuron_labels) %>%
             `[[`,
             pattern = map(neuron_labels_uni)
             ),
  tar_target(xe_obj_cca_td_neuron_labels_polar_label_markers_uni,
             transfer_data_cca_00_polar_label_unimodal_v02(xe_obj_cca_td_neuron_labels, obj_d5_neuron_01, xenium_genes, neuron_labels) %>%
             resolve_find_all_markers(idents='polar_label'),
             pattern = map(neuron_labels_uni)
             )
)
stage_05_other_uni = list(
  tar_target(other_labels_uni,
             obj_d5_other_01 %>% `[[` %>% distinct(labels) %>% pull(labels)),
  tar_target(xe_obj_cca_td_other_labels_polar_label_metadata_uni,
             transfer_data_cca_00_polar_label_unimodal_v02(xe_obj_cca_td_other_labels, obj_d5_other_01, xenium_genes, other_labels) %>%
             `[[`,
             pattern = map(other_labels_uni)
             ),
  tar_target(xe_obj_cca_td_other_labels_polar_label_markers_uni,
             transfer_data_cca_00_polar_label_unimodal_v02(xe_obj_cca_td_other_labels, obj_d5_other_01, xenium_genes, other_labels) %>%
             resolve_find_all_markers(idents='polar_label'),
             pattern = map(other_labels_uni)
             )
)
stage_05_both_uni = list(
  tar_combine(neuron_polar_label_metadata_uni,
              stage_05_neuron_uni %>% tar_select_targets(starts_with("xe_obj_cca_td_neuron_labels_polar_label_metadata_uni"))),
  tar_combine(neuron_polar_label_markers_uni,
              stage_05_neuron_uni %>% tar_select_targets(starts_with("xe_obj_cca_td_neuron_labels_polar_label_markers_uni"))),
  tar_combine(other_polar_label_metadata_uni,
              stage_05_other_uni %>% tar_select_targets(starts_with("xe_obj_cca_td_other_labels_polar_label_metadata_uni"))),
  tar_combine(other_polar_label_markers_uni,
              stage_05_other_uni %>% tar_select_targets(starts_with("xe_obj_cca_td_other_labels_polar_label_markers_uni"))),            
  tar_target(polar_label_meta_uni,
             bind_rows(neuron_polar_label_metadata_uni,
                       other_polar_label_metadata_uni)
             ),
  tar_target(polar_label_markers_uni,
            bind_rows(neuron_polar_label_markers_uni,
                      other_polar_label_markers_uni)
             ),           
  tar_target(xe_obj_cca_td_3s_uni,
             merge(xe_obj_cca_td_neuron_labels_uni, xe_obj_cca_td_other_labels_uni) %>%
             AddMetaData(metadata = polar_label_meta[c('polar_label', 'polar_label_prediction.score.max')]) %>%
             process_xenium(xenium_genes=xenium_genes)
             ),
  tar_target(xe_obj_cca_td_2s_uni,
            merge(xe_obj_cca_td_neuron_2s_uni, xe_obj_cca_td_other_2s_uni) %>%
             AddMetaData(metadata = polar_label_meta[c('polar_label', 'polar_label_prediction.score.max')]) %>%
             process_xenium(xenium_genes=xenium_genes)
            )
)
stage_05_uni = list(stage_05_neuron_uni, stage_05_other_uni, stage_05_both_uni)



stage_01 = list(stage_01a, stage_01b, stage_01c)

stage_05 = list(stage_05_neuron,
                stage_05_other, 
                stage_05_both
               )
stage_05_uni = list(stage_05_neuron_uni, stage_05_other_uni, stage_05_both_uni)

run_list = list(
  stage_01,
#   stage_02,
#   stage_03,
  stage_04,
  stage_05,
  stage_04_uni,
  stage_05_uni,
  list()
)

run_list
