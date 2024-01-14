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
  packages = c("Seurat", "tidyverse"), # packages that your targets need to run
  format = "qs", # default storage format,
  error = "null",
  retrieval = "worker",
  storage = "worker",
    workspace_on_error = TRUE # Save a workspace file for a target that errors out.
  # Set other options as needed.
)

PROJECT_DIR = Sys.getenv('PROJECT_DIR')
source(paste0(PROJECT_DIR, "/00_preprocessing/preprocessing.R"))
source(paste0(PROJECT_DIR, "/code/resolve2xe/LoadResolveBaysor.R"))
source(paste0("resolve.R"))

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore",
        clustermq.ssh.timeout=36000,
        clustermq.worker.timeout=36000,
        clustermq.error.timeout=36000
        )

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
# tar_source()
# source("other_functions.R") # Source other scripts as needed. # nolint

ingest_tibble = qs::qread('00_ingest_tibble.qs')
stage_01 = tar_map(
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
             priority=1))

stage_02 = list(
  tar_target(obj_merged, 
             merge(obj_A1, list(obj_A2, obj_B1, obj_B2, obj_C1, obj_C2, obj_D1, obj_D2))
             ),
  tar_target(obob_d5_path,
             '../genebasisr_obob5v5/_targets/objects/exp_all_n1xo',
             format='file'),
  tar_target(obj_fgf1,
             qs::qread(obob_d5_path)),
  tar_target(xenium_genes_all,
            get_xe_genes_all(obj_merged)
            ),
  tar_target(xenium_genes,
             get_xe_genes(obj_merged, obj_fgf1)
            ),
  tar_target(obj_fgf1_sct_xeg,
             obj_fgf1 %>% 
             Seurat::RunUMAP(assay='SCT', slot='counts', features = xenium_genes, return.model = TRUE)
             ),
  tar_target(obj_merged_sct,
             obj_merged %>%
             filter_down_cells %>%
             sc_transform_resolve %>%
             Seurat::RunUMAP(assay='SCT', slot='counts', features = xenium_genes, return.model = TRUE)
            )
  )

  stage_03 = list(
      tar_target(obob_d5_neurons_path,
                '../genebasisr_obob5v5/_targets/objects/exp_neuron_obob5v5',
                format='file'
                ),
      tar_target(obob_d5_neurons_obj,
                qs::qread(obob_d5_neurons_path)
                ),
      tar_target(obj_fgf1_neurons_sct_xeg,
                obob_d5_neurons_obj %>% 
                Seurat::RunUMAP(assay='SCT', slot='counts', features = xenium_genes, return.model = TRUE)
                ),
      tar_target(obob_d5_other_path,
                '../genebasisr_obob5v5/_targets/objects/exp_other_obob5v5',
                format='file'
                ),
      tar_target(obob_d5_other_obj,
                qs::qread(obob_d5_other_path)
                ),
      tar_target(obj_fgf1_other_sct_xeg,
                obob_d5_other_obj %>% 
                Seurat::RunUMAP(assay='SCT', slot='counts', features = xenium_genes, return.model = TRUE)
                )
)

stage_04 = list(
    tar_target(xe_obj_cca_td_class,
               transfer_data_cca_00(obj_merged_sct, obj_fgf1, 'cell_class') %>%
              #  reclass_by_gene('Agrp', 1, 'neuron') %>%
              #  reclass_by_gene('Pomc', 1, 'neuron')
              reclass_by_gene_hilo('Agrp', 1, 10, 'other', 'neuron') %>%
              reclass_by_gene_hilo('Pomc', 1, 10, 'other', 'neuron')
            ),
    tar_target(xe_obj_cca_td_neuron,
               xe_obj_cca_td_class %>% split_cell_class('neuron') %>% sc_transform_resolve %>%
             Seurat::RunUMAP(assay='Xenium',
                           slot='counts',
                           features = .@assays$SCT@meta.features %>% 
                                      rownames_to_column(var = 'genes') %>% 
                                      filter(genes != 'Lmx1a') %>% 
                                      pull(genes),
                           return.model = TRUE)
               ),
    tar_target(xe_obj_cca_td_neuron_2s,
            xe_obj_cca_td_neuron %>%
            drop_blah_class %>% 
            sc_transform_resolve %>%
            transfer_data_cca_00(obj_fgf1_neurons_sct_xeg, 'polar_label') %>%
            labels_from_polar_label %>%
            relabel_by_gene_hilo('Agrp', 1, 1000, 'Agrp') %>%
            drop_blah_labels %>% 
            sc_transform_resolve %>%
             Seurat::RunUMAP(assay='Xenium',
                           slot='counts',
                           features = .@assays$SCT@meta.features %>% 
                                      rownames_to_column(var = 'genes') %>% 
                                      filter(genes != 'Lmx1a') %>% 
                                      pull(genes),
                           return.model = TRUE)
        ),
    tar_target(xe_obj_cca_td_neuron_labels,
               transfer_data_cca_00(xe_obj_cca_td_neuron, obj_fgf1_neurons_sct_xeg, "labels") %>%
               relabel_by_gene_hilo('Agrp', 1, 10, 'Agrp')
               ),
    tar_target(xe_obj_cca_td_other,
               xe_obj_cca_td_class %>% split_cell_class('other') %>% sc_transform_resolve %>%
             Seurat::RunUMAP(assay='Xenium',
                           slot='counts',
                           features = .@assays$SCT@meta.features %>% 
                                      rownames_to_column(var = 'genes') %>% 
                                      filter(genes != 'Lmx1a') %>% 
                                      pull(genes),
                           return.model = TRUE)
               ),
    tar_target(xe_obj_cca_td_other_2s,
            xe_obj_cca_td_other %>%
            drop_blah_class %>% 
            sc_transform_resolve %>%
            transfer_data_cca_00(obj_fgf1_other_sct_xeg, 'polar_label') %>%
            labels_from_polar_label %>%
             Seurat::RunUMAP(assay='Xenium',
                           slot='counts',
                           features = .@assays$SCT@meta.features %>% 
                                      rownames_to_column(var = 'genes') %>% 
                                      filter(genes != 'Lmx1a') %>% 
                                      pull(genes),
                           return.model = TRUE)
        ),
    tar_target(xe_obj_cca_td_other_labels,
               transfer_data_cca_00(xe_obj_cca_td_other, obj_fgf1_other_sct_xeg, 'labels')
               )
)

stage_05_neuron = list(
  tar_target(neuron_labels,
             obj_fgf1_neurons_sct_xeg %>% `[[` %>% distinct(labels) %>% pull(labels)),
  tar_target(xe_obj_cca_td_neuron_labels_polar_label_metadata,
             transfer_data_cca_00_polar_label(xe_obj_cca_td_neuron_labels, obj_fgf1_neurons_sct_xeg, neuron_labels) %>%
             `[[`,
             pattern = map(neuron_labels)
             ),
  tar_target(xe_obj_cca_td_neuron_labels_polar_label_markers,
             transfer_data_cca_00_polar_label(xe_obj_cca_td_neuron_labels, obj_fgf1_neurons_sct_xeg, neuron_labels) %>%
             resolve_find_all_markers(idents='polar_label'),
             pattern = map(neuron_labels)
             ),
  tar_target(xe_obj_cca_td_neuron_2s_labels_polar_label_markers,
            transfer_data_cca_00_polar_label(xe_obj_cca_td_neuron_2s, obj_fgf1_neurons_sct_xeg, neuron_labels) %>%
            resolve_find_all_markers(idents='polar_label'),
            pattern = map(neuron_labels)
            )
)
stage_05_other = list(
  tar_target(other_labels,
             obj_fgf1_other_sct_xeg %>% `[[` %>% distinct(labels) %>% pull(labels)),
  tar_target(xe_obj_cca_td_other_labels_polar_label_metadata,
             transfer_data_cca_00_polar_label(xe_obj_cca_td_other_labels, obj_fgf1_other_sct_xeg, other_labels) %>%
             `[[`,
             pattern = map(other_labels)
             ),
  tar_target(xe_obj_cca_td_other_labels_polar_label_markers,
             transfer_data_cca_00_polar_label(xe_obj_cca_td_other_labels, obj_fgf1_other_sct_xeg, other_labels) %>%
             resolve_find_all_markers(idents='polar_label'),
             pattern = map(other_labels)
             ),
  tar_target(xe_obj_cca_td_other_2s_labels_polar_label_markers,
             transfer_data_cca_00_polar_label(xe_obj_cca_td_other_2s, obj_fgf1_other_sct_xeg, other_labels) %>%
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
             sc_transform_resolve %>%
             Seurat::RunUMAP(assay='Xenium', slot='counts', features = xenium_genes_all, return.model = TRUE)
             ),
  tar_target(xe_obj_cca_td_2s,
            merge(xe_obj_cca_td_neuron_2s, xe_obj_cca_td_other_2s) %>%
            sc_transform_resolve %>%
             Seurat::RunUMAP(assay='Xenium',
                           slot='counts',
                           features = .@assays$SCT@meta.features %>% 
                                      rownames_to_column(var = 'genes') %>% 
                                      filter(genes != 'Lmx1a') %>% 
                                      pull(genes),
                           return.model = TRUE)
            )
)
stage_05 = list(stage_05_neuron, stage_05_other, stage_05_both)


####
stage_04_uni = list(
    tar_target(xe_obj_cca_td_class_uni,
               transfer_data_cca_00_unimodal(obj_merged_sct, obj_fgf1, 'cell_class') %>%
              #  reclass_by_gene('Agrp', 1, 'neuron') %>%
              #  reclass_by_gene('Pomc', 1, 'neuron')
              reclass_by_gene_hilo('Agrp', 1, 10, 'other', 'neuron') %>%
              reclass_by_gene_hilo('Pomc', 1, 10, 'other', 'neuron')
            ),
     tar_target(xe_obj_cca_td_neuron_2s_uni,
            xe_obj_cca_td_neuron_uni %>%
            drop_blah_class %>% 
            sc_transform_resolve %>%
            transfer_data_cca_00_unimodal(obj_fgf1_neurons_sct_xeg, 'polar_label') %>%
            labels_from_polar_label %>%
            relabel_by_gene_hilo('Agrp', 1, 1000, 'Agrp') %>%
            drop_blah_labels %>% 
            sc_transform_resolve %>%
             Seurat::RunUMAP(assay='Xenium',
                           slot='counts',
                           features = .@assays$SCT@meta.features %>% 
                                      rownames_to_column(var = 'genes') %>% 
                                      filter(genes != 'Lmx1a') %>% 
                                      pull(genes),
                           return.model = TRUE)
        ),
    tar_target(xe_obj_cca_td_neuron_uni,
               xe_obj_cca_td_class_uni %>% split_cell_class('neuron') %>% sc_transform_resolve
               ),
    tar_target(xe_obj_cca_td_neuron_labels_uni,
               transfer_data_cca_00_unimodal(xe_obj_cca_td_neuron_uni, obj_fgf1_neurons_sct_xeg, "labels") %>%
               relabel_by_gene_hilo('Agrp', 1, 10, 'Agrp')
               ),
    tar_target(xe_obj_cca_td_other_uni,
               xe_obj_cca_td_class_uni %>% split_cell_class('other') %>% sc_transform_resolve
               ),
    tar_target(xe_obj_cca_td_other_2s_uni,
            xe_obj_cca_td_other_uni %>%
            drop_blah_class %>% 
            sc_transform_resolve %>%
            transfer_data_cca_00_unimodal(obj_fgf1_other_sct_xeg, 'polar_label') %>%
            labels_from_polar_label %>%
             Seurat::RunUMAP(assay='Xenium',
                           slot='counts',
                           features = .@assays$SCT@meta.features %>% 
                                      rownames_to_column(var = 'genes') %>% 
                                      filter(genes != 'Lmx1a') %>% 
                                      pull(genes),
                           return.model = TRUE)
        ),
    tar_target(xe_obj_cca_td_other_labels_uni,
               transfer_data_cca_00_unimodal(xe_obj_cca_td_other_uni, obj_fgf1_other_sct_xeg, 'labels')
               )
)

stage_05_neuron_uni = list(
  tar_target(neuron_labels_uni,
             obj_fgf1_neurons_sct_xeg %>% `[[` %>% distinct(labels) %>% pull(labels)),
  tar_target(xe_obj_cca_td_neuron_labels_polar_label_metadata_uni,
             transfer_data_cca_00_polar_label_unimodal(xe_obj_cca_td_neuron_labels_uni, obj_fgf1_neurons_sct_xeg, neuron_labels_uni) %>%
             `[[`,
             pattern = map(neuron_labels_uni)
             ),
  tar_target(xe_obj_cca_td_neuron_labels_polar_label_markers_uni,
             transfer_data_cca_00_polar_label_unimodal(xe_obj_cca_td_neuron_labels_uni, obj_fgf1_neurons_sct_xeg, neuron_labels_uni) %>%
             resolve_find_all_markers(idents='polar_label'),
             pattern = map(neuron_labels_uni)
             )
)
stage_05_other_uni = list(
  tar_target(other_labels_uni,
             obj_fgf1_other_sct_xeg %>% `[[` %>% distinct(labels) %>% pull(labels)),
  tar_target(xe_obj_cca_td_other_labels_polar_label_metadata_uni,
             transfer_data_cca_00_polar_label_unimodal(xe_obj_cca_td_other_labels_uni, obj_fgf1_other_sct_xeg, other_labels_uni) %>%
             `[[`,
             pattern = map(other_labels_uni)
             ),
  tar_target(xe_obj_cca_td_other_labels_polar_label_markers_uni,
             transfer_data_cca_00_polar_label_unimodal(xe_obj_cca_td_other_labels_uni, obj_fgf1_other_sct_xeg, other_labels_uni) %>%
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
             AddMetaData(metadata = polar_label_meta_uni[c('polar_label', 'polar_label_prediction.score.max')]) %>%
             sc_transform_resolve %>%
             Seurat::RunUMAP(assay='Xenium',
                           slot='counts',
                           features = .@assays$SCT@meta.features %>% 
                                      rownames_to_column(var = 'genes') %>% 
                                      filter(genes != 'Lmx1a') %>% 
                                      pull(genes),
                           return.model = TRUE)
             ),
  tar_target(xe_obj_cca_td_2s_uni,
            merge(xe_obj_cca_td_neuron_2s_uni, xe_obj_cca_td_other_2s_uni) %>%
            sc_transform_resolve %>%
             Seurat::RunUMAP(assay='Xenium',
                           slot='counts',
                           features = .@assays$SCT@meta.features %>% 
                                      rownames_to_column(var = 'genes') %>% 
                                      filter(genes != 'Lmx1a') %>% 
                                      pull(genes),
                           return.model = TRUE)
            )
)
stage_05_uni = list(stage_05_neuron_uni, stage_05_other_uni, stage_05_both_uni)





run_list = list(
  stage_01,
  stage_02,
  stage_03,
  stage_04,
  stage_05,
  stage_04_uni,
  stage_05_uni
)

run_list
