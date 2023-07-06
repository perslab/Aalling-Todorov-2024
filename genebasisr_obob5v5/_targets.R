# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed. # nolint
library(future)
library(future.callr)
plan(callr)

# Set target options:
tar_option_set(
#   packages = c("tibble"), # packages that your targets need to run
  format = "qs", # default storage format,
  error = "continue",
  retrieval = "worker",
  storage = "worker"
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore")
options(future.globals.maxSize= (200 * 1000)*1024^2) #200GB

source(paste0("../00_preprocessing/preprocessing.R"))
source(paste0("../01_milo/milo.R"))
source(paste0("../genebasisr/genebasisr.R"))
source(paste0("../genebasisr/prep_probe_selection.R"))

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
# tar_source()
# source("other_functions.R") # Source other scripts as needed. # nolint

load_inputs_list = list(
  tar_target(exp_labelled_other_path, paste0(Sys.getenv("PROJECT_DIR"), '/00_preprocessing/_targets/objects/exp_labelled_other'), format = "file"),
  tar_target(exp_labelled_neuron_path, paste0(Sys.getenv("PROJECT_DIR"), '/00_preprocessing/_targets/objects/exp_RElabelled_neuron'), format = "file"),
  tar_target(pldf_obob5v5_path, paste0(Sys.getenv("PROJECT_DIR"), '/01_milo/_targets/objects/combined_pldf_obob5v5'), format = "file"),
    tar_target(pldf_obob5v5, 
             qs::qread(pldf_obob5v5_path)),
  tar_target(exp_labelled_other, 
             qs::qread(exp_labelled_other_path) %>%
             add_metadata_label_to_seurat("other", "cell_class") %>% 
             prep_obj(pldf_obob5v5)),
  tar_target(exp_labelled_neuron,
             qs::qread(exp_labelled_neuron_path) %>%
             add_metadata_label_to_seurat("neuron", "cell_class") %>%
             prep_obj(pldf_obob5v5))
)

transform_inputs_tibble <- tibble(
  input_obj = rlang::syms(c("exp_labelled_other", "exp_labelled_neuron")),
  names = c("other_obob5v5", "neuron_obob5v5")
)
transform_inputs_list = tar_map(
    values = transform_inputs_tibble,
    names = "names",
    tar_target(exp,
               input_obj %>%
               subset_obob5v5_from_obj %>%
               sc_transform_fgf1)
)


# n other cells = 13594
resampled_selection_inputs_list = list(
  tar_target(exp_neuron_1xo,
             exp_neuron_obob5v5 %>%
             sample_n_cells_from_obj(13594) %>%
             sc_transform_fgf1),
  tar_target(exp_neuron_2xo,
             exp_neuron_obob5v5 %>%
             sample_n_cells_from_obj(13594*2) %>%
             sc_transform_fgf1),
  tar_target(exp_all_n1xo,
             merge_neuron_and_other(exp_neuron_1xo, exp_other_obob5v5) %>%
             sc_transform_fgf1),
  tar_target(exp_all_n2xo,
             merge_neuron_and_other(exp_neuron_2xo, exp_other_obob5v5) %>%
             sc_transform_fgf1)
)

sce_inputs_tibble = tibble(
  exp_obj = rlang::syms(c("exp_other_obob5v5", "exp_neuron_1xo", "exp_neuron_2xo", "exp_all_n1xo", "exp_all_n2xo")),
  names = c("other", "neuron_1xo", "neuron_2xo", "all_n1xo", "all_n2xo")
)
sce_inputs_list = tar_map(
  values = sce_inputs_tibble,
  names = "names",
  tar_target(sce_00,
             exp_obj %>%
             convert_obj_to_sce) #note that this is not re-sctransformed, just filtered
)


batched_inputs_list = list(
  tar_target(literature_genes_df,
             get_literature_genes_df()),
  tar_target(other_var_genes,
             sce_00_other %>% gbr_retain_informative_genes %>% get_sce_genes),
  tar_target(neuron_1xo_var_genes,
             sce_00_neuron_1xo %>% gbr_retain_informative_genes %>% get_sce_genes),
  tar_target(all_n1xo_var_genes,
             sce_00_all_n1xo %>% gbr_retain_informative_genes %>% get_sce_genes),
  tar_target(informative_genes_combo,
             literature_genes_df$gene %>% 
             union(neuron_1xo_var_genes) %>% 
             union(all_n1xo_var_genes) %>%
             union(other_var_genes)),
  tar_target(exp_all_n1xo_nosct,
             merge_neuron_and_other(exp_neuron_1xo, exp_other_obob5v5)),
  tar_target(sce_00_all_n1xo_nosct,
             exp_all_n1xo_nosct %>%
             convert_obj_to_sce)
)


batched_selection_tibble = qs::qread("batched_selection_tibble.qs")
batched_selection_list = tar_map(
  values = batched_selection_tibble,
  names = "selection_combined_name",
  tar_target(sce_01_batched__TEMP,
             selection_obj %>%
             gbr_retain_informative_genes_preselected(n=20e3, preselected_genes=literature_genes_df$gene),
             cue=tar_cue(file=FALSE)),
  tar_target(genes_stat_100_batched,
             gbr_genes_stat_preselected(sce_01_batched__TEMP, base_genes=NULL, n_genes=100, batch=selection_batch, drop_nng=TRUE)),
  tar_target(genes_100_batched,
             get_selected_genes(genes_stat_100_batched)),
  tar_target(genes_stat_100_batched_litgenes,
             gbr_genes_stat_preselected(sce_01_batched__TEMP, base_genes=literature_genes_df$gene, n_genes=100, batch=selection_batch, drop_nng=TRUE)),
  tar_target(genes_stat_11.to.100_batched_litgenes,
             genes_stat_100_batched_litgenes %>% `[`(11:100,)),
  tar_target(genes_100_batched_litgenes,
             get_selected_genes(genes_stat_100_batched_litgenes))
)


batched_selection_eval_tibble = qs::qread("batched_selection_eval_tibble.qs")
batched_selection_eval_list = tar_map(
  values = batched_selection_eval_tibble,
  names = "combined_name",
  tar_target(sce_01_batched__TEMP,
             gbr_retain_informative_genes_preselected(evaluation_obj, n=20e3, preselected_genes=genes_stat_obj$gene),
             cue=tar_cue(file=FALSE)),
  tar_target(ctm_100_batched,
             gbr_get_celltype_mapping(sce=sce_01_batched__TEMP,
                                      selected_genes=genes_stat_obj,
                                      celltype_id=celltype_id,
                                      batch=evaluation_batch)),
  tar_target(lib_stat_batched,
             gbr_evaluate_lib(sce=sce_01_batched__TEMP,
                              selected_genes=genes_stat_obj,
                              celltype_id=celltype_id,
                              library.size_type = "series",
                              step=10,
                              batch=evaluation_batch)),
  tar_target(redundancy_stat_batched,
             gbr_calc_redundancy_stat(sce=sce_01_batched__TEMP,
                                      selected_genes=genes_stat_obj,
                                      celltype_id=celltype_id,
                                      batch=evaluation_batch)),
  tar_target(plotof_ctm_100_batched,
             plot_ctm(ctm_100_batched, name=combined_name))
)

reranked_eval_tibble = qs::qread("reranked_eval_tibble.qs")
reranked_eval_list = tar_map(
  values = reranked_eval_tibble,
  names = "combined_name",
  tar_target(genes_stat_reranked,
              rerank_genes_stat(genes_stat_obj, lib_stat_obj)),
  tar_target(ctm_100_batched_reranked,
             gbr_get_celltype_mapping(sce=sce_obj,
                                      selected_genes=genes_stat_reranked,
                                      celltype_id=celltype_id,
                                      batch=evaluation_batch)),
  tar_target(lib_stat_batched_reranked,
             gbr_evaluate_lib(sce=sce_obj,
                              selected_genes=genes_stat_reranked,
                              celltype_id=celltype_id,
                              library.size_type = "series",
                              step=10,
                              batch=evaluation_batch)),
  tar_target(redundancy_stat_batched_reranked,
             gbr_calc_redundancy_stat(sce=sce_obj,
                                      selected_genes=genes_stat_reranked,
                                      celltype_id=celltype_id,
                                      batch=evaluation_batch)),
  tar_target(plotof_ctm_100_batched_reranked,
             plot_ctm(ctm_100_batched_reranked, name=combined_name))
)


load_genes_stat_85_list = list(
  tar_target(genes_stat_85_path, "genes_stat_85.qs", format = "file"),
  tar_target(genes_stat_85, qs::qread(genes_stat_85_path))
)
gs85_tibble = qs::qread("gs85_tibble.qs")
gs85_tibble_list = tar_map(
  values = gs85_tibble,
  names = "combined_name",
  tar_target(sce_01_gs85_select__TEMP,
             selection_obj %>%
             gbr_retain_informative_genes_preselected(n=20e3, preselected_genes=genes_stat_85$gene),
             cue=tar_cue(file=FALSE)),
  tar_target(genes_stat_100_gs85,
             gbr_genes_stat_preselected(sce_01_gs85_select__TEMP, base_genes=genes_stat_85$gene, n_genes=100, batch=selection_batch, drop_nng=TRUE)),
  tar_target(sce_01_gs85_evaluate__TEMP,
             evaluation_obj %>%
             gbr_retain_informative_genes_preselected(n=20e3, preselected_genes=genes_stat_100_gs85$gene),
             cue=tar_cue(file=FALSE)),
  tar_target(ctm_100_gs85,
             gbr_get_celltype_mapping(sce=sce_01_gs85_evaluate__TEMP,
                                      selected_genes=genes_stat_100_gs85,
                                      celltype_id=celltype_id,
                                      batch=evaluation_batch)),
  tar_target(lib_stat_gs85,
             gbr_evaluate_lib(sce=sce_01_gs85_evaluate__TEMP,
                              selected_genes=genes_stat_100_gs85,
                              celltype_id=celltype_id,
                              library.size_type = "series",
                              step=10,
                              batch=evaluation_batch)),
  tar_target(redundancy_stat_gs85,
             gbr_calc_redundancy_stat(sce=sce_01_gs85_evaluate__TEMP,
                                      selected_genes=genes_stat_100_gs85,
                                      celltype_id=celltype_id,
                                      batch=evaluation_batch)),
  tar_target(plotof_ctm_100_gs85,
             plot_ctm(ctm_100_gs85, name=combined_name))
)



load_genes_stat_97_list = list(
  tar_target(genes_stat_97_path, "genes_stat_97.qs", format = "file"),
  tar_target(genes_stat_97, qs::qread(genes_stat_97_path))
)
gs97_tibble = qs::qread("gs97_tibble.qs")
gs97_tibble_list = tar_map(
  values = gs97_tibble,
  names = "combined_name",
  tar_target(sce_01_gs97_select__TEMP,
             selection_obj %>%
             gbr_retain_informative_genes_preselected(n=20e3, preselected_genes=genes_stat_97$gene),
             cue=tar_cue(file=FALSE)),
  tar_target(genes_stat_100_gs97,
             gbr_genes_stat_preselected_discard(sce_01_gs97_select__TEMP, base_genes=genes_stat_97$gene, n_genes=100, batch=selection_batch, drop_nng=TRUE, discard_genes=c('Sprr2a2'))),
  tar_target(sce_01_gs97_evaluate__TEMP,
             evaluation_obj %>%
             gbr_retain_informative_genes_preselected(n=20e3, preselected_genes=genes_stat_100_gs97$gene),
             cue=tar_cue(file=FALSE)),
  tar_target(ctm_100_gs97,
             gbr_get_celltype_mapping(sce=sce_01_gs97_evaluate__TEMP,
                                      selected_genes=genes_stat_100_gs97,
                                      celltype_id=celltype_id,
                                      batch=evaluation_batch)),
  tar_target(lib_stat_gs97,
             gbr_evaluate_lib(sce=sce_01_gs97_evaluate__TEMP,
                              selected_genes=genes_stat_100_gs97,
                              celltype_id=celltype_id,
                              library.size_type = "series",
                              step=10,
                              batch=evaluation_batch)),
  tar_target(redundancy_stat_gs97,
             gbr_calc_redundancy_stat(sce=sce_01_gs97_evaluate__TEMP,
                                      selected_genes=genes_stat_100_gs97,
                                      celltype_id=celltype_id,
                                      batch=evaluation_batch)),
  tar_target(plotof_ctm_100_gs97,
             plot_ctm(ctm_100_gs97, name=combined_name))
)



# combine_tibble = qs::qread('combine_tibble.qs')
# combine_targets_list = tar_map(
#   values = combine_tibble,
#   names = "label",
#   tar_target(genes_stat_pl_100_00,
#              genes_stat_obj %>%
#              mutate(label = label))
# )


# sw_selection_tibble = qs::qread('sw_selection_tibble.qs')
# sw_selection_targets_list = tar_map(
#   values = sw_selection_tibble,
#   names = "name",
#   tar_target(genes_stat_1st,
#              do_first_sw_selection(first_selection_obj, n_first_selection)),
#   tar_target(genes_stat_2nd,
#              do_next_sw_selection(second_selection_obj, n_second_selection, genes_stat_1st)),
#   tar_target(genes_stat_3rd,
#              do_next_sw_selection(third_selection_obj, n_third_selection, genes_stat_2nd))
# )


# sw_evaluation_tibble = qs::qread("sw_evaluation_tibble.qs")
# sw_evaluation_targets_list = tar_map(
#   values = sw_evaluation_tibble,
#   names = "name",
#   tar_target(ctm_3rd,
#              make_ctm_sw_selection(evaluation_obj,
#                                    celltype_id,
#                                    genes_stat_3rd_obj,
#                                    batch=NULL),
#              priority=0.5),
#   tar_target(lib_stat_3rd,
#              make_lib_stat_sw_selection(evaluation_obj,
#                                         genes_stat_3rd_obj,
#                                         celltype_id,
#                                         library.size_type = "single",
#                                         step=10,
#                                         batch=NULL),
#              priority=0.99),
#   tar_target(cell_score_stat_3rd,
#              make_cell_score_stat_sw_selection(lib_stat_3rd, name, evaluation_obj)),
#   tar_target(csss_sw,
#              summarise_cell_score_stat_sw_selection(cell_score_stat_3rd),
#              deployment = "main")
# )

list(load_inputs_list,
transform_inputs_list,
resampled_selection_inputs_list,
sce_inputs_list,
batched_inputs_list,
batched_selection_list,
batched_selection_eval_list,
reranked_eval_list,
load_genes_stat_85_list,
gs85_tibble_list,
load_genes_stat_97_list,
gs97_tibble_list
)


## why don't these genes_stat_pl_100_00 show up in tar_meta?
# turns out meta only updated when target created
# tar_combines should be wrapped in a list
# return empty if required targets don't exist in tar_meta
# combine_tibble = qs::qread('combine_tibble.qs')
# combine_targets_list = tar_map(
#   values = combine_tibble,
#   names = "label",
#   tar_target(genes_stat_pl_100_00,
#              genes_stat_obj %>%
#              mutate(label = label))
# )
  # combine_targets_list,
  # tarchetypes::tar_combine(combined_genes_stat, 
  #                            make_split_objs %>% 
  #                            tar_select_targets(contains(tar_meta(fields = NULL) %>% 
  #                             filter(is.na(error)) %>% 
  #                             filter(str_detect(name, 'genes_stat_pl_100')) %>% 
  #                             pull(name))) %>%
  #                             tar_select_targets(!contains(tar_errored()))
  # )
  #   tarchetypes::tar_combine(combined_genes_stat, 
  #                            combine_targets_list %>% 
  #                             tar_select_targets(contains("genes_stat_pl_100_00")) %>%
  #                             tar_select_targets(!contains(c('`match` must be a character vector of non empty strings.', tar_errored())))
  # ),

# single_selection_100g_tibble = qs::qread("single_selection_100g_tibble.qs")
# single_selection_100g_list = tar_map(
#   values = single_selection_100g_tibble,
#   names = "selection_name",
#   tar_target(sce_01__TEMP,
#              selection_obj %>%
#              gbr_retain_informative_genes),
#   tar_target(genes_stat_100,
#              gbr_genes_stat(sce_01__TEMP, 100)),
#   tar_target(genes_100,
#              get_selected_genes(genes_stat_100))  
# )


# single_selection_100g_eval_tibble = qs::qread("single_selection_100g_eval_tibble.qs")
# single_selection_100g_eval_list = tar_map(
#   values = single_selection_100g_eval_tibble,
#   names = "combined_name",
#   tar_target(sce_01__TEMP,
#              gbr_retain_informative_genes_preselected(selection_obj, n=NULL, preselected_genes=genes_stat_obj$gene)),
#   tar_target(ctm_100,
#              gbr_get_celltype_mapping(sce=sce_01__TEMP, selected_genes=genes_stat_obj, celltype_id=celltype_id)),
#   tar_target(lib_stat,
#              gbr_evaluate_lib(sce=sce_01__TEMP, selected_genes=genes_stat_obj, celltype_id=celltype_id,  library.size_type = "series", step=10)),
#   tar_target(redundancy_stat,
#              gbr_calc_redundancy_stat(sce=sce_01__TEMP, selected_genes=genes_stat_obj, celltype_id=celltype_id))
# )