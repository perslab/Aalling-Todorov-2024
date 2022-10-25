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
source(paste0("genebasisr.R"))
source(paste0("prep_probe_selection.R"))

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
# tar_source()
# source("other_functions.R") # Source other scripts as needed. # nolint

# Replace the target list below with your own:
clusters_tibble = qs::qread('clusters_tibble.qs')
make_split_objs = tar_map(
  values = clusters_tibble,
  names = label,
  tar_target(obj_ngo_pl,
             prep_obj_for_milo(object) %>%
             single_split(label) %>%
             prep_obj(pldf_obob5v5) %>%
             make_new_seurat_ngo %>%
             sc_transform_fgf1),
  tar_target(sce_ngo_pl_00, convert_obj_to_sce(obj_ngo_pl)),
  tar_target(sce_ngo_pl_01,
             gbr_retain_informative_genes(sce_ngo_pl_00)),
  tar_target(genes_stat_pl_100,
             gbr_genes_stat(sce_ngo_pl_01, 100)),
  tar_target(genes_50_vec,
             get_selected_genes(genes_stat_pl_100)),
  tar_target(ctm_pl,
             gbr_get_celltype_mapping(sce=sce_ngo_pl_01, selected_genes=genes_stat_pl_100, celltype_id="polar_label")),
  tar_target(lib_stat_pl,
             gbr_evaluate_lib(sce=sce_ngo_pl_01, selected_genes=genes_stat_pl_100, celltype_id="polar_label", library.size_type = "series", step=1)),
  tar_target(redundancy_stat_pl,
             gbr_calc_redundancy_stat(sce=sce_ngo_pl_01, selected_genes=genes_stat_pl_100, celltype_id="polar_label"))
)



big_ngo_targets = list(
  tar_target(exp_labelled_other_path, paste0(Sys.getenv("PROJECT_DIR"), '/00_preprocessing/_targets/objects/exp_labelled_other'), format = "file"),
  tar_target(exp_labelled_neuron_path, paste0(Sys.getenv("PROJECT_DIR"), '/00_preprocessing/_targets/objects/exp_RElabelled_neuron'), format = "file"),
  tar_target(pldf_obob5v5_path, paste0(Sys.getenv("PROJECT_DIR"), '/01_milo/_targets/objects/combined_pldf_obob5v5'), format = "file"),
  tar_target(exp_labelled_other, qs::qread(exp_labelled_other_path)),
  tar_target(exp_labelled_neuron, qs::qread(exp_labelled_neuron_path)),
  tar_target(pldf_obob5v5, 
             qs::qread(pldf_obob5v5_path)),
  tar_target(exp_labelled_other_ngo,
             exp_labelled_other %>%
             prep_obj(pldf_obob5v5) %>%
             make_new_seurat_ngo %>%
             sc_transform_fgf1),
  tar_target(exp_labelled_neuron_ngo,
             exp_labelled_neuron %>%
               prep_obj(pldf_obob5v5) %>%
               make_new_seurat_ngo %>%
               sc_transform_fgf1),
  tar_target(exp_ngo_all,
             merge_neuron_and_other(exp_labelled_neuron_ngo, exp_labelled_other_ngo) %>%
             sc_transform_fgf1),
  tar_target(exp_ngo_all_40k,
             sample_n_cells_from_obj(exp_ngo_all, 40000) %>%
             sc_transform_fgf1),
    tar_target(exp_ngo_all_70k,
             sample_n_cells_from_obj(exp_ngo_all, 70000) %>%
             sc_transform_fgf1),
  tar_target(exp_labelled_neuron_ngo_40k,
             sample_n_cells_from_obj(exp_labelled_neuron_ngo, 40000) %>%
             sc_transform_fgf1)
  
)


big_targets_100g_tibble <- tibble(
  exp_ngo_obj = rlang::syms(c("exp_ngo_all_70k", "exp_ngo_all_40k", "exp_labelled_neuron_ngo_40k", "exp_labelled_other_ngo")),
  names = c("all_70k", "all_40k", "neuron_40k", "other")
)
big_targets_100g = tar_map(
  values = big_targets_100g_tibble,
  names = "names",
  tar_target(sce_ngo_00,
             convert_obj_to_sce(exp_ngo_obj)),
  tar_target(sce_ngo_01, gbr_retain_informative_genes(sce_ngo_00)),
    tar_target(genes_stat_100,
               gbr_genes_stat(sce_ngo_01, 100)),
  tar_target(genes_100,
             get_selected_genes(genes_stat_100)),
  tar_target(ctm_00_100,
             gbr_get_celltype_mapping(sce=sce_ngo_01, selected_genes=genes_stat_100, celltype_id="polar_label")),
  tar_target(lib_stat_00,
             gbr_evaluate_lib(sce=sce_ngo_01, selected_genes=genes_stat_100, celltype_id="polar_label",  library.size_type = "series", step=10)),
  tar_target(redundancy_stat,
             gbr_calc_redundancy_stat(sce=sce_ngo_01, selected_genes=genes_stat_100, celltype_id="polar_label"))
  
)


preselected_genes_from_milo = list(
  tar_target(ngo_panel_5ea_path, paste0(Sys.getenv("PROJECT_DIR"), '/01_milo/_targets/objects/ngo_panel_5ea'), format = "file"),
  tar_target(ngo_panel_5ea, qs::qread(ngo_panel_5ea_path))
)


preselected_genes_targets_tibble <- tibble(
  exp_ngo_obj = rlang::syms(c("exp_ngo_all_40k", "exp_ngo_all_70k")),
  names = c("all_40k_pre5ea", "all_70k_pre5ea")
)
preselected_genes_targets = tar_map(
  values = preselected_genes_targets_tibble,
  names = "names",
  tar_target(sce_ngo_00,
             convert_obj_to_sce(exp_ngo_obj)),
  tar_target(sce_ngo_01,
             gbr_retain_informative_genes_preselected(sce_ngo_00, n=NULL, preselected_genes=ngo_panel_5ea)),
    tar_target(genes_stat_100,
               gbr_genes_stat_preselected(sce_ngo_01, ngo_panel_5ea, n_genes=100, batch=NULL)),
  tar_target(genes_100,
             get_selected_genes(genes_stat_100)),
  tar_target(ctm_00_100,
             gbr_get_celltype_mapping(sce=sce_ngo_01, selected_genes=genes_stat_100, celltype_id="polar_label")),
    tar_target(ctm_01_100,
             gbr_get_celltype_mapping(sce=sce_ngo_01, selected_genes=genes_stat_100, celltype_id="labels")),
  tar_target(lib_stat_00,
             gbr_evaluate_lib(sce=sce_ngo_01, selected_genes=genes_stat_100, celltype_id="polar_label",  library.size_type = "series", step=10)),
  tar_target(lib_stat_01,
             gbr_evaluate_lib(sce=sce_ngo_01, selected_genes=genes_stat_100, celltype_id="labels",  library.size_type = "series", step=10)),
  tar_target(redundancy_stat,
             gbr_calc_redundancy_stat(sce=sce_ngo_01, selected_genes=genes_stat_100, celltype_id="polar_label"))
)  


combine_tibble = qs::qread('combine_tibble.qs')
combine_targets_list = tar_map(
  values = combine_tibble,
  names = "label",
  tar_target(genes_stat_pl_100_00,
             genes_stat_obj %>%
             mutate(label = label))
)


sw_selection_tibble = qs::qread('sw_selection_tibble.qs')
sw_selection_targets_list = tar_map(
  values = sw_selection_tibble,
  names = "name",
  tar_target(genes_stat_1st,
             do_first_sw_selection(first_selection_obj, n_first_selection)),
  tar_target(genes_stat_2nd,
             do_next_sw_selection(second_selection_obj, n_second_selection, genes_stat_1st)),
  tar_target(genes_stat_3rd,
             do_next_sw_selection(third_selection_obj, n_third_selection, genes_stat_2nd))
)


sw_evaluation_tibble = qs::qread("sw_evaluation_tibble.qs")
sw_evaluation_targets_list = tar_map(
  values = sw_evaluation_tibble,
  names = "name",
  tar_target(ctm_3rd,
             make_ctm_sw_selection(evaluation_obj,
                                   celltype_id,
                                   genes_stat_3rd_obj,
                                   batch=NULL),
             priority=0.5),
  tar_target(lib_stat_3rd,
             make_lib_stat_sw_selection(evaluation_obj,
                                        genes_stat_3rd_obj,
                                        celltype_id,
                                        library.size_type = "single",
                                        step=10,
                                        batch=NULL),
             priority=0.99),
  tar_target(cell_score_stat_3rd,
             make_cell_score_stat_sw_selection(lib_stat_3rd, name, evaluation_obj)),
  tar_target(csss_sw,
             summarise_cell_score_stat_sw_selection(cell_score_stat_3rd),
             deployment = "main")
)

list(big_ngo_targets,
  make_split_objs,
  big_targets_100g,
  preselected_genes_from_milo,
  preselected_genes_targets,
  combine_targets_list,
  tarchetypes::tar_combine(combined_genes_stat, 
                            combine_targets_list %>% 
                            tar_select_targets(contains("genes_stat_pl_100_00")) %>%
                            tar_select_targets(!contains(c('`match` must be a character vector of non empty strings.', tar_errored()))),
                           deployment="main"),
  sw_selection_targets_list,
  sw_evaluation_targets_list,
  tarchetypes::tar_combine(combined_csss_sw, 
                            sw_evaluation_targets_list %>% 
                            tar_select_targets(contains("csss_sw")) %>%
                            tar_select_targets(!contains(c('`match` must be a character vector of non empty strings.', tar_errored()))),
                          deployment="main")

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