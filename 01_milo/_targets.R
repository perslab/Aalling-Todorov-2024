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
  error = "continue"
  # Set other options as needed.
)

source("../00_preprocessing/splitwrapper.R")
source(paste0("../00_preprocessing/preprocessing.R"))
source(paste0("milo.R"))

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
# tar_source()
# source("other_functions.R") # Source other scripts as needed. # nolint

clusters_tibble = qs::qread('clusters_tibble.qs')
make_split_objs = tar_map(
  values = clusters_tibble,
  names = label,
  tar_target(obj,
             prep_obj_for_milo(object) %>%
             single_split(label) %>%
             sc_transform_fgf1 %>%
             run_sct_chaser),
  tar_target(milo, 
             convert_obj_to_sce(obj) %>%
             make_milo),
  tar_target(design_df,
             make_design_df(milo)),
  tar_target(mm,
             make_model_matrix(design_df))
  
)


da_recipe = qs::qread('da_recipe.qs')
make_da_results = tar_map(
  values = da_recipe,
  names = label,
  tar_target(da_results_00,
             get_da_results(milo_obj=milo_obj,
                            model_matrix=mm, design_df=design_df, model_contrasts_str=contrast)),
  tar_target(da_results_01,
             add_polarity_to_da_results(da_results_00, spatial_fdr_cutoff=0.1)),
  tar_target(da_results_02,
             annotate_nhoods(da_results_01, milo_obj, col="labels")),
  tar_target(pldf,
             make_polar_labels_df(milo_obj, da_results_01, top_frac=0.1))
)



list(    
  tar_target(exp_labelled_other_path, paste0(Sys.getenv("PROJECT_DIR"), '/00_preprocessing/_targets/objects/exp_labelled_other'), format = "file"),
  tar_target(exp_labelled_neuron_path, paste0(Sys.getenv("PROJECT_DIR"), '/00_preprocessing/_targets/objects/exp_RElabelled_neuron'), format = "file"),
  tar_target(exp_labelled_other, qs::qread(exp_labelled_other_path)),
  tar_target(exp_labelled_neuron, qs::qread(exp_labelled_neuron_path)),
#   tar_target(milo_obj,
#              make_milo_from_sce(sce_obj) %>%
#              .milo_buildGraph() %>%
#              .milo_makeNhoods() %>%
#              .milo_countCells() %>%
#              .milo_calcNhoodDistance()),
  make_split_objs,
  make_da_results,
  # tarchetypes::tar_combine(combined, 
  #                          make_da_results %>% 
  #                          tar_select_names(starts_with("da_results_01"),
  #                                           ends_with("obob5v5")))
    # tarchetypes::tar_combine(combined_da_results_obob5v5, 
    #                          make_da_results %>% 
    #                          tar_select_targets(starts_with("da_results_01")) %>% 
    #                          tar_select_targets(ends_with("obob5v5")),
    #                          command = rbind(!!!.x)
    # )
    tarchetypes::tar_combine(combined_da_results_obob5v5, 
                             make_da_results %>% 
                             tar_select_targets(starts_with("da_results_02")) %>% 
                             tar_select_targets(ends_with("obob5v5")) %>%
                             tar_select_targets(!contains(tar_errored()))
    ),
    tarchetypes::tar_combine(combined_pldf_obob5v5, 
                             make_da_results %>% 
                             tar_select_targets(contains(tar_meta(fields = NULL) %>% 
                              filter(is.na(error)) %>% 
                              filter(str_detect(name, 'pldf')) %>% 
                              filter(str_detect(name, 'obob5v5')) %>% 
                              pull(name)))
    )
)

