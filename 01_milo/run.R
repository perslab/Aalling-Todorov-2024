# This is a helper script to run the pipeline.
# Choose how to execute the pipeline below.
# See https://books.ropensci.org/targets/hpc.html
# to learn about your options.
library(tidyverse)

STAGE = 5
SHORTCUT = FALSE
n_workers = 80

do_staging = function(){
    if (STAGE == 1){
        targets::tar_make_future(workers=n_workers)
    }
    if (STAGE == 2){
        combination_recipe = qs::qread("combination_recipe.qs")
        target_names = combination_recipe$output_name
        # targets::tar_invalidate(target_names)
        targets::tar_make(names = unlist(target_names), shortcut=SHORTCUT)
    }
    if (STAGE == 3){
        restored_recipe = qs::qread("restored_recipe.qs")
        target_names = c(paste0('restored_df_', c(restored_recipe$restored_suffix)),
                         paste0('restored_summary_', restored_recipe$restored_suffix),
                         c("all_restored_summary"))
        # targets::tar_invalidate(target_names)
        targets::tar_make(names = unlist(target_names), shortcut=SHORTCUT)
    }
    if (STAGE == 4){
        annotate_nh_grouping_recipe = qs::qread("annotate_nh_grouping_recipe.qs")
        target_names = tidyr::crossing(target_base = c("da_results_nhg_", "nhgc_"),
            output_suffix = annotate_nh_grouping_recipe$da_results_nhg_output_suffix) %>%
            rowwise %>% 
            mutate(target_names = paste0(target_base, output_suffix)) %>% 
            ungroup %>% pull(target_names)
        # targets::tar_invalidate(target_names)
        targets::tar_make(names = unlist(target_names), shortcut=SHORTCUT)
        # targets::tar_make_future(names = unlist(target_names), shortcut=SHORTCUT, workers=10)
    }
    if (STAGE == 5){
        # stage_05_meta = qs::qread("stage_05_meta.qs")
        # target_names = stage_05_meta %>%
        #     filter(is.na(error)) %>%
        #     pull(name)
        # targets::tar_invalidate(target_names)
        # targets::tar_make(names = unlist(target_names))
        # targets::tar_make_future(names = unlist(target_names), shortcut=SHORTCUT, workers=n_workers)
        targets::tar_make_future(workers=n_workers)
    }
    if (STAGE == 6){
        combination_recipe_degs = qs::qread("combined_recipe_degs.qs")
        target_names = combination_recipe_degs$output_name
        targets::tar_make(names = unlist(target_names), shortcut=SHORTCUT)
    }
    if (STAGE == 7){
        target_names = c("combined_deg_seurat_formatted", "combined_gsea_seurat_formatted")
        targets::tar_make(names = unlist(target_names), shortcut=SHORTCUT)
    }
}

do_staging()

# targets::tar_make_clustermq(workers = 2) # nolint
# targets::tar_make_future(workers = 2) # nolint
