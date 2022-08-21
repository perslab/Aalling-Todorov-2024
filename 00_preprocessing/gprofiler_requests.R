library(reticulate)
library(tidyverse)

source_python('gprofiler_requests.py')

.get_gosts_df = function(top_tags_item, fdr_cutoff=0.05){
    item_name = paste0(top_tags_item$result$cluster_name,
                       '___',
                       top_tags_item$result$comparison)
    gene_names = top_tags_item$result$table %>%
        filter(FDR < fdr_cutoff) %>%
        rownames
    module_members = list(gene_names) %>% rlang::set_names(item_name)
    df = make_gost_df(module_members)
    df
}


get_gosts_df = function(toptags_list){
    get_gosts_df_p = purrr::safely(.get_gosts_df, quiet=FALSE)
    gpdf_list = toptags_list %>%
      map(~get_gosts_df_p(.x))
    names_vec = gpdf_list %>% 
        map(~paste0(.x$result$cluster_name, '___', .x$result$comparison))
    gpdf_list = rlang::set_names(gpdf_list, names_vec)
    gpdf_list
}

bind_gosts_dfs = function(gpdf_list){
    gpdfs_bound = gpdf_list %>%
     keep( ~ !is.null(.x$result) ) %>%
     map(~.x$result) %>%
     data.table::rbindlist()
    gpdfs_bound
}

