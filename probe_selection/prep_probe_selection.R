library(tidyverse)
library(Seurat)

prep_pldf = function(pldf){
    pldf = pldf %>% 
    mutate(polar_label = case_when(
        stringr::str_detect(polar_label, 'NFOL') ~ 'NFOL.none',
        TRUE ~ polar_label
    ))
    pldf
}


prep_obj = function(obj, pldf){
    obj = AddMetaData(obj, metadata = pldf)
    obj@meta.data = obj@meta.data %>%
      mutate(polar_label = ifelse(is.na(polar_label), stringr::str_replace(labels, '-', '__') %>% paste0('.none'), polar_label))
    obj
}
