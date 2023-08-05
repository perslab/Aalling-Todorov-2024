filter_down_cells = function(xenium.obj){
    keep_cells = xenium.obj %>% `[[` %>%
        filter(cell_area < 7000) %>%
        filter(cell_area > 100) %>%
        filter(avg_confidence >= 0.99) %>%
        filter(nCount_Xenium >= 10) %>%
        filter(nCount_Xenium <= 600) %>%
        rownames
    xenium.obj = subset(xenium.obj, cells = keep_cells)
    xenium.obj
}

sc_transform_resolve = function(xenium.obj, keep_cells){
    xenium.obj %>%
    Seurat::SCTransform(assay='Xenium',
                        method="glmGamPoi",
                        vst.flavor="v2",
                        verbose=TRUE)
}

resolve_find_all_markers = function(xe_obj, idents='predicted.label'){
    Idents(xe_obj) = idents
    all_markers = xe_obj %>% FindAllMarkers(assay='SCT')
    all_markers
}


subset_xe_obj_cell_class = function(xe_obj, cell_class_score_thr = 0.99){
    barcodes_good_class_score = xe_obj %>%
        `[[` %>%
        filter(predicted.cell_class.score >= cell_class_score_thr) %>%
        rownames
    
    xe_obj = xe_obj %>% subset(cells = barcodes_good_class_score)
    xe_obj
}

get_xe_genes_all = function(xe_obj){
    xenium_genes_all = xe_obj %>% 
             `@`('assays') %>% 
             `$`('Xenium') %>% 
             `@`('meta.features') %>% 
             rownames_to_column %>%
             filter(rowname != 'Lmx1a') %>%
             pull(rowname)
    xenium_genes_all
}

get_xe_genes = function(xe_obj, obj_fgf1){
    ref_genes_all= obj_fgf1 %>% 
                `@`('assays') %>% 
                `$`('SCT') %>% 
                `@`('meta.features') %>% 
                rownames

        xenium_genes = xe_obj %>% 
                `@`('assays') %>% 
                `$`('Xenium') %>% 
                `@`('meta.features') %>% 
                rownames %>%
                intersect(ref_genes_all) %>%
                tibble %>%
                filter(. != 'Agrp') %>% #this is not common between both datasets
                pull
        xenium_genes
}

transfer_data_cca_00 = function(xe_obj, obj_fgf1, refdata_column){
    xenium_genes = get_xe_genes(xe_obj, obj_fgf1)
    xenium_genes_all = get_xe_genes_all(xe_obj)
    obj_fgf1 = obj_fgf1 %>% Seurat::RunUMAP(assay='SCT', slot='counts', features = xenium_genes, return.model = TRUE)
    refdata = obj_fgf1 %>% `[[` %>% pull(refdata_column)
    anch <- FindTransferAnchors(reference = obj_fgf1,
                                query = xe_obj,
                                features = xenium_genes,
                                normalization.method = "SCT",
                                recompute.residuals = F)
    predictions <- TransferData(anchorset = anch, refdata = refdata) %>% 
        dplyr::select(predicted.id, prediction.score.max)
    predictions[[refdata_column]] = predictions$predicted.id
    predictions[[paste0(refdata_column, '_prediction.score.max')]] = predictions$prediction.score.max
    predictions = predictions %>% select(-predicted.id, -prediction.score.max)
    xe_obj <- AddMetaData(xe_obj, metadata = predictions)
    Misc(object = xe_obj, slot = paste0("predictions_", refdata_column)) <- predictions
    xe_obj
}

split_cell_class = function(xe_obj, selected_cell_class){
    selected_cells = xe_obj %>% 
        `[[` %>% 
        filter(cell_class == selected_cell_class) %>%
        rownames
    xe_obj = xe_obj %>% subset(cells = selected_cells)
    xe_obj
}

split_cell_labels = function(xe_obj, selected_labels){
    selected_cells = xe_obj %>% 
        `[[` %>% 
        filter(labels == selected_labels) %>%
        rownames
    xe_obj = xe_obj %>% subset(cells = selected_cells)
    xe_obj
}

transfer_data_cca_00_polar_label = function(xe_obj, obj_fgf1, selected_label) {
  tryCatch(
    {
      xe_obj = xe_obj %>%
        split_cell_labels(selected_label) %>%
        sc_transform_resolve
      obj_fgf1 = obj_fgf1 %>%
        split_cell_labels(selected_label) %>%
        sc_transform_fgf1
      xe_obj = transfer_data_cca_00(xe_obj, obj_fgf1, "polar_label")
      xe_obj
    },
    error = function(e) {
      message("An error occurred: ", e$message)
      return(NA)
    }
  )
}
