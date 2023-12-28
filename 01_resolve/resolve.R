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
                # filter(. != 'Agrp') %>% #this is not common between both datasets
                pull
        xenium_genes
}

transfer_data_cca_00 = function(xe_obj, obj_fgf1, refdata_column){
    xenium_genes = get_xe_genes(xe_obj, obj_fgf1)
    xenium_genes_all = get_xe_genes_all(xe_obj)
    obj_fgf1 = obj_fgf1 %>% Seurat::RunUMAP(assay='SCT', slot='counts', features = xenium_genes, return.model = TRUE)
    refdata = obj_fgf1 %>% `[[` %>% pull(refdata_column)
    # anch <- FindTransferAnchors(reference = obj_fgf1,
    #                             query = xe_obj,
    #                             features = xenium_genes,
    #                             normalization.method = "SCT",
    #                             reduction='cca',
    #                             recompute.residuals = T,
    #                             dims = 1:20)
    anch = FindTransferAnchors(reference = obj_fgf1, query = xe_obj, features = xenium_genes,  
                               normalization.method = 'SCT',
                               reduction = 'cca',
                               query.assay = 'SCT',
                               reference.assay = 'SCT',
                               recompute.residuals = F,
                               verbose=TRUE)
    predictions <- TransferData(anchorset = anch, refdata = refdata, weight.reduction='cca') %>% 
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


transfer_data_cca_00_unimodal = function(xe_obj, obj_fgf1, refdata_column){
    xenium_genes = get_xe_genes(xe_obj, obj_fgf1)
    xenium_genes_all = get_xe_genes_all(xe_obj)
    obj_fgf1 = obj_fgf1 %>% 
        Seurat::RunUMAP(assay='SCT', slot='counts', features = xenium_genes, return.model = TRUE)
        # RunCCA(object1 = ., object2 = xe_obj, features=xenium_genes, num.cc = 30)
    refdata = obj_fgf1 %>% `[[` %>% pull(refdata_column)
    # anch <- FindTransferAnchors(reference = obj_fgf1,
    #                             query = xe_obj,
    #                             features = xenium_genes,
    #                             normalization.method = "SCT",
    #                             reduction='cca',
    #                             recompute.residuals = T,
    #                             dims=1:20)
    anch = FindTransferAnchors(reference = obj_fgf1, query = xe_obj, features = xenium_genes,  
                               normalization.method = 'SCT',
                               reduction = 'cca',
                               query.assay = 'SCT',
                               reference.assay = 'SCT',
                               recompute.residuals = F,
                               verbose=TRUE)
    # predictions <- TransferData(anchorset = anch, refdata = refdata, weight.reduction = "cca", dims=1:30) %>% 
    #     dplyr::select(predicted.id, prediction.score.max)
    xe_obj = MapQuery(anchorset = anch,
                           reference = obj_fgf1,
                           query = xe_obj,
                           refdata = refdata,
                           reduction.model = "umap")
    predictions = xe_obj %>% 
                           `[[` %>%
                  dplyr::select(predicted.id, predicted.id.score)
    predictions[[refdata_column]] = predictions$predicted.id
    predictions[[paste0(refdata_column, '_prediction.score.max')]] = predictions$predicted.id.score
    predictions = predictions %>% select(-predicted.id, -predicted.id.score)
    xe_obj <- AddMetaData(xe_obj, metadata = predictions)
    Misc(object = xe_obj, slot = paste0("predictions_", refdata_column)) <- predictions
    xe_obj
}

reclass_by_gene = function(xe_obj, gene, gene_thr, corrected_class){
    meta = xe_obj %>% `[[`
    meta = meta %>%
        mutate(gene_count = xe_obj[["Xenium"]]@counts[gene, ]) %>%
        mutate(cell_class = case_when(gene_count >= gene_thr ~ corrected_class,
                                   TRUE ~ cell_class)) %>%
        mutate(cell_class_prediction.score.max = case_when(gene_count >= gene_thr ~ 1,
                                                           TRUE ~ cell_class_prediction.score.max)) %>%
        select(-gene_count)
    xe_obj@meta.data = meta
    xe_obj
}

reclass_by_gene_hilo = function(xe_obj, gene, gene_thr_lo, gene_thr_hi, wrong_class, correct_class){
    meta = xe_obj %>% `[[`
    meta = meta %>%
        mutate(gene_count = xe_obj[["Xenium"]]@counts[gene, ]) %>%
        mutate(cell_class = case_when((gene_count >= gene_thr_hi &
                                      gene_count >= gene_thr_lo &
                                      cell_class == wrong_class) ~ correct_class,
                                      (gene_count <= gene_thr_hi &
                                      gene_count >= gene_thr_lo &
                                      cell_class == wrong_class) ~ 'blah',
                                      TRUE ~ cell_class)) %>%
        mutate(cell_class_prediction.score.max = case_when((gene_count >= gene_thr_hi &
                                                            gene_count >= gene_thr_lo &
                                                            cell_class == wrong_class) ~ 1,
                                                            (gene_count <= gene_thr_hi &
                                                            gene_count >= gene_thr_lo &
                                                            cell_class == wrong_class) ~ 1,
                                                            TRUE ~ cell_class_prediction.score.max)) %>%
        select(-gene_count)
    xe_obj@meta.data = meta
    xe_obj
}

relabel_by_gene_hilo = function(xe_obj, gene, gene_thr_lo, gene_thr_hi, correct_label){
    meta = xe_obj %>% `[[`
    meta = meta %>%
        mutate(gene_count = xe_obj[["Xenium"]]@counts[gene, ]) %>%
        mutate(labels = case_when((gene_count >= gene_thr_hi &
                                      gene_count >= gene_thr_lo &
                                      !str_detect(labels, gene)) ~ correct_label,
                                      (gene_count <= gene_thr_hi &
                                      gene_count >= gene_thr_lo &
                                      !str_detect(labels, gene)) ~ 'blah',
                                      TRUE ~ labels)) %>%
        mutate(labels_prediction.score.max = case_when((gene_count >= gene_thr_hi &
                                                            gene_count >= gene_thr_lo &
                                                            !str_detect(labels, gene)) ~ 1,
                                                            (gene_count <= gene_thr_hi &
                                                            gene_count >= gene_thr_lo &
                                                            !str_detect(labels, gene)) ~ 1,
                                                            TRUE ~ labels_prediction.score.max)) %>%
        select(-gene_count)
    xe_obj@meta.data = meta
    xe_obj
}

labels_from_polar_label = function(xe_obj){
  meta = xe_obj %>% `[[` %>%
    mutate(labels = polar_label %>% 
                    str_replace(fixed('.neg'), '') %>% 
                    str_replace(fixed('.none'), '') %>% 
                    str_replace(fixed('.pos'), '')) %>%
    mutate(labels_prediction.score.max = 2)
  xe_obj@meta.data = meta
  xe_obj
}

drop_blah_class = function(xe_obj){
  meta = xe_obj %>% `[[`
  good_cells = meta %>%
    filter(cell_class != 'blah') %>%
    rownames
  xe_obj = xe_obj %>% subset(cells = good_cells)
  xe_obj
}

drop_blah_labels = function(xe_obj){
  meta = xe_obj %>% `[[`
  good_cells = meta %>%
    filter(labels != 'blah') %>%
    rownames
  xe_obj = xe_obj %>% subset(cells = good_cells)
  xe_obj
}

transfer_data_cca_00_polar_label_unimodal = function(xe_obj, obj_fgf1, selected_label) {
  tryCatch(
    {
      xe_obj = xe_obj %>%
        split_cell_labels(selected_label) %>%
        sc_transform_resolve
      obj_fgf1 = obj_fgf1 %>%
        split_cell_labels(selected_label) %>%
        sc_transform_fgf1
      xe_obj = transfer_data_cca_00_unimodal(xe_obj, obj_fgf1, "polar_label")
      xe_obj
    },
    error = function(e) {
      message("An error occurred: ", e$message)
      return(NA)
    }
  )
}