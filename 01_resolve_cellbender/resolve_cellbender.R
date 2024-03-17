make_obj_d5 = function(nhgc_d5_path, neuron_path, other_path, tany_path){
    nhgc_d5 = qs::qread(nhgc_d5_path) 
    neuron_meta = qs::qread(neuron_path) %>% `[[`
    other_meta = qs::qread(other_path) %>% `[[`
    tany_meta = qs::qread(tany_path) %>% `[[`
    d5_meta = bind_rows((neuron_meta %>% rownames_to_column), 
                        (other_meta %>% rownames_to_column),
                        (tany_meta %>% rownames_to_column)) %>%
        filter(!(labels_lvl1 == 'Tanycytes' & is.na(labels_lvl2))) %>%
        filter(!str_detect(labels_chunk, 'drop')) %>%
        right_join(nhgc_d5, by = 'rowname') %>%
        filter(time == 'Day5') %>%
        filter(strain == 'obob')  %>% 
        mutate(labels_lvl2 = case_when(is.na(labels_lvl2) ~ labels_lvl1,
                                       TRUE ~ labels_lvl2))
    d5_cells = d5_meta %>% pull(rowname)
    neuron_obj = qs::qread(neuron_path) %>% subset(cells = d5_cells)
    DefaultAssay(neuron_obj) = 'RNA'
    other_obj = qs::qread(other_path) %>% subset(cells = d5_cells)
    DefaultAssay(other_obj) = 'RNA'
    obj_d5 = merge(neuron_obj, other_obj) %>%
        AddMetaData(d5_meta %>% column_to_rownames)
    obj_d5[['RNA']] = JoinLayers(obj_d5[['RNA']])
    obj_d5
}


filter_down_cells_v02 = function(xenium.obj){
    keep_cells = xenium.obj %>% `[[` %>%
        filter(cell_area < 10000) %>%
        filter(cell_area > 100) %>%
        filter(avg_confidence >= 0.99) %>%
        filter(nCount_Xenium >= 100) %>%
        filter(nCount_Xenium <= 10000) %>%
        rownames
    xenium.obj = subset(xenium.obj, cells = keep_cells)
    xenium.obj
}


process_xenium = function(obj, xenium_genes=NULL){
    DefaultAssay(obj) = 'Xenium'
    
    if (is.null(xenium_genes)){
        xenium_genes = get_xe_genes_all_v02(obj)
    } else {
        xenium_genes = xenium_genes
    }
    

    obj = obj %>% JoinLayers

    obj <- NormalizeData(obj)
    VariableFeatures(obj) = xenium_genes
    obj <- ScaleData(obj, features = xenium_genes)
    obj <- RunPCA(obj, features = xenium_genes, npcs = 80)
    obj <- FindNeighbors(obj, dims = 1:80)
    obj <- FindClusters(obj, resolution = 2, cluster.name = "seurat_clusters")
    obj = obj %>% RunUMAP(assay='Xenium', slot='counts', features = xenium_genes, reduction.name = "umap", return.model=TRUE)
    obj
}


get_high_conf_cells = function(xe_obj,     cutoff_score = .90){
    high_conf_cells = xe_obj %>%
        AddMetaData(xe_obj@misc$ann_by_level$scores$l1) %>%
        `[[` %>% 
        mutate(low_conf_score = case_when((neuron < cutoff_score) & (neuron > 1-cutoff_score) ~ TRUE,
                                          TRUE ~ FALSE)) %>%
        filter(!low_conf_score) %>%
        rownames

    high_conf_cells
}


set_default_assay = function(obj, assay){
    DefaultAssay(obj) = assay
    obj
}

annotate_by_level_v02 = function(obj, 
                             counts_assay='RNA', 
                             graph_name='RNA_snn', 
                             classification_data_path=classification_data_path, 
                             annotation_col=annotation_col){
    #uses conos, cellannotater
    if (is.null(counts_assay)) {
        counts_assay = obj@active.assay
    }
    
    cm <- obj@assays[[counts_assay]] %>% 
#         JoinLayers %>% 
        `$`('counts')
    cm_norm <- Matrix::t(obj@assays[[counts_assay]]$data)
    
    
    if (!is.null(graph_name)){
        if (graph_name == 'active_snn'){
            graph_name = paste0(counts_assay, '_snn')
            graph <- obj@graphs[[graph_name]]
        } else {
            graph <- obj@graphs[[graph_name]]
        }
    } else {
        graph = NULL
    }
    emb <- obj@reductions$umap@cell.embeddings
    clusters <- setNames(obj@meta.data$seurat_clusters, rownames(obj@meta.data))
    clf_data <- getClassificationData(cm, classification_data_path)
    ann_by_level <- assignCellsByScores(graph=graph, clf_data, clusters=clusters)
    Idents(obj) <- ann_by_level$annotation$l1
#     p2 = DimPlot(obj, reduction = "umap", label=T, repel = T)+NoLegend() 
    obj[[annotation_col]] <- obj@active.ident
    obj@misc$clf_data = clf_data
    obj@misc$ann_by_level = ann_by_level
    obj
}


get_xe_genes_all_v02 = function(xe_obj){
    DefaultAssay(xe_obj) = 'Xenium'
    xenium_genes_all = xe_obj %>% 
             Features %>%
             tibble(rowname = .) %>%
             filter(rowname != 'Lmx1a') %>%
             pull(rowname)
    xenium_genes_all
}

get_xe_genes_v02 = function(xe_obj, obj_fgf1){
#     DefaultAssay(obj_fgf1) = 'integrated'
    ref_genes_all= obj_fgf1 %>% 
                Features

        xenium_genes = xe_obj %>% 
                get_xe_genes_all_v02 %>%
                intersect(ref_genes_all) %>%
                tibble %>%
                # filter(. != 'Agrp') %>% #this is not common between both datasets
                pull
        xenium_genes
}


reclass_by_gene_hilo_v02 = function(xe_obj, gene, gene_thr_lo, gene_thr_hi, wrong_class, correct_class){
    meta = xe_obj %>% `[[`
    meta = meta %>%
        mutate(gene_count = xe_obj[["Xenium"]]$counts[gene, ]) %>%
        mutate(cell_class = case_when((gene_count >= gene_thr_hi &
                                      gene_count >= gene_thr_lo &
                                      cell_class == wrong_class) ~ correct_class,
                                      (gene_count <= gene_thr_hi &
                                      gene_count >= gene_thr_lo &
                                      cell_class == wrong_class) ~ 'blah',
                                      TRUE ~ cell_class)) %>%
        select(-gene_count)
    xe_obj@meta.data = meta
    xe_obj
}


relabel_by_gene_hilo_v02 = function(xe_obj, gene, gene_thr_lo, gene_thr_hi, correct_label){
    meta = xe_obj %>% `[[`
    meta = meta %>%
        mutate(gene_count = xe_obj[["Xenium"]]$counts[gene, ]) %>%
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


add_cell_class_and_polar_label = function(obj_fgf1){
    meta = obj_fgf1 %>%
        `[[` %>% 
        mutate(cell_class = class) %>%
        mutate(polar_label = paste0(labels, '.', fgf1_grouping))
    obj_fgf1 = obj_fgf1 %>% AddMetaData(meta)
    obj_fgf1
}



transfer_data_cca_00_v02 = function(xe_obj, obj_fgf1, refdata_column){
    xenium_genes = get_xe_genes_v02(xe_obj, obj_fgf1)
    xenium_genes_all = get_xe_genes_all_v02(xe_obj)
    obj_fgf1 = obj_fgf1 %>% Seurat::RunUMAP(assay='RNA', slot='counts', features = xenium_genes, return.model = TRUE)
    refdata = obj_fgf1 %>% `[[` %>% pull(refdata_column)
    anch = FindTransferAnchors(reference = obj_fgf1, query = xe_obj, features = xenium_genes,  
                               normalization.method = 'LogNormalize',
                               reduction = 'cca',
                               query.assay = 'Xenium',
                               reference.assay = 'integrated',
                               recompute.residuals = T,
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


transfer_data_cca_00_polar_label_v02 = function(xe_obj, obj_fgf1, xenium_genes, selected_label) {
  tryCatch(
    {
      xe_obj = xe_obj %>%
        split_cell_labels(selected_label) %>%
        process_xenium(xenium_genes=xenium_genes)
      obj_fgf1 = obj_fgf1 %>%
        split_cell_labels(selected_label) %>%
        milo_prep_and_proc
      xe_obj = transfer_data_cca_00_v02(xe_obj, obj_fgf1, "polar_label")
      xe_obj
    },
    error = function(e) {
      message("An error occurred: ", e$message)
      return(NA)
    }
  )
}


transfer_data_cca_00_unimodal_v02 = function(xe_obj, obj_fgf1, refdata_column){
    xenium_genes = get_xe_genes_v02(xe_obj, obj_fgf1)
    xenium_genes_all = get_xe_genes_all_v02(xe_obj)
    obj_fgf1 = obj_fgf1 %>% Seurat::RunUMAP(assay='RNA', slot='counts', features = xenium_genes, return.model = TRUE)
    refdata = obj_fgf1 %>% `[[` %>% pull(refdata_column)
    # anch <- FindTransferAnchors(reference = obj_fgf1,
    #                             query = xe_obj,
    #                             features = xenium_genes,
    #                             normalization.method = "SCT",
    #                             reduction='cca',
    #                             recompute.residuals = T,
    #                             dims=1:20)
    anch = FindTransferAnchors(reference = obj_fgf1, query = xe_obj, features = xenium_genes,  
                               normalization.method = 'LogNormalize',
                               reduction = 'cca',
                               query.assay = 'Xenium',
                               reference.assay = 'integrated',
                               recompute.residuals = T,
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


transfer_data_cca_00_polar_label_unimodal_v02 = function(xe_obj, obj_fgf1, xenium_genes, selected_label) {
  tryCatch(
    {
      xe_obj = xe_obj %>%
        split_cell_labels(selected_label) %>%
        process_xenium(xenium_genes=xenium_genes)
      obj_fgf1 = obj_fgf1 %>%
        split_cell_labels(selected_label) %>%
        milo_prep_and_proc
      xe_obj = transfer_data_cca_00_unimodal_v02(xe_obj, obj_fgf1, "polar_label")
      xe_obj
    },
    error = function(e) {
      message("An error occurred: ", e$message)
      return(NA)
    }
  )
}


make_RCTD_ref = function(obj, ident){
    Idents(obj) <- ident
    allowed_idents = obj %>% `[[` %>% group_by(.data[[ident]]) %>% summarise(n=n()) %>% filter(n >=25) %>% select(-n) %>% pull
    selected_cells = obj %>% `[[` %>% filter(.data[[ident]] %in% allowed_idents) %>% rownames
    obj = obj %>% subset(cells = selected_cells )
    counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
    cluster <- obj[[ident]] %>% rownames_to_column %>% deframe %>% as.factor
    nUMI <- obj[['nCount_RNA']] %>% rownames_to_column %>% deframe
    nUMI <- colSums(counts)
    reference <- Reference(counts, cluster, nUMI)
    reference
}


make_RCTD_query = function(xe_obj, fov){
    query.counts <- GetAssayData(xe_obj, assay = "Xenium", slot = "counts")[, Cells(xe_obj[[fov]])]
    coords <- GetTissueCoordinates(xe_obj[[fov]], which = "centroids")
    rownames(coords) <- coords$cell
    coords$cell <- NULL
    query <- SpatialRNA(coords, query.counts, colSums(query.counts))
    query
}


get_rctd_annotations_df = function(rctd){
    rctd %>% `@`('results') %>% `$`('results_df')
}


label_transfer_to_fov = function(xe_obj, rctd_ref, fov, confidence_threshold = 5, doublet_threshold=20){
    fov_annot = xe_obj %>% 
        make_RCTD_query(fov) %>%
        create.RCTD(rctd_ref, max_cores = 1, CONFIDENCE_THRESHOLD = confidence_threshold, DOUBLET_THRESHOLD=doublet_threshold) %>%
        run.RCTD(doublet_mode = "doublet") %>%
        get_rctd_annotations_df
    fov_annot
}


run_rctd = function(xe_obj, ref_obj, transfer_col, confidence_threshold = 5, doublet_threshold=20, return_obj=FALSE, misc_slot_name='rctd'){
    rctd_ref = make_RCTD_ref(ref_obj, transfer_col)
    fov.1_annot = label_transfer_to_fov(xe_obj, rctd_ref, 'fov', confidence_threshold = confidence_threshold, doublet_threshold=doublet_threshold)
    fov.2_annot = label_transfer_to_fov(xe_obj, rctd_ref, 'fov.2', confidence_threshold = confidence_threshold, doublet_threshold=doublet_threshold)
    fov.3_annot = label_transfer_to_fov(xe_obj, rctd_ref, 'fov.3', confidence_threshold = confidence_threshold, doublet_threshold=doublet_threshold)
    fov.4_annot = label_transfer_to_fov(xe_obj, rctd_ref, 'fov.4', confidence_threshold = confidence_threshold, doublet_threshold=doublet_threshold)
    fov.5_annot = label_transfer_to_fov(xe_obj, rctd_ref, 'fov.5', confidence_threshold = confidence_threshold, doublet_threshold=doublet_threshold)
    fov.6_annot = label_transfer_to_fov(xe_obj, rctd_ref, 'fov.6', confidence_threshold = confidence_threshold, doublet_threshold=doublet_threshold)
    fov.7_annot = label_transfer_to_fov(xe_obj, rctd_ref, 'fov.7', confidence_threshold = confidence_threshold, doublet_threshold=doublet_threshold)
    fov.8_annot = label_transfer_to_fov(xe_obj, rctd_ref, 'fov.8', confidence_threshold = confidence_threshold, doublet_threshold=doublet_threshold)
    annot_df = bind_rows(fov.1_annot,
                         fov.2_annot,
                         fov.3_annot,
                         fov.4_annot,
                         fov.5_annot,
                         fov.6_annot,
                         fov.7_annot,
                         fov.8_annot)
    if (return_obj){
        xe_obj = xe_obj@misc[[misc_slot_name]] = annot_df
        xe_obj
    } else {
        annot_df
    }
}


store_in_misc = function(obj, slot, payload){
    Misc(object = obj, slot = slot) <- payload
    obj
}



milo_prep_and_proc = function(object){
object %>%
#  set_labels_to_lvl1 %>%
#  prep_obj_for_milo_cb_v01 %>%
#  set_batch_to_lane %>% # do not set batch to lane for cluster splits, will error on design or model matrix. reset later
#  prep_obj_for_milo_cb_v01(set_orig.batch = FALSE) %>%
# #  subset_exp_by_time(day) %>%
# #  single_split(cluster) %>%
 reconsitute_rna_seurat %>%
#              process_seurat(method = "integrate", 
#                             batch ="batch", # batch to batch
#                             dims = 30, res = 0.8,
#                             k.weight = 40)
 run_until_success(function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                  batch = "batch",
                                  dims = 15, res = 0.8,
                                  k.weight = 100, k.anchor=25),
                   function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                  batch = "batch",
                                  dims = 15, res = 0.8,
                                  k.weight = 50, k.anchor=25),
                   function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                  batch = "batch",
                                  dims = 15, res = 0.8,
                                  k.weight = 20, k.anchor=25),
                   function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                  batch = "batch",
                                  dims = 15, res = 0.8,
                                  k.weight = 10, k.anchor=25),
                   function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                  batch = "orig.batch",
                                  dims = 15, res = 0.8,
                                  k.weight = 100, k.anchor=25),
                   function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                  batch = "orig.batch",
                                  dims = 15, res = 0.8,
                                  k.weight = 50, k.anchor=25),
                   function() subset(eval(.)) %>% process_seurat(method = "integrate", 
                                  batch = "orig.batch",
                                  dims = 15, res = 0.8,
                                  k.weight = 20, k.anchor=25),
                   function() subset(eval(.)) %>% Seurat::SCTransform(assay='RNA',
                                   method="glmGamPoi",
                                   vars.to.regress= 'batch',
                                   vst.flavor="v2",
                                   verbose=TRUE) %>% run_sct_chaser,
                   function() subset(eval(.)) %>% Seurat::SCTransform(assay='RNA',
                                   method="glmGamPoi",
                                   vars.to.regress= 'orig.batch',
                                   vst.flavor="v2",
                                   verbose=TRUE) %>% run_sct_chaser
                  ) # %>% #reset batch for downstream applications
#  reset_orig.batch %>% 
#  prep_obj_for_milo_cb_v01
}