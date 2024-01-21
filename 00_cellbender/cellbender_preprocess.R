make_mat = function(path_to_cb){
    mat = Read_CellBender_h5_Multi_File(data_dir = path_to_cb, 
                                        custom_name = '.h5', 
                                        merge = TRUE, 
                                        parallel = TRUE, 
                                        num_cores = 20)
    mat
}


fix_mat_colnames = function(mat){
    new_colnames = mat %>%
        colnames %>%
        str_replace('_out_', '_') %>% 
        str_replace('-1', '')
    colnames(mat) = new_colnames
    mat
}


fix_mat_colnames_ranger = function(mat){
    new_colnames = mat %>%
        colnames %>%
        str_replace('_out_', '_') %>% 
        str_replace('-1', '') %>%
        str_replace('_raw_', '_')
    colnames(mat) = new_colnames
    mat
}


make_obj_from_mat = function(mat){
    obj = CreateSeuratObject(counts = mat, names.field = 1, names.delim = "_",
                             min.features = 500, min.cells = 10)
    good_cells = obj %>%
        `[[` %>%
        filter(nCount_RNA > 1100) %>%
        rownames
    obj = obj %>% subset(cells = good_cells)
    obj
}


read_meta = function(file_name){
    seurat_obj = readRDS(file_name)
    meta = seurat_obj@meta.data
    meta
}

make_metadata = function(obj, path_to_cb, path_to_sample_meta){
    obj_metadata = obj %>% 
        `[[` %>% 
        mutate(rowname = rownames(.)) %>%
        separate(rowname, into = c("Index.10x", "barcode"), sep = "_", remove=FALSE) %>%
        select(-orig.ident)

    #get meta from all seurats from SCOP
    # List all .rds files in the folder
    file_names <- list.files(path = path_to_cb, pattern = "\\.rds$", full.names = TRUE)
    # Function to extract the name from the file name
    extract_name <- function(file_name) {
      # Extract the part of the file name right before "_seurat"
      sub(".*_([[:alnum:]]+)_seurat\\.rds", "\\1", basename(file_name))
    }
    # Read all RDS files and store them in a list with extracted names
    combined_metadata <- setNames(lapply(file_names, read_meta), sapply(file_names, extract_name)) %>%
        bind_rows %>%
        mutate(barcode = rownames(.)) %>%
        select(-orig.ident) %>%
        mutate(rowname = paste0(Index.10x, '_', barcode)) %>%
        tibble %>% data.frame %>%
        column_to_rownames %>%
        rename_with(~ paste0(., "_SCOP")) %>%
        rownames_to_column
    
    # join with metadata from the cellbender obj
    combined_metadata = combined_metadata %>%
        left_join(obj_metadata, by='rowname')

    # bring in meta for each sample
    more_meta = read_csv(path_to_sample_meta) %>%
        rename(hash.mcl.ID_SCOP = sample_name)
    # add to the combined metadata
    combined_metadata = combined_metadata %>% 
        left_join(more_meta, by = 'hash.mcl.ID_SCOP')
    # restore rownames as rownames
    combined_metadata = combined_metadata %>%
        column_to_rownames
    #return
    combined_metadata
}

filter_down_cells = function(obj){
    selected_cells = obj %>%
        `[[` %>%
        filter(HTO_mcl_classification.global_SCOP == 'Singlet') %>%
#         filter(HTO_mcl_margin_SCOP > 0.1) %>%
        rownames
    obj = obj %>% subset(cells = selected_cells)
    obj
}


annotate_by_level = function(obj, 
                             counts_assay='RNA', 
                             graph_name='RNA_snn', 
                             classification_data_path=classification_data_path, 
                             annotation_col=classification_data_path){
    #uses conos, cellannotater
    if (is.null(counts_assay)) {
        counts_assay = obj@active.assay
    }
    
    cm <- obj@assays[[counts_assay]]$counts
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
    obj
}


prepend_obj_labels = function(obj, label_col, prefix){
meta = obj %>%
    `[[`
meta$labels = meta[[label_col]]
meta = meta %>%
    rownames_to_column %>%
    rowwise %>%
    mutate(labels = paste0(prefix, labels, collapse = '_')) %>%
    ungroup %>%
    column_to_rownames
obj = obj %>% AddMetaData(meta)
obj
}


rpca_sct_integrate = function(obj=obj,
                              batch=NULL,
                              nfeats=3000,
                              dims=30,
                              k.anchor=5,
                              k.weight=100){
    obj.list <- SplitObject(obj, split.by = batch)
    obj.list <- lapply(X = obj.list, FUN = SCTransform, method = "glmGamPoi")
    features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = nfeats)
    obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)
    obj.list <- lapply(X = obj.list, FUN = RunPCA, features = features)
    
    obj.anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT",
    anchor.features = features, dims = seq(dims), reduction = "rpca", k.anchor = k.anchor)
    
    obj.combined.sct <- IntegrateData(anchorset = obj.anchors, normalization.method = "SCT", dims = seq(dims), k.weight = k.weight)
    obj.combined.sct <- RunPCA(obj.combined.sct, verbose = FALSE)
    obj.combined.sct <- RunUMAP(obj.combined.sct, reduction = "pca", dims = seq(dims))
    obj
}


transfer_labels_via_map_query = function(obj = obj, 
                                         dims = 30, 
                                         ref_path = ref_path){
    ref = qs::qread(ref_path)
    ref <- subset(ref, cells = sample(Cells(ref), 30000))
    anchors <- FindTransferAnchors(
      reference = ref, query = obj,
      reference.reduction = 'pca',
      dims = 1:30,
      normalization.method = "LogNormalize",
      reference.assay = "integrated",
      query.assay = "RNA", 
      k.filter = 300, max.features = 100
      )
    query <- MapQuery(
        reference = ref,
        query = obj,
        anchorset = anchors,
        refdata = list(celltype = "labels_202310"),
        reference.reduction = "pca",
        reduction.model = "umap"
      )
    query_labels = query %>% 
    `[[`
    obj = obj %>% AddMetaData(query_labels)
    obj
}


map_camp_seurat_cluster = function(query, ref_path=''){
  dims=30
  ref_obj <- qs::qread(ref_path)
  ref_sub <- subset(ref_obj, cells = sample(Cells(ref_obj), 20000))
  ref_sub$labels_202310 <- ifelse(grepl("Agrp",ref_sub$labels_202310), "Agrp", ref_sub$labels_202310)
  anchors <- FindTransferAnchors(reference = ref_sub, query = query, dims = seq(dims), reference.reduction = "pca")
  predictions <- TransferData(anchorset = anchors, refdata = ref_sub$labels_202310, dims = seq(dims))
  Misc(object = query, slot = paste0("predictions_", "labels_202310ß")) <- predictions
  predictions = predictions %>% dplyr::select(predicted.id, prediction.score.max)

  query <- AddMetaData(query, metadata = predictions)
  #add predictions object to dataset
  
  
  # get proportion of predicted.ids that map to each cluster
  clustering_res_name = .DollarNames(query) %>% #default is SCT_snn_res.0.8
    str_detect('res.') %>% 
    which %>% 
    .DollarNames(query)[.] %>% 
    .[order(., decreasing=TRUE)] %>% 
    `[[`(1)# take the highest resolution
  cluster_names <- 
    prop.table(table(query@meta.data[[clustering_res_name]], query$predicted.id), margin = 1) %>% 
    data.frame()

  # if more than 25% of cells in a cluster map to 2 predicted.ids, label cluster as a combination
  cluster_split <- 
    cluster_names %>% 
    group_by(Var1) %>% 
    slice_max(n=2, order_by = Freq) %>% 
    group_split() %>% 
    map_chr(~if_else(.x$Freq[2]<0.25, 
                     as.character(.x$Var2[1]), 
                     paste0( as.character(.x$Var2[1]))))
  
  # add labels to object  
  labels <- unique(cluster_names$Var1) %>%  as.character()
  names(labels) <- cluster_split
  query$labels <- fct_recode(query@meta.data[[clustering_res_name]], !!!labels)
  
  return(query)
#       predictions <- TransferData(anchorset = anchors, refdata = ref_sub$labels_202310, dims = seq(dims)) 
#     predictions
}




