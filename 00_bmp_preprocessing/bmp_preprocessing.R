.munge_bmp_meta = function(bmp_meta){
    colnames(bmp_meta) = bmp_meta %>% colnames %>% str_replace(' ', '_') %>% str_replace('/', '_') %>% str_to_lower
    bmp_meta = bmp_meta %>% select(sample_name, condition_1,  collection_date)  %>%
        mutate(sample_name = factor(sample_name)) %>%
        mutate(sample_name = paste0("B", sample_name)) %>%
        mutate(condition_1 = factor(condition_1)) %>%
        mutate(collection_date = factor(collection_date)) %>%
        rename(hash.mcl.ID = sample_name) %>%
        rename(treatment = condition_1) %>%
        mutate(batch = 'Batch 4') %>%
        mutate(strain = 'obob') %>%
        mutate(time = "Day5")
    bmp_meta
}

add_meta_to_bmp_obj = function(bmp_obj, bmp_meta){
    bmp_meta = .munge_bmp_meta(bmp_meta)
    md = bmp_obj %>% `[[` %>% 
        rownames_to_column %>%  
        mutate(hash.mcl.ID = paste0("B", hash.mcl.ID)) %>%
        left_join(bmp_meta, by = "hash.mcl.ID") %>% 
        column_to_rownames
    bmp_obj = AddMetaData(bmp_obj, metadata=md)
    bmp_obj
}


bmp_combine_astrocytes = function(obj){
    md = obj %>% 
        `[[` %>% 
        rownames_to_column("rownames") %>%
        rowwise() %>% 
        mutate(seurat_clusters = ifelse(seurat_clusters %in% c(5, 0, 1), 5, seurat_clusters)) %>%
        mutate(seurat_clusters = factor(seurat_clusters)) %>%
        column_to_rownames("rownames") %>%
        select(seurat_clusters)

    obj = AddMetaData(obj, metadata = md)

    obj@meta.data = obj@meta.data %>%
        mutate(SCT_snn_res.0.8	= seurat_clusters)
    obj
}


add_umap_model = function(exp){
    exp <- RunUMAP(exp, dims = 1:30, reduction = "pca", return.model = TRUE)
    exp
}


find_fgf1_anchors = function(exp, fgf1, recompute_residuals=FALSE){
    future::plan(future::sequential)
    transfer_anchors <- Seurat::FindTransferAnchors(reference = fgf1,
                                            query = exp,
                                            dims = 1:30,
                                            reference.reduction = "pca",
                                            normalization.method="SCT",
                                            recompute.residuals=recompute_residuals)
    transfer_anchors
}


map_bmp1_to_fgf1_ref = function(bmp1_obj, bmp1_anchors, fgf1_ref){
    future::plan(future::sequential)
    bmp1_obj = MapQuery(anchorset = bmp1_anchors, reference = fgf1_ref, query = bmp1_obj,
                        refdata = list(labels = "labels"), reference.reduction = "pca", reduction.model = "umap")
    bmp1_obj@meta.data = bmp1_obj@meta.data %>%
        mutate(labels = predicted.labels) %>%
        select(-predicted.labels)
    bmp1_obj
}


merge_bmp1_mp_fgf1_ref = function(bmp1_obj, fgf1_ref){
    emb = Embeddings(bmp1_obj, reduction = "ref.umap")
    dro = CreateDimReducObject(embeddings = emb, key = "UMAP_", assay = DefaultAssay(bmp1_obj))
    dro_2 = fgf1_ref[['umap']]
    umap_merged = merge(dro, dro_2)
    emb = Embeddings(bmp1_obj, reduction = "ref.pca")
    dro = CreateDimReducObject(embeddings = emb, key = "PC_", assay = DefaultAssay(bmp1_obj))
    dro_2 = fgf1_ref[['pca']]
    pca_merged = merge(dro, dro_2)
    
    obj_merged = merge(bmp1_obj, y=fgf1_ref) %>%
        sc_transform_fgf1

    obj_merged[['mergedumap']] = umap_merged
    obj_merged[['mergedpca']] = pca_merged

    md = obj_merged %>%
        `[[` %>%
        select(all_of(c("treatment", "time", "strain"))) %>%
        mutate(group = interaction(.) %>% fct_drop())
    
    obj_merged = AddMetaData(obj_merged, md)

    obj_merged
}


project_bmp1_into_fgf1_unimod = function(both_obj){
    fgf1_ref = subset(x = both_obj, subset = batch != 'Batch 4') %>% sc_transform_fgf1
    bmp1_obj = subset(x = both_obj, subset = batch == 'Batch 4') %>% sc_transform_fgf1_nobatch

    fgf1_ref = add_umap_model(fgf1_ref)
    bmp1_transfer_anchors = find_fgf1_anchors(bmp1_obj, fgf1_ref)
    bmp1_obj = map_bmp1_to_fgf1_ref(bmp1_obj, bmp1_transfer_anchors, fgf1_ref)
    obj = merge_bmp1_mp_fgf1_ref(bmp1_obj, fgf1_ref)
    obj
}


set_unimod_embedding = function(obj){
    unimod_umap = Embeddings(obj, reduction = "mergedumap") %>%
        CreateDimReducObject(embeddings = ., key = "UMAP_", assay = DefaultAssay(obj))
    sc_umap = Embeddings(obj, reduction = "umap") %>%
        CreateDimReducObject(embeddings = ., key = "UMAP_", assay = DefaultAssay(obj))
    obj[['umap']] = unimod_umap
    obj[['sctumap']] = sc_umap
    unimod_pca = Embeddings(obj, reduction = "mergedpca") %>%
        CreateDimReducObject(embeddings = ., key = "PC_", assay = DefaultAssay(obj))
    sc_pca = Embeddings(obj, reduction = "pca") %>%
        CreateDimReducObject(embeddings = ., key = "PC_", assay = DefaultAssay(obj))
    obj[['pca']] = unimod_pca
    obj[['sctpca']] = sc_pca
    obj
}

