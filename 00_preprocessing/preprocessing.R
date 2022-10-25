library(tidyverse)
library(SingleCellExperiment)
library(Seurat)
library(scRNAseq)

PROJECT_DIR = Sys.getenv('PROJECT_DIR')

read_sc = function(path_to_rds){
    exp = readRDS(path_to_rds)
    exp
}

open_obj = function(fp){
    obj = readRDS(fp)
    obj
}

save_obj = function(obj, fn){
    saveRDS(obj, fn)
    fn
}

make_exp_list = function(...){
    exp_list = tibble::lst(...)
    exp_list
}

select_integration_featuers = function(exp_list){
    exp_list = Seurat::SelectIntegrationFeatures(object.list = exp_list, nfeatures = 5000)
    exp_list
}

prep_sc_transform = function(exp_list, features){
    exp_list = Seurat::PrepSCTIntegration(object.list = exp_list, anchor.features = features)
}

find_integration_anchors = function(exp_list, features){
    anchors = Seurat::FindIntegrationAnchors(object.list = exp_list,
                                 normalization.method = 'SCT',
                                 reduction = 'rpca',
                                 anchor.features = features)
    anchors
}

integrate_data = function(anchors){
    exp = Seurat::IntegrateData(anchorset = anchors, normalization.method = 'SCT')
    exp
}

run_pca = function(exp){
    exp = Seurat::RunPCA(exp, verbose = FALSE)
    exp
}

run_umap = function(exp){
    exp = Seurat::RunUMAP(exp, reduction = "pca", dims = 1:30)
    exp
}

read_metadata = function(path_to_meta){
    meta = readr::read_csv(path_to_meta)
    meta
}

add_meta_to_exp = function(exp, meta){
    meta = meta %>% mutate_all(factor)
    rn = rownames(exp@meta.data)
    meta = dplyr::left_join(exp@meta.data, meta, by=c("hash.mcl.ID" = "sample_name"))
    rownames(meta) = rn
    exp@meta.data = meta
    exp
}

load_campbell = function(){
    campbell <- scRNAseq::CampbellBrainData()
    campbell <- Seurat::as.Seurat(campbell, counts = "counts", data="counts")
    campbell
}

sc_transform_campbell = function(campbell){
    campbell <- Seurat::SCTransform(campbell,
                        assay='originalexp',
                        method="glmGamPoi",
                        vars.to.regress="batches",
                        verbose=TRUE)
    campbell
}

run_find_neighbors = function(exp){
    exp = Seurat::FindNeighbors(exp, dims = 1:30, verbose = FALSE)
    exp
}

run_find_clusters = function(exp, resolution=0.8){
    exp = Seurat::FindClusters(exp, verbose = FALSE, resolution=resolution)
    exp@meta.data$seurat_clusters = exp$seurat_clusters
    exp
}

run_sct_chaser = function(exp, resolution=0.8){
    exp = run_pca(exp)
    exp = run_umap(exp)
    exp = run_find_neighbors(exp)
    exp = run_find_clusters(exp, resolution=resolution)
    exp
}

sc_transform_fgf1 = function(fgf1){
    fgf1 <- Seurat::SCTransform(fgf1,
                        assay='RNA',
                        method="glmGamPoi",
                        vars.to.regress="batch",
                        vst.flavor="v2",
                        verbose=TRUE)
    fgf1 = run_sct_chaser(fgf1, resolution=0.8)
    fgf1
}

downsample_seurat_obj = function(obj, size){
    obj <- obj[, sample(colnames(obj), size = size, replace=F)]
}

find_anchors_campbell = function(exp, campbell){
    transfer_anchors <- Seurat::FindTransferAnchors(reference = campbell,
                                            query = exp,
                                            dims = 1:30,
                                            normalization.method="SCT",
                                            recompute.residuals=FALSE)
    transfer_anchors
}

transfer_campbell_labels = function(exp, campbell, transfer_anchors){
    for (cluster_col in c('clust_all', 'clust_all_neurons', 'clust_all_micro', 'clust_neurons')) {
        predictions <- Seurat::TransferData(anchorset = transfer_anchors,
                                refdata = dplyr::pull(campbell@meta.data, cluster_col),
                                dims = 1:30,
                                k.weight = 30)
        colnames(predictions)[colnames(predictions) == 'predicted.id'] <- cluster_col
        predictions_f <- predictions[c(cluster_col)]
        exp = Seurat::AddMetaData(exp, metadata = predictions_f)
    }
    exp
}


subset_exp = function(exp, class="all"){
    if (class == 'all') {
        exp = exp
    } else if (class=='neuron'){
        exp_neuron = subset(x = exp, subset = predicted.id == 'neuron')
        exp_miss = subset(x = exp, subset = predicted.id == 'miss')
        exp = merge(exp_neuron, exp_miss)
        exp
    } else if (class=='other'){
        exp = subset(x = exp, subset = predicted.id != 'neuron')
    }
    exp
}

subset_exp_by_strain = function(exp, strain_id){
    exp = subset(x = exp,
                 subset = strain == strain_id)
    exp
}

do_neuron_cluster_surgery = function(neuron){
    md = neuron %>% 
        `[[` %>% 
        rownames_to_column("rownames") %>%
        rowwise() %>% 
        mutate(labels = ifelse(labels %in% c("Sst_Pthlh"), paste0(labels, '.sc', seurat_clusters), paste0(labels))) %>%
        mutate(labels = ifelse(labels %in% c("Unassigned1"), paste0(labels, '.sc', seurat_clusters), paste0(labels))) %>% 
        mutate(labels = ifelse(labels %in% c("Unassigned1.sc28", "Unassigned1.sc17"), "Unassigned1.sc17.sc28", paste0(labels))) %>%
        mutate(labels = ifelse(labels %in% c("Unassigned1.sc49"), "Unassigned2", paste0(labels))) %>%  # this looks like its split off from sc29, an Unassigned2 clust
        column_to_rownames("rownames") %>%
        select(labels)
    neuron = AddMetaData(neuron, metadata = md)
    neuron
}


.get_only_named_genes = function(exp){
    not_a_rik = exp@assays$RNA@meta.features %>% rownames %>% stringr::str_ends("Rik") %>% `!`
    not_an_ensmusg = exp@assays$RNA@meta.features %>% rownames %>% stringr::str_starts("ENSMUSG") %>% `!`
    not_a_Gm = exp@assays$RNA@meta.features %>% rownames %>% stringr::str_starts("Gm\\d") %>% `!`
    not_4_or_more_digits = exp@assays$RNA@meta.features %>% rownames %>% stringr::str_ends("\\d{4,}") %>% `!`
    only_named_genes = (not_a_rik & not_an_ensmusg & not_a_Gm & not_4_or_more_digits) %>% which
    only_named_genes
}

get_only_named_genes_vec = function(vec_of_genes){
    not_a_rik = vec_of_genes %>% stringr::str_ends("Rik") %>% `!`
    not_an_ensmusg = vec_of_genes %>% stringr::str_starts("ENSMUSG") %>% `!`
    not_a_Gm = vec_of_genes %>% stringr::str_starts("Gm\\d") %>% `!`
    not_4_or_more_digits = vec_of_genes %>% stringr::str_ends("\\d{4,}") %>% `!`
    only_named_genes = (not_a_rik & not_an_ensmusg & not_a_Gm & not_4_or_more_digits) %>% which
    vec_of_genes[only_named_genes]
}


sce_retain_ngo_genes = function(sce){
    sce_genes = rownames(logcounts(sce))
    only_named_genes_vec = get_only_named_genes_vec(sce_genes)
    sce_filtered = sce[only_named_genes_vec,]
    sce_filtered
}


make_new_seurat_ngo = function(exp){
    counts_exp = GetAssayData(exp, slot="counts", assay="RNA") 
    counts_exp = counts_exp[.get_only_named_genes(exp),]
    new_exp = CreateSeuratObject(counts=counts_exp, meta.data = exp@meta.data)
    new_exp
}

merge_neuron_and_other = function(exp_neuron, exp_other){
    exp = merge(exp_neuron, y=exp_other)
    exp
}

subset_seurat_obj_genes = function(exp, genes){
    counts_exp = GetAssayData(exp, slot="counts", assay="RNA") 
    counts_exp = counts_exp[genes,]
    new_exp = CreateSeuratObject(counts=counts_exp, meta.data = exp@meta.data)
    new_exp
}

sample_n_cells_from_obj = function(obj, n_cells){
    subsampled_obj = obj[, sample(colnames(obj), size = n_cells, replace=F)]
    subsampled_obj
}


subset_obob5v5_from_obj = function(obj){
    cell_ids = obj@meta.data %>% filter(strain == 'obob' & time == 'Day5') %>% rownames
    obj = subset(obj, cells=cell_ids)
    obj
}


