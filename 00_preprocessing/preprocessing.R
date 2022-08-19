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

sc_transform_fgf1 = function(fgf1){
    fgf1 <- Seurat::SCTransform(fgf1,
                        assay='RNA',
                        method="glmGamPoi",
                        vars.to.regress="comments",
                        verbose=TRUE)
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


subset_neurons = function(exp){
    exp = subset(x = exp, subset = clust_all_neurons == 'neuron')
    exp
}

subset_other = function(exp){
    exp = subset(x = exp, subset = clust_all_neurons != 'neuron')
    exp
}