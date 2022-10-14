library(tidyverse)
library(Seurat)
library(zellkonverter)
library(basilisk)
library(scRNAseq)

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


make_basilisk_env = function(){
    bask_env = zellkonverter::zellkonverterAnnDataEnv()
    bask_env
}


run_basilisk_convert = function(sce, path, basilisk_env){
    roundtrip = basiliskRun(fun = function(sce) {
        # Convert SCE to AnnData:
    #      adata <- SCE2AnnData(sce)

        # Maybe do some work in Python on 'adata':
        # BLAH BLAH BLAH

        # Convert back to an SCE:
    #     writeH5AD(sce, '/scratch/nmq407/exp_04_labelled_h5ad')
        writeH5AD(sce, path)
    }, env = basilisk_env, sce = sce)
    roundtrip
} 


get_sct_variable_features = function(obj, nfeats){
    obj <- Seurat::SCTransform(obj ,
                        assay='RNA',
                        method="glmGamPoi",
                        vars.to.regress="batch",
                        vst.flavor="v2",
                        variable.features.n = nfeats,          
                        verbose=TRUE)
    var_feats = obj@assays$SCT@var.features
    var_feats
}