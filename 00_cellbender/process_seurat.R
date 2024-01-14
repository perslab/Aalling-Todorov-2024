#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param .
#' @param method
#' @param res
#' @param features
#' @param dims
#' @return
#' @author dylanmr
#' @export


process_seurat <- function(obj, method, ref_datasets = NULL, k.anchor=5,k.weight=100,
                           res=NULL, features=NULL, dims=NULL, 
                           batch=NULL, return_model = F, cluster=T, 
                           type = "seur", nfeats = 5000, neighbor=F) {
  
  if(type == "sce") {
    obj <- CreateSeuratObject(counts = counts(obj))
  } else {
    obj <- obj
  }

  
  if(!is.null(features)) {
    features <- rownames(obj)[!grepl(features, rownames(obj))]
    obj <- subset(obj, features = features)
  } 
  
  if(method == "integrate") {
    obj <- .integrate_seurat(obj, split = batch, nfeats = nfeats, ref_datasets=ref_datasets, k.anchor = k.anchor, k.weight = k.weight)
    DefaultAssay(obj) <- "integrated"
    obj <- ScaleData(obj) %>% RunPCA(.)
  } else if (method == "log") {
    
    if (inherits(obj, "list")) {
      obj <- obj[[1]]
    }
    
    DefaultAssay(obj) <- "RNA"
    obj <-
      NormalizeData(obj) %>%
      FindVariableFeatures(., selection.method = "vst", nfeatures = nfeats) %>%
      ScaleData(., vars.to.regress=batch) %>%
      RunPCA(.)
  } else if (method == "glm") {
    obj <-
      SCTransform(obj, method = "glmGamPoi", batch_var=batch, variable.features.n = nfeats) %>%
      RunPCA(.)
  } else if (method == "qpoisson") {
    obj <-
      SCTransform(obj, method = "qpoisson", variable.features.n = nfeats, vars.to.regress = batch) %>%
      RunPCA(.)
  }
  
  if(cluster==T & method == "integrate"){
    obj <-
      obj %>%
      RunUMAP(., dims = seq(dims), return.model=return_model) %>%
      FindNeighbors(., dims = seq(dims), return.neighbor=neighbor) %>%
      FindClusters(., resolution = res)
  } else if(cluster==T & method != "integrate") {
    obj <-
      obj %>%
      RunUMAP(., dims = seq(dims), return.model=return_model) %>%
      FindNeighbors(., dims = seq(dims), return.neighbor=neighbor) %>%
      FindClusters(., resolution = res)
  }
  
  return(obj)
}

.integrate_seurat <- function(obj, split, nfeats, ref_datasets, k.anchor, k.weight) {
  
  if(is.list(obj)) {
    list <- obj
  } else {
    DefaultAssay(obj) <- "RNA"
    list <- SplitObject(obj, split.by = split)
  }
  
  list <- lapply(X = list, FUN = function(x) {
    if(sum("SCT" %in% names(x@assays))>0) {
      x[["SCT"]] <- NULL
    }
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeats)
  })
  
  features <- SelectIntegrationFeatures(object.list = list)
  
  list <- lapply(X = list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  
  anchors <- FindIntegrationAnchors(object.list = list, reference = ref_datasets, k.anchor = k.anchor, anchor.features = features, reduction = "rpca")
  
  integrated <- IntegrateData(anchorset = anchors, k.weight=k.weight)
  
  DefaultAssay(integrated) <- "integrated"
  integrated[["RNA"]] = JoinLayers(JoinLayers(integrated@assays$RNA))
  return(integrated)
  
}


reconsitute_rna_seurat = function(obj){
    mat = obj@assays$RNA$counts
    meta = obj %>% `[[`
    obj = CreateSeuratObject(counts = mat,
                             meta.data = meta,
                             names.field = 1, names.delim = "_",
                             min.features = 500, min.cells = 10)
    obj
}


set_empty_labels_lvl2_to_lvl1 = function(obj){
    meta = obj %>% `[[` %>%
        mutate(labels_lvl2 = case_when(is.na(labels_lvl2) ~ labels_lvl1,
                                       TRUE ~ labels_lvl2,))
    obj = obj %>%
        AddMetaData(meta)
    obj
}

