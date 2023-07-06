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

process_seurat <- function(obj, method, res=NULL, features=NULL, dims=NULL, batch=NULL, return_model = F, cluster=T, type = "seur", nfeats = 3000) {
  
  if(type == "sce") {
    obj <- Seurat::CreateSeuratObject(counts = counts(obj))
  } else {
    obj <- obj
  }
  
  if(!is.null(features)) {
    features <- rownames(obj)[!grepl(features, rownames(obj))]
    obj <- subset(obj, features = features)
  } 
  
  
  if(method == "log") {
    obj <-
      Seurat::NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000) %>%
      Seurat::FindVariableFeatures(., selection.method = "vst", nfeatures = nfeats) %>%
      Seurat::ScaleData(.) %>%
      Seurat::RunPCA(.)
  } else if (method == "glm") {
    #plan("sequential")
    obj <-
      Seurat::SCTransform(obj, method = "glmGamPoi", batch_var=batch, variable.features.n = nfeats) %>%
      Seurat::RunPCA(.)
  } else if (method == "qpoisson") {
    #plan("sequential")
    obj <-
      Seurat::SCTransform(obj, method = "qpoisson", variable.features.n = nfeats) %>%
      Seurat::RunPCA(.)
  }
  
  if(cluster==T){
    obj <-
      obj %>%
      Seurat::RunUMAP(., dims = seq(dims), return.model=return_model) %>%
      Seurat::FindNeighbors(., dims = seq(dims)) %>%
      Seurat::FindClusters(., resolution = res)
  }
  
  return(obj)
}
