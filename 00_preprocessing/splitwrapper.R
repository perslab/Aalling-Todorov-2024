##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param neur
##' @param glia
##' @return
##' @author dylanmr
##' @export
splitwrapper <- function(seur, split.by = "predicted.id") {
  
  future::plan(future::sequential)
  Seurat::DefaultAssay(seur) <- "SCT"
  split <- Seurat::SplitObject(seur, split.by = split.by)
  return(split)

}
