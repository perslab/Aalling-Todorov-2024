#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param obj
#' @param type
#' @param dims
#' @return
#' @author dylanmr
#' @export
map_ref <- function(obj = obj, dims = 30, 
                    ref = "/projects/amj/my_projects/reference_ARCDVC/20231010_full_integration/1_dat/1_seurat_obj/20231027_ARC_filtered.qs", 
                    column_name = "final_labels") {

  ref_obj <- qs::qread(ref)
  ref_sub <- subset(ref_obj, cells = sample(Cells(ref_obj), 20000))
  ref_sub$labels_202310 <- ifelse(grepl("Agrp",ref_sub$labels_202310), "Agrp", ref_sub$labels_202310)
  anchors <- FindTransferAnchors(reference = ref_sub, query = obj, dims = seq(dims), reference.reduction = "pca")
  predictions <- TransferData(anchorset = anchors, refdata = ref_sub$labels_202310, dims = seq(dims))
  obj[[column_name]] <- predictions$predicted.id
  obj[[column_name]] <- ifelse(grepl("Agrp", obj[[column_name]][,1]), "Agrp", obj[[column_name]][,1])
  obj[["prediction.score.max"]] <- predictions$prediction.score.max
  return(obj)

}



map_ref_merged_clusters <- function(obj = obj, dims = 30, 
                    ref = "/projects/amj/my_projects/reference_ARCDVC/20231010_full_integration/1_dat/1_seurat_obj/20231027_ARC_filtered.qs", 
                    column_name = "final_labels") {

  ref_obj <- qs::qread(ref)
  ref_obj 
  ref_sub <- subset(ref_obj, cells = sample(Cells(ref_obj), 20000))
  ref_sub$labels_202310 <- ifelse(grepl("Agrp",ref_sub$labels_202310), "Agrp", ref_sub$labels_202310)
  ref_sub@meta.data = ref_sub %>% # this merges a bunch of very specific labels under one
    `[[` %>%
    mutate(labels_202310 = labels_202310 %>% 
               str_replace_all('^Esr.*', 'Esr1') %>% # these have bad specificty to clusters, bad prediction scores
               str_replace_all('^Fez.*', 'Fez1') %>% # bad predictions scores
               str_replace_all('^Fez.*', 'Fez1') %>% # bad prediction scores
               str_replace_all('^Lepr_1|Lepr_2|Lepr_3|Lepr_4|Lepr_5', 'Lepr15') %>% # as directed, competing for specificity
               str_replace_all('^Lepr_6|Lepr_7', 'Lepr67') %>% # as directed
               str_replace_all('^Lepr_8|Lepr_9', 'Lepr89') %>% 
               str_replace_all('^Dlk1.*', 'Dlk1') %>% # bad specificity
               str_replace_all('^Sim1.*', 'Sim1') # bad specificity
          )
  ref_sub = ref_sub %>% subset(subset = labels_202310 != 'Foxp2_3') #waaay too small of a cluster, just adds noise
  anchors <- FindTransferAnchors(reference = ref_sub, query = obj, dims = seq(dims), reference.reduction = "pca")
  predictions <- TransferData(anchorset = anchors, refdata = ref_sub$labels_202310, dims = seq(dims))
  obj[[column_name]] <- predictions$predicted.id
  obj[[column_name]] <- ifelse(grepl("Agrp", obj[[column_name]][,1]), "Agrp", obj[[column_name]][,1])
  obj[["prediction.score.max"]] <- predictions$prediction.score.max
  return(obj)

}
