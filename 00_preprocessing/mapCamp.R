##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param sc_obj
##' @return
##' @author dylanmr
##' @export

mapCamp <- function(query, dims=40, class = "clust_all_neurons") {
  
  camp <- scRNAseq::CampbellBrainData()

  neurlabs <- 
    read.table("campbell_meta.txt", header=T) %>% 
    distinct(NeuronsOnlyClusters, .keep_all=T) %>% 
    filter(NeuronsOnlyClusters != "Non-neuronal") %>% 
    mutate(cluster_labs = stringr::str_split_fixed(NeuronsOnlyClusters, pattern = "[.]", n = 2)[,2] %>%  
             str_replace_all(pattern = "/", replacement = "_"),
           numeric_labs = stringr::str_split_fixed(NeuronsOnlyClusters, pattern = "[.]", n = 2)[,1],
           cluster_labs = case_when(grepl("Agrp", cluster_labs) ~ "Agrp",
                                    grepl("Ttr|Anxa2", cluster_labs) ~ "Pomc_Lepr",
                                    grepl("Kiss", cluster_labs) ~ "Kisspeptin",
                                    T ~ cluster_labs)
    )
  
  glialabs <- 
    read.table("campbell_meta.txt", header=T)  %>% 
    distinct(AllCellSubclusters, .keep_all=T) %>% 
    filter(NeuronsOnlyClusters == "Non-neuronal") %>% 
    mutate(cluster_labs = stringr::str_split_fixed(AllCellSubclusters, pattern = "[.]", n = 2)[,2] %>%  str_replace_all(pattern = "/", replacement = "_"),
           numeric_labs = stringr::str_split_fixed(AllCellSubclusters, pattern = "[.]", n = 2)[,1],
           cluster_labs = case_when(grepl("b2", cluster_labs) ~ "b2_Tany",
                                    grepl("a1", cluster_labs) ~ "a1_Tany",
                                    grepl("Epend", cluster_labs) ~ "Epend",
                                    grepl("Tuber1", cluster_labs) ~ "ParsTuber1",
                                    grepl("Tuber2", cluster_labs) ~ "ParsTuber2",
                                    grepl("Olig.*[5|6]", cluster_labs) ~ "NFOL",
                                    grepl("Olig.*[1|2|3|4]", cluster_labs) ~ "MOL",
                                    T ~ cluster_labs))
  
  if(class == 'all') {
    refdata <- camp$clust_all_neurons
  } else if(class == "neuron") {
    camp <- camp[,camp$clust_neurons%in%neurlabs$numeric_labs]
    rename <- neurlabs$numeric_labs
    names(rename) <- neurlabs$cluster_labs
    refdata <- fct_recode(camp$clust_neurons, !!!rename)
  } else if(class == "other"){
    camp <- camp[,camp$clust_all_micro%in%glialabs$numeric_labs]
    rename <- glialabs$numeric_labs
    names(rename) <- glialabs$cluster_labs
    refdata <- fct_recode(camp$clust_all_micro, !!!rename)
  }

  future::plan(future::sequential)
  
  camp <- CreateSeuratObject(counts = counts(camp)) %>% 
    SCTransform(vst.flavor = "v2")
#   camp <- scRNAseq::CampbellBrainData()
#   camp <- Seurat::as.Seurat(camp, counts = "counts", data="counts")
#   camp <- Seurat::SCTransform(camp,
#                               assay='originalexp',
#                               method="glmGamPoi",
#                               vars.to.regress="batches",
#                               vst.flavor="v2",
#                               verbose=TRUE)
  anch <- FindTransferAnchors(reference = camp, query = query, normalization.method = "SCT", recompute.residuals = F)
  predictions <- TransferData(anchorset = anch, refdata = refdata) %>% 
    dplyr::select(predicted.id, prediction.score.max)
  query <- AddMetaData(query, metadata = predictions)
  #add predictions object to dataset
  Misc(object = query, slot = paste0("predictions_", class)) <- predictions
  
  # get proportion of predicted.ids that map to each cluster
  clustering_res_name = .DollarNames(query) %>% #default is SCT_snn_res.0.8
    str_detect('res.') %>% 
    which %>% 
    .DollarNames(query)[.] %>% 
    .[order(., decreasing=TRUE)] # take the highest resolution
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
                     paste0( as.character(.x$Var2[1]),"-", as.character(.x$Var2[2]))))
  
  # add labels to object  
  labels <- unique(cluster_names$Var1) %>%  as.character()
  names(labels) <- cluster_split
  query$labels <- fct_recode(query@meta.data[[clustering_res_name]], !!!labels)
  
  return(query)

}
