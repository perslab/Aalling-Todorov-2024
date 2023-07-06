build_wgcna <- function(obj, gamma) {
  # Compute metacells using SuperCell package
  MC <- SuperCell::SCimplify(
    X = Seurat::GetAssayData(obj), # single-cell log-normalized gene expression data
    genes.use = Seurat::VariableFeatures(obj), 
    gamma = gamma,
    n.pc = 10
  )
  
  MC.counts <- SuperCell::supercell_GE(
    ge = Seurat::GetAssayData(obj, slot = "counts", assay="SCT"),
    mode = "sum", # summing counts instead of the default averaging
    groups = MC$membership
  )
  
  keep <- Matrix::rowSums(MC.counts > 5) > 10
  MC.sub <- MC.counts[keep,]
  MC.norm <- as.matrix(Seurat::NormalizeData(MC.sub))
  datExpr <- t(MC.norm) %>%  as.matrix()
  return(datExpr)
}


.plot_threshold <- function(sft, outname, powers) {
  cex1 <- 0.9
  pdf(paste0(here::here("output_files/figures/"),outname,".pdf"))
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", main = paste("Scale independence"))
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, cex = cex1, col = "red")
  abline(h = 0.80, col = "red")
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity"))
  text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
  dev.off()
}

find_threshold <- function(datExpr, outname, corfnc = "cor") {
  powers <- c(c(1:10), seq(from = 12, to = 40, by = 2))
  if(corfnc == "cor") {
    sft <- WGCNA::pickSoftThreshold(datExpr, corFnc = corfnc,
                             dataIsExpr = TRUE, powerVector = powers, corOptions = list(use = "p"),
                             networkType = "signed")
  } else {
    sft <- WGCNA::pickSoftThreshold(datExpr, corFnc = corfnc,
                             dataIsExpr = TRUE, powerVector = powers, corOptions = list(use = "p", maxPOutliers = 0.05),
                             networkType = "signed")
  }
  # .plot_threshold(sft, outname, powers)
  print(sft$powerEstimate)
  return(sft$powerEstimate)
} 


generate_tom <- function(datExpr, sftpwr, method = "average", cortype) {

  if(cortype == 1) {
    adj <- WGCNA::adjacency(datExpr, type = "signed", power = sftpwr, corFnc = "cor", corOptions = list(use = 'p'))
    diag(adj) <- 0
    TOM <- WGCNA::TOMsimilarityFromExpr(datExpr, networkType = "signed", 
                                 TOMType = "signed",maxPOutliers = 0.1,
                                 power = sftpwr)
  } else {
    adj <- WGCNA::adjacency(datExpr, type = "signed", power = sftpwr, corFnc = "bicor", corOptions = list(use = 'p', maxPOutliers = 0.1))
    diag(adj) <- 0
    TOM <- WGCNA::TOMsimilarityFromExpr(datExpr, networkType = "signed", corType = "bicor",
                                 TOMType = "signed",maxPOutliers = 0.1,
                                 power = sftpwr)
  }
  
  
  colnames(TOM) <- rownames(TOM) <- colnames(datExpr)
  dissTOM <- 1 - TOM
  # use complete for method rather than average (gives better results)
  geneTree <- hclust(as.dist(dissTOM), method = method) 
  return(list(gt = geneTree, dt = dissTOM))
}


find_modules <- function(datExpr, tom, modsize = 15, deepsplit = 1) {
  # cluster genes
  dynamicMods <- dynamicTreeCut::cutreeDynamic(
    cutHeight = 0.9999,
    dendro = tom$gt, distM = as.matrix(tom$dt),
    method = "hybrid", pamStage = F, deepSplit = deepsplit,
    minClusterSize = modsize
  )
  # name modules
  # dynamicColors <- WGCNA::labels2colors(dynamicMods) 
  return(dynamicMods)
}


merge_modules <- function(datExpr, dynamicColors, thresh = 0.2) {
  MEs <- WGCNA::moduleEigengenes(datExpr, dynamicColors)$eigengenes
  MEDiss <- 1 - cor(MEs)
  METree <- hclust(as.dist(MEDiss), method = "average")
  merge <- WGCNA::mergeCloseModules(datExpr, dynamicColors, cutHeight = thresh, verbose = 3)
  return(merge)
}

get_hubgenes <- function(merged, datExpr, label) {
  
  mergedMEs <- merged$newMEs
  modulekME <- WGCNA::signedKME(datExpr, mergedMEs)
  
  hubgenes <- 
    dplyr::bind_cols(modulekME, cols = merged$colors) %>% 
    rownames_to_column("gene") %>% 
    pivot_longer(c(-cols,-gene)) %>% 
    mutate(cols = paste0("kME", cols)) %>% 
    filter(cols == name) %>%
    filter(!grepl("^kME0$", cols)) %>% 
    dplyr::select(-cols) %>%
    mutate(label = label)
   
   hubgenes
}


.me_2_mod = function(vec, prefix){
    vec %>%
    stringr::str_replace('ME', '') %>%
    stringr::str_pad(2, 'left', pad='0') %>%
    paste0(prefix, '_M', .)
}


calc_eigengenes <- function(obj, merged, datExpr, label) {
  
  obj <- Seurat::NormalizeData(obj)
  ME <- WGCNA::moduleEigengenes(t(as.matrix(obj@assays$SCT@data[match(colnames(datExpr),rownames(obj@assays$SCT@data)),])), merged$colors)$eigengenes
  ME <- dplyr::bind_cols(ME, obj@meta.data)
  colnames(ME) = colnames(ME) %>%
    ifelse(stringr::str_detect(., 'ME'), .me_2_mod(., label), .)
  ME
}


make_gost = function(hubgenes, label){
    hubgenes = hubgenes %>% 
      mutate(new_name = name) %>%
      mutate(new_name = stringr::str_replace(new_name, 'kME', '')) %>% 
      mutate(new_name = stringr::str_pad(new_name, 2, 'left', pad='0')) %>%
      mutate(new_name = paste0(label, '_M', new_name))
    module_members = split(hubgenes$gene,hubgenes$new_name)
    gost_results = gprofiler2::gost(query = module_members,
                                    organism = "mmusculus",
                                    ordered_query = FALSE,
                                    exclude_iea = FALSE, 
                                    evcodes = TRUE)
    gost_results$result = gost_results$result %>%
      rowwise %>%
      mutate(gene_ids = str_split(intersection, fixed(','))) %>%
      mutate(gene_ids = {hubgenes %>% 
                          filter(new_name %in% query) %>%
                          filter(gene %in% gene_ids) %>% 
                          pull(gene) %>% 
                          paste0(collapse=',')}) %>%
        ungroup
    gost_results
}


get_module_cols = function(ME){
    module_cols = ME %>%
    colnames %>%
    stringr::str_subset('_M\\d{2,}$') %>%
    stringr::str_subset('_M00', negate = TRUE)
    module_cols
}


filter_contrasts = function(ME, contrasts){
    pattern = paste0({ME %>% pull(group) %>% unique}, collapse='|')
    contrasts = contrasts %>%
        mutate(expression_replaced = str_replace_all(expression, pattern, '')) %>%
        mutate(contains_all_terms = !str_detect(expression_replaced, '[a-zA-Z0-9]')) %>%
        filter(contains_all_terms == TRUE) %>%
        select(name, expression)
    contrasts
}


fit_limma_me = function(ME, formula_str){
    contrasts= tribble(
            ~"name", ~"expression",
            "Day5.obob.LDNvsVehicle", "LDN.Day5.obob - Vehicle.Day5.obob",
            "Day5.obob.VehiclevsVehPF", "Vehicle.Day5.obob - Veh_PF.Day5.obob",
            "Day5.obob.LDNvsFGF1", "LDN.Day5.obob - FGF1.Day5.obob",
            "Day5.obob.VehiclevsFGF1", "Vehicle.Day5.obob - FGF1.Day5.obob",
            "Day5.obob.FGF1vsVehPF", "FGF1.Day5.obob - Veh_PF.Day5.obob",
            "Day14.obob.FGF1vsVehPF", "FGF1.Day14.obob - Veh_PF.Day14.obob",
            "Day5vsDay14.BL6.VehPF", "Veh_PF.Day5.BL6 - Veh_PF.Day14.BL6",
            "Day5vsDay14.obob.VehPF", "Veh_PF.Day5.obob - Veh_PF.Day14.obob")
    contrasts = filter_contrasts(ME, contrasts)
    contrasts = contrasts %>%
        deframe
    module_cols = get_module_cols(ME)
    eset = ME[module_cols] %>% t
    design = model.matrix(formula(formula_str), data = ME)
    colnames(design) = colnames(design) %>%
        str_replace('group', '') %>%
        str_replace('batch', '')
    corfit <- limma::duplicateCorrelation(eset,design,block=ME$hash.mcl.ID)
    fit <- limma::lmFit(eset,design,block=ME$hash.mcl.ID,correlation=corfit$consensus)
    cm = limma::makeContrasts(
        contrasts=contrasts,
        levels=design)
    colnames(cm) = names(contrasts)
    fit2 <- limma::contrasts.fit(fit, cm)
    fit2 <- limma::eBayes(fit2)
    fit2
}



.get_sig_mod_single = function(fit2, coef){
    limma::topTable(fit2, coef = coef, number=Inf, sort.by="none") %>%
    mutate(significant = case_when(adj.P.Val < 0.05 ~ TRUE,
                                   adj.P.Val >= 0.05 ~ FALSE)) %>%
    mutate(comparison = coef) %>%
    rownames_to_column(var = 'module') %>%
    relocate(comparison, module)
}


get_sig_mods = function(fit2){
    top_table = colnames(fit2$contrasts) %>% 
        purrr::map(\(x) .get_sig_mod_single(fit2, x)) %>% 
        do.call(rbind, .)
    top_table
}


make_gost_tt_summary = function(gost_results, top_table, max_term_size=500){
    result = gost_results$result = gost_results$result %>%
        filter(term_size <= max_term_size) %>%
        rename(module = query)
    #join
    result = left_join(x=top_table, y=result, by = 'module')
    #make it smaller
    drop_cols = c("AveExpr", "t", "P.Value", "B", "precision", "recall", 
                  "source_order", "parents", "evidence_codes", "effective_domain_size")
    result = result %>%
    select(-any_of(drop_cols)) %>%
    select(-significant.y) %>% #gprofiler returns only significant
    rename(lm_logFC = logFC) %>%
    rename(lm_padj = adj.P.Val) %>%
    rename(lm_significant = significant.x) %>%
    rename(gp_padj = p_value) %>%
    relocate(term_id, .after=gp_padj) %>%
    relocate(term_name, .after=term_id)
    result
}


make_hubgenes_summary = function(hubgenes){
    hubgenes_summary = hubgenes %>% 
    mutate(new_name = name) %>%
    mutate(new_name = stringr::str_replace(new_name, 'kME', '')) %>% 
    mutate(new_name = stringr::str_pad(new_name, 2, 'left', pad='0')) %>%
    mutate(new_name = paste0(label, '_M', new_name)) %>%
    group_by(new_name) %>%
    arrange(desc(value)) %>%
    select(gene, new_name) %>%
    distinct %>%
    mutate(genes = paste0(gene, collapse=', ')) %>%
    mutate(module_size = n()) %>%
    select(new_name, module_size, genes) %>%
    distinct %>%
    arrange(new_name) %>%
    rowwise %>% 
    mutate(label = str_replace(new_name, "_M\\d{2,}$", '')) %>%
    rename(module_name = new_name) %>%
    relocate(label) %>%
    ungroup
    
    hubgenes_summary
}


make_hubgenes_tt_summary = function(hubgenes, top_table){
    hubgenes_summary = make_hubgenes_summary(hubgenes)
    
    hubgenes_tt_summary = top_table %>%
        rename(module_name = module) %>%
        select(module_name, comparison, logFC, adj.P.Val, significant) %>%
        left_join(hubgenes_summary, ., by = 'module_name') %>%
        relocate(genes, .after = significant) %>%
        arrange(label, module_name, comparison)
    
    hubgenes_tt_summary    
}

#summarize the make_hubgenes_tt_summary with directional chrs
make_combined_hubgenes_summary_tt_all_summary = function(combined_hubgenes_summary_tt_all){
    comparison_levels = c("Day5.obob.LDNvsVehicle",
                          "Day5.obob.VehiclevsVehPF",
                          "Day5.obob.LDNvsFGF1",
                          "Day5.obob.VehiclevsFGF1",
                          "Day5.obob.FGF1vsVehPF",
                          "Day14.obob.FGF1vsVehPF",
                          "Day5vsDay14.obob.VehPF",
                          "Day5vsDay14.BL6.VehPF")
    combined_hubgenes_summary_tt_all %>% 
    filter(!is.na(comparison)) %>% 
    mutate(adj_logFC = case_when(significant == FALSE ~ 0,
                             TRUE ~ logFC)) %>%
    mutate(direction = case_when(adj_logFC == 0 ~ '-',
                                 adj_logFC > 0 ~ 'UP',
                                 adj_logFC < 0 ~ 'DOWN')) %>%
    select(label, module_name, module_size, comparison, direction) %>% 
    pivot_wider(names_from = comparison, values_from = direction) %>%
    relocate(any_of(comparison_levels), .after=module_size) %>%
    mutate(FGF1_Day5 = case_when((Day5.obob.FGF1vsVehPF != '-') & !(Day5vsDay14.obob.VehPF != '-' | Day5vsDay14.BL6.VehPF != '-') ~ Day5.obob.FGF1vsVehPF,
                                 TRUE ~ '-')) %>%
    mutate(FGF1_Day14 = case_when((Day14.obob.FGF1vsVehPF != '-') & !(Day5vsDay14.obob.VehPF != '-' | Day5vsDay14.BL6.VehPF != '-') ~ Day14.obob.FGF1vsVehPF,
                                  TRUE ~ '-')) %>%
    mutate(FGF1_both_days_same = case_when( ((FGF1_Day5 == FGF1_Day14) & (FGF1_Day5 != '-') & (FGF1_Day14 != '-')) ~ FGF1_Day5,
                                           TRUE ~ '-')) %>%
    mutate(FGF1_both_days_diff = case_when( ((FGF1_Day5 != FGF1_Day14) & (FGF1_Day5 != '-') & (FGF1_Day14 != '-')) ~ FGF1_Day5,
                                           TRUE ~ '-')) %>%
    mutate(LDN = case_when((Day5.obob.LDNvsVehicle != '-') ~ Day5.obob.LDNvsVehicle,
                                 TRUE ~ '-')) %>%
    mutate(LDN_FGF1_same = case_when( ((LDN != '-') & (LDN == FGF1_Day5) | (LDN == FGF1_Day14)) ~ LDN,
                                      TRUE ~ '-')) %>%
    mutate(LDN_FGF1_diff = case_when( ((LDN != '-') & ((LDN != FGF1_Day5) & (FGF1_Day5 != '-')) | ((LDN == FGF1_Day14) & (FGF1_Day14 != '-'))) ~ LDN,
                                      TRUE ~ '-')) %>%
    relocate(any_of(c("LDN", "LDN_FGF1_same", "FGF1_Day5", "FGF1_Day14", "FGF1_both_days_same", "FGF1_both_days_diff")), .after=module_size)
}


#summarize the make_hubgenes_tt_summary with directional chrs
make_combined_hubgenes_summary_tt_all_summary_bmp1_only = function(combined_hubgenes_summary_tt_all){
    comparison_levels = c("Day5.obob.LDNvsVehicle")
    combined_hubgenes_summary_tt_all %>% 
    filter(!is.na(comparison)) %>% 
    mutate(adj_logFC = case_when(significant == FALSE ~ 0,
                             TRUE ~ logFC)) %>%
    mutate(direction = case_when(adj_logFC == 0 ~ '-',
                                 adj_logFC > 0 ~ 'UP',
                                 adj_logFC < 0 ~ 'DOWN')) %>%
    select(label, module_name, module_size, comparison, direction) %>% 
    pivot_wider(names_from = comparison, values_from = direction) %>%
    relocate(any_of(comparison_levels), .after=module_size) %>%
    mutate(LDN = case_when((Day5.obob.LDNvsVehicle != '-') ~ Day5.obob.LDNvsVehicle,
                                 TRUE ~ '-')) %>%
    relocate(any_of(c("LDN")), .after=module_size)
}

