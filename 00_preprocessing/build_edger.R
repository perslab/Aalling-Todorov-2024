#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param ncell min number of cells per sample
#' @param seur seurat object of a single cluster
#' @param design metadata from seurat object to include in design matrix
#' @param batch any batch variables to include in design matrix
#' @return edger object
#' @author dylanmr
#' @export
#' 
#' design = c(list(c("geno", "diet", "age")), list(c("geno", "time", "treatment")), list('treatment'), list('geno'), list(c("geno", "diet", "age"))),
#' batch = c("seq_pool", "seq_pool", "hash_pool", "hash_pool", "seq_pool"),
library(edgeR)
.build_edger <-  function(obj, ncell, design, batch) {
  
  groupvars <- unlist(design)
  
  meta <- 
    obj[[]] %>%  
    plyr::rename(c('hash.mcl.ID'='hash')) %>%
    mutate(hash = str_replace_all(hash, "-", '__')) %>%
    # mutate(Group = str_replace_all(Group, "-", '__')) %>%
    # mutate(Genotype = str_replace_all(Genotype, "-", '__')) %>%
    mutate(batch = str_replace_all(batch, fixed(" "), '__')) %>%
    rownames_to_column("cell") %>% 
    group_by(hash) %>% 
    mutate(n = n()) %>% 
    filter(n > ncell) %>% 
    ungroup()
  
  meta$group <- interaction(meta %>% dplyr::select(all_of(groupvars))) %>% fct_drop()
  
  animal <- factor(meta$hash)
  mm <- model.matrix(~ 0 + animal)
  mat.summary.mm <- obj@assays$SCT@counts[,meta$cell] %*% mm

  meta <- 
    colnames(mat.summary.mm) %>% 
    enframe() %>% 
    mutate(hash = factor(str_remove_all(value, "animal"))) %>% 
    left_join(meta %>% distinct(hash, .keep_all = T))
  
  meta <- meta[Matrix::colSums(mat.summary.mm)>0,]
  mat <- mat.summary.mm[,Matrix::colSums(mat.summary.mm)>0]

  y <- DGEList(counts=mat, samples = meta)
  keep <- filterByExpr(y, group = meta$group)
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  mds_plot <- plotMDS(y, plot=F, labels=y$samples$group)

  design <- formula(paste("~0 +", paste(c("group", batch), collapse = "+"), collapse = " "))
  mm <- model.matrix(design, data=y$samples)
  y <- estimateDisp(y, mm, robust=TRUE)
  y$symbols <- rownames(y)
  res = list(y=y, modelmat = mm, mds = mds_plot)
  return(res)
}


add_cluster_names = function(edger){
    edger = map2(edger, names(edger), function(a, b){
    a$result$name = b
    a
    })
    names_vec = edger %>% 
        map(~.x$name) %>%
        names
    
    edger
}

.fit_qlf = function(edger_cluster, contrast_str_vec){
    lev <- edger_cluster$result$modelmat
    mycontrasts = contrast_str_vec
    con = makeContrasts(contrasts=contrast_str_vec, levels=lev)
    # cmd <- paste("con <- makeContrasts(", mycontrasts, ", levels = lev)", sep ='"')
    # eval(parse(text = cmd))
    y = edger_cluster$result$y
    fit <- glmQLFit(y, lev, robust=TRUE)
    qlf <- glmQLFTest(fit, contrast=con)
    qlf
}

.fit_qlf = function(edger_cluster, contrast_str_vec){
    lev <- edger_cluster$result$modelmat
    mycontrasts = contrast_str_vec
    con = makeContrasts(contrasts=contrast_str_vec, levels=lev)
    # cmd <- paste("con <- makeContrasts(", mycontrasts, ", levels = lev)", sep ='"')
    # eval(parse(text = cmd))
    y = edger_cluster$result$y
    fit <- glmQLFit(y, lev, robust=TRUE)
    qlf <- glmQLFTest(fit, contrast=con)
    qlf$cluster_name = edger_cluster$result$name
    qlf
}

.fit_qlf_p = purrr::safely(.fit_qlf, quiet=FALSE)

get_qlf <- function(edger_list, contrasts_list){
  edger_contrasts_list = purrr::cross(list(edger=edger_list, contrasts_str=contrasts_list))
  qlf_list = edger_contrasts_list %>%
    map(~.fit_qlf_p(.x$edger, c(.x$contrasts_str))) 
  names_vec = qlf_list %>% 
    map(~paste0(.x$result$cluster_name, '___', .x$result$comparison))
  qlf_list = rlang::set_names(qlf_list, names_vec)
  qlf_list
}

.qlf2toptags = function(qlf){
    toptags = edgeR::topTags(qlf$result, p.value=1, n=1e6)
    toptags$cluster_name = qlf$result$cluster_name
    toptags
}
.qlf2toptags_p = purrr::safely(.qlf2toptags, quiet=FALSE)
get_toptags <- function(qlf_list){
    toptags_list = qlf_list %>%
      map(~.qlf2toptags_p(.x))
    toptags_list
}

build_edger <- purrr::safely(.build_edger, quiet = FALSE)
# get_qlf_p <- purrr::possibly(.get_qlf, otherwise = NULL)
safe_toptags <- purrr::possibly(edgeR::topTags, otherwise = NULL)

