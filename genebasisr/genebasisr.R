library(geneBasisR)
library(SingleCellExperiment)
library(tibble) 
library(ggplot2)
library(ggpubr)
library(Seurat)


gbr_retain_informative_genes = function(sce, n=NULL){
    sce = geneBasisR::retain_informative_genes(sce, n=n)
    sce
}

get_not_named_genes_vec = function(vec_of_genes){
    a_rik = vec_of_genes %>% stringr::str_ends("Rik")
    an_ensmusg = vec_of_genes %>% stringr::str_starts("ENSMUSG")
    a_Gm = vec_of_genes %>% stringr::str_starts("Gm\\d")
    starts_w4_or_more_digits = vec_of_genes %>% stringr::str_ends("\\d{4,}")
    not_named_genes = (a_rik | an_ensmusg | a_Gm | starts_w4_or_more_digits) %>% which
    vec_of_genes[not_named_genes]
}


gbr_retain_informative_genes_preselected = function(sce, n=NULL, preselected_genes=preselected_genes){
    sce_genes_excluding_panel = rownames(logcounts(sce)) %>% setdiff(preselected_genes)
    sce_excluding_panel = sce[sce_genes_excluding_panel,]
    sce_excluding_panel = geneBasisR::retain_informative_genes(sce, n=n)
    infomrative_genes_excluding_panel = rownames(logcounts(sce_excluding_panel))
    selected_genes = union(preselected_genes, infomrative_genes_excluding_panel)
    sce_filtered = sce[selected_genes,]
    sce_filtered
}


    genes_stat = geneBasisR::gene_search(sce,
                             n_genes_total = n_genes,
                             verbose = T)
    genes_stat
}

get_selected_genes = function(genes_stat){
    genes = genes_stat$gene
    genes
}

gbr_get_celltype_mapping = function(sce, genes, celltype_id){
    celltype_mapping = get_celltype_mapping(sce , genes.selection = genes, celltype.id = celltype_id, return.stat = F)
    celltype_mapping
}