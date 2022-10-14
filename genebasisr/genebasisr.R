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

gbr_genes_stat = function(sce, n_genes){
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