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


gbr_genes_stat = function(sce, n_genes, batch=NULL){
    genes_stat = geneBasisR::gene_search(sce,
                             n_genes_total = n_genes,
                             verbose = T)
    genes_stat
}


gbr_genes_stat_preselected_discard = function(sce, base_genes=NULL, n_genes=100, batch=NULL, drop_nng=TRUE, discard_genes=NULL){
    if (drop_nng){
        nng_vec = sce %>% get_sce_genes %>% get_not_named_genes_vec
        nng_vec = c(discard_genes, nng_vec)
    } else {
        nng_vec = discard_genes
    }
    genes_stat = geneBasisR::gene_search(sce,
                                         genes_base=base_genes,
                             n_genes_total = n_genes,
                                         genes.discard = nng_vec,
                                         verbose = T)
    genes_stat
}


gbr_genes_stat_preselected = function(sce, base_genes=NULL, n_genes=100, batch=NULL, drop_nng=TRUE){
    if (drop_nng){
        nng_vec = sce %>% get_sce_genes %>% get_not_named_genes_vec
    } else {
        nng_vec = NULL
    }
    genes_stat = geneBasisR::gene_search(sce,
                                         genes_base=base_genes,
                                         n_genes_total = n_genes,
                                         genes.discard = nng_vec,
                             verbose = T)
    genes_stat
}

get_selected_genes = function(genes_stat){
    genes = genes_stat$gene
    genes
}

gbr_get_celltype_mapping = function(sce, selected_genes, celltype_id, batch=NULL){
    celltype_mapping = geneBasisR::get_celltype_mapping(sce, 
                                                        genes.selection = selected_genes$gene,
                                                        celltype.id = celltype_id,
                                                        batch = batch,
                                                        return.stat = T)
    celltype_mapping
}


gbr_evaluate_lib = function(sce, selected_genes, celltype_id, library.size_type = "single", step=10, batch=NULL){
    lib_stat = geneBasisR::evaluate_library(sce, 
                                            selected_genes$gene, 
                                            genes.all = rownames(sce), 
                                            library.size_type = library.size_type,
                                            celltype.id = celltype_id,
                                            batch = batch,
                                            n_genes.step = step,
                                            return.cell_score_stat = T, return.gene_score_stat = T, return.celltype_stat = T,
                                            verbose = T)
    lib_stat
}


gbr_calc_redundancy_stat = function(sce, selected_genes, celltype_id, batch=NULL){
    redundancy_stat = geneBasisR::get_redundancy_stat(sce,
                                                      selected_genes$gene,
                                                      genes_to_assess = selected_genes$gene,
                                                      batch = batch,
                                                      celltype.id = celltype_id
                                                      ) 
}


do_first_sw_selection = function(genes_stat, n_first_selection){
    genes_stat_first_selection = genes_stat[1:n_first_selection,]
    genes_stat_first_selection
}

do_next_sw_selection = function(next_selection_obj, n_next_selection, gene_stat_prior){
    if (n_next_selection != 0){
        n_preselected = dim(gene_stat_prior)[1]
        genes_vec_prior_selection = gene_stat_prior %>% pull(gene)
        n_genes = n_preselected + n_next_selection
        sce_f = gbr_retain_informative_genes_preselected(next_selection_obj,
                                                        n=NULL,
                                                        preselected_genes=genes_vec_prior_selection)
        genes_stat_next_selection = gbr_genes_stat_preselected(sce_f,
                                                               genes_vec_prior_selection,
                                                               n_genes=n_genes,
                                                               batch=NULL)
    } else {
        genes_stat_next_selection = gene_stat_prior
    }
    genes_stat_next_selection
}


make_ctm_sw_selection = function(sce, celltype_id, gene_stat_prior, batch=NULL){
    genes_vec_prior_selection = gene_stat_prior %>% pull(gene)
    sce_f = gbr_retain_informative_genes_preselected(sce,
                                                     n=NULL,
                                                     preselected_genes=genes_vec_prior_selection)
    ctm = gbr_get_celltype_mapping(sce_f, gene_stat_prior, celltype_id, batch=batch)
    ctm
}

make_lib_stat_sw_selection = function(sce, gene_stat_prior, celltype_id, library.size_type = "single", step=10, batch=NULL){
    genes_vec_prior_selection = gene_stat_prior %>% pull(gene)
    sce_f = gbr_retain_informative_genes_preselected(sce,
                                                     n=NULL,
                                                     preselected_genes=genes_vec_prior_selection)
    lib_stat = gbr_evaluate_lib(sce_f, gene_stat_prior, celltype_id, library.size_type=library.size_type, step=step, batch=batch)
    lib_stat
}


make_cell_score_stat_sw_selection = function(lib_stat, name, sce){
    cell_score_stat = lib_stat %>% 
        `[[`("cell_score_stat") %>%
        mutate(name = name)
    meta = as.data.frame(colData(sce)) %>% 
        select(labels, polar_label) %>%
        rownames_to_column(var = 'cell')
    cell_score_stat = left_join(cell_score_stat, meta) %>% tibble
    cell_score_stat
}


summarise_cell_score_stat_sw_selection = function(cell_score_stat){
    cell_score_stat_summary = cell_score_stat %>%
        group_by(labels) %>%
        mutate(mean_cell_score.labels = mean(cell_score)) %>%
        ungroup %>% group_by(polar_label) %>%
        mutate(mean_cell_score.polar_label = mean(cell_score)) %>%
        select(name, labels, polar_label, mean_cell_score.labels, mean_cell_score.polar_label) %>%
        distinct %>% 
        arrange(polar_label)
    cell_score_stat_summary
}


get_sce_genes = function(sce){
    gene_names = rownames(logcounts(sce))
    gene_names
}

get_literature_genes_df = function(){
    lgdf = data.frame(gene = c("Pdgfra", # (OPC)
                               "Bmp4", # (NFOL)
                               "Plp1", # (MOL)
                               "Aqp4", # (Astro)
                               "Rax", # (tany)
                               "Agrp",
                               "Pomc",
                               "Htr3b",
                               "Lef1",
                               "Lmx1a"))
    lgdf
}


plot_ctm = function(ctm, name='', col_order=NULL){
    if (is.null(col_order)){
        col_order = ctm$mapping %>% pull(celltype) %>% unique %>% sort
    } 
    ctm$mapping = ctm$mapping %>%
    mutate(celltype = factor(celltype, levels=col_order)) %>%
    arrange(celltype)
    p = geneBasisR::plot_mapping_heatmap(ctm$mapping, title = paste0("Cell type confusion matrix    ", name))
    p
}


plot_cell_preservation_score = function(sce, lib_stat, n_genes_total, column){
    meta = colData(sce) %>% data.frame %>% select(labels, polar_label, cell_class)
    cell_score_stat = lib_stat$cell_score_stat[lib_stat$cell_score_stat$n_genes == n_genes_total , ]
    rownames(cell_score_stat) = cell_score_stat$cell
    cell_score_stat = left_join(rownames_to_column(cell_score_stat), rownames_to_column(meta)) %>% column_to_rownames
    column = ensym(column)
    p = ggplot(cell_score_stat , aes(x = !!column, y = cell_score, fill = !!column)) + 
      geom_boxplot() + 
      theme_classic() + 
      labs(y = "Cell neighborhood preservation score" , x = "# genes") + 
      theme(legend.position = "none") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ylim(0,1)
    p
}
