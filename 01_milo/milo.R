pyget = function(some_list, some_item, default_value){
    some_list %>%
    `[[`(some_item) %>%
    ifelse(is.null(.), default_value, .)
}


prep_obj_for_milo = function(obj){
    obj@meta.data = obj@meta.data %>% mutate(batch = stringr::str_replace_all(batch, stringr::fixed(" "), '__'))
    obj@meta.data = obj@meta.data %>% mutate(labels = stringr::str_replace_all(labels, stringr::fixed("-"), '__'))
    obj@meta.data$group = interaction(obj@meta.data$treatment, obj@meta.data$time, obj@meta.data$strain, drop = TRUE)
    obj
}


convert_obj_to_sce = function(obj){
    sce = Seurat::as.SingleCellExperiment(obj)
    sce
}


make_milo = function(sce){
    milo_obj = miloR::Milo(sce) %>%
    .milo_buildGraph() %>%
    .milo_makeNhoods() %>%
    .milo_countCells() %>%
    .milo_calcNhoodDistance() %>%
    .milo_buildNhoodGraph()
    milo_obj
}


.milo_buildGraph = function(milo_obj, ...){
    milo_obj <- miloR::buildGraph(milo_obj, k=40, d=30, reduced.dim = 'PCA')
    milo_obj
}

.milo_makeNhoods = function(milo_obj){
    milo_obj <- miloR::makeNhoods(milo_obj, prop = 0.1, k= 40, d=30, refined = T, reduced_dims = 'PCA')
    milo_obj
}


.milo_countCells = function(milo_obj, sample='hash.mcl.ID'){
    milo_obj = miloR::countCells(milo_obj, meta.data = as.data.frame(SingleCellExperiment::colData(milo_obj)), sample = sample)
    milo_obj
}

.milo_calcNhoodDistance = function(milo_obj){
    milo_obj <- miloR::calcNhoodDistance(milo_obj, d=30, reduced.dim = 'PCA')
    milo_obj
}

.milo_buildNhoodGraph = function(milo_obj){
    milo_obj = miloR::buildNhoodGraph(milo_obj)
    milo_obj
}


make_design_df = function(milo_obj){
    design_df <- data.frame(SingleCellExperiment::colData(milo_obj))[, c('hash.mcl.ID', 'group', 'batch', "strain")]
    #convert seq-pool to factor
    design_df$batch <- as.factor(design_df$batch)
    #keep unique rows
    design_df <- dplyr::distinct(design_df)
    #change rownames
    rownames(design_df) <- design_df$hash.mcl.ID
    design_df
}

make_design_df_bmp1 = function(milo_obj){
    design_df <- data.frame(SingleCellExperiment::colData(milo_obj))[, c('hash.mcl.ID', 'group', "strain")]
    #convert seq-pool to factor
    #keep unique rows
    design_df <- dplyr::distinct(design_df)
    #change rownames
    rownames(design_df) <- design_df$hash.mcl.ID
    design_df
}


make_model_matrix = function(design_df){
    model <- model.matrix(~0 + group + batch, data=design_df)
    model
}

make_model_matrix_bmp1 = function(design_df){
    model <- model.matrix(~0 + group, data=design_df)
    model
}


get_da_results = function(milo_obj=milo_obj, model_matrix=model_matrix, design_df=design_df, model_contrasts_str=model_contrasts_str){
    da_results <- miloR::testNhoods(milo_obj,
                                    design = model_matrix,
                                    design.df = design_df,
                                    model.contrasts = c(model_contrasts_str))
    da_results = da_results %>%
        dplyr::arrange(SpatialFDR) 
    da_results
}

single_split = function(obj, label) {
    obj = subset(x = obj, subset = labels == label)
    obj
}


add_polarity_to_da_results = function(da_results, spatial_fdr_cutoff=0.1){
    da_results = da_results %>% 
     mutate(
        polarity = case_when(
          logFC > 0 & SpatialFDR < spatial_fdr_cutoff ~ "pos",
          logFC < 0 & SpatialFDR < spatial_fdr_cutoff ~ "neg",
          TRUE ~ "none"
        )
      )
    da_results
}


.get_polar_cells = function(milo, da_results, direction, frac){
    barcodes = c()
    dgc = milo@nhoods
    polar_hoods = da_results %>% filter(polarity == direction) %>% pull(Nhood)
    if (length(polar_hoods) > 0){
        ranks_df = dgc[, polar_hoods] %>% rowSums %>% rank(ties='min') %>% data.frame()
        colnames(ranks_df) = 'rank'
        if (length(unique(ranks_df$rank)) != 1){
            how_many = round(dim(ranks_df)[1]*frac)
            barcodes = ranks_df %>% arrange(desc(rank)) %>% top_n(how_many, rank) %>% rownames
        }
    }
    barcodes
}


.assign_cell_polarity = function(milo, da_results, frac=0.1){
    pos_cells = .get_polar_cells(milo, da_results, "pos", frac)
    neg_cells = .get_polar_cells(milo, da_results, "neg", frac)
    none_cells = rownames(milo@nhoods)
    pos_and_neg_cells = intersect(pos_cells, neg_cells)
    pos_cells = setdiff(pos_cells, pos_and_neg_cells)
    neg_cells = setdiff(neg_cells, pos_and_neg_cells)
    none_cells = setdiff(none_cells, union(pos_cells, neg_cells))
    none_cells = union(none_cells, pos_and_neg_cells)
    if (!is.null(pos_cells)){
        pos_tib = tibble(barcode = pos_cells, polarity='pos')
    } else {
        pos_tib = tibble()
    }
    if (!is.null(neg_cells)){
        neg_tib = tibble(barcode = neg_cells, polarity='neg')
    } else {
        neg_tib = tibble()
    }
    none_tib = tibble(barcode = none_cells, polarity='none')
    all_tib = rbind(pos_tib, neg_tib, none_tib)
    polarity_df = data.frame(all_tib)
    rownames(polarity_df) = polarity_df$barcode
    polarity_df
}


make_polar_labels_df = function(milo, da_results, top_frac=0.1){
    labels_df = milo %>% `@`('colData') %>% data.frame %>% select(labels)
    polarity_df = .assign_cell_polarity(milo, da_results, frac=top_frac)
    labels_df = merge(polarity_df, labels_df, by=0) %>%
      tibble::column_to_rownames("Row.names") %>%
      mutate(polar_label = paste0(labels, '.', polarity)) %>%
      select(-c(labels, polarity)) %>%
      mutate(barcode =  factor(barcode, levels = rownames(milo %>% `@`('colData')))) %>%
      arrange(barcode) %>%
      select(-barcode)
    labels_df
}

plot_polarity_beeswarm = function(da_results){
    options(repr.plot.width=10, repr.plot.height=5)
    pdab = plotDAbeeswarm(da_results, group.by = "polarity", )
    pdab
}


annotate_nhoods = function(da_results, milo_obj, col){
    if (milo_obj@colData %>% 
        as.data.frame %>% 
        pull(col) %>% 
        unique %>% 
        length %>% 
        `==`(1)){
        da_results[col] = milo_obj@colData %>% as.data.frame %>% pull(col) %>% unique
        da_results[paste0(col, '_fraction')] = as.double(1)
        # this checks if there is only one unique value
        # because some milo internals use a factor, this gets a dim error
    } else {
        da_results <- miloR::annotateNhoods(milo_obj, da_results, coldata_col = col)
    }
    da_results
}


#### new annotation methods

make_logFC_df = function(da_results){
    logFC_df = da_results %>% 
        arrange(Nhood) %>% select(logFC) %>% t %>% as.data.frame
    logFC_df 
}


make_nhm = function(milo_obj){
    nh = milo_obj@nhoods
    nh@Dimnames[[2]] = seq(1:length(nh@Dimnames[[2]]))
    nhm = nh %>% as.matrix %>% as.data.frame
    nhm
}


annotate_nhood_counts = function(da_results, nhm){
    nhood_counts = nhm %>% colSums %>% enframe %>%
      transmute(Nhood=name, n_cells=value) %>%
      mutate(Nhood = as.double(Nhood))
    da_results = dplyr::left_join(x=da_results, y=nhood_counts)
    da_results
}


make_nhood_summary = function(da_results, nhm){
    nhood_summary = da_results %>% 
        # filter(SpatialFDR < 0.1) %>%
        group_by(polarity) %>%
        summarise(n_cells = sum(n_cells),
                  n_nhoods = n()) %>%
        mutate(cell_freq = n_cells/sum(n_cells),
               nhood_freq = n_nhoods/sum(n_nhoods)) %>%
        mutate(n_cells_nhood = n_cells)
    nhood_summary = nhood_summary %>%
        mutate(n_cells = round(dim(nhm)[1] * cell_freq))
    nhood_summary
}


check_split_validity = function(nhood_summary, min_cluster_frac=0.025, min_cells=200){
    nhood_summary = nhood_summary %>%
        mutate(split_valid = case_when((cell_freq >= min_cluster_frac & n_cells >= min_cells) ~ TRUE,
                                       TRUE ~ FALSE)) %>%
        mutate(adj_split = NA)
    if (dim(nhood_summary)[1] > 1){
        nhood_summary = nhood_summary %>% 
        mutate(adj_split = case_when(polarity == 'neg' ~ "none",
                                     polarity == 'pos' ~ "none",
                                     polarity == 'none' ~ 'huh'))
        if (nhood_summary %>% pull(adj_split) %>% `==`('huh') %>% any()){
            # this gets the smallest adj split to none, if 'none' is too small
            pos_or_neg = nhood_summary %>% filter(polarity != 'none') %>% arrange(n_cells) %>% pull(polarity) %>% `[[`(1)
            nhood_summary = nhood_summary %>% 
              mutate(adj_split = case_when(adj_split == "huh" ~ pos_or_neg,
                                           TRUE ~ adj_split))
        }
    }
    nhood_summary
}


combine_splits = function(da_results, nhm, min_cluster_frac=0.025, min_cells=100){
    while (TRUE){
        nhood_summary = make_nhood_summary(da_results, nhm) %>%
          check_split_validity(min_cluster_frac=0.025, min_cells=200)
        if (dim(nhood_summary)[1] == 1){
            # already combined!
            return(da_results)
        } else if (nhood_summary %>% pull(split_valid) %>% all()) {
            # all splits valid
            return(da_results)
        } else {
            # merges the smallest cluster
            split_to_merge = nhood_summary %>% 
                filter(split_valid == FALSE) %>%
                arrange(n_cells) %>%
                slice(1) 
            old_polarity = split_to_merge$polarity
            new_polarity = split_to_merge$adj_split
            da_results = da_results %>%
                mutate(polarity = case_when(polarity == old_polarity ~ new_polarity,
                                            TRUE ~ polarity))
        }
    }
}


make_logFC_cells_df = function(nhm, logFC_df){
    logFC_cells = mapply(`*`, nhm, (logFC_df)) %>% as.data.frame
    rownames(logFC_cells) = rownames(nhm)
    logFC_cells
}


nonzero_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  x = x[x != 0]
  if (length(x) > 0){
   nzm = mean(x)
  } else {
  nzm = 0
  }
  nzm
}


make_nzdf = function(logFC_cells){
    nonzero_mean_df = logFC_cells %>% rowwise() %>%
    summarise(nonzero_mean = nonzero_mean(c_across(everything())),
              nonzero_sum = sum(c_across(everything())))
    nonzero_mean_df$rownames = rownames(logFC_cells)
    nonzero_mean_df = nonzero_mean_df %>% as.data.frame %>% column_to_rownames(var = 'rownames')
    nonzero_mean_df = nonzero_mean_df %>% arrange(nonzero_sum)
    nonzero_mean_df
}


get_cell_freq = function(nhood_summary, split_name){
  cell_freq = ifelse(nhood_summary %>% filter(polarity == split_name) %>% nrow %>% `>`(0),
                     nhood_summary %>% filter(polarity == split_name) %>% pull(cell_freq),
                     as.double(0))
  cell_freq
}


make_barcode_polarity_df = function(nonzero_mean_df, nhood_summary){
    if (nhood_summary %>% filter(polarity == 'neg') %>% dim %>% `[`(1) %>% `>`(0)){
        neg_cells = top_frac(nonzero_mean_df, -get_cell_freq(nhood_summary, 'neg'), wt=nonzero_sum) %>% 
                      rownames %>% 
                      data.frame(barcodes=., polarity=rep('neg', length(.)))
    } else {
        neg_cells = data.frame(barcodes=character(), polarity=character())
    }
    if (nhood_summary %>% filter(polarity == 'pos') %>% dim %>% `[`(1) %>% `>`(0)){
        pos_cells = top_frac(nonzero_mean_df, get_cell_freq(nhood_summary, 'pos'), wt=nonzero_sum) %>% 
                      rownames %>% 
                      data.frame(barcodes=., polarity=rep('pos', length(.)))
    } else {
        pos_cells = data.frame(barcodes=character(), polarity=character())
    }
    none_cells = setdiff(rownames(nonzero_mean_df), union(neg_cells$barcodes, pos_cells$barcodes)) %>%
                 data.frame(barcodes=., polarity=rep('none', length(.)))
    barcode_df = rbind(neg_cells, none_cells, pos_cells)
    barcode_df
}


make_pldf = function(barcode_df, milo_obj){
    labels_df = milo_obj %>% `@`('colData') %>% data.frame %>% select(labels) %>% rownames_to_column(var = "barcodes")
    pldf = left_join(labels_df, barcode_df) %>% 
        mutate(polar_label = paste0(labels, '.', polarity)) %>%
        select(-c(labels, polarity)) %>%
        mutate(barcodes =  factor(barcodes, levels = labels_df$barcodes)) %>%
        column_to_rownames("barcodes")
    pldf
}


.get_only_named_genes_vec = function(vec_of_genes){
    not_a_rik = vec_of_genes %>% stringr::str_ends("Rik") %>% `!`
    not_an_ensmusg = vec_of_genes %>% stringr::str_starts("ENSMUSG") %>% `!`
    not_a_Gm = vec_of_genes %>% stringr::str_starts("Gm\\d") %>% `!`
    not_4_or_more_digits = vec_of_genes %>%  stringr::str_ends("\\d{4,}") %>% `!`
    only_named_genes = (not_a_rik & not_an_ensmusg & not_a_Gm & not_4_or_more_digits) %>% which
    only_named_genes
}


get_milo_ngo_genes_vec = function(milo_assay){
    #milo assay such as milo_obj@assays@data$logcounts
    ngo_genes_idx = milo_assay@Dimnames[[1]] %>% .get_only_named_genes_vec
    ngo_genes_vec = milo_assay@Dimnames[[1]][ngo_genes_idx]
    ngo_genes_vec
}


set_polarity_da_results_NhoodGroup = function(da_results){
    da_results$NhoodGroup = da_results$polarity
    da_results
}


get_nhood_markers = function(milo_obj, da_results, sample_col, 
                             grouping_col="NhoodGroup",
                             gene_subset=NULL,
                             subset_groups=NULL,
                             tag=''){
    subset_nhoods = NULL
    if (!is.null(subset_groups)){
        subset_groups = stringr::str_split(subset_groups, pattern=fixed('.')) %>% unlist
        if (length(subset_groups) > 1){
            subset_nhoods = da_results$NhoodGroup %in% subset_groups
        }
    }
    if (grouping_col != "NhoodGroup"){
        da_results['NhoodGroup'] = da_results[grouping_col]
    }
    nhood_markers <- miloR::findNhoodGroupMarkers(milo_obj, da_results, subset.row = gene_subset,
                                                aggregate.samples = TRUE, sample_col = sample_col,
                                                subset.groups=subset_groups,
                                                subset.nhoods = subset_nhoods)
    nhood_markers = nhood_markers %>% 
        relocate(GeneID) %>%
        remove_rownames
    nhood_markers = nhood_markers %>% arrange(.[[3]])
    if (length(tag) > 0){
        nhood_markers$tag = tag
    }
    nhood_markers
}


get_top_milo_degs = function(nhood_markers, subset_groups, n_degs=5, tag){
    subset_groups = stringr::str_split(subset_groups, pattern=fixed('.')) %>% unlist
    p_adj_1 = paste0("adj.P.Val_", subset_groups[1])
    p_adj_2 = paste0("adj.P.Val_", subset_groups[2])
    top_degs = nhood_markers %>% 
        filter(((!!rlang::sym(p_adj_1)) < 0.05) & ((!!rlang::sym(p_adj_2)) < 0.05)) %>%
        mutate(p_product = (!!rlang::sym(p_adj_1)) * (!!rlang::sym(p_adj_2))) %>%
        arrange(p_product) %>%
        head(n_degs) %>%
        pull(GeneID) %>%
        enframe %>%
        rename(gene=value, rank=name) %>%
        mutate(tag = tag)
    top_degs
}


get_top_deg_panel = function(combined_tmd_ngo, n_genes){
    panel = combined_tmd_ngo %>%
        filter(rank <= n_genes) %>%
        select(-c(tag, rank)) %>%
        distinct %>%
        pull(gene)
    panel
}


make_nh_restored_tibble = function(results_fgf1, results_bl6, comparison=""){
    # find_restored_nh
    results_fgf1 = results_fgf1 %>% select(labels, polarity, Nhood, SpatialFDR)
    results_bl6 = results_bl6 %>% select(labels, polarity, Nhood, SpatialFDR)
    avb = left_join(results_fgf1, results_bl6, by = c("labels", "Nhood"),
                    suffix = c(".fgf1", ".BL6")) %>%
    dplyr::select(order(colnames(.))) %>%
    relocate(labels, Nhood) %>%
    mutate(restored = case_when(((SpatialFDR.fgf1 < 0.1) & # treated - veh
                                 (SpatialFDR.BL6 < 0.1) & # condition - obob
                                 (polarity.fgf1 != 'none') &
                                 (polarity.BL6 != 'none') &
                                 (polarity.fgf1 != polarity.BL6)) ~ TRUE,
                                TRUE ~ FALSE))
    if (comparison != ""){
        avb['comparison'] = comparison
    }
    avb
}


make_nh_restored_summary_tibble = function(avb, comparison=""){
    nh_restored = avb %>%
        filter(restored == TRUE) %>%
        group_by(labels) %>%
        summarise(n_restored = n()) %>% arrange(desc(n_restored))
    nh_restored.pos = avb %>%
        filter(restored == TRUE) %>%
        group_by(labels) %>%
        filter(polarity.fgf1 == 'pos') %>%
        summarise(n_restored.pos = n())
    nh_restored.neg = avb %>%
        filter(restored == TRUE) %>%
        group_by(labels) %>%
        filter(polarity.fgf1 == 'neg') %>%
        summarise(n_restored.neg = n())


    total_summary = avb %>%
        group_by(labels) %>%
        summarise(n_total = n())

    fgf1_summary = avb %>%
        group_by(labels) %>%
        filter(SpatialFDR.fgf1 < 0.1) %>%
        summarise(n_fgf1 = n())
    fgf1_summary.pos = avb %>%
        group_by(labels) %>%
        filter(SpatialFDR.fgf1 < 0.1) %>%
        filter(polarity.fgf1 == 'pos') %>%
        summarise(n_fgf1.pos = n())
    fgf1_summary.neg = avb %>%
        group_by(labels) %>%
        filter(SpatialFDR.fgf1 < 0.1) %>%
        filter(polarity.fgf1 == 'neg') %>%
        summarise(n_fgf1.neg = n())

    bl6_summary = avb %>%
        group_by(labels) %>%
        filter(SpatialFDR.BL6 < 0.1) %>%
        summarise(n_BL6 = n())

    total_summary = total_summary %>% 
        full_join(fgf1_summary) %>%
        full_join(fgf1_summary.pos) %>%
        full_join(fgf1_summary.neg) %>%
        mutate(prop_fgf1.pos = n_fgf1.pos/n_fgf1) %>%
        full_join(bl6_summary) %>%
        mutate(n_fgf1 = replace_na(n_fgf1, 0)) %>%
        mutate(n_BL6 = replace_na(n_BL6, 0))

    total_summary = left_join(total_summary, nh_restored) %>%
        left_join(nh_restored.pos) %>%
        left_join(nh_restored.neg) %>%
        mutate(n_restored = replace_na(n_restored, 0)) %>%
        mutate(n_restored.pos = replace_na(n_restored.pos, 0)) %>%
        mutate(n_restored.neg = replace_na(n_restored.neg, 0)) %>%
        mutate(prop_restored_BL6 = n_restored/n_BL6) %>%
        mutate(prop_restored_FGF1 = n_restored/n_fgf1) %>%
        arrange(desc(n_fgf1))

    if (comparison != ""){
        total_summary["comparison"] = comparison
    }

    total_summary

}


annotate_summary_groupings = function(da_results, restored_df){
    da_results = left_join(da_results, restored_df, by =c("labels", "Nhood")) %>%
        mutate(exact_grouping = case_when( (polarity.fgf1 == 'pos') & (polarity.BL6 == 'neg') ~ 'pos_restored',
                                           (polarity.fgf1 == 'neg') & (polarity.BL6 == 'pos') ~ 'neg_restored',
                                           (polarity.fgf1 == 'pos') & (polarity.BL6 == 'none') ~ 'pos_FGF1',
                                           (polarity.fgf1 == 'neg') & (polarity.BL6 == 'none') ~ 'neg_FGF1',
                                           (polarity.fgf1 == 'none') & (polarity.BL6 == 'none') ~ 'none',
                                           (polarity.fgf1 == 'none') & (polarity.BL6 == 'pos') ~ 'pos_BL6',
                                           (polarity.fgf1 == 'none') & (polarity.BL6 == 'neg') ~ 'neg_BL6',
                                           (polarity.fgf1 == 'pos') & (polarity.BL6 == 'pos') ~ 'pos_away',
                                           (polarity.fgf1 == 'neg') & (polarity.BL6 == 'neg') ~ 'neg_away'))%>%
        mutate(restored_grouping = case_when(SpatialFDR.fgf1 < 0.1 ~ exact_grouping,
                                             TRUE ~ 'none')) %>%
        mutate(fgf1_grouping = polarity.fgf1) %>%
        mutate(bl6_grouping = polarity.BL6) %>%
        select(all_of(c(colnames(da_results), c("comparison", "restored", 'exact_grouping', 'restored_grouping', 'fgf1_grouping', 'bl6_grouping'))))
    da_results
}


add_gsea_cols_to_deg_results = function(deg_results){
    converted = gprofiler2::gconvert(query = deg_results$GeneID,
                                     organism = "mmusculus",
                                     target = "ENSG",
                                     mthreshold = 1,
                                     filter_na = FALSE) %>%
                    mutate(ensmusg = target) %>%
                    mutate(GeneID = input) %>%
                    select(GeneID, ensmusg)
    deg_results = deg_results %>%
        left_join(converted) %>%
        distinct(ensmusg, .keep_all=TRUE) %>%
        mutate(gsea_sort_score = (-log10(.[[3]]))) %>%
        mutate(gsea_sort_score = .[[2]] * gsea_sort_score) %>%
        arrange(desc(gsea_sort_score))
    deg_results
}

deg2vec = function(deg_results){
    deg_gsea_vec = deg_results %>%
        filter(str_detect(ensmusg, "ENSMUSG")) %>%
        dplyr::select(ensmusg, gsea_sort_score) %>%
        distinct(gsea_sort_score, .keep_all = TRUE) %>%
        deframe
    deg_gsea_vec
}

do_gseGO = function(deg_results){
    deg_gsea_vec = deg2vec(deg_results)
    gsego_result <- clusterProfiler::gseGO(geneList     = deg_gsea_vec,
                                            OrgDb        = org.Mm.eg.db::org.Mm.eg.db,
                                            ont          = "ALL",
                                            keyType      = "ENSEMBL",
                                            minGSSize    = 15,
                                            maxGSSize    = 500,
                                            pvalueCutoff = 0.05,
                                            verbose      = FALSE)
    gsego_result
}

run_all_fgsea = function(deg_results, padj_cutoff=0.10, tag=''){
    deg_gsea_vec = deg2vec(deg_results)
    all_gene_sets = msigdbr::msigdbr(species = "Mus musculus")
    gset_desc = all_gene_sets %>%
        dplyr::select(gs_cat, gs_subcat, gs_name, gs_description) %>%
        mutate(pathway = gs_name) %>%
        dplyr::select(-gs_name) %>%
        distinct
    test_list = all_gene_sets %>%
        group_by(gs_subcat) %>%
        group_split %>%
        purrr::map(~split(x = .x$ensembl_gene, f = .x$gs_name))
    fgseaResults = test_list %>%
        purrr::map(~fgsea::fgsea(pathways = .x, 
                                  stats    = deg_gsea_vec,
                                  minSize  = 15,
                                  maxSize  = 500))
    merged_results = do.call("rbind", fgseaResults) %>%
        left_join(gset_desc) %>%
        relocate(gs_cat, gs_subcat, pathway, gs_description) %>%
        filter(padj < padj_cutoff)
    if (length(tag) > 0){
        merged_results$tag = tag
    }
    merged_results
}


annotate_nhg = function(da_results_nhg){
    nhg_annotation = da_results_nhg %>%
        select(Nhood, restored_grouping, bl6_grouping) %>%
        distinct %>%
        mutate(Nhood = as.character(Nhood))
    nhg_annotation
}


summarise_nhg_annotation = function(nhg_annotation){
    grouping_summary = nhg_annotation %>% 
        group_by(restored_grouping) %>% 
        summarise(n = n()) %>%
        ungroup() %>%
        mutate(frac_nhoods = n/sum(n))
    grouping_summary
}


nhg2cell = function(nhm, da_results_nhg) {
    nhg_annotation = annotate_nhg(da_results_nhg)
    grouping_summary = summarise_nhg_annotation(nhg_annotation)
    nhg_tib = nhm %>% 
        rownames_to_column %>% 
        pivot_longer(cols = !contains("row")) %>%
        dplyr::rename(Nhood = name) %>%
        filter(value != 0) %>%
        select(-value) %>%
        left_join(nhg_annotation, by = "Nhood") %>%
        left_join(grouping_summary, by="restored_grouping") %>%
        select(-n) %>%
        select(-Nhood) %>%
        mutate(restored_grouping = as.factor(restored_grouping)) %>%
        group_by(rowname, restored_grouping) %>%
        mutate(group_weight = n() * (1-frac_nhoods)) %>%
        ungroup() %>%
        group_by(rowname) %>%
        mutate(total_count = n()) %>%
        ungroup() %>%
        distinct %>%
        mutate(weight = group_weight/total_count) %>%
        arrange(desc(weight)) %>%
        distinct(rowname, .keep_all = TRUE) %>%
        select(-group_weight, -total_count, -weight, -frac_nhoods) %>%
        mutate(fgf1_grouping = case_when((str_detect(restored_grouping, "pos") & 
                                          !str_detect(restored_grouping, "BL6")) ~ "pos",
                                         (str_detect(restored_grouping, "neg") & 
                                          !str_detect(restored_grouping, "BL6")) ~ "neg",
                                         TRUE ~ "none")
              )
    nhg_tib
}


get_seurat_nhg_markers = function(seurat_obj, nhgc, grouping_col, group_a, group_b='', tag=''){
    nhgc['grouping'] = nhgc[grouping_col]
    group_a = stringr::str_split(group_a, pattern=fixed('.')) %>% unlist
    cells_a = nhgc %>%
        filter(grouping %in% group_a) %>%
        pull(rowname)
    if (group_b == ''){
        group_b = nhgc %>%
            filter(!(grouping %in% group_a)) %>%
            pull(grouping) %>%
            as.character %>%
            unique %>%
            paste0(collapse='.')
    }
    group_b = stringr::str_split(group_b, pattern=fixed('.')) %>% unlist
    cells_b = nhgc %>%
        filter(grouping %in% group_b) %>%
        pull(rowname)
    markers = Seurat::FindMarkers(seurat_obj, ident.1=cells_a, ident.2=cells_b, slot="data", assay="SCT", verbose=TRUE,
                                  min.cells.group = 10, 
                                  min.cells.feature = 10,
                                  min.pct = 0.01,
                                  logfc.threshold = 0,
                                  only.pos = FALSE) 
    markers['tag'] = tag
    markers
}


add_gsea_cols_to_seurat_marker_results = function(sm_results){
    sm_results = sm_results %>% 
        rownames_to_column(var = "GeneID")
    converted = gprofiler2::gconvert(query = sm_results$GeneID,
                                     organism = "mmusculus",
                                     target = "ENSG",
                                     mthreshold = 1,
                                     filter_na = FALSE) %>%
                    mutate(ensmusg = target) %>%
                    mutate(GeneID = input) %>%
                    select(GeneID, ensmusg)
    sm_results = sm_results %>%
        left_join(converted) %>%
        distinct(ensmusg, .keep_all=TRUE) %>%
        mutate(gsea_sort_score = -log10(p_val_adj) * avg_log2FC) %>%
        arrange(desc(gsea_sort_score))
    sm_results
}


rehydrate_deg_tag = function(df){
    df %>%
        rowwise %>%
        mutate(tag_split = str_split(tag, '___')) %>%
        mutate(cluster = tag_split[[1]]) %>%
        mutate(data_day = case_when(str_detect(cluster, "Day5") ~ 'Day5',
                                    str_detect(cluster, "Day14") ~ 'Day14',
                                    TRUE ~ 'all')) %>%
        mutate(cluster = str_replace(cluster, fixed('Day5.'), '')) %>%
        mutate(cluster = str_replace(cluster, fixed('Day14.'), '')) %>%
        mutate(comparison = tag_split[[2]]) %>%
        mutate(grouping = tag_split[[3]]) %>%
        select(-tag_split) %>%
        mutate(comparison_split = str_split(comparison, '__v__')) %>%
        mutate(fgf1_conditions = comparison_split[[1]]) %>%
        mutate(bl6_conditions = comparison_split[[2]]) %>%
        select(-comparison, -comparison_split) %>%
        mutate(fgf1_split = str_split(fgf1_conditions, fixed('.'))) %>%
        mutate(fgf1_day = fgf1_split[[1]]) %>%
        mutate(fgf1_comparison = fgf1_split[[2]]) %>%
        select(-fgf1_conditions, -fgf1_split) %>%
        mutate(bl6_split = str_split(bl6_conditions, fixed('.'))) %>%
        mutate(bl6_day = bl6_split[[1]]) %>%
        mutate(bl6_comparison = bl6_split[[2]]) %>%
        select(-bl6_conditions, -bl6_split) %>%
        mutate(grouping_split = str_split(grouping, fixed("."))) %>%
        mutate(grouping = grouping_split[[1]]) %>%
        mutate(cell_comparison = grouping_split[[2]]) %>%
        select(-grouping_split) %>%
        mutate(cell_split = str_split(cell_comparison, '_vs_')) %>%
        mutate(cells_a = cell_split[[1]]) %>%
        mutate(cells_b = cell_split[[2]]) %>%
        select(-cell_split, -cell_comparison) %>%
        relocate(tag, data_day, cluster, fgf1_day, fgf1_comparison, bl6_day, bl6_comparison, grouping, cells_a, cells_b) %>%
        ungroup
}


make_gost = function(degs){
    degs = degs %>% filter(p_val_adj < 0.05)
    gene_list_all = degs %>% pull(ensmusg)
    gene_list_pos = degs %>% filter(avg_log2FC > 0) %>% pull(ensmusg)
    gene_list_neg = degs %>% filter(avg_log2FC < 0) %>% arrange(gsea_sort_score) %>% pull(ensmusg)
    gene_list = list(all = gene_list_all,
                     up = gene_list_pos,
                     down = gene_list_neg)
    gost_results = gprofiler2::gost(query = gene_list,
                                    organism = "mmusculus",
                                    ordered_query = TRUE,
                                    exclude_iea = FALSE, 
                                    evcodes = TRUE)
    gost_results$result = gost_results$result %>%
    rowwise %>%
    mutate(gene_ids = str_split(intersection, fixed(','))) %>%
    mutate(gene_ids = {degs %>% 
                       filter(ensmusg %in% gene_ids) %>% 
                       pull(GeneID) %>% 
                       paste0(collapse=',')}) %>%
    ungroup
    gost_results
}


make_milo_index_tibble = function(milo_obj, label){
    index_hash = milo_obj %>%
        `@`('nhoodIndex') %>% unlist %>%
        colnames(milo_obj)[.] 
    nhood_ncells = milo_obj@nhoodCounts %>% as.matrix %>% rowSums    
    index_tibble = tibble(labels=label, Nhood=1:length(index_hash), hash.mcl.ID=index_hash, n_cells = nhood_ncells)
    index_tibble
}