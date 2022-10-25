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


.milo_buildGraph = function(milo_obj){
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


make_model_matrix = function(design_df){
    model <- model.matrix(~0 + group + batch, data=design_df)
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


get_nhood_markers = function(milo_obj, da_results, sample_col, gene_subset=NULL, subset_groups=NULL){
    if (!is.null(subset_groups)){
        subset_groups = stringr::str_split(subset_groups, pattern=fixed('.')) %>% unlist
    }
    nhood_markers <- miloR::findNhoodGroupMarkers(milo_obj, da_results, subset.row = gene_subset,
                                                  aggregate.samples = TRUE, sample_col = sample_col,
                                                  subset.groups=subset_groups)
    nhood_markers
}


