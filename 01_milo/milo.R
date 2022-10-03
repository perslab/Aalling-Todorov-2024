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