prep_obj_for_milo_cb_v01 = function(obj, set_orig.batch=TRUE){
    obj@meta.data = obj@meta.data %>% mutate(hash.mcl.ID = hash.mcl.ID_SCOP)
    if (set_orig.batch){
        obj@meta.data = obj@meta.data %>% mutate(orig.batch = batch)
    }
    obj@meta.data = obj@meta.data %>% mutate(batch = stringr::str_replace_all(batch, stringr::fixed(" "), '__'))
    obj@meta.data = obj@meta.data %>% 
        mutate(labels = stringr::str_replace_all(labels, stringr::fixed("-"), '__')) %>%
        mutate(labels = stringr::str_replace_all(labels, stringr::fixed("/"), '__')) %>%
        mutate(labels = stringr::str_replace_all(labels, stringr::fixed("("), '__')) %>%
        mutate(labels = stringr::str_replace_all(labels, stringr::fixed(")"), '__'))
    obj@meta.data$group = interaction(obj@meta.data$treatment, 
                                      obj@meta.data$time, 
                                      obj@meta.data$strain, 
                                      drop = TRUE)
    obj
}


reset_orig.batch = function(obj){
    obj@meta.data = obj@meta.data %>% mutate(batch = orig.batch)
    obj
}


set_batch_to_lane = function(obj){
    obj@meta.data = obj@meta.data %>% mutate(batch = Index.10x_SCOP)
    obj
}

set_labels_to_class = function(obj){
    obj@meta.data = obj@meta.data %>% mutate(labels = class)
    obj
}

set_labels_to_lvl1 = function(obj){
    obj@meta.data = obj@meta.data %>% mutate(labels = labels_lvl1)
    obj
}

set_labels_to_lvl2 = function(obj){
    obj@meta.data = obj@meta.data %>% mutate(labels = labels_lvl2)
    obj
}

set_labels_to_labels_chunk = function(obj){
    obj@meta.data = obj@meta.data %>% mutate(labels = labels_chunk)
    obj
}



process_seurat_with_batch_tryCatch = function(obj, batch_to_correct){
    tryCatch({
                # Attempt to process seurat k.weight 100
                obj %>% process_seurat(method = "integrate", 
                            batch = batch_to_correct,
                            dims = 30, res = 0.8,
                            k.weight = 100)
              }, error = function(e) {
                    obj %>% Seurat::SCTransform(assay='RNA',
                                               method="glmGamPoi",
                                               vars.to.regress= batch_to_correct,
                                               vst.flavor="v2",
                                               verbose=TRUE) %>%
                           run_sct_chaser
              })
}

run_until_success <- function(...){
  funcs <- list(...)
  for (func in funcs) {
    result <- tryCatch({
      func()
    }, error = function(e){
      return(NULL)
    })
    if (!is.null(result)) {
      return(result)
    }
  }
  stop("All functions resulted in an error")
}





process_seurat_with_batch_tryCatch = function(obj, batch){
    tryCatch({
                # Attempt to process seurat k.weight 100
                obj %>% process_seurat(method = "integrate", 
                            batch = batch,
                            dims = 30, res = 0.8,
                            k.weight = 100)
              }, error = function(e) {
                #
                tryCatch({
                    obj %>% process_seurat(method = "integrate", 
                            batch = batch,
                            dims = 30, res = 0.8,
                            k.weight = 40)
                }, error = function(e){
                   obj %>% Seurat::SCTransform(assay='RNA',
                                               method="glmGamPoi",
                                               vars.to.regress= batch,
                                               vst.flavor="v2",
                                               verbose=TRUE) %>%
                           run_sct_chaser
                })
                
                
              })
}




chi_squared_test <- function(batches, all_batches) {
  # Calculate the table of observed frequencies for batches
  observed <- table(batches)

  # Calculate the table of frequencies for all_batches
  all_batches_freq <- table(all_batches)

  # Calculate the expected frequencies
  expected <- all_batches_freq * length(batches) / length(all_batches)

  # Match the names in observed and expected
  # Assign 0 to categories in expected that are not in observed
  expected <- expected[names(all_batches_freq)]
  names(expected) <- names(all_batches_freq)
  observed <- observed[names(all_batches_freq)]
  observed[is.na(observed)] <- 0

  # Perform the Chi-Squared test
  test_result <- chisq.test(x = observed, p = expected / sum(expected))

  # Return the test result
  return(test_result)
}


calc_batchy_score = function(obj=obj,
                        batch = batch,
                        dims = dims,
                        reduction = 'pca',
                        k.param = 40){

    obj = obj %>% FindNeighbors(reduction = reduction, 
                            dims = seq(dims), 
                            return.neighbor = TRUE, 
                            k.param = k.param, 
                            graph.name = c('nn_graph'))
    all_batches = obj[[batch]]

    df = Cells(obj) %>% 
        tibble(barcode=.) %>%
        rowwise %>%
        mutate(top_neighbors = list(TopNeighbors(obj@neighbors$nn_graph, barcode, k.param))) %>%
        mutate(batches = obj %>% 
                       `[[` %>%
                       filter(row.names(.) %in% unlist(top_neighbors)) %>%
                       `[[`(batch) %>%
                       list) %>%
        mutate(batch_xsqp = chi_squared_test(batches %>% unlist, all_batches) %>%
                            `$`('p.value')) %>%
        ungroup
    
    df_meta = df %>%
        select(barcode, batch_xsqp) %>%
        column_to_rownames(var = 'barcode') %>%
        rowwise %>%
        mutate(logp_bxsqp = -log10(batch_xsqp)) %>%
        ungroup

    obj = obj %>% AddMetaData(df_meta)
    obj
}