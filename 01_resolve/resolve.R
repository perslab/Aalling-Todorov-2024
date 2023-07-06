filter_down_cells = function(xenium.obj){
    keep_cells = xenium.obj %>% `[[` %>%
        filter(cell_area < 7000) %>%
        filter(cell_area > 100) %>%
        filter(avg_confidence >= 0.99) %>%
        filter(nCount_Xenium >= 10) %>%
        filter(nCount_Xenium <= 600) %>%
        rownames
    xenium.obj = subset(xenium.obj, cells = keep_cells)
    xenium.obj
}

sc_transform_resolve = function(xenium.obj, keep_cells){
    Seurat::SCTransform(assay='Xenium',
                        method="glmGamPoi",
                        verbose=TRUE)
    xenium.obj
}

