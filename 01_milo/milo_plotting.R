make_kde2d = function(nhgc_coords){
    range_factor = 0.25
    x = nhgc_coords$UMAP_1
    y = nhgc_coords$UMAP_2
    x_min = min(x)
    x_max = max(x)
    x_range = x_max - x_min
    x_min = x_min - (range_factor * x_range)
    x_max = x_max + (range_factor * x_range)
    y_min = min(y)
    y_max = max(y)
    y_range = y_max - y_min
    y_min = y_min - (range_factor * y_range)
    y_max = y_max + (range_factor * y_range)
    kde = MASS::kde2d(x, y, n=100, lims = c(x_min, x_max, y_min, y_max))
    kde$labels = nhgc_coords$labels
    kde    
}


enframe_kde = function(kde){
    kde_df = data.frame(expand.grid(x=kde$x, y=kde$y, 
                                    labels=kde$labels), 
                                    z=as.vector(kde$z)) %>% distinct
    kde_df
}


filter_obj_to_nhgc = function(nhgc, obj){
    subset_cells = nhgc %>% pull(rowname)
    obj = subset(obj, cells = subset_cells)
    obj
}


assign_colors = function(tib, categories_in){
    tib %>%
        mutate(group_colors = case_when(labels == 'pos_restored' ~ "#2166ac",
                                        labels == 'neg_restored' ~ "#b2182b",
                                        labels == 'pos_FGF1' ~ "#4393c3",
                                        labels == 'neg_FGF1' ~ "#d6604d",
                                        labels == 'none' ~ "#707070",
                                        labels == 'pos_BL6' ~ "#d1e5f0",
                                        labels == 'neg_BL6' ~ "#f4a582",
                                        labels == 'pos_away' ~ "#92c5de",
                                        labels == 'neg_away' ~ "#fddbc7",
                                        labels == 'pos' ~ "#2166ac",
                                        labels == 'neg' ~ "#b2182b")) %>%
        mutate(group_colors = case_when(labels %in% categories_in ~ group_colors,
                                        TRUE ~ "#707070")) %>%
        mutate(labels = factor(labels, #level and order for plotting
                               levels=c("none", 
                                        "pos_away", "neg_away",
                                        "pos_BL6", "neg_BL6", 
                                        "pos_FGF1", "neg_FGF1", 
                                        "pos", "neg", 
                                        "pos_restored", "neg_restored"))) %>%
        arrange(labels)
}


make_nhgc_coords = function(nhgc, obj, grouping_col){
    nhgc$labels = nhgc[[grouping_col]]
    umap_coords = obj@reductions$umap@cell.embeddings %>%
        as.data.frame %>%
        rownames_to_column
    nhgc_coords = nhgc %>%
        left_join(umap_coords, by = "rowname")  %>%
        assign_colors(c('pos_restored', 'neg_restored', 'none', "pos", "neg"))
    nhgc_coords
}

make_kde_df = function(nhgc_coords){
    kdes = nhgc_coords %>%
        group_by(labels) %>%
        group_map(~ make_kde2d(.x), .keep=TRUE) %>%
        map(~ enframe_kde(.x))
    kde_df = do.call(rbind, kdes) %>%
        relocate(x, y, z) %>%
        group_by(labels) %>%
        mutate(break_min_val = quantile(z, 0.75)) %>%
        ungroup %>%
        filter(labels %in% c('pos_restored', 'neg_restored', 'none', "pos", "neg")) %>%
        assign_colors(c('pos_restored', 'neg_restored', 'none', "pos", "neg")) %>%
        dplyr::rename(contour_colors = group_colors)
    kde_df
}


add_contours = function(gg, kde_df, line_size=1, label_size=5){
    for (tag in unique(kde_df$labels)){
        kde_df_f = kde_df %>% filter(labels == tag)
        contour_color = kde_df_f %>% pull(contour_colors) %>% unique
        break_min_val = kde_df_f$break_min_val %>% unique
        gg = gg + metR::geom_contour2(aes(x=x, y=y, z=z, label=labels), data= kde_df_f, 
                                      breaks=c(break_min_val), 
                                      colour=contour_color,
                                      size=line_size, label_size=label_size, 
                                      fontface='bold')
    }
    gg
}


add_more_contours = function(gg, kde_df, line_size=1, label_size=5, drop_none=TRUE){
    if (drop_none == TRUE){
        kde_df = kde_df %>% filter(labels != 'none')
    }
    for (tag in unique(kde_df$labels)){
        kde_df_f = kde_df %>% filter(labels == tag)
        contour_color = kde_df_f %>% pull(contour_colors) %>% unique
        contour_breaks = quantile(kde_df$z , c(0.90, 0.95, 0.99))
        gg = gg + metR::geom_contour2(aes(x=x, y=y, z=z), data= kde_df_f, 
                                      breaks=contour_breaks, 
                                      colour=contour_color,
                                      size=line_size)
    }
    gg
}


add_contours_no_label = function(gg, kde_df, line_size=1, label_size=5){
    for (tag in unique(kde_df$labels)){
        kde_df_f = kde_df %>% filter(labels == tag)
        contour_color = kde_df_f %>% pull(contour_colors) %>% unique
        break_min_val = kde_df_f$break_min_val %>% unique
        gg = gg + metR::geom_contour2(aes(x=x, y=y, z=z), data= kde_df_f, 
                                      breaks=c(break_min_val), 
                                      colour=contour_color,
                                      size=line_size)
    }
    gg
}


plot_group_overview = function(obj, nhgc, grouping_col, plot_title){
    nhgc_coords = make_nhgc_coords(nhgc, obj, grouping_col)
    color_scale = nhgc_coords %>% select(labels, group_colors) %>% distinct %>% deframe
    kde_df = make_kde_df(nhgc_coords)
    gg = ggplot() + 
        labs(title=plot_title) +
        geom_point(data=nhgc_coords, aes(x = UMAP_1, y = UMAP_2, color = labels), alpha=0.2) + 
        scale_color_manual(values=color_scale) +
        theme_void() + 
        theme(legend.position="none") + 
        theme(plot.title = element_text(size=20,face = "bold",
              margin = margin(b = -70, t=70), hjust = 0.03)) +
        theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
    gg = add_contours(gg, kde_df)
    gg = add_more_contours(gg, kde_df)
    gg
}


make_fp_plot = function(obj, value){
    list(FeaturePlot(obj,
                slot='counts',
               cols = c("#eaeaea", "#323232"),
               features = c(value),
               pt.size=8,
               order=TRUE,
               min.cutoff="q01",
               max.cutoff="q99",
               raster=TRUE,
               raster.dpi=c(1024, 1024),
               ncol=1
               ) + 
                 xlim(xl[1], xl[2]) +
                 ylim(yl[1], yl[2]) +
                 theme_void() + 
                 theme(legend.position="none") + 
                 theme(plot.title = element_text(size=25,face = "bold",
                     margin = margin(b = -30, t=30), hjust = 0.03)) +
                 theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
}


make_plot_data = function(obj, nhgc_coords, gene){
    gene_vals = obj %>% `@`('assays') %>% `$`('SCT') %>% `@`('counts') %>% 
        `[`(gene,) %>% 
        enframe(name = 'rowname')
    max_cutoff_q = gene_vals %>% pull(value) %>% quantile(0.99)
    min_cutoff_q = gene_vals %>% pull(value) %>% quantile(0.01)
    plot_data = nhgc_coords %>% 
    left_join(gene_vals, by = 'rowname') %>%
    mutate(value = case_when(value > max_cutoff_q ~ max_cutoff_q,
                             value < min_cutoff_q ~ min_cutoff_q,
                             TRUE ~ value)) %>%
    arrange(value)
}

make_featureplot = function(obj, nhgc, grouping_col, gene){
    nhgc_coords = make_nhgc_coords(nhgc, obj, grouping_col)
    plot_data = make_plot_data(obj, nhgc_coords, gene)
    kde_df = make_kde_df(nhgc_coords)
    gg = ggplot() + 
        geom_point(data=plot_data, aes(x = UMAP_1, y = UMAP_2, colour=value)) +
        scale_colour_gradient(low = "#cce7d7", high = "#006f2d")
    gg = add_contours_no_label(gg, kde_df)
    gg = gg + labs(title=gene)
    gg = gg + 
         theme_void() + 
         theme(legend.position="none") + 
         theme(plot.title = element_text(size=25,face = "bold",
               margin = margin(b = -30, t=30), hjust = 0.03)) +
         theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
    gg
}


make_summary_layout = function(){
    layout <- c(
        area(103, 0, 647, 544),
        ###top 4.1
        area(26, 576, 207, 757),
        area(26, 757, 207, 938),
        area(26, 938, 207, 1119),
        area(26, 1119, 207, 1300),
        ## top 4.2
        area(207, 576, 388, 757),
        area(207, 757, 388, 938),
        area(207, 938, 388, 1119),
        area(207, 1119, 388, 1300),
        ### bottom 4.1
        area(420, 576, 601, 757),
        area(420, 757, 601, 938),
        area(420, 938, 601, 1119),
        area(420, 1119, 601, 1300),
        area(601, 576, 782, 757),
        area(601, 757, 782, 938),
        area(601, 938, 782, 1119),
        area(601, 1119, 782, 1300))
    layout
}


make_summary_deg_plot = function(obj, nhgc, degs, grouping_col, title){
    obj = filter_obj_to_nhgc(nhgc, obj)
    grouping_overview = plot_group_overview(obj, nhgc, grouping_col, title)
    up_markers = degs %>% pull(GeneID) %>% head(8)
    down_markers = degs %>% arrange(gsea_sort_score) %>% pull(GeneID) %>% head(8)
    markers = c(up_markers, down_markers)
    plot_list = markers %>% 
        enframe(value = 'GeneID') %>%
        rowwise() %>%
        mutate(fp = list(make_featureplot(obj, nhgc, grouping_col, GeneID))) %>%
        ungroup %>%
        pull(fp)
    plot_list = append(list(grouping_overview), plot_list)
    layout = make_summary_layout()
    gg_wrap = wrap_plots(plot_list, design = layout)
    gg_wrap
}


save_summary_plot = function(summary_plot, name_components){
    dirpath = 'outputs/plots/deg_summary/'
    dir.create(dirpath, recursive = TRUE, showWarnings=FALSE)
    name_components = c(dirpath, name_components, '.pdf')
    plot_path = paste0(name_components, collapse = '')
    ggsave(plot_path,
       plot=summary_plot,
       device='pdf',
       width=26.66,
       height=15,
       units=c("in"))
    plot_path
}


make_and_save_summary_plot = function(obj, nhgc, degs, grouping_col, title, name_components, sct_type, deg_output_suffix){
    result = tryCatch({
        make_summary_deg_plot(obj, nhgc, degs, grouping_col, title) %>%
        save_summary_plot(., c(sct_type, deg_output_suffix))
    },
    error = function(e) {
        return(NULL)
    })
}