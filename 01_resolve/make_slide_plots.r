library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(forcats)
library(Polychrome)
library(grid)

# this looks for a string in your clusters. something like "Glu" will find all Glu neurons for example
create_grouping_col <- function(seurat_obj, prefix) {
  new_metadata <- case_when(str_detect(seurat_obj@meta.data$predicted.labels, fixed(prefix)) ~ seurat_obj@meta.data$predicted.labels,
                            TRUE ~ NA_character_) %>% as.factor
  new_metadata = new_metadata %>% fct_relevel(sort(na.omit(levels(new_metadata))), after = Inf)
  return(new_metadata)
}

# this plots all cells in a given fov which contain cell_str 
# (e.g. something like "Glu" or "Htr3b" to hit all clusters with that in the name)
plot_celltype_idp = function(seurat_obj, cell_str, fov='fov'){
    new_grouping = seurat_obj %>%
        create_grouping_col(cell_str)
    seurat_obj = AddMetaData(seurat_obj, new_grouping, col.name = 'grouping_col')

    # Get the levels of the lab_Pomc column, including the NA_character_ level
    grouping_levels <- levels(seurat_obj@meta.data$grouping_col)
    # Generate colors using the polychrome palette
    how_many_colors = length(grouping_levels)
    colors <- Polychrome::sky.colors(how_many_colors) %>% as.character %>% `[`(1:how_many_colors) # Excluding the NA_character_ level
    # Add 'gray10' for the NA_character_ level
    colors <- c(colors, "gray10")
    # Create a named vector of colors, with levels as names
    named_colors <- setNames(colors, grouping_levels)
    # Map the lab_Pomc column values to the corresponding colors
    # color_vector <- mapvalues(xenium.obj@meta.data$lab_Pomc, from = lab_Pomc_levels, to = named_colors)
    
    ggp = ImageDimPlot(seurat_obj,
                       group.by='grouping_col',
                       boundaries = 'segmentation',
                       fov = fov,
                       border.size =0.1, 
                       border.color =NA,
                       na.value='gray20', 
                       cols=named_colors, 
                       axes = FALSE,
                       size=0)
    ggp
}

# this saves all fovs in the order fov, fov.1, ..., fov.7 as a pdf
# shaped like this
# fov    fov.4
# fov.1  fov.5
# fov.2  fov.6
# fov.3  fov.7
# make sure to create pdf_folder ahead of time
#
make_idp_pdf = function(seurat_obj, cell_type, pdf_folder=""){
    a1 = plot_celltype_idp(seurat_obj, cell_type, fov='fov')  + theme(legend.position = "none")
    a2 = plot_celltype_idp(seurat_obj, cell_type, fov='fov.1') + theme(legend.position = "none")
    b1 = plot_celltype_idp(seurat_obj, cell_type, fov='fov.2') + theme(legend.position = "none")
    b2 = plot_celltype_idp(seurat_obj, cell_type, fov='fov.3') + theme(legend.position = "none")
    c1 = plot_celltype_idp(seurat_obj, cell_type, fov='fov.4') + theme(legend.position = "none")
    c2 = plot_celltype_idp(seurat_obj, cell_type, fov='fov.5') + theme(legend.position = "none")
    d1 = plot_celltype_idp(seurat_obj, cell_type, fov='fov.6') + theme(legend.position = "none")
    d2 = plot_celltype_idp(seurat_obj, cell_type, fov='fov.7') + theme(legend.position = "none")
    title_theme <- theme(plot.title = element_text(color = "white", size = 8, hjust = 0.5))
    # Set black background theme for each plot
    legend_theme <- theme(legend.position = "bottom",
                          legend.direction = "horizontal",
                          legend.box = "horizontal",
                          legend.text = element_text(size = 4),
                          legend.key.size = unit(2, "mm"))
    # black_background_theme <- theme(plot.background = element_rect(fill = "black"),
    #                                 panel.background = element_rect(fill = "black"))
    black_background_theme <- theme(plot.background = element_rect(fill = "black", color = 'black'),
                          panel.background = element_rect(fill = "black", color = 'black'),
                          plot.margin = margin(0, 0, 0, 0),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank())

    a1_with_title <- a1 + theme(legend.position = "right") + ggtitle("a1_Veh_PF") + title_theme + legend_theme + black_background_theme
    a2_no_legend_with_title <- a2 + guides(scale="none") + ggtitle("a2_Veh_PF") + title_theme + black_background_theme
    b1_no_legend_with_title <- b1 + guides(scale="none") + ggtitle("b1_Veh_PF") + title_theme + black_background_theme
    b2_no_legend_with_title <- b2 + guides(scale="none") + ggtitle("b2_Veh_PF") + title_theme + black_background_theme
    c1_no_legend_with_title <- c1 + guides(scale="none") + ggtitle("c1_FGF1") + title_theme + black_background_theme
    c2_no_legend_with_title <- c2 + guides(scale="none") + ggtitle("c2_FGF1") + title_theme + black_background_theme
    d1_no_legend_with_title <- d1 + guides(scale="none") + ggtitle("d1_FGF1") + title_theme + black_background_theme
    d2_no_legend_with_title <- d2 + guides(scale="none") + ggtitle("d2_FGF1") + title_theme + black_background_theme

    # Combine the plots using patchwork with a black background
    combined_plots <- a1_with_title + c1_no_legend_with_title +
                      a2_no_legend_with_title + c2_no_legend_with_title +
                      b1_no_legend_with_title + d1_no_legend_with_title +
                      b2_no_legend_with_title + d2_no_legend_with_title +
                      plot_layout(guides = "collect", ncol = 2, heights = c(1, 1, 1, 1), widths = c(1, 1))

    # combined_plots = combined_plots & black_background_theme & legend_theme 
    combined_plots = combined_plots + plot_annotation(title=cell_type, theme=black_background_theme) + plot_annotation(theme=legend_theme) + plot_annotation(theme=title_theme)
    # Save the combined_plots to a PDF with a black background
    pdf_name = paste0(pdf_folder, cell_type, '_ct_mapping.pdf')
    pdf(pdf_name, width = 8.27, height = 11.69, bg = "black")
    # grid.newpage()
    grid.draw(combined_plots)
    dev.off()
}

cell_types = xenium.obj %>% `[[` %>% group_by(predicted.labels) %>% summarise(n = n()) %>% arrange(desc(n)) %>% pull(predicted.labels) %>% unique
cell_types

for (cell_type in cell_types){
    print(cell_type)
    make_idp_pdf(xenium.obj, cell_type, pdf_folder="ct_mapping_230502_xenium_sct_unimapped_cca/")
}