library(tidyverse)

# grab the matching Seurat object
obj = qs::qread('_targets/objects/obj_Day5.Agrp')
obj

# (n)eighbor(h)ood (ma)atrix is an object which tells you which cells belong to which neighborhoods
nhm = qs::qread('_targets/objects/nhm_Day5.Agrp___obob5v5')
nhm %>%  `[`(1:10, 1:20) %>% print

library(Seurat)
library(Matrix)
library(tidyverse)

# Extract Counts Matrix (ensure it's a sparse matrix)
counts_matrix <- GetAssayData(obj, slot = "counts")
# Ensure the matrix is a dgCMatrix for efficient operations
counts_matrix <- as(counts_matrix, "dgCMatrix")

# convert nhm to a dgc matrix so we can multiply it with the seurat counts dgc matrix
nhm_matrix <- as(as.matrix(nhm), "dgCMatrix")

# Matrix Multiplication
# This effectively sums the counts for each nh
summed_counts <- counts_matrix %*% nhm_matrix

# Calculate the Number of Cells per Neighborhood
cells_per_neighborhood <- colSums(nhm_matrix)

# Average Calculation
average_expression <- summed_counts / cells_per_neighborhood

# Convert to Data Frame for easier handling
average_expression_df <- as.data.frame(as.matrix(average_expression))
rownames(average_expression_df) <- rownames(counts_matrix)

# Viewing the result
head(average_expression_df)


degs = qs::qread('_targets/objects/combined_deg_seurat_formatted')

degs %>%
filter(cluster == 'Agrp') %>%
filter(fgf1_comparison == 'obob14v14') %>%
filter(bl6_comparison == 'obobBL6d14') %>%
filter(grouping == 'restored_grouping') %>%
filter(cells_a == 'neg_restored') %>%
filter(cells_b == 'none') %>%
filter(fgf1_day == 'Day14') %>%
filter(bl6_day == 'Day14') %>%
filter(avg_log2FC > 0) %>% # I only counted genes with avg_log2FC > 0 in the other charts %>%
head(10)

# filter to the NMGs for Agrp, Day14, restored_grouping, neg_restored, avg_log2FC >0

these_degs = degs %>%
filter(cluster == 'Agrp') %>%
filter(fgf1_comparison == 'obob14v14') %>%
filter(bl6_comparison == 'obobBL6d14') %>%
filter(grouping == 'restored_grouping') %>%
filter(cells_a == 'neg_restored') %>%
filter(cells_b == 'none') %>%
filter(fgf1_day == 'Day14') %>%
filter(bl6_day == 'Day14') %>%
filter(avg_log2FC > 0) %>% # I only counted genes with avg_log2FC > 0 in the other charts
pull(GeneID)

these_degs %>% head

# da_results_nhg - da results with neighborhood grouping 
# this object gives you what "grouping" neighborhoods fall under depending on the logic conditions met
# there is a restored_grouping, fgf1_grouping, and bl6 grouping column one can select from
da_results_nhg = qs::qread('_targets/objects/da_results_nhg_Agrp___Day5.obob5v5__v__Day5.obobBL6d5')
da_results_nhg %>% head

# I pull all the neg_restored nh from restored_grouping
restored_nhoods = da_results_nhg %>%
filter(restored_grouping == 'neg_restored') %>%
mutate(Nhood = as.character(as.integer(Nhood))) %>% 
pull(Nhood)
restored_nhoods

# here I also pull some restored nh to compare against
none_nhoods = da_results_nhg %>%
filter(restored_grouping == 'none') %>%
mutate(Nhood = as.character(as.integer(Nhood))) %>% 
# sample_n(length(restored_nhoods)) %>% # optionally you can randomly sample an equal amount of nones, or pull all of them by commenting this line out
pull(Nhood) 

none_nhoods

# I select restored nhoods and none_nhoods to become columns in a dataframe later
selected_nhoods = c(restored_nhoods, none_nhoods)

# this is the NMGs x selected_nhoods dataframe 
# on visual inspection, the columns to the left (neg_restored nh) have higher expression of these genes than the columns to the right (none nh)
heatmap_df = average_expression_df[these_degs, selected_nhoods]
heatmap_df %>% dim
heatmap_df %>% head

# we can also average these out to compare them at a glance

neg_restored_exp = average_expression_df[these_degs,restored_nhoods] %>%
rowMeans

none_exp = average_expression_df[these_degs,none_nhoods] %>%
rowMeans

# here is the comparison. It mostly holds - we expect all NMGs to have a positive diff
data.frame(neg_restored_exp, none_exp) %>%
mutate(diff = neg_restored_exp - none_exp) %>%
head(20)

# not all NMGs have a positive difference because this was computed on cells, not neighborhoods
# at one point, I use a voting strategy described in the manuscript to assign cells to neighborhoods exclusively
# that's kinda necessary to make seurat's FindMarkers work
# but we can also examine this on the level of cells that got the neg_restored and pos_restored labels

# nhgc - (n)eighbor(h)ood (g)rouping by (c)ell
# gives the vote outcome for each grouping col, now applied to cells
nhgc = qs::qread('_targets/objects/nhgc_Agrp___Day5.obob5v5__v__Day5.obobBL6d5')
nhgc %>% head

# select cells_a (neg_restored in this case) and cells_b (none)
cells_a = nhgc %>%
filter(restored_grouping == 'neg_restored') %>%
pull(rowname)

cells_b = nhgc %>% 
filter(restored_grouping == 'none') %>%
pull(rowname)



# select cells_a and cells_b
selected_cells = c(cells_a, cells_b)
# make cell_heatmap_df
cell_heatmap_df = counts_matrix[these_degs, selected_cells] %>%
as.data.frame



# let's look at means
exp_a = counts_matrix[these_degs,cells_a] %>% rowMeans
exp_b = counts_matrix[these_degs,cells_b] %>% rowMeans

data.frame(exp_a, exp_b) %>%
mutate(diff = exp_a - exp_b) %>% 
head(20)


