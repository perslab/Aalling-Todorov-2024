import scanpy as sc
import spapros as sp
import numpy as np
import pandas as pd
from matplotlib.pyplot import figure
import pickle
import os
import scipy

figure(figsize=(15, 10), dpi=300)
adata = sc.read_h5ad("/scratch/nmq407/ds_obj_sce_h5ad_from_notebook")


adata_00 = adata[adata.obs_names, adata.var_names]
adata_00.obs['celltype'] = adata.obs['polar_label']
adata_00.obsm['X_umap'] = adata.obsm['UMAP'].to_numpy()
adata_00.X = scipy.sparse.csr_matrix(adata_00.X)
sc.pp.filter_genes(adata_00, min_cells=1)
sc.pp.highly_variable_genes(adata_00, flavor="seurat_v3",n_top_genes=5000)


selector = sp.se.ProbesetSelector(adata_00,
                                  n=100,
                                  celltype_key="celltype",
                                  verbosity=1,
                                  save_dir='run_spapros_100/',
                                  n_jobs=40)

selector.select_probeset()

selected_set = selector.probeset.index[selector.probeset.selection].tolist()