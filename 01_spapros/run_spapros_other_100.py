import scanpy as sc
import spapros as sp
import numpy as np
import pandas as pd
from matplotlib.pyplot import figure
import pickle
import os

#sc.settings.verbosity = 0
sc.logging.print_header()
print(f"spapros=={sp.__version__}")

figure(figsize=(15, 10), dpi=300)
adata = sc.read_h5ad("/scratch/nmq407/ds_obj_sce_other_h5ad_from_notebook")


adata_00 = adata[adata.obs_names, adata.var_names]
adata_00.obs['celltype'] = adata.obs['polar_label']
adata_00.obsm['X_umap'] = adata.obsm['UMAP'].to_numpy()
sc.pp.highly_variable_genes(adata_00,flavor="cell_ranger",n_top_genes=5000)


selector = sp.se.ProbesetSelector(adata_00,
                                  n=100,
                                  celltype_key="polar_label",
                                  verbosity=1,
                                  save_dir='run_spapros_100_other/',
                                  n_jobs=40)

selector.select_probeset()

selected_set = selector.probeset.index[selector.probeset.selection].tolist()

evaluator = sp.ev.ProbesetEvaluator(adata_00,
                                    metrics=['cluster_similarity', 'knn_overlap', 'forest_clfs', 'marker_corr', 'gene_corr'],
                                    verbosity=2, 
                                    results_dir='run_spapros_100_other_eval/')
evaluator.evaluate_probeset(selected_set, set_id="Spapros")
reference_sets = sp.se.select_reference_probesets(adata_00, n=100)

for set_id, df in reference_sets.items():
    gene_set = df[df["selection"]].index.to_list()
    evaluator.evaluate_probeset(gene_set, set_id=set_id)

