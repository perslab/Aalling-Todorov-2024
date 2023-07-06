import pandas as pd
from collections import defaultdict
from scipy.io import mmwrite
from scipy.sparse import coo_matrix
import gzip
import requests

# Read the input CSV file
df = pd.read_csv('/data/petar/fgf1/resolve/baysor_segmentation/32810-1377-slide3_A1-1_results_baysor_results/baysor_results.csv')

# Create a dictionary to store gene names and assign them unique indices
genes = defaultdict(lambda: len(genes))

# Assign unique indices to cell names
df['cell_idx'] = df['cell'].astype('category').cat.codes

# Assign unique indices to gene names
df['gene_idx'] = df['gene'].apply(lambda x: genes[x])

# Group the data by cell and gene, then count the molecules for each group
grouped_data = df.groupby(['cell_idx', 'gene_idx']).size().reset_index(name='counts')

# Save the matrix to a Matrix Market format file
rows, cols, data = grouped_data['cell_idx'].values, grouped_data['gene_idx'].values, grouped_data['counts'].values
# sparse_matrix = coo_matrix((data, (rows, cols)), dtype=int)
sparse_matrix = coo_matrix((data, (rows, cols)), shape=(max(rows) + 1, max(cols) + 1)).T
with gzip.open('cell_feature_matrix/matrix.mtx.gz', 'wb') as f:
    mmwrite(f, sparse_matrix, field='integer', precision=None, symmetry='general')

# Query gprofiler2 API to convert gene names to Ensembl gene IDs (ENSG)
unique_genes = list(genes.keys())
r = requests.post(
    url='https://biit.cs.ut.ee/gprofiler/api/convert/convert/',
    json={
        'organism': 'mmusculus',
        'target': 'ENSG',
        'query': unique_genes,
    }
)

converted_genes = {gene: result['converted'] for gene, result in zip(unique_genes, r.json()['result'])}

# Save the gene and cell names to separate files
with gzip.open('cell_feature_matrix/features.tsv.gz', 'wt') as f:
    for gene, idx in sorted(genes.items(), key=lambda x: x[1]):
        ensembl_id = converted_genes[gene]
        f.write(f"{ensembl_id}\t{gene}\tGene Expression\n")

with gzip.open('cell_feature_matrix/barcodes.tsv.gz', 'wt') as f:
    for cell in df['cell'].astype('category').cat.categories:
        f.write(f"{cell}\n")

# Check the dimensions of the output files
print(f"Matrix dimensions: {sparse_matrix.shape}")
print(f"Number of genes: {len(genes)}")
print(f"Number of cells: {len(df['cell'].astype('category').cat.categories)}")
