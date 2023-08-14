#!/usr/bin/env python
# coding: utf-8

# In[1]:


import anndata as ad
import scanpy as sc

import numpy as np
import pandas as pd
import sklearn as sk
import matplotlib.pyplot as plt
import torch

from persist import PERSIST, ExpressionDataset


# In[ ]:





# In[2]:


adata = ad.read_h5ad('/scratch/nmq407/dvc_neurons_rat_20kc_3kvf.h5Seurat.h5ad')
adata


# In[ ]:





# In[ ]:





# In[ ]:





# In[3]:


sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=3000, inplace=True)


# In[4]:


adata


# In[5]:


adata.var['highly_variable']


# In[6]:


adata = adata[:,adata.var['highly_variable']]


# In[7]:


adata.obs['ludwig.predictions']


# In[8]:


adata.obs['ludwig.predictions_codes'] = pd.Categorical(adata.obs['ludwig.predictions']).codes


# In[9]:


adata


# In[10]:


adata.layers['bin'] = (adata.X>0).astype(np.float32)


# In[11]:


print(adata)


# In[12]:


# Choose training and validation splits. 
# You may want to use a different strategy to choose these - see https://scikit-learn.org/stable/modules/classes.html#module-sklearn.model_selection
train_ind, val_ind = sk.model_selection.train_test_split(np.arange(adata.shape[0]), train_size=0.8)

print(f'{adata.shape[0]} total samples')
print(f'{np.size(train_ind)} in training set')
print(f'{np.size(val_ind)} in validation set')

# These are views, so they do not take up memory
adata_train = adata[train_ind,:]
adata_val = adata[val_ind,:]


# In[ ]:


import time
# Get the start time
start_time = time.time()
print(start_time)

# Initialize the dataset for PERSIST
# Note: Here, data_train.layers['bin'] is a sparse array
# data_train.layers['bin'].A converts it to a dense array
train_dataset = ExpressionDataset(adata_train.layers['bin'], adata_train.obs['ludwig.predictions_codes'])
val_dataset = ExpressionDataset(adata_val.layers['bin'], adata_val.obs['ludwig.predictions_codes'])


# Use GPU device if available -- we highly recommend using a GPU!
device = torch.device(torch.cuda.current_device() if torch.cuda.is_available() else 'cpu')

# Number of genes to select within the current selection process.
num_genes = (32, 64, 100)
persist_results = {}

# Set up the PERSIST selector
selector = PERSIST(train_dataset,
                   val_dataset,
                   loss_fn=torch.nn.CrossEntropyLoss(),
                   device=device)
print(device)

# Coarse removal of genes
print('Starting initial elimination...')
candidates, model = selector.eliminate(target=500, max_nepochs=250)
print('Completed initial elimination.')

print('Selecting specific number of genes...')
for num in num_genes:
    inds, model = selector.select(num_genes=num, max_nepochs=250)
    persist_results[num] = inds
print('Done')

# Get the end time
end_time = time.time()
print(time.localtime(end_time))
# Calculate the execution time
execution_time = end_time - start_time

# Format the execution time in a human-readable format
minutes, seconds = divmod(execution_time, 60)
hours, minutes = divmod(minutes, 60)
formatted_time = f"{int(hours)} hours, {int(minutes)} minutes, {int(seconds)} seconds"
print("Execution time:", formatted_time)


# In[19]:


persist_results


# In[ ]:




