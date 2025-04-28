import scanpy as sc
import anndata as ad
import pandas as pd
import os
import sys

# # Load the datasetthe dataset

adata = ad.read_h5ad('D:/OneDrive/dhaynessimmons/OneDrive - UMC Utrecht/Documenten/projects/Ageing flies/Data/adata_body_S_v1.0.h5ad')

# # View basic structure and metadata
print(adata.shape, '\n')

# # Explore expression matrix
# print(adata.X, '\n')


# # Explore metadata for genes (variables)
# print(ada, '\n')

#Print the number of unique values in the metadata
print(adata.obs.nunique(), '\n')

# the unique afca_annotation categories
afca_types = adata.obs['afca_annotation'].unique().to_list()
# print(afca_types, '\n')


# list of selected afca_annotation categories
selected_cell_types = [
    'adult differentiating enterocyte', 'adult midgut enterocyte', 'enterocyte of anterior adult midgut epithelium', 
    'enterocyte of posterior adult midgut epithelium', 'enterocyte-like', 'adult fat body body', 'intestinal stem cell'
    ]

# # Filter the dataset to only include cells from the selected list
adata = adata[adata.obs['afca_annotation'].isin(selected_cell_types)].copy()
print(adata.shape, '\n')

# # Explore metadata for cells (observations)
print(adata.obs.sex_age.unique().to_list(), '\n')
print(adata.obs.sex_age.value_counts(), '\n')

# Remove the samples with age_sex = 'mix_5'
adata = adata[adata.obs['sex_age'] != 'mix_5']
print(adata.shape, '\n')