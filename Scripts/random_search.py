import pandas as pd
import os
import sys
import anndata as ad


directory = "/hpc/shared/onco_janssen/dhaynessimmons/projects/ageing_flies/data/ct_specific"

for file in os.listdir(directory):
    if file.endswith(".h5ad"):
        print(file)
        adata = ad.read_h5ad(os.path.join(directory, file))
        print(adata.obs.iloc[-3:].head(), '\n\n')