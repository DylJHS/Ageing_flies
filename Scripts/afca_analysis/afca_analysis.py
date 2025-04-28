# Import the required library before running this script.
import re 
import os
import anndata as ad
import scanpy as sc
import diffxpy.api as de
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------
# Load the filtered dataset from the specified h5ad file.
# ---------------------------------------------------------------------
adata = ad.read_h5ad('D:\OneDrive\dhaynessimmons\OneDrive - UMC Utrecht\Documenten\projects\Ageing flies\Data\adata_body_filtered_2.h5ad')

# ---------------------------------------------------------------------
# Loop over each of the cell types and batches and assess the expression profiles.
# ---------------------------------------------------------------------
# for cell_type in adata.obs.afca_annotation.unique():
#     adata_ctype = adata[adata.obs.afca_annotation == cell_type].copy()
#     print(str(cell_type), '\n',
#           adata_ctype.obs[['batch', 'sex_age']].value_counts().sort_index(),
#           '\n\n\n')
    
    # # Recalculate the umap for the cell type
    # print(adata_ctype.obsm['X_umap'])
    # sc.pp.normalize_total(adata_ctype)
    # sc.pp.log1p(adata_ctype)

    # sc.pp.neighbors(adata_ctype)
    # sc.tl.umap(adata_ctype)

    # # Plot the umap
    # sc.pl.umap(adata_ctype, color="batch", title=f"{cell_type} - Batch")

    # # Save the figure
    # folder = f"/hpc/shared/onco_janssen/dhaynessimmons/Figures/ageing_flies/{cell_type}"
    # os.makedirs(folder, exist_ok=True)
    # plt.savefig(f"{folder}/umap_batch.png")
    # plt.close()


#     for batch in adata_ctype.obs.batch.unique():
#         print(f"Processing batch {batch}")
#         adata_batch = adata_ctype[adata_ctype.obs.batch == batch].copy()

#         # Make violin plots of n_genes_by_counts, total_counts, and total_counts_mt using sns
#         ax = sc.pl.violin(
#             adata_batch,
#             keys=['n_genes_by_counts', 'total_counts', 'total_counts_mt'],
#             jitter=0.4,
#             multi_panel=True,
#             show=False,
#             xlabel=f'Batch: {batch}',
#         )
        
#         # Save the figure
#         plt.savefig(f"{cell_type}_{batch}_violin.png")
#         plt.close()