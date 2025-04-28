import os
import scanpy as sc
import anndata as ad

# Define directories
main_dir = "/hpc/shared/onco_janssen/dhaynessimmons/data/ageing_flies/afca_data"
sample_dir = os.path.join(main_dir, "afcaCellRangerSummary/allMatrixRaw_filter")
temp_dir = os.path.join(main_dir, "afcaCellRanger_processed/allMatrixRawFilter_proccessed_heads")
output_file = os.path.join(main_dir, "afcaCellRanger_processed/allMatrixRaw_filter_combined.h5ad")
os.makedirs(temp_dir, exist_ok=True)
os.makedirs(output_file, exist_ok=True)
os.makedirs(sample_dir, exist_ok=True)


# List subfolders (one per individual)
sample_list = [f for f in os.listdir(sample_dir) if os.path.isdir(os.path.join(sample_dir, f))]

for sample in sample_list:
    print(f"\nProcessing sample: {sample}")

    if "body" not in sample:
        print(f"\t\tSkipping non-body sample: {sample}")
        continue

    sample_path = os.path.join(main_dir, sample)
    output_path = os.path.join(temp_dir, f"{sample}.h5ad")

    # Skip if already processed
    if os.path.exists(output_path):
        print(f"\t\t\tAlready exists: {output_path}")
        continue

    try:
        adata = sc.read_10x_mtx(sample_path, var_names='gene_symbols', cache=False)
        adata.obs['sample_id'] = sample
        adata.write(output_path)
        print(f"\tSaved: {output_path}")
        del adata
    except Exception as e:
        print(f"\tError processing {sample}: {e}")


# Get list of saved individual .h5ad files
processed_files = [os.path.join(temp_dir, f) for f in os.listdir(temp_dir) if f.endswith(".h5ad")]
processed_files.sort()  # optional: keep order stable

# Check if files are present
if not processed_files:
    raise ValueError("No processed .h5ad files found in: " + temp_dir)

# Load and concatenate with backing
print(f"\nConcatenating {len(processed_files)} samples...")
adatas = [sc.read_h5ad(f, backed='r') for f in processed_files]

adata_all = ad.concat(
    adatas,
    label='sample_id',
    keys=[os.path.basename(f).replace('.h5ad', '') for f in processed_files],
    index_unique='-'
)

# Save the final combined file
adata_all.write(output_file)
print(f"\n Combined AnnData saved to: {output_file}")