# Ageing Fly Cell Atlas – Differential Ageing Analysis

## Table of contents

1. [Project overview](#project-overview)

2. [Dataset](#dataset)

3. [Directory layout](#directory-layout)

4. [Analysis pipeline](#analysis-pipeline)

   * 4.1 [Environment set-up](#environment-set-up)
   * 4.2 [Data pre-processing](#data-pre-processing)
   * 4.3 [Differential gene-expression (DGE)](#differential-gene-expression-dge)
   * 4.4 [Visualisation](#visualisation)

5. [Running the pipeline](#running-the-pipeline)

6. [Expected outputs](#expected-outputs)

7. [References](#references)

---

## Project overview

This repository contains all code and documentation for the **Ageing Fly Cell Atlas (AFCA) differential-ageing analysis** described in *AFCA analysis report v3*.
We quantify age-related transcriptional changes across 30 adult *Drosophila melanogaster* cell types in the publicly available AFCA single-nucleus RNA-seq dataset (Lu *et al.*, 2023).  Analyses are performed for the whole body and for an aggregated **Enterocyte** lineage with an emphasis on a curated chromatin-associated gene panel:

```
Su(var)205 Su(var)3-9 G9a HP1b HP1c HP4 HP5 HP6 ADD1 Su(var)2-HP2 Su(var)3-7 Lam LamC LBR Kdm4A Kdm4B His2Av His3.3A His3.3B
```

---

## Dataset

| File                     | Description                                                             |
| ------------------------ | ----------------------------------------------------------------------- |
| `adata_body_S_v1.0.h5ad` | AFCA body atlas (5 d/30 d/50 d/70 d; ⟶ 276 k nuclei)  |
| `FCA.h5ad`               | 5-day Fly Cell Atlas reference (harmonised with AFCA) |

All raw files are downloaded from the **Hongjie Li lab** portal ([https://hongjielilab.org/afca](https://hongjielilab.org/afca)).  
Concerning the metadata, only the following columns are required by the downstream scripts: `age`, `sex`, `afca_annotation`, `dataset`, `total_counts`, `pct_counts_mt`.

---

## Directory layout

```
├── Data
│   ├── raw/                         ← .h5ad files (AFCA, FCA)
│   └── ct_specific/                 ← pre-filtered per-cell-type .h5ad files
├── Scripts
│   ├── preprocessing/
│   │   ├── afca_analysis_dgea_set_data_prep_280425.py
│   │   └── afca_analysis_entero_geneset_data_prep_160525.py
│   ├── dgea/
│   │   ├── run_dgea_single_hpc.r
│   │   ├── run_dgea_single_hpc.sh
│   │   └── afca_analyis_MAST_entero_dgea_160525.R
│   └── visualisation/
│       ├── afca_analysis_viz.R
│       └── afca_analysis_entero_vis_170525.R
├── Results
│   └── age_scdgea/YYYYMMDD/         ← csv outputs (`markers_<cell-type>.csv`)
├── Figures
│   └── afca_analysis/YYYYMMDD/      ← volcano plots (.png)
└── README.md                        ← *you are here*
```

---

## Analysis pipeline

### Environment set-up

* **Python 3.10** – `scanpy`, `anndata`, `pandas`, `numpy`, `matplotlib`, `scanpy-recipes`
* **R 4.3**       – `Seurat (≥5)`, `zellkonverter`, `MAST`, `EnhancedVolcano`, `patchwork`, `cowplot`

Create the conda & renv environments via:

### Data pre-processing

`afca_analysis_dgea_set_data_prep_280425.py` loops over **all** AFCA cell-types, applies QC, collapses each into a cell-type-specific `.h5ad` and stores it under `Data/ct_specific/`  .  Key filters:

* ≥200 nuclei overall **and** ≥100 nuclei per age group
* Genes expressed in ≥3 cells **or** in *genes of interest* list
* Add unique sample identifier (`indiv`) extracted from cell barcodes

`afca_analysis_entero_geneset_data_prep_160525.py` repeats the above but merges all *enterocyte* sub-annotations into a single **Enterocyte** object to boost statistical power  .

### Differential gene-expression (DGE)

`run_dgea_single_hpc.r` converts each `.h5ad` to a Seurat object and runs **MAST** via `FindMarkers()` comparing **5-day** (young) nuclei versus 30 d / 50 d / 70 d, with `indiv` as a latent variable.
`run_dgea_single_hpc.sh` dispatches the job across a SLURM array (max 5 simultaneous tasks).

For the **aggregated Enterocyte** lineage, a separate script `afca_analyis_MAST_entero_dgea_160525.R` is used to run MAST locally (not on HPC), starting from the `.h5ad` prepared in step 2  .

### Visualisation

Volcano plots are generated with `afca_analysis_viz.R` (all cell-types) and `afca_analysis_entero_vis_170525.R` (enterocytes only).  Genes of interest are highlighted in **red**; housekeeping controls in **green**.

---

## Running the pipeline

```bash
# 1. Pre-filter all cell-types (~40 min, 16 GB RAM)
python Scripts/preprocessing/afca_analysis_dgea_set_data_prep_280425.py

# 2. Submit DGE jobs (HPC)
cd Scripts/dgea
sbatch run_dgea_single_hpc.sh   # ~1 h per array chunk

# 3. Run enterocyte DGE locally
Rscript afca_analyis_MAST_entero_dgea_160525.R

# 4. Combine individual marker files (optional helper script)
python Scripts/dgea/combine_markers.py

# 5. Generate plots
Rscript Scripts/visualisation/afca_analysis_viz.R
Rscript Scripts/visualisation/afca_analysis_entero_vis_170525.R
```

Local execution without SLURM is also possible:

```bash
Rscript run_dgea_single_hpc.r Data/ct_specific/follicle\ cell.h5ad Results/age_scdgea/20250520
```

---

## Expected outputs

| Path                                                     | Content                                                        |
| -------------------------------------------------------- | -------------------------------------------------------------- |
| `Results/age_scdgea/YYYYMMDD/markers_<cell-type>.csv`    | MAST statistics (avg\_log2FC, pct.1/2, p\_val, p\_val\_adj, …) |
| `Figures/afca_analysis/YYYYMMDD/volcano_<cell-type>.png` | Volcano plot per comparison                                    |
| `Figures/afca_analysis/YYYYMMDD/volcano_enterocyte.png`  | Aggregated enterocyte plot                                     |

The main biological findings are summarised in *AFCA analysis report v3* (Section 3).  Briefly, adult **oenocytes** show the strongest age-associated down-regulation of nuclear lamina and HP1 family genes, whereas **indirect flight muscle** displays progressive up-regulation of Kdm4 histone demethylases with age.

---

## References

* T.-C. Lu *et al.* **The Aging Fly Cell Atlas identifies exhaustive ageing features at cellular resolution**. *Science* 380 e-adg0934 (2023).
* H. Li *et al.* **Fly Cell Atlas: A single-nucleus transcriptomic atlas of the adult fruit fly**. *Science* 375 e-abk2432 (2022).
* **AFCA analysis report v3** (2025-05-19).  Internal project report, Ageing Fly project.
