# 📦 Section 3: AnnData - Annotated Data for Single-Cell Analysis

Two notebooks exploring the **AnnData** data structure: first building one from scratch, then dissecting a real preprocessed PBMC dataset to understand every component in depth.

> 📚 **Tutorial Sources:**
> - [AnnData Getting Started - ReadTheDocs](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html)
> - [scverse AnnData Tutorial](https://scverse-tutorials.readthedocs.io/en/latest/notebooks/anndata_getting_started.html)

---

## 📌 Overview

`AnnData` is the core data container of the scverse ecosystem. It stores the expression matrix alongside rich cell and gene metadata, embeddings, graphs, and arbitrary annotations - all in one object that can be saved to and loaded from disk efficiently.

---

## 📦 Dependencies

```bash
pip install anndata scanpy scipy pooch
```

| Library | Purpose |
|---------|---------|
| `anndata` | Core AnnData data structure |
| `scanpy` | Used in Notebook 2 for preprocessing helpers & plots |
| `scipy` | Sparse matrix support (`csr_matrix`) |
| `numpy` | Array operations |
| `pandas` | DataFrame-based metadata |
| `pooch` | Reproducible dataset downloads |
| `matplotlib` | Visualizations in Notebook 2 |

---

## 📓 Notebooks

| File | Dataset | Focus |
|------|---------|-------|
| `anndata_intro_tutorial.ipynb` | Synthetic (100 cells × 2000 genes, Poisson counts) | Building AnnData from scratch, all slots, I/O |
| `anndata_tutorial.ipynb` | Real PBMC3k preprocessed dataset (~100 MB `.h5ad`) | Exploring a real-world AnnData in depth |

---

## 🏗️ AnnData Object Structure

```
AnnData  n_obs × n_vars  (cells × genes)
│
├── .X            Main data matrix              [sparse/dense, n_obs × n_vars]
│
├── .obs          Cell-level metadata           [pandas DataFrame, n_obs rows]
│                 e.g. cell_type, leiden_cluster, pct_mito, is_low_quality
│
├── .var          Gene-level metadata           [pandas DataFrame, n_vars rows]
│                 e.g. gene_ids, gene_names, highly_variable, mean
│
├── .obsm         Multi-dim cell embeddings     [dict of arrays]
│                 e.g. X_pca, X_umap, X_tsne
│
├── .varm         Multi-dim gene matrices       [dict of arrays]
│                 e.g. gene_stuff (PCA loadings)
│
├── .obsp         Pairwise cell matrices        [dict of sparse matrices]
│                 e.g. distances_all, connectivities
│
├── .layers       Alternative data matrices     [dict of arrays]
│                 e.g. counts (raw), log_transformed, counts_per_million
│
└── .uns          Unstructured metadata         [dict]
                  e.g. louvain_colors, pca variance, neighbor params
```

---

## 📖 Notebook 1 - `anndata_intro_tutorial.ipynb`

**Build an AnnData from scratch** using synthetic Poisson count data.

### Topics Covered

| Step | What is demonstrated |
|------|----------------------|
| Initialize | `ad.AnnData(csr_matrix(...))` — sparse matrix as `.X` |
| Name cells/genes | `adata.obs_names`, `adata.var_names` |
| Subset | By name, boolean mask, and integer index |
| `.obs` metadata | Add `cell_type` as a pandas Categorical |
| `.obsm` / `.varm` | Store UMAP coordinates and gene-level matrices |
| `.uns` | Store free-form metadata (lists, dicts) |
| `.layers` | Add `log_transformed` layer alongside raw `.X` |
| Save to disk | `adata.write('my_results.h5ad', compression='gzip')` |
| Views vs. copies | Subsetting returns a memory-efficient **view**; `.copy()` makes it independent |
| Backed mode | `ad.read_h5ad(..., backed='r')` for partial reading of large files |

### Key Code Snippets

```python
import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

# Create AnnData with sparse count matrix
counts = csr_matrix(np.random.poisson(1, size=(100, 2000)), dtype=np.float32)
adata = ad.AnnData(counts)
adata.obs_names = [f'Cell_{i:d}' for i in range(adata.n_obs)]
adata.var_names = [f'Gene_{i:d}' for i in range(adata.n_vars)]

# Add cell type metadata
adata.obs['cell_type'] = pd.Categorical(
    np.random.choice(['B', 'T', 'Monocyte'], size=adata.n_obs)
)

# Add layers and embeddings
adata.layers['log_transformed'] = np.log1p(adata.X)
adata.obsm['X_umap'] = np.random.normal(0, 1, size=(adata.n_obs, 2))

# Save
adata.write('my_results.h5ad', compression='gzip')

# View vs. copy
view  = adata[:5, ['Gene_1', 'Gene_3']]        # lightweight view
copy  = adata[:5, ['Gene_1', 'Gene_3']].copy() # independent copy
```

---

## 📖 Notebook 2 - `anndata_tutorial.ipynb`

**Explore a real preprocessed PBMC3k dataset** — all AnnData slots with real biological data.

### Dataset
Preprocessed PBMC3k dataset downloaded via `pooch`:
```
https://exampledata.scverse.org/tutorials/scverse-getting-started-anndata-pbmc3k_processed.h5ad
```

### Topics Covered

| Step | What is demonstrated |
|------|----------------------|
| `.X` inspection | Sparse matrix, non-zero fraction, raw data values |
| `.layers` | Accessing `raw` layer; adding a new `counts_per_million` (CPM) layer |
| `.obs` | Cell annotations: `louvain_cell_types`, `percent_mito`; adding `is_low_quality` flag |
| `.var` | Gene annotations: `gene_ids`, `gene_names`; switching index between Ensembl IDs and symbols |
| `obs_names` / `var_names` | Shared index between matrix and metadata |
| Subsetting | By name list, boolean mask (high-quality cells only) |
| `.obsm` | PCA, t-SNE, and UMAP embeddings; scatter plot of all three |
| `.obsp` | `distances_all` cell-cell distance matrix; visualised with/without cell type ordering |
| `.uns` | Louvain params, cluster colours, PCA variance explained |
| Views vs. copies | Demonstrated on real data — parent modification propagates to view but not to copy |

### Key Code Snippets

```python
import anndata, scanpy as sc, numpy as np

# Load preprocessed PBMC
adata = anndata.read_h5ad('pbmc3k_processed.h5ad')

# Add CPM layer
adata.layers['counts_per_million'] = adata.layers['raw'].copy()
sc.pp.normalize_total(adata, target_sum=1e6, layer='counts_per_million')

# Flag low-quality cells
adata.obs['is_low_quality'] = adata.obs['percent_mito'] > 0.03
adata_hq = adata[~adata.obs['is_low_quality'], :]

# Visualise PCA / t-SNE / UMAP
import matplotlib.pyplot as plt
for i, (key, title) in enumerate([('X_pca','PCA'),('X_tsne','t-SNE'),('X_umap','UMAP')]):
    plt.subplot(1, 3, i+1)
    plt.scatter(adata.obsm[key][:,0], adata.obsm[key][:,1],
                c=adata.obs['louvain_cell_types']=='B cells', s=3, cmap='coolwarm')
    plt.title(title); plt.axis('off')
```

---

## 📂 Output Files

```
section3-anndata/
├── README.md
├── anndata_intro_tutorial.ipynb    ← Synthetic data, build from scratch
├── anndata_tutorial.ipynb          ← Real PBMC3k data, deep exploration
└── data/
    └── my_results.h5ad             ← Written by Notebook 1
```

---

## 🔗 References

- [AnnData Official Documentation](https://anndata.readthedocs.io/en/latest/)
- [scverse Tutorials](https://scverse-tutorials.readthedocs.io/)
- [AnnData GitHub Repository](https://github.com/scverse/anndata)
- [Virshup et al. 2021 — AnnData paper](https://doi.org/10.1101/2021.12.16.473007)
