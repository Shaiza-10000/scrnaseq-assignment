# 🧬 Section 2: Basic scRNA-seq Analysis Tutorial (Updated)

End-to-end single-cell RNA-seq analysis of **bone marrow mononuclear cells** from the NeurIPS 2021 benchmarking dataset - from raw count matrices to annotated cell types.

> 📚 **Tutorial Source:** [Basic scRNA-seq Tutorial (Updated) - scverse](https://scverse.org/training)

---

## 📌 Overview

This notebook takes the count matrices produced in Section 1 and performs a complete downstream analysis: quality control, normalization, dimensionality reduction, clustering, marker gene identification, and automated cell type annotation using CellTypist.

**Dataset:** NeurIPS 2021 Benchmarking Dataset = samples `s1d1` and `s1d3`  
**Cell type:** Bone marrow mononuclear cells (immune cell types)  
**Output file:** `pbmc_processed.h5ad`

---

## 📦 Dependencies

```bash
pip install anndata scanpy pooch celltypist omnipath leidenalg
```

| Library | Purpose |
|---------|---------|
| `scanpy` | Core single-cell analysis framework |
| `anndata` | AnnData data structure |
| `pooch` | Reproducible data downloads with hash verification |
| `celltypist` | Automated immune cell type annotation |
| `omnipath` | PanglaoDB marker gene database access |
| `leidenalg` | Leiden community-detection clustering |

---

## 📓 Notebook

| File | Description |
|------|-------------|
| `basic_scrna_tutorial.ipynb` | Full pipeline: QC → normalization → UMAP → clustering → cell type annotation |

---

## 🔄 Analysis Pipeline

```
10X .h5 files (s1d1, s1d3)
        ↓
1.  Load & Concatenate samples
    → ad.concat(adatas, label='sample')
        ↓
2.  Quality Control
    - Flag mitochondrial (MT-), ribosomal (RPS/RPL), haemoglobin (HB) genes
    - sc.pp.calculate_qc_metrics(qc_vars=['mt','ribo','hb'])
    - Filter: min_genes=100, pct_counts_mt < 20
        ↓
3.  Doublet Detection
    - sc.external.pp.scrublet(batch_key='sample')
        ↓
4.  Normalization
    - Save raw counts → adata.layers['counts']
    - sc.pp.normalize_total()  ← normalize to median depth
    - sc.pp.log1p()            ← log-transform
        ↓
5.  Feature Selection
    - Top 2000 Highly Variable Genes (HVGs)
    - sc.pp.highly_variable_genes(n_top_genes=2000, batch_key='sample')
        ↓
6.  PCA → Neighbors Graph → UMAP
    - sc.tl.pca()
    - sc.pp.neighbors()
    - sc.tl.umap()
        ↓
7.  Leiden Clustering (3 resolutions)
    - Resolution 0.5, 1.0, 2.0
        ↓
8.  Remove Doublets
    - Filter predicted_doublet == True
        ↓
9.  Marker Gene Analysis
    - Wilcoxon rank-sum test per Leiden cluster
    - Manual marker lookup (CD14, CD4, CD8A, NKG7, MS4A1 …)
        ↓
10. Automated Annotation — CellTypist
    - Model: Immune_All_Low.pkl
    - majority_voting over leiden_res0_5
        ↓
11. PanglaoDB Marker Lookup (OmniPath)
    - Cross-reference gene markers with known immune cell types
        ↓
Output: pbmc_processed.h5ad
```

---

## 🦠 Cell Types Identified

| Cluster Marker Genes | Cell Type |
|----------------------|-----------|
| `FCN1`, `CD14` | CD14+ Monocyte |
| `TCF7L2`, `FCGR3A`, `LYN` | CD16+ Monocyte |
| `SOX4`, `IL3RA` | Plasmacytoid Dendritic Cell (pDC) |
| `CLEC9A`, `CADM1` | cDC1 |
| `CST3`, `FCER1A` | cDC2 |
| `SPINK2`, `GATA2`, `CD34` | Haematopoietic Stem Cell (HSC) |
| `CD4`, `TCF7`, `CCR7` | CD4+ T cell |
| `CD8A`, `GZMK`, `GZMB` | CD8+ T cell |
| `GNLY`, `NKG7` | NK cell |
| `MS4A1`, `IGHM` | B cell |
| `MZB1`, `IGHG1` | Plasma cell |

---

## 📈 Key Outputs & Plots

| Output | Description |
|--------|-------------|
| Violin plots | Distribution of `n_genes_by_counts`, `total_counts`, `pct_counts_mt` |
| Scatter plot | `total_counts` vs `n_genes_by_counts` coloured by MT% |
| HVG plot | Top 2000 highly variable genes highlighted |
| PCA variance ratio | Scree plot for choosing number of PCs |
| UMAP (sample) | Cells coloured by sample of origin |
| UMAP (doublets) | Predicted doublets and doublet scores |
| UMAP (clusters) | Leiden clusters at 3 resolutions |
| Dot plot | Top 5 marker genes per cluster |
| CellTypist UMAP | `majority_voting` and `predicted_labels` overlaid |

---

## 🔗 References

- [Scanpy Documentation](https://scanpy.readthedocs.io/)
- [CellTypist](https://www.celltypist.org/)
- [OmniPath / PanglaoDB](https://omnipathdb.org/)
- [NeurIPS 2021 scRNA-seq Benchmarking Dataset](https://openproblems.bio/neurips_2021/)
- [Luecken & Theis 2019 — Best Practices for scRNA-seq](https://www.embopress.org/doi/full/10.15252/msb.20188746)
