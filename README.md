# 🧫 Single-Cell RNA-seq Analysis

Individual assignment covering the complete scRNA-seq workflow across three sections: raw data pre-processing on Galaxy, downstream analysis with Scanpy, and deep exploration of the AnnData data structure.

---

## 📁 Repository Structure

```
📦 scrnaseq-assignment/
│
├── 📂 section1-galaxy/
│   └── README.md                        ← Galaxy pre-processing walkthrough
│
├── 📂 section2-scrnaseq/
│   ├── README.md                        ← Analysis pipeline documentation
│   └── basic_scrna_tutorial.ipynb       ← NeurIPS 2021 bone marrow dataset
│
├── 📂 section3-anndata/
│   ├── README.md                        ← AnnData structure documentation
│   ├── anndata_intro_tutorial.ipynb     ← Build AnnData from scratch
│   └── anndata_tutorial.ipynb           ← Explore real PBMC3k dataset
│
└── README.md                            ← This file
```

---

## 📋 Sections

| # | Section | Description | Tools |
|---|---------|-------------|-------|
| 1 | [Pre-processing (Galaxy)](./section1-galaxy/README.md) | QC, alignment, and count matrix generation from raw 10X FASTQ reads | Galaxy, STARsolo, FastQC |
| 2 | [Basic scRNA-seq Tutorial](./section2-scrnaseq/README.md) | QC filtering, normalization, UMAP, Leiden clustering, CellTypist annotation | Scanpy, CellTypist, OmniPath |
| 3 | [AnnData](./section3-anndata/README.md) | AnnData object structure, all slots, views, layers, I/O | anndata, scanpy |

---

## 🔗 Tutorial Links

- [Hands-on: Pre-processing of 10X Single-Cell RNA Datasets (Galaxy)](https://training.galaxyproject.org/training-material/topics/single-cell/)
- [Basic scRNA-seq Tutorial — Updated (scverse)](https://scverse.org/training)
- [AnnData Getting Started](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html)
- [scverse AnnData Tutorial](https://scverse-tutorials.readthedocs.io/en/latest/notebooks/anndata_getting_started.html)

---

## ⚙️ Environment Setup

```bash
conda create -n scrna python=3.10
conda activate scrna
pip install anndata scanpy pooch celltypist omnipath leidenalg scipy matplotlib seaborn
```

