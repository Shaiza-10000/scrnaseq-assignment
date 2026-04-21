# 🔬 Section 1: Pre-processing of 10X Single-Cell RNA Datasets (Galaxy)

This section covers the full pre-processing pipeline for 10X Genomics Single-Cell RNA sequencing data using the **Galaxy** bioinformatics platform, following the hands-on tutorial from the Galaxy Training Network.

> 📚 **Tutorial Source:** [Hands-on: Pre-processing of 10X Single-Cell RNA Datasets — Single Cell / Galaxy Training Network](https://training.galaxyproject.org/training-material/topics/single-cell/)

---

## 📌 Overview

Single-cell RNA sequencing (scRNA-seq) data from 10X Genomics must go through a rigorous pre-processing pipeline before any downstream analysis. Starting from raw FASTQ reads, this workflow produces a clean, filtered gene-expression count matrix that can be imported directly into tools like Scanpy or Seurat.

---

## 🧰 Tools Used (Galaxy Platform)

| Tool | Purpose |
|------|---------|
| **Galaxy Platform** | Browser-based bioinformatics workflow environment - no local installation required |
| **FastQC** | Per-read quality assessment of raw FASTQ files |
| **MultiQC** | Aggregated QC report across all samples |
| **STARsolo / Cell Ranger** | Alignment of reads to reference genome + cell barcode & UMI demultiplexing |
| **Knee Plot Tool** | Distinguish true cells from empty droplets |

---

## 🔄 Pre-processing Workflow

```
Raw FASTQ Files (R1: barcodes+UMI, R2: cDNA)
            ↓
1. Quality Control
   - FastQC on raw reads
   - MultiQC summary report
            ↓
2. Alignment & Quantification (STARsolo)
   - Map reads to reference genome
   - Demultiplex by cell barcode (16 bp)
   - Deduplicate by UMI (10 bp)
            ↓
3. Cell Filtering
   - Knee / elbow plot to identify real cells
   - Separate filtered vs. raw matrix
            ↓
Output: filtered_feature_bc_matrix/
  ├── barcodes.tsv.gz    ← valid cell barcodes
  ├── features.tsv.gz    ← gene names + IDs
  └── matrix.mtx.gz      ← sparse UMI count matrix
```

---

## 📂 Output File Structure

```
output/
├── filtered_feature_bc_matrix/
│   ├── barcodes.tsv.gz      # Cell barcodes that passed filtering
│   ├── features.tsv.gz      # Ensembl gene IDs and gene symbols
│   └── matrix.mtx.gz        # Sparse count matrix (Market Exchange Format)
├── raw_feature_bc_matrix/
│   └── ...                  # All barcodes before cell-calling
└── web_summary.html         # Run summary with alignment statistics
```

> 📁 This output is used as input in **Section 2** (loaded via `sc.read_10x_h5` or `sc.read_10x_mtx`).

---

## 🔑 Key Concepts

| Concept | Explanation |
|---------|-------------|
| **Cell Barcode** | 16 bp sequence identifying which droplet/cell each read came from |
| **UMI** | Unique Molecular Identifier — 10 bp tag that removes PCR duplicate bias |
| **Knee Plot** | Plot of barcode rank vs. UMI count; the "knee" separates cells from empty droplets |
| **Ambient RNA** | Background mRNA from lysed cells that contaminates all droplets |
| **Filtered Matrix** | Only barcodes called as real cells (above the knee threshold) |

---

## 📊 Quality Metrics to Evaluate

| Metric | Typical Acceptable Range |
|--------|--------------------------|
| Reads Mapped to Genome | > 70% |
| Reads Mapped Confidently to Transcriptome | > 30% |
| Valid Barcodes | > 75% |
| Estimated Number of Cells | Close to expected cell count loaded |
| Median Genes Detected per Cell | 500 – 5,000 (tissue-dependent) |
| Sequencing Saturation | > 50% preferred |

---

## 🔗 References

- [Galaxy Training Network — Single Cell Tutorials](https://training.galaxyproject.org/training-material/topics/single-cell/)
- [10X Genomics: What is Cell Ranger?](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)
- [STARsolo Documentation](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md)
