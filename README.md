# Single-Cell RNA-Seq Standard Pipeline in R or Scanpy 🧬

This repository provides a comprehensive pipeline for processing and analyzing single-cell RNA-sequencing data using either R or Python/Scanpy.  Each workflow is fully scripted (RMarkdown for R, Jupyter Notebook for Scanpy), modular, and checkpointed so you can save intermediate objects and resume without re-running everything.


## Pipeline Overview 🚀

### Section 1: QC & Standard Data Preparation 📊

| Step                                      | R / Seurat                                             | Python / Scanpy                                              |
|-------------------------------------------|--------------------------------------------------------|--------------------------------------------------------------|
| **Load packages & initialize object**     | `library(Seurat)`<br>`CreateSeuratObject()`            | `import scanpy as sc`<br>`adata = sc.read_10x_mtx(...)`      |
| **Assess mitochondrial content**          | `PercentageFeatureSet(..., pattern = "^MT-")`          | `adata.obs['percent_mt'] = sc.pp.calculate_qc_metrics(adata, percent_top=NULL)['pct_counts_mt']` |
| **Visualize metadata distributions**      | `VlnPlot()`<br>`FeatureScatter()`                      | `sc.pl.violin(adata, ...)`<br>`sc.pl.scatter(adata, ...)`    |
| **Filter low-quality cells & features**   | `seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > X & percent.mt < Y)` | `adata = adata[adata.obs.n_genes > X & adata.obs.pct_counts_mt < Y]` |
| **Normalize data**                        | `NormalizeData(method = "LogNormalize")`               | `sc.pp.normalize_total(adata)`<br>`sc.pp.log1p(adata)`       |
| **Cell cycle scoring**                    | `CellCycleScoring()`                                   | `sc.tl.score_genes_cell_cycle(adata)`                        |
| **Identify variable features**            | `FindVariableFeatures()`                               | `sc.pp.highly_variable_genes(adata)`                         |

### Section 2: Doublet Detection 🔬

| Step                                           | R / Seurat + DoubletFinder                            | Python / Scanpy + Scrublet                                   |
|------------------------------------------------|--------------------------------------------------------|--------------------------------------------------------------|
| **Load doublet-detection tools**               | `library(DoubletFinder)`                               | `import scrublet as scrb`                                    |
| **Import processed data**                      | *use existing Seurat object with QC’d & normalized data* | `adata = sc.read_h5ad("filtered_norm.h5ad")`                 |
| **Parameter sweep (pN/pK vs. expected rate)**   | `paramSweep_v3(); find.pK()`                           | `scrub = scrb.Scrublet(adata.X, expected_doublet_rate=0.06)` |
| **Classify & annotate singlets/doublets**      | `doubletFinder_v3()`                                   | `doublet_scores, predicted = scrub.scrub_doublets()`         |
| **Re-cluster & validate singlet populations**  | `RunPCA(); RunUMAP(); FindNeighbors(); FindClusters()`  | `sc.tl.pca(adata); sc.pp.neighbors(adata); sc.tl.umap(adata); sc.tl.leiden(adata)` |
| **Visualize doublet labels**                   | `DimPlot(seurat_obj, group.by = "DF.classifications")` | `sc.pl.umap(adata, color='predicted_doublet', size=10)`      |


All code lives in parallel directories:
	•	R/ with .Rmd files, plots and saved Seurat objects
	•	Scanpy/ with Jupyter notebooks, figures and saved .h5ad files

Feel free to pick either implementation depending on your language preference!
