# Single-Cell RNA-Seq Standard Pipeline in R or Scanpy 🧬

This repository provides a comprehensive pipeline for processing and analyzing single‑cell RNA‑sequencing data using either **R** or **Python/Scanpy**.  
Each workflow is fully scripted (R Markdown for R, Jupyter Notebook for Scanpy), modular, and checkpointed so you can pick up from any stage without re‑running everything.  

---

## Pipeline Overview 🚀

### Section 1 – QC & Standard Data Preparation 📊

| Step                                      | R / Seurat                                             | Python / Scanpy                                              |
|-------------------------------------------|--------------------------------------------------------|--------------------------------------------------------------|
| **Load packages & initialize object**     | `library(Seurat)`<br>`CreateSeuratObject()`            | `import scanpy as sc`<br>`adata = sc.read_10x_mtx(...)`      |
| **Assess mitochondrial content**          | `PercentageFeatureSet(..., pattern = "^MT-")`          | `adata.obs['percent_mt'] = sc.pp.calculate_qc_metrics(adata, percent_top=None)['pct_counts_mt']` |
| **Visualize metadata distributions**      | `VlnPlot()` · `FeatureScatter()`                       | `sc.pl.violin(adata, ...)` · `sc.pl.scatter(adata, ...)`     |
| **Filter low‑quality cells & features**   | `subset(seurat_obj, subset = nFeature_RNA > X & percent.mt < Y)` | `adata = adata[adata.obs.n_genes > X & adata.obs.pct_counts_mt < Y]` |

### Section 2 – Doublet Detection 🔬

| Step                                           | R / Seurat + DoubletFinder                            | Python / Scanpy + Scrublet                                   |
|------------------------------------------------|-------------------------------------------------------|--------------------------------------------------------------|
| **Load doublet‑detection tools**               | `library(DoubletFinder)`                              | `import scrublet as scrb`                                    |
| **Parameter sweep & optimal pK**               | `paramSweep_v3(); find.pK()`                          | `scrub = scrb.Scrublet(adata.X, expected_doublet_rate=0.06)` |
| **Classify singlets/doublets**                 | `doubletFinder_v3()`                                  | `doublet_scores, predicted = scrub.scrub_doublets()`         |
| **Re‑cluster & validate singlets**             | `RunPCA(); RunUMAP(); FindClusters()`                 | `sc.tl.pca(); sc.pp.neighbors(); sc.tl.umap(); sc.tl.leiden()` |
| **Visualize doublet labels**                   | `DimPlot(..., group.by="DF.classifications")`         | `sc.pl.umap(adata, color='predicted_doublet')`               |

---

### Section 3 – Normalize & Log‑Transform the Data 🔄

| Step                           | R / Seurat                                   | Python / Scanpy                                   |
|--------------------------------|----------------------------------------------|---------------------------------------------------|
| **Total‑count normalisation**  | `NormalizeData(seurat_obj, method="LogNormalize")` | `sc.pp.normalize_total(adata)`                    |
| **Log‑transform counts**       | Performed automatically by `NormalizeData()` | `sc.pp.log1p(adata)`                              |
| **Store result**               | `seurat_obj[["RNA"]]@data`                   | Transformed matrix is in `adata.X`                |

---

### Section 4 – Mouse‑Compatible Cell‑Cycle Gene Lists 🧬

| Step                                   | R / Seurat                                          | Python / Scanpy                                             |
|----------------------------------------|-----------------------------------------------------|--------------------------------------------------------------|
| **Prepare mouse CC genes**             | `cc.genes <- ConvertHumanGeneList(cc.genes.updated.2019)` | `mouse_s, mouse_g2m = sc.queries.mgi_gene_sets()` |
| **Score S & G2M phases**               | `CellCycleScoring(seurat_obj, s.features=cc.genes$S, g2m.features=cc.genes$G2M)` | `sc.tl.score_genes_cell_cycle(adata, mouse_s, mouse_g2m)` |
| **Regress out CC effect (optional)**   | `ScaleData(vars.to.regress = c("S.Score","G2M.Score"))` | `sc.pp.regress_out(adata, ['S_score','G2M_score'])`         |

---

### Section 5 – Feature Selection & Dimensionality Reduction 📉

| Step                              | R / Seurat                                               | Python / Scanpy                                                     |
|-----------------------------------|----------------------------------------------------------|---------------------------------------------------------------------|
| **Identify highly‑variable genes**| `FindVariableFeatures(seurat_obj, selection.method="vst")` | `sc.pp.highly_variable_genes(adata, flavor='seurat_v3')`            |
| **Scale data (variable genes)**   | `ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))` | `sc.pp.scale(adata, max_value=10)`                                  |
| **Principal component analysis**  | `RunPCA(seurat_obj, npcs=50)`                            | `sc.tl.pca(adata, svd_solver='arpack')`                             |
| **Neighborhood graph & UMAP/tSNE**| `FindNeighbors(); RunUMAP()`                             | `sc.pp.neighbors(adata); sc.tl.umap(adata)`                         |

---

All code lives in parallel directories:

* `R/` – RMarkdown notebooks, plots and saved **Seurat** objects  
* `Scanpy/` – Jupyter notebooks, figures and saved **.h5ad** files  

Feel free to use either implementation depending on your language preference!

---

*Last updated 2025‑05‑09*
