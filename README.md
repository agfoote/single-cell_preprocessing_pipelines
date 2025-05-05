# Single-Cell RNA-Seq Standard Pipeline in R 🧬

This repository provides a comprehensive pipeline for processing and analyzing single-cell RNA-sequencing data using R. 

The workflow leverages several R packages—including Seurat, tidyverse, and DoubletFinder (among others)—to facilitate quality control, normalization, doublet detection, downstream analysis, and initial visualization. 
Each folder contains a Markdown (.Rmd) file, including corresponding code, visualizations (plots), and detailed outputs, providing a clear record of the analysis workflow and results.
The pipeline is designed to be modular, allowing users to save intermediate results (e.g., normalized data objects) and resume analysis without re-running the entire workflow.

## _Pipeline Overview 🚀_ 
## Section 1: QC Standard Data Preparation 📊
* Load required packages and create the initial Seurat object
* Assess mitochondrial gene expression
* Visualize distribution of metadata values
* Generate scatter plots for quality control
* Filter out low-quality cells and features
* Normalize data
* Assign cell cycle scores
* Identify variable genes

## Section 2: Doublet Detection 🔬
* Load required packages
* Import processed data for DoubletFinder analysis
* Select optimal parameters (pN and pK)
* Recluster data to identify and validate singlets
* Confirm compartments and feature counts

## Section 3: Initial Refinement 🧹
* Load necessary packages
* Import previously filtered and normalized data
* Perform additional quality checks
* Refine clusters based on updated parameters
* Validate refined datasets
  
## Section 4: Integration & Additional Refinement 🛠️
* Load packages and import refined data
* Integrate datasets using RPCA anchor-based method
* Add and validate metadata
* Visualize integrated data and metadata
* Conduct clustering analysis on integrated data
* Identify marker genes within clusters
* Analyze and refine epithelial populations specifically
* Generate final marker gene lists for clusters
* Produce visualizations and save final integrated datasets
