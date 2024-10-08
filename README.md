# scRNA-seq Data Analysis Scripts

## Description
This repository contains a set of Python scripts tailored for the analysis of single-cell RNA sequencing (scRNA-seq) data. Each script is designed to address specific stages of the analysis pipeline, including preprocessing, quality control, cell type annotation, data integration, and sample selection. These scripts were developed to streamline scRNA-seq analysis workflows and ensure reproducibility. The repository is suitable for researchers aiming to analyze complex cellular datasets, identify novel cell populations, and achieve robust data integration.

## Key Scripts
- **`sc_pp.py`**: Performs preprocessing and quality control of scRNA-seq data using `scanpy` and `doubletdetection`. This script filters out low-quality cells and doublets, ensuring high-quality data for downstream analyses.
- **`celltypist_model.py`**: Utilizes `celltypist` for automated cell type annotation, allowing for rapid and accurate identification of cell populations across diverse scRNA-seq datasets.
- **`integration.py`**: Integrates multiple scRNA-seq datasets using `scvi` for deep generative modeling and `ray` for hyperparameter optimization. This script ensures that datasets from different batches are harmonized, enabling cross-sample comparisons.
- **`select_samples.py`**: Randomly selects subsets of samples for validation and testing, supporting the generation of balanced datasets for downstream analyses.

## Contents
- `sc_pp.py`: Preprocessing and QC using scanpy and doubletdetection.
- `celltypist_model.py`: Automated cell type annotation with celltypist.
- `integration.py`: Data integration with scvi and hyperparameter tuning using ray.
- `select_samples.py`: Randomly selects samples for analysis.

## Keywords
scRNA-seq, single-cell RNA sequencing, data analysis, bioinformatics, Python, scanpy, doubletdetection, celltypist, scvi, ray
