# scRNA
This repository contains a collection of Python scripts for the anlaysis of single-cell RNA sequencing (scRNA-seq) data. Each script is designed to address specific stages of the analysis pipeline, including preprocessing, quality control, cell type annotation, data integration, and sample selection. These scripts were developed to streamline scRNA-seq analysis workflows and ensure reproducibility. The repository is suitable for researchers aiming to analyze complex cellular datasets, identify novel cell populations, and achieve robust data integration.

Contents:
sc_pp.py: Preprocessing and QC using scanpy and doubletdetection.
celltypist_model.py: Automated cell type annotation with celltypist.
integration.py: Data integration with scvi and hyperparameter tuning using ray.
select_samples.py: Randomly selects samples for analysis.


Keywords: scRNA-seq, single-cell RNA sequencing, data analysis, bioinformatics, Python, scanpy, doubletdetection, celltypist, scvi, ray
