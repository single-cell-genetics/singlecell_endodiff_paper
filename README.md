# Single-cell RNA-sequencing of differentiating iPS cells reveals dynamic genetic effects on gene expression

This repository contains scripts for data processing, analysis and figure generation using scRNA-Seq, bulk RNA-seq and ChIP-seq data for our paper:

[Single-cell RNA-sequencing of differentiating iPS cells reveals dynamic genetic effects on gene expression](https://www.nature.com/articles/s41467-020-14457-z)

## Analysis scripts

The following folders contain scripts for data processing and analysis.
A short description can be found below:

* [Preprocessing steps](../master/scrnaseq_preprocessing/) contains snakemake files to process the sequencing data (including alignment, donor assignment etc.).

* [QC and merging steps](../master/merging_and_qc/) contains jupyter notebooks to merge experiment-level SCE objects and perform QC and normalization steps to obtain the final SCE object used for all following analyses.

* [Plotting Notebooks](../master/plotting_notebooks/) contains all jupyter notebooks to reproduce the individual figures (main and supplements).


## Data availability 
All HipSci data can be accessed from http://www.hipsci.org.

### Bulk RNA-seq
Bulk RNA-seq data are available under accession numbers: ERP007111 (ENA project) and EGAS00001001137, EGAS00001000593 (EGA projects). 
### Single cell RNA-seq
Single cell RNA-seq data are available under the accession numbers ERP016000 (ENA project) and EGAS00001002278, EGAD00001005741(EGA project: study ID, dataset ID). 
### ChIP-seq
All Chip-seq data used is available at PRJNA593217. 

### Processed Data 
Processed single cell count data and metadata are available from Zenodo: https://zenodo.org/record/3625024#.Xil-0y2cZ0s



