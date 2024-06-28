# Single-cell RNA-sequencing of differentiating iPS cells reveals dynamic genetic effects on gene expression

This repository contains scripts for data processing, analysis and figure generation using scRNA-Seq, bulk RNA-seq and ChIP-seq data for our paper:

Cuomo*, Seaton*, McCarthy* et al. [Single-cell RNA-sequencing of differentiating iPS cells reveals dynamic genetic effects on gene expression](https://www.nature.com/articles/s41467-020-14457-z), Nature Communications, 2020.

## Analysis scripts

The following folders contain scripts for data processing and analysis.

A short description can be found below:

* [Preprocessing steps](scrnaseq_preprocessing/) contains snakemake files to process the sequencing data (including alignment, donor assignment etc.).

* [QC and merging steps](merging_and_qc/) contains jupyter notebooks to merge experiment-level SCE objects and perform QC and normalization steps to obtain the final SCE object used for all following analyses.

* [Plotting Notebooks](plotting_notebooks/) contains all jupyter notebooks to reproduce the individual figures (main and supplements).


## Data availability 

All HipSci data can be accessed from http://www.hipsci.org.

### Bulk RNA-seq

Bulk RNA-seq data are available under accession numbers: ERP007111 (ENA project) and EGAS00001001137, EGAS00001000593 (EGA projects). 

### Single cell RNA-seq

Single cell RNA-seq data are available under the accession numbers ERP016000 (ENA project) and EGAS00001002278, EGAD00001005741(EGA project: study ID, dataset ID). 

### ChIP-seq

All Chip-seq data used is available at PRJNA593217. 

### Processed Data 

Processed and raw single cell count data, metadata, as well as donor-level allele-specific expression (ASE) data are available at [this Zenodo link](https://zenodo.org/record/3625024). 

Note that raw counts are not integer numbers due to feature quantification being performed using [Salmon](https://www.nature.com/articles/nmeth.4197). 

See details of preprocessing steps in [the Snakefile provided](../main/scrnaseq_preprocessing/Snakefile).



