---
title: "scNT-seq_README"
author: "Peng Hu"
date: "2020/8/31"
output: html_document
---

Source code of the manuscript **Massively parallel and time-resolved RNA sequencing in single cells with scNT-seq**. Nature Methods (in press), https://www.biorxiv.org/content/10.1101/2019.12.19.882050v2).

# Content

### scNT-seq TC calling pipeline includes following steps:

- Step1_alignment.sh: utilize the Drop-seq computational pipeline (James Nemesh, McCarroll Lab, version 1.12; Macosko et al., 2015) to map the reads to the genome and tag the reads with cell barcode, UMI barcode and gene annotation in bam files. Next, we extracted intronic reads in bam file because the legacy Drop-seq computational pipeline (version 1.12) only consider exonic reads.

- Step2_extract_alignment_info.sh: sam2tsv (https://github.com/lindenb/jvarkit/; version ec2c2364) is used to extract detailed alignment information from bam files and then T-to-C substitutions are identified in both experimental and control samples (without Timelapse chemical conversion reaction, as a control for background mutations).

- Step3_substract_background_locus.sh: exclude the genomic sites with background T-to-C substitutions from the downstream analysis.

- Step4_genetare_TC_matrix.sh: generate labeled and unlabeled gene expression matrix.

### Raw and processed data files for this study:

- Raw data files are available at NCBI Gene Expression Omnibus (GEO) (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141851).

- The folder "notebook_for_figures" contains the R code to reproduce the main figures. The input files can be downloaded from [here](https://drive.google.com/drive/folders/1CTdrLUpzye_nlZXWJH9ggS7BRzM-VSqQ?usp=sharing). Additional data analysis related information will be available upon request.

- The `neuron_revision_figures_n_s_velocity.ipynb` and `neuron_revision_figures.ipynb` files from the "notebook_for_figures" folder provide the time-resolved RNA velocity analysis and the conventional scRNA-seq RNA velocity analysis with [Dynamo](https://github.com/aristoteleo/dynamo-release). 

- To reproduce the exact figures for time-resolved RNA velocity analysis in Figure 3a, please ensure installing the [Dynamo](https://github.com/aristoteleo/dynamo-release) version as printed out in the corresponding notebooks. Make sure also that `anndata==0.7.1` and `umap-learn==0.3.9`. Tutorials on using the newest dynamo for the scNT-seq dataset and more can be found [here](https://dynamo-release.readthedocs.io/en/latest/scNT_seq.html)


# Contact
* Peng Hu
* penghu@upenn.edu
* [Wu lab](https://www.wulabupenn.org), Department of Genetics, University of Pennsylvania


