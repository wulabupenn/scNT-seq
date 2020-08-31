---
title: "scNT-seq_README"
author: "Peng Hu"
date: "2020/8/31"
output: html_document
---

Source code of the manuscript **Massively parallel and time-resolved RNA sequencing in single cells with scNT-seq**. Nature Methods (in press) https://www.biorxiv.org/content/10.1101/2019.12.19.882050v2) 

# Content

### scNT-seq TC calling pipeline includes following steps:

Step1_alignment.sh: utilize the Drop-seq computational pipeline (James Nemesh, McCarroll Lab, version 1.12; Macosko et al., 2015) to map the reads to the genome and tag the reads with cell barcode, UMI barcode and gene annotation in bam files. Next, we extracted intronic reads in bam file because the legacy Drop-seq computational pipeline (version 1.12) only consider exonic reads.

Step2_extract_alignment_info.sh: sam2tsv (https://github.com/lindenb/jvarkit/; version ec2c2364) is used to extract detailed alignment information from bam files and then T-to-C substitutions were identified in both experimental and control samples (without Timelapse chemical conversion reaction, as a control for background mutations).

Step3_substract_background_locus.sh: exclude the genomic sites with background T-to-C substitutions from the downstream analysis.

Step4_genetare_TC_matrix.sh: generate labeled and unlabeled gene expression matrix.

### The folder "notebook_for_figures" contains the code to reproduce the main figures

The input files can be downloaded from [here] (https://drive.google.com/drive/folders/1CTdrLUpzye_nlZXWJH9ggS7BRzM-VSqQ?usp=sharing).
Raw data is available at NCBI Gene Expression Omnibus (GEO) (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141851).

For time-resolved RNA velocity analysis in Figure 3a, we use Dynamo (https://github.com/aristoteleo/dynamo-release, developed by Dr. Xiaojie Qiu), an inclusive model of expression dynamics in metabolic labeling based scRNA-seq. 


# Contact
* Peng Hu
* penghu@upenn.edu
* Wu lab, Department of Genetics, University of Pennsylvania (https://www.wulabupenn.org)


