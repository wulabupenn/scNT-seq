---
title: "scNTseq_README"
author: "Peng Hu"
date: "2020/8/26"
output: html_document
---

Source code of the manuscript **Massively parallel, time-resolved single-cell RNA sequencing with scNT-Seq**. https://www.biorxiv.org/content/10.1101/2019.12.19.882050v2) 

# Content

## scNT-seq TC calling pipeline will including following steps:

1. First, ulitize the Drop-seq computational pipeline (James Nemesh, McCarroll Lab, version 1.0.1; Macosko et al., 2015) to map the reads into the genome and tag the reads with cell barcode, UMI barcode and Gene annotation in bam files. Second, annotate intron mapped reads in bam file since default Drop-seq computational pipeline (version 1.0.1) don't consider intronic reads.

2. sam2tsv (https://github.com/lindenb/jvarkit/; version ec2c2364) was used to extract aligment details from bam files and then T->C mutations were identified in both with and without TFEA chemical reaction (served as control) libraries.

3. Exclude the locus which has T->C mutations in control sample

4. Generate labeled and unlabeled matrix

## The folder "notebook_for_figures" contains the code to reproduce the main figures

## The input files could be downloaded from [here] (https://drive.google.com/drive/folders/1CTdrLUpzye_nlZXWJH9ggS7BRzM-VSqQ?usp=sharing)

# Contact
* Peng Hu
* penghu@upenn.edu
* Wu lab, Department of Genetics, University of Pennsylvania


