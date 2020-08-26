#!/bin/bash

## This script ulitizes sam2tsv (https://github.com/lindenb/jvarkit/; version ec2c2364) to extract aligment details from bam files and then T->C mutations were identified in both with and without TFEA chemical reaction (served as control) libraries.

## Pre-install requirment: sam2tsv, samtools, perl

bam_dir=$1 # bam file folder including bam file which was generated in step1
i=$2 # sample names: including TFEA treated and untreated samples (control)
num_barcode=$3 # number of barcodes shown in the bam file
num_core_barcode=$4 # half of $3
genome_fa=$5 # genome fasta file (mm10 or hg38)

analysis_dir=${bam_dir}/TC_calling
mkdir -p ${analysis_dir}

output_dir=${analysis_dir} 

TMP_dir=${analysis_dir}/tmp/
mkdir -p ${TMP_dir}


############ 0. check the avalibality of bam file

bam_file_prex=${i}_starAligned.sorted.merged.GeneExonTagged.TagIntronic.clean.${num_barcode}

if [ ! -f ${bam_dir}/${bam_file_prex}.bam ]; then
   echo "Please verify that bam file ${bam_dir}/${bam_file_prex}.bam exists!"
   exit
fi

############ 1. bam to sam
samtools view -@ 10 -h ${bam_dir}/${bam_file_prex}.bam > ${output_dir}/${i}.sam
  
############ 2. extract read ID with GE:Z annotation (get read ID with gene annotation)
grep "GE:Z" ${i}.sam | awk '{print $1}' > ${i}_ann_read.txt
  
############ 3. sam to tsv by sam2tsv
cd ${output_dir}
sam2tsv --reference ${genome_fa} ${i}.sam | awk '{if ($9 ~/M|=|X/ ) print $0}'  > ${i}_both_strand_all.tsv

############ 4. Extract T->C mutation from _all.tsv considering both strands
awk '{if ($2 ==0 && $5=="C" && $8=="T") print $0; else if ( $2==16 && $5=="G" && $8=="A" ) print $0}' ${i}_both_strand_all.tsv > ${i}_both_strand_all_TC.tsv
  
############ 5. Retain T-to-C substitutions with a base Phred quality score of > 27
perl extrac_refT_readC.pl -tsv ${i}_both_strand_all_TC.tsv -qual 27

############ 6. QC output: summarise the all types of mutation rate
perl extrac_conversion_frequency_gene_annotate.pl -read ${i}_ann_read.txt -tsv ${i}_both_strand_all.tsv -qual 27

### In the end you will get 2 tsv files of read ID with T->C mutations. 
### One is from the sample with TFEA treatment: #{sample}_both_strand_all_TC.tsv_q27.tsv
### The other is from the Control sample without TFEA treatment: #{control}_both_strand_all_TC.tsv_q27.tsv


