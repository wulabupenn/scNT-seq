#!/bin/bash

## This script ulitizes the Drop-seq computational pipeline (James Nemesh, McCarroll Lab, version 1.0.1, related documents: http://mccarrolllab.org/wp-content/uploads/2016/03/Drop-seqAlignmentCookbookv1.2Jan2016.pdf; Macosko et al., 2015) to map the reads into the genome and tag the reads with cell barcode, UMI barcode and Gene annotation in bam files

## Pre-install requirment: Drop-seq computational pipeline(version 1.0.1), Java, Picardtools, samtools, STAR, perl

sample=$1 # sample name
index=$2 # sequencing index (the order in the samplesheet for demultliplexing: "_S1","_S2")
analysis_dir=$3  # analysis directory 

fastq_dir=$4 # fastq file directory after bcl2fastq

dropseq_root=$5 ## the dropseq computational pipeline directory
star_index_dir=$6 # reference indexed by STAR
fasta=$7 # reference sequence
gtf=$8 # reference gtf

num_barcodes=$9 # NUM_BARCODES= <roughly 2x the number of cells> 
num_core_barcodes=$10 # the number of cells from scNT-seq run

TMP_dir=$11 # directory for temporary
picard_root=$12 # picard tools directory


############ 0. create the folder

output_dir=${analysis_dir}/${sample}_${index}_${genome}
mkdir -p ${output_dir}

output_dir_STAR_modified=${analysis_dir}/${sample}_${index}_${genome}/output_STAR_modified
mkdir -p ${output_dir_STAR_modified}

star_result_dir=${output_dir}/STAR_results_modified0.3
mkdir -p ${star_result_dir}


############ 1. convert Fastq to unaligned BAM (picard)
java -Xmx64g -Djava.io.tmpdir=${TMP_dir} -jar ${picard_root}/picard.jar FastqToSam \
FASTQ=${fastq_dir}/${sample}_${index}_R1_001.fastq.gz \
FASTQ2=${fastq_dir}/${sample}_${index}_R2_001.fastq.gz \
QUALITY_FORMAT=Standard \
OUTPUT=${output_dir}/${sample}_unaligned.bam \
SAMPLE_NAME=${sample} \
SORT_ORDER=queryname


############ 2. Tag cell barcode 
java -Xmx64g -Djava.io.tmpdir=${TMP_dir} -jar ${dropseq_root}/jar/dropseq.jar TagBamWithReadSequenceExtended \
SUMMARY=${output_dir}/${sample}_unaligned_tagged_Cellular.bam_summary.txt \
BASE_RANGE=1-12 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=False \
TAG_NAME=XC \
NUM_BASES_BELOW_QUALITY=1 \
INPUT=${output_dir}/${sample}_unaligned.bam OUTPUT=${output_dir}/${sample}_unaligned_tagged_Cell.bam

############ 3. Tag molecular barcode
java -Xmx64g -Djava.io.tmpdir=${TMP_dir} -jar ${dropseq_root}/jar/dropseq.jar TagBamWithReadSequenceExtended \
SUMMARY=${output_dir}/${sample}_unaligned_tagged_Molecular.bam_summary.txt \
BASE_RANGE=13-20 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=True \
TAG_NAME=XM \
NUM_BASES_BELOW_QUALITY=1 \
INPUT=${output_dir}/${sample}_unaligned_tagged_Cell.bam \
OUTPUT=${output_dir}/${sample}_unaligned_tagged_CellMolecular.bam

############ 4. Filter tagged bam files
java -Xmx64g -Djava.io.tmpdir=${TMP_dir} -jar ${dropseq_root}/jar/dropseq.jar FilterBAM \
TAG_REJECT=XQ INPUT=${output_dir}/${sample}_unaligned_tagged_CellMolecular.bam \
OUTPUT=${output_dir}/${sample}_unaligned_barcode.tagged_filtered.bam

############ 5. Trim 5' SMART adaptor (TSO) sequence from read 2
java -Xmx64g -Djava.io.tmpdir=${TMP_dir} -jar ${dropseq_root}/jar/dropseq.jar TrimStartingSequence \
OUTPUT_SUMMARY=${output_dir}/${sample}_adaptor_trimming_report.txt \
SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
MISMATCHES=0 NUM_BASES=5 \
INPUT=${output_dir}/${sample}_unaligned_barcode.tagged_filtered.bam \
OUTPUT=${output_dir}/${sample}_unaligned_barcode.tagged_filtered_smart.trimmed.bam

############ 6. Trim 3' polyA tails from read 2
java -Xmx64g -Djava.io.tmpdir=${TMP_dir} -jar ${dropseq_root}/jar/dropseq.jar PolyATrimmer \
OUTPUT_SUMMARY=${output_dir}/${sample}_polyA_trimming_report.txt \
MISMATCHES=0 NUM_BASES=6 \
INPUT=${output_dir}/${sample}_unaligned_barcode.tagged_filtered_smart.trimmed.bam \
OUTPUT=${output_dir}/${sample}_unaligned_barcode.tagged_filtered_smart.polyA.trimmed.bam

############ 7. Convert unaligned BAM back to Fastq (picard)
java -Xmx64g -Djava.io.tmpdir=${TMP_dir} -jar ${picard_root}/picard.jar SamToFastq \
INPUT=${output_dir}/${sample}_unaligned_barcode.tagged_filtered_smart.polyA.trimmed.bam \
FASTQ=${output_dir}/${sample}_unaligned_barcode.tagged_filtered_smart.polyA.trimmed.fastq

############ 8. Aligment to genome
gzip ${output_dir}/${sample}_unaligned_barcode.tagged_filtered_smart.polyA.trimmed.fastq
STAR --runThreadN 10 --genomeDir ${star_index_dir} --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --readFilesCommand zcat --readFilesIn ${output_dir}/${sample}_unaligned_barcode.tagged_filtered_smart.polyA.trimmed.fastq --outFileNamePrefix ${star_result_dir}/${sample}_star

############ 9. sort aligned sam file
java -Xmx64g -Djava.io.tmpdir=${TMP_dir} -jar ${picard_root}/picard.jar SortSam \
INPUT=${star_result_dir}/${sample}_starAligned.out.sam \
OUTPUT=${output_dir_STAR_modified}/${sample}_starAligned.sorted.bam \
SORT_ORDER=queryname

############ 10. merge the sorted alignment output from STAR with unaligned/barcode.tagged BAM file
java -Xmx64g -Djava.io.tmpdir=${TMP_dir} -jar ${picard_root}/picard.jar MergeBamAlignment \
REFERENCE_SEQUENCE=$fasta \
UNMAPPED_BAM=${output_dir}/${sample}_unaligned_barcode.tagged_filtered_smart.polyA.trimmed.bam \
ALIGNED_BAM=${output_dir_STAR_modified}/${sample}_starAligned.sorted.bam \
OUTPUT=${output_dir_STAR_modified}/${sample}_starAligned.sorted.merged.bam \
INCLUDE_SECONDARY_ALIGNMENTS=false \
PAIRED_RUN=false

############ 11. Adds a BAM tag "GE" onto reads when the read overlaps the exon of a gene
java -Xmx64g -Djava.io.tmpdir=${TMP_dir} -jar ${dropseq_root}/jar/dropseq.jar TagReadWithGeneExon \
INPUT=${output_dir_STAR_modified}/${sample}_starAligned.sorted.merged.bam \
OUTPUT=${output_dir_STAR_modified}/${sample}_starAligned.sorted.merged.GeneExonTagged.bam \
ANNOTATIONS_FILE=$gtf TAG=GE

############ 12. reterive the intronic reads
perl TagIntronicRead_V3.pl -gtf $gtf -bam ${output_dir_STAR_modified}/${sample}_starAligned.sorted.merged.GeneExonTagged.bam

############ 13. Detect and repair barcode synthesis error
java -Xmx64g -Djava.io.tmpdir=${TMP_dir} -jar ${dropseq_root}/jar/dropseq.jar DetectBeadSynthesisErrors \
  INPUT=${output_dir_STAR_modified}/${sample}_starAligned.sorted.merged.GeneExonTagged.TagIntronic.bam \
  OUTPUT=${output_dir_STAR_modified}/${sample}_starAligned.sorted.merged.GeneExonTagged.TagIntronic.clean.${num_barcodes}.bam \
  OUTPUT_STATS=${output_dir_STAR_modified}/${sample}_exonic.intronic_synthesis_stats_num.barcode.${num_barcodes}.txt \
  SUMMARY=${output_dir_STAR_modified}/${sample}_exonic.intronic_synthesis_stats.summary_num.barcode.${num_barcodes}.txt \
  NUM_BARCODES=${num_barcodes} PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC

### In the end you will get a bam file with GE:Z:(genename) annotation tag. The file name is: #{sample}_starAligned.sorted.merged.GeneExonTagged.TagIntronic.clean.1200.bam



