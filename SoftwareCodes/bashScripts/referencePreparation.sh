#!/bin/bash

## Define paths to references
GENOMEDIR=/path_to_genome_reference
TRANSCDIR=/path_to_transcriptome_reference
ANNOTDIR=/path_to_annotation_files
DEXSeqPyScripts=/path_to_DEXSeq_python_scripts

## Building STAR indexes
STAR --runMode genomeGenerate --genomeDir $GENOMEDIR/STAR --genomeFastaFiles $GENOMEDIR \
--runThreadN 20

## Building Bowtie indexes

bowtie-build $TRANSCDIR/transcriptome_reference.fasta $TRANSCDIR/transcriptome_reference

## Building KALLISTO index
kallisto index -i $TRANSCDIR/transcripts.idx $TRANSCDIR/transcriptome_reference.fasta
## Formatting GTF annotation file 

python $DEXSeqPyScripts/dexseq_prepare_annotation.py --aggregate='no' \
$ANNOTDIR/annotation_file.gtf $ANNOTDIR/flattened.dexseq.gtf
