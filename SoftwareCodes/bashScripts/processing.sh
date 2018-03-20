#!/bin/bash

## Define paths to files

# Reference files
TRANSCDIR=/path_to_transcriptome_reference
GENOMEDIR=/path_to_genome_reference
ANNOTDIR=/path_to_annotation_files

#file containing the samples names
SAMPLES=/path_to_samplesNames_file/samplesNames

# Path to simulated RNA-seq reads
ISOSIM=/path_to_simulated_samples

# Path to STAR alignment
GENOMEALIGN=/path_to_save_genomeAlign

# Path to Kallisto quantification
KALLISTOEX=/path_to_quantification_Kallisto

# DEXSeq expression
DEXSeqPyScripts=/path_to_DEXSeq_python_scripts
DEXSEQEXP=/path_to_quantification_exonLev_DEXSeq

for file in `cat $SAMPLES`
do
## STAR alignment
	cd $GENOMEALIGN
	mkdir ${file}_STAR
	cd ${file}_STAR
	STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --genomeDir $GENOMEDIR/STAR \
	--runThreadN 20 --outFileNamePrefix ${file} --readFilesIn $ISOSIM/${file}_1.fq $ISOSIM/${file}_2.fq 
## Exon level quantification
	cd $DEXSEQEXP
	mkdir ${file}_dxsq
	cd ${file}_dxsq
	python $DEXSeqPyScripts/dexseq_count.py -p yes -f bam -a 20 -r pos \
	$ANNOTDIR/flattened.dexseq.gtf $GENOMEALIGN/${file}_STAR/${file}.bam \
	$file.htseq.counts
## Kallisto 
	cd $KALLISTOEX
	mkdir ${file}_kallisto
	cd ${file}_kallisto
	kallisto quant -i $TRANSCDIR/transcripts.idx -b 100 -t 10 $ISOSIM/${file}_1.fq $ISOSIM/${file}_2.fq 
done
