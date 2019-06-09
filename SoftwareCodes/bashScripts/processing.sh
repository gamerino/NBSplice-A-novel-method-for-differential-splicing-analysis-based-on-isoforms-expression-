#!/bin/bash

## Define paths to files

# Reference files
TRANSCDIR=/path_to_transcriptome_reference

#file containing the samples names
SAMPLES=/path_to_samplesNames_file/samplesNames

# Path to simulated RNA-seq reads
ISOSIM=/path_to_simulated_samples

# Path to Kallisto quantification
KALLISTOEX=/path_to_quantification_Kallisto

for file in `cat $SAMPLES`
do
## Kallisto 
	cd $KALLISTOEX
	mkdir ${file}_kallisto
	cd ${file}_kallisto
	kallisto quant -i $TRANSCDIR/transcripts.idx -b 100 -t 10 $ISOSIM/${file}_1.fq $ISOSIM/${file}_2.fq 
done
