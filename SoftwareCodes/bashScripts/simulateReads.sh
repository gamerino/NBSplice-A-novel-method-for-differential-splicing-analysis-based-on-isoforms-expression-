#!/bin/bash
## Define paths to files

# Reference files
TRANSCDIR=/path_to_transcriptome_reference

# Simulated reads of an specific scenario
SIMREADS=/path_to_simulated_reads

# Real expression profiles
ISOQUANT=/path_to_save_isoform_quant

# Simulated expression profiles

ISOSIM=/path_to_simulated_expression

sampleName=$1

cd $SIMREADS

rsem-simulate-reads $TRANSCDIR/transcriptome_reference \
$ISOQUANT/$sampleName/${sampleName}.stat/${sampleName}.model \
$ISOSIM/${sampleName}.sim.iso.results 0.1 40000000 $sampleName



