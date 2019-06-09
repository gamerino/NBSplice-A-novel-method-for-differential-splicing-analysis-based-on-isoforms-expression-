This file provides the analysis steps performed to simulate and analyze RNA-seq experiments using the provided scripts. Note that all the scripts have definitions about files location that you need to change BEFORE their use. Also, if the program files are installed in other directories than your $PATH, you will need to specify the programs' directory on the script files. 

1-	Samples downloading: Samples mus be obtained from the NCBI ftp and saved those together, for example in a folder called *sraFiles*. 

2-	Extracting fastq files from SRA files running the `fastqExtraction.sh` script. ATTENTION: Before you begin, you need to create a text file (*samplesNames*) listing the names of the sra files, one per line. 

For instance, if you have three sra files called "file1.sra", "file2.sra" and "file3.sra", your *samplesNames* file should look as 

        file1
        file2
        file3

3-	Formatting and indexing reference files using the `referencePreparation.sh` script. 

4-	Obtaining information from real samples: Each real RNA-seq sample is aligned against the reference transcriptome (HG19 here) using `Bowtie`, and alignments are quantified to obtain the expression profile at the isoform level using the `initialProcessing.sh` script.

5-	Modification of the expression profiles: The samples' expression profiles are used to build the expression matrix where differential splicing changes will be simulated. This step is performed running the `R` script: `profilesSimulation.R` which is contained in the Simulation/scripts/ directory.

6-	Samples simulation: Each sample of each realization of the simulated experiment is generated using the `simulateReads.sh` script. The use of this script requires the specification of two ordered parameters, the first one is the number of sequencing reads to simulate (*N*) and the second one is the name of the sample to be simulated (*sampleName*).

7-	Pseudoalignment and quantification of isoforms: This step is performed running the `processing.sh` script, which analyzes each sample in the *samplesNames* file with `Kallisto`. Note that you must run this step for each realization of the simulated experiment or you could put into a `for` cycle.

8-	Differential splicing analysis: It is performed running the `DSAnalysis.R` script stored in the `R_scripts` directory. 

In order to simulate an experimental scenario replicated *N* times, we suggest creating a file structure like the following

    * SIM/
       * REPLICATION_I/
       * ...
       * REPLICATION_J/
       * ...
       * REPLICATION_N/
    
The pipeline explained above allows the generation of one replication (*J*) of one scenario. 


