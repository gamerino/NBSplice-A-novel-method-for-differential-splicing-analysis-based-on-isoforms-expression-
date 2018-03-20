# NBSplice: A novel method for differential splicing analysis based on isoforms expression.
Code for generation of simulated RNA-seq datasets, NBSplice evaluation and comparison with commonly-used R packages for differential splicing detection.

Merino, G.A. & Fernández, E.A. (2018). NBSplice: A novel method for differential splicing analysis based on isoforms expression.

The aim of our work was to develop an R package able to predict isoform relative expression and its change to infer differential gene alternative splicing. Our study presented the evaluation of the developed tool, based on simulated RNA-seq dataset generated from a real experiment. 

The structure of this repository is as follows:

- SoftwareCodes
  - bash_scripts: Directory containing the scripts used for processing RNA-seq data 
  - R_scripts: Directory having the R scripts used to perform differential splicing analysis

- AnalysisCode
  - R_scripts: Directory containing the code used to perform the comparative analysis

Each directory contain a README and a HOWTO file. 

The fastq files corresponding to the human samples can be downloaded from ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP002/SRP002628

- The normal samples used are: 
    - SRR057649
    - SRR057650
    - SRR057651
    - SRR057652
    
- The tumor samples used are:
    - SRR057631
    - SRR057643
    - SRR057645
    - SRR057648

The human reference files can be donloaded from: 

   * Genome FASTA file,  ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
    
   * Transcriptome FASTA file, ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz
    
   * Annotation GTF file, ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz

Installation of the following software is necessary:

- SRA toolkit (http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software)

- STAR (https://github.com/alexdobin/STAR)

- KALLISTO (https://github.com/pachterlab/kallisto)

Installation of the following R packages is required:

- NBSplice (https://github.com/gamerino/NBSplice)
- Limma (http://bioconductor.org/packages/Limma/)
- DEXSeq (http://bioconductor.org/packages/DEXSeq/)
- edgeR (https://bioconductor.org/packages/edgeR/)
- BiocParallel (http://bioconductor.org/packages/BiocParallel/)








