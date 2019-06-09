###############################################
## Simulated counts analysis and exploration ##
###############################################
stop()
# run one time for each simulation
library("ggplot2")
library("gridExtra")
library(BiocParallel)
library(tximport)
library("NBSplice")
library("DEXSeq")
library("DRIMSeq")
library(iCOBRA)
library("parallel")
# setting the amount of simulated replications
n<-10
setwd("path_to_simulated_samples")
# file having the gene-transcript relations.
g2t<-read.delim("gene2transc.txt", header=FALSE)
colnames(g2t)<-c("gene_id", "transcript_id")
# for each simulated replication do
mclapply(1:n, function(n.sim){
    setwd(paste("sim", n.sim,"/", sep=""))
    # load the quantification files for all samples
    iso_files<-paste("kallisto_SRR0576", c(49:52,31,43,45,48),"/abundance.h5",sep="")
    txi <- tximport(iso_files, type="kallisto", txOut=TRUE,countsFromAbundance="no")
    iso_cm<-txi$counts
    colnames(iso_cm)<-c("C1R1","C1R2","C1R3","C1R4","C2R1","C2R2","C2R3","C2R4")
    # removing zero counts transcripts
    indexZero<-which(rowSums(iso_cm) ==0)
    iso_cm<-iso_cm[-indexZero,]
    g2t_ok<-g2t[(g2t$transcript_id%in% rownames(iso_cm)),]
    g2t_ok<-g2t_ok[order(g2t_ok$gene_id, g2t_ok$transcript_id),]
    genes<-as.character(unique(g2t_ok$gene_id))    
    iso_cm<-iso_cm[as.character(g2t_ok$transcript_id),]
    # loading the simulated expression matrix
    iso_info_sim<-read.delim(paste("../preparingCM/sim", n.sim,"/iso_info.tab", sep=""), header=T, stringsAsFactor=FALSE)
    gene_info_original<-read.delim(paste("../preparingCM/sim", n.sim,"/gene_info.tab", sep=""), header=T, stringsAsFactor=F)
    # conserving the samples of interest and the non-filtered transcripts
    iso_info<-iso_info_sim[iso_info_sim$transcript_id %in% rownames(iso_cm),c(1:4,27,19:24,11,12)]
    iso_cm<-iso_cm[which(rownames(iso_cm) %in% iso_info_sim$transcript_id ),]
    iso_cm<-iso_cm[as.character(iso_info$transcript_id ),]
    g2t_ok<-g2t_ok[as.character(g2t_ok$transcript_id) %in% rownames(iso_cm),]
    # computing average expression for each isoform along replicates
    iso_info$iso_meanC1<-rowMeans(iso_cm[,1:4])
    iso_info$iso_meanC2<-rowMeans(iso_cm[,5:8])
    # identifiyng the genes simulated as differentially spliced
    DSgenes<-unique(iso_info[ iso_info$DS  | iso_info$DIEDS , "gene_id"])
    id_mayorDS<-as.character(g2t_ok$transcript_id[g2t_ok$gene_id %in% DSgenes & iso_info$mayor])
    # saving the full-matrix containing information about simulated status of each gene and its isoforms of replication n.sim
    save(iso_info, file=paste("iso_info_final_sim", n.sim,".RData", sep=""), compress="xz")
    # saving the simulated expression matrix of replication n.sim
    save(iso_cm, file=paste("iso_cm_sim", n.sim,".RData", sep=""), compress="xz")
    # saving the gene-isoform relation matrix of replication n.sim
    save(g2t_ok, file=paste("gene_iso_sim", n.sim,".RData", sep=""), compress="xz") 
}, mc.cores=8) # change 'mc.cores' value to set another number of CPUs to be used.
