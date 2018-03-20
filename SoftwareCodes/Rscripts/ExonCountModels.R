#############################################################################
## Differential splicing analysis with packages based on exon-count models ##
#############################################################################

# run one time for each simulation

library("NBSplice");library("DEXSeq");library("edgeR");library("limma")
library("BiocParallel")
BPPARAM<-MulticoreParam(20)
#setting path to files
setwd("/path_to_quantification_exonLev_DEXSeq/")
# 
g2t<-read.delim("/path2geneIsoFile/gene2transc.txt", header=FALSE)
colnames(g2t)<-c("gene_id", "transcript_id")

## Data preparation
# load exon counts from DEXSeq python script
exon_files<-paste("/htseq_sim_SRR0576",
    c(49:52, 31, 43, 45, 48), "/SRR0576", c(49:52, 31, 43, 45, 48), 
    ".htseq.counts", sep="")

exon_res<-lapply(exon_files, function(iso_f){
  return(read.delim(iso_f, header=FALSE))
  })
names(exon_res)<-sapply(exon_files,function(iso_f){
  return(strsplit(iso_f, split="[.]")[[1]][1])
  })

#remove last quantification result rows
exon_cm<-as.data.frame(do.call(cbind, lapply(1:length(exon_res), function(x){
 return(exon_res[[x]][,2][-(554628:554632)])

})))

rownames(exon_cm)<-exon_res[[1]][-(554628:554632),1]
colnames(exon_cm)<-c("C1R1","C1R2","C1R3","C1R4","C2R1","C2R2","C2R3","C2R4")

# We have ten realizations of the RNA-seq experiment, set the working directory
# to one of those
setwd("/path_to_i_th_experiment_realization")

# filtering exons with 0 counts in all samples
indexZero<-which(rowSums(exon_cm) ==0)
length(indexZero)
# [1] 184404

exon_cm<-exon_cm[-indexZero,]
dim(exon_cm)
# [1] 370223     8

gene<-do.call(rbind,strsplit(rownames(exon_cm), split="[:]"))[,1]
geneExon<-data.frame(gene_id=gene, exon=rownames(exon_cm))

# Number of expressed genes
genes<-as.character(unique(geneExon$gene_id))
genes<-sort(genes)
length(genes)
# [1] 23910

# Only considering simulated genes
iso_info_sim<-read.delim("/path_to_simulation_results/iso_info.tab", header=T,
    stringsAsFactor=FALSE)
gene_info_original<-read.delim("/path_to_simulation_results/gene_info.tab", 
    header=T, stringsAsFactor=F)

iso_info<-iso_info_sim[iso_info_sim$gene_id %in% genes,c(1:4,27,19:24,11,12)]

## Differential splicing analysis

# Design Matrix
designMatrix<-data.frame(samples=colnames(exon_cm), condition=factor(c(rep(
    "Normal",4), rep("Tumor",4)), levels=c("Normal", "Tumor")))
rownames(designMatrix)<-designMatrix$samples

exon_cm<-round(exon_cm)


## DEXSeq
cm<-DEXSeqDataSet( countData=exon_cm, sampleData=designMatrix, design= ~ sample + exon + condition:exon, 
    featureID=geneExon$exon, groupID=geneExon$gene_id)
cm<-estimateSizeFactors(cm)
cm<-estimateDispersions(cm, maxit=500, BPPARAM=BPPARAM)
fullModel<- ~ sample + exon + condition:exon
reducedModel<- ~ sample + exon
cm<-testForDEU(cm, fullModel=fullModel, reducedModel=reducedModel, BPPARAM=BPPARAM)

cm<-estimateExonFoldChanges(cm, fitExpToVar="condition", BPPARAM=BPPARAM, denominator="Normal")
myresults<- DEXSeqResults( cm )
myresultsDF<-as.data.frame(myresults)
#
myresultsDF<-myresultsDF[!is.na(myresultsDF$padj),]
perGeneQ<-perGeneQValue(myresults)

myresultsDF$qvalGene<-do.call(c, bplapply(1:nrow(myresultsDF), function(i){

return(perGeneQ[names(perGeneQ) == myresultsDF$groupID[i]])

}, BPPARAM=BPPARAM))

GenesDEX<-(unique(myresultsDF$groupID[myresultsDF$qvalGene < 0.05]))

## Limma
# Filter low-expressed exons using NBSplice filter method
myExonDataSet<-IsoDataSet(exon_cm, designMatrix, colName="condition", geneExon, BPPARAM)

ratioThres<-0.05
countThres<-1

myExonDataSet<-buildLowExpIdx(myExonDataSet, "condition", ratioThres, countThres, BPPARAM)

exon_cm_express<-counts(myExonDataSet)

geneExon2<-geneExon[rownames(exon_cm_express),]

y.all <- DGEList(counts=exon_cm_express, genes=geneExon2)
y.all<- calcNormFactors(y.all)
condition<-factor(c(rep("Normal",4), rep("Tumor",4)), levels=c("Normal", "Tumor"))
design <- model.matrix(~ condition)
v <- voom(y.all,design,plot=FALSE)
fit <- lmFit(v, design)
ex <- diffSplice(fit[,"conditionTumor"], geneid = "gene_id", exonid = "exon")
DSRes<-topSplice(ex, test="simes", n=length(ex))
ASGenes<-DSRes[DSRes$FDR < 0.05,]
LimmaGenes<-ASGenes$gene_id

##EDGER

y <- estimateDisp(y.all, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)
qlf <- glmQLFTest(fit, coef=2)
qlf<- diffSpliceDGE(fit, coef=2, geneid="gene_id", exonid="exon")
DSResEdgeR<-topSpliceDGE(qlf,  n=length(ex))
ASGenesEdgeR<-DSResEdgeR[DSResEdgeR$FDR < 0.05,]
edgeRGenes<-ASGenesEdgeR$gene_id

#setting path to files
setwd("/path_to_quantification_Kallisto/")
# 
geneIso<-read.delim("/path2geneIsoFile/gene2transc.txt", header=FALSE)
colnames(geneIso)<-c("gene_id", "isoform_id")
###################################################
### Differential splicing analysis with NBSplice ##
###################################################

## Data preparation

# isoforms
iso_files<-c(paste("kallisto_SRR0576", c(49:52,31,43,45,48), "/abundance.tsv",
    sep=""))
    
iso_res<-lapply(iso_files, function(iso_f){
  return(read.delim(iso_f))
  })
names(iso_res)<-sapply(iso_files,function(iso_f){
    return(strsplit(strsplit(iso_f[1], split="_")[[1]][2], split="/")[[1]][1])
})

iso_cm<-as.data.frame(do.call(cbind, lapply(1:length(iso_res), function(x){
  iso_sample<-iso_res[[x]]
  return(iso_sample$est_counts)

})))
rownames(iso_cm)<-iso_res[[1]]$target_id
colnames(iso_cm)<-c("C1R1","C1R2","C1R3","C1R4","C2R1","C2R2","C2R3","C2R4")
geneIso<-g2t[(g2t$transcript_id%in% rownames(iso_cm)),]
colnames(geneIso)[2]<-"isoform_id"

setwd("/path_to_i_th_experiment_realization")

# Expressed genes
genes<-as.character(unique(geneIso$gene_id))
genes<-sort(genes)

## ONLY CONSIDER SIMULATED GENES## 
iso_info_sim<-read.delim("/path_to_simulation_results/iso_info.tab", header=T,
    stringsAsFactor=FALSE)

iso_info<-iso_info_sim[iso_info_sim$transcript_id %in% rownames(iso_cm),
    c(1:4,27,19:24,11,12)]
iso_cm<-iso_cm[which(rownames(iso_cm) %in% iso_info_sim$transcript_id ),]
iso_cm<-iso_cm[as.character(iso_info$transcript_id ),]
geneIso<-geneIso[(geneIso$isoform_id%in% rownames(iso_cm)),]


## NBSplice
designMatrix<-data.frame(samples=colnames(iso_cm), condition=factor(c(rep(
    "Normal",4), rep("Tumor",4)), levels=c("Normal", "Tumor")))
rownames(designMatrix)<-designMatrix$samples

myIsoDataSet<-IsoDataSet(iso_cm,designMatrix, colName="condition", geneIso, BPPARAM)
myIsoDataSet<-buildLowExpIdx(myIsoDataSet, BPPARAM)

myNBResults<-NBTest(myIsoDataSet, test="F", BPPARAM=BPPARAM)

NBGenes<-GetDSGenes(myNBResults)


