###############################################
## Simulated counts analysis and exploration ##
###############################################
stop()
# run one time for each simulation
library("ggplot2")
library("gridExtra")
library(BiocParallel)
#setting path to files
setwd("~/prostateDS2/simulation/ExonCounts/")
# # gene2transcript
g2t<-read.delim("../preparingCM/gene2transc.txt", header=FALSE)
colnames(g2t)<-c("gene_id", "transcript_id")
# # gene
## change sim number

# isoforms
exon_files<-paste("../",c("sim1/htseq_sim_SRR057649/SRR057649.htseq.counts","sim1/htseq_sim_SRR057650/SRR057650.htseq.counts","sim1/htseq_sim_SRR057651/SRR057651.htseq.counts","sim1/htseq_sim_SRR057652/SRR057652.htseq.counts","sim1/htseq_sim_SRR057631/SRR057631.htseq.counts","sim1/htseq_sim_SRR057643/SRR057643.htseq.counts","sim1/htseq_sim_SRR057645/SRR057645.htseq.counts","sim1/htseq_sim_SRR057648/SRR057648.htseq.counts"), sep="")

exon_res<-lapply(exon_files, function(iso_f){
  return(read.delim(iso_f, header=FALSE))
  })
names(exon_res)<-sapply(exon_files,function(iso_f){
  return(strsplit(iso_f, split="[.]")[[1]][1])
  })


exon_cm<-as.data.frame(do.call(cbind, lapply(1:length(exon_res), function(x){
 return(exon_res[[x]][,2][-(554628:554632)])

})))

rownames(exon_cm)<-exon_res[[1]][-(554628:554632),1]
colnames(exon_cm)<-c("C1R1","C1R2","C1R3","C1R4","C2R1","C2R2","C2R3","C2R4")

setwd("./sim1")

# PCA plot
pr_comp_y<-prcomp(t((exon_cm)))
resumen <- summary(pr_comp_y)
labX <- signif((resumen$importance[2,1])*100, 3)
labY <- signif((resumen$importance[2,2])*100, 3)
PCAdata<-cbind(as.data.frame(pr_comp_y$x), condition=factor(c(rep("Normal",4), rep("Tumor",4)), levels=c("Normal", "Tumor")))
# df_plot<-data.frame(counts=rbind(exon_cm[,1],exon_cm[,2],exon_cm[,3],exon_cm[,4],exon_cm[,5],exon_cm[,6],exon_cm[,7],exon_cm[,8]), sample=c(rep(colnames(exon_cm)[1], nrow(exon_cm)),rep(colnames(exon_cm)[2], nrow(exon_cm)),rep(colnames(exon_cm)[3], nrow(exon_cm)),rep(colnames(exon_cm)[4], nrow(exon_cm)),rep(colnames(exon_cm)[5], nrow(exon_cm)),rep(colnames(exon_cm)[6], nrow(exon_cm)),rep(colnames(exon_cm)[7], nrow(exon_cm)),rep(colnames(exon_cm)[8], nrow(exon_cm))))

x11(type="cairo")
ggplot(PCAdata, aes(x=PC1, y=PC2, fill=condition, color=condition))+labs(title="PCA of isoform counts", x=paste("PC1 (",labX,"%)", sep=""), y=paste("PC2 (",labY,"%)", sep=""))+geom_point(aes(fill = condition), size = 10)

ggplot(exon_cm)+geom_boxplot(aes(x="C1R1",y=log2(C1R1+1), fill=1))+geom_boxplot(aes(x="C1R2",y=log2(C1R2+1), fill=2))+geom_boxplot(aes(x="C1R3",y=log2(C1R3+1), fill=3))+ geom_boxplot(aes(x="C1R4",y=log2(C1R4+1), fill=4))+geom_boxplot(aes(x="C2R1",y=log2(C2R1+1), fill=100))+geom_boxplot(aes(x="C2R2",y=log2(C2R2+1), fill=101))+ geom_boxplot(aes(x="C2R3",y=log2(C2R3+1), fill=102))+geom_boxplot(aes(x="C2R4",y=log2(C2R4+1), fill=103))+scale_y_continuous(limits = c(0,8))+xlab("samples")+ ylab("log2(counts+1)")+scale_fill_continuous()+guides(fill=FALSE) 
# ggsave(g,file="comparisson_of_boxplot_distributions.png", height=6, width=10)

# filtering isoforms with 0 counts in all samples
indexZero<-which(rowSums(exon_cm) ==0)
length(indexZero)
# [1] 184404

# # filter
exon_cm<-exon_cm[-indexZero,]
dim(exon_cm)
# [1] 370223     8

gene<-do.call(rbind,strsplit(rownames(exon_cm), split="[:]"))[,1]
geneExon<-data.frame(gene_id=gene, exon=rownames(exon_cm))

pr_comp_y<-prcomp(t((exon_cm)))
resumen <- summary(pr_comp_y)
labX <- signif((resumen$importance[2,1])*100, 3)
labY <- signif((resumen$importance[2,2])*100, 3)
PCAdata<-cbind(as.data.frame(pr_comp_y$x), condition=factor(c(rep("Normal",4), rep("Tumor",4)), levels=c("Normal", "Tumor")))

x11(type="cairo")
ggplot(PCAdata, aes(x=PC1, y=PC2, fill=condition, color=condition))+labs(title="PCA of isoform counts", x=paste("PC1 (",labX,"%)", sep=""), y=paste("PC2 (",labY,"%)", sep=""))+geom_point(aes(fill = condition), size = 10)

g<-ggplot(exon_cm)+geom_boxplot(aes(x="C1R1",y=log2(C1R1), fill=1))+geom_boxplot(aes(x="C1R2",y=log2(C1R2), fill=2))+geom_boxplot(aes(x="C1R3",y=log2(C1R3), fill=3))+ geom_boxplot(aes(x="C1R4",y=log2(C1R4), fill=4))+geom_boxplot(aes(x="C2R1",y=log2(C2R1), fill=100))+geom_boxplot(aes(x="C2R2",y=log2(C2R2), fill=101))+ geom_boxplot(aes(x="C2R3",y=log2(C2R3), fill=102))+geom_boxplot(aes(x="C2R4",y=log2(C2R4), fill=103))+scale_y_continuous(limits = c(0,8))+xlab("samples")+ ylab("log2(counts)")+scale_fill_continuous()+guides(fill=FALSE) 
x11(type="cairo");g

# Number of expressed genes
genes<-as.character(unique(geneExon$gene_id))
genes<-sort(genes)
length(genes)
# [1] 23910

## ONLY CONSIDER SIMULATED GENES## 
iso_info_sim<-read.delim("../../preparingCM/sim1/iso_info.tab", header=T, stringsAsFactor=FALSE)
gene_info_original<-read.delim("../../preparingCM/sim1/gene_info.tab", header=T, stringsAsFactor=F)

iso_info<-iso_info_sim[iso_info_sim$gene_id %in% genes,c(1:4,27,19:24,11,12)]

##  exploring target genes
DSgenes<-unique(iso_info_sim[ iso_info_sim$DS  | iso_info_sim$DIEDS , "gene_id"])
length(DSgenes)
# [1] 1680

# deberian ser 1680
table(unique(iso_info$gene_id) %in% DSgenes)
# FALSE  TRUE 
# 14085  1637 


############################################################################
library(BiocParallel)
library(MASS)
library(car)
library(mppa)

sampleData<-data.frame(samples=colnames(exon_cm), condition=factor(c(rep("Normal",4), rep("Tumor",4)), levels=c("Normal", "Tumor")))


##DEXSeq
library(DEXSeq)
cm<-DEXSeqDataSet( countData=round(exon_cm), sampleData, design= ~ sample + exon + condition:exon, featureID=geneExon$exon, groupID=geneExon$gene_id)
cm<-estimateSizeFactors(cm)
cm<-estimateDispersions(cm, maxit=500, BPPARAM=MulticoreParam(workers=18))
fullModel<- ~ sample + exon + condition:exon
reducedModel<- ~ sample + exon
cm<-testForDEU(cm, fullModel=fullModel, reducedModel=reducedModel, BPPARAM=MulticoreParam(workers=18))

cm<-estimateExonFoldChanges(cm, fitExpToVar="condition", BPPARAM=MulticoreParam(workers=18), denominator="Normal")
myresults<- DEXSeqResults( cm )
myresultsDF<-as.data.frame(myresults)
#
myresultsDF<-myresultsDF[!is.na(myresultsDF$padj),]
table(myresultsDF$padj < 0.05)
# FALSE  TRUE 
#330973  13290 

perGeneQ<-perGeneQValue(myresults)

myresultsDF$qvalGene<-do.call(c, bplapply(1:nrow(myresultsDF), function(i){

return(perGeneQ[names(perGeneQ) == myresultsDF$groupID[i]])

}, BPPARAM=MulticoreParam(10)))

GenesDEX<-(unique(myresultsDF$groupID[myresultsDF$qvalGene < 0.05]))

table(GenesDEX %in% unique(iso_info$gene_id[ iso_info$DIEDS | iso_info$DS]))
# FALSE  TRUE 
#     681  1194 
###################
exon_cm<-round(exon_cm)
designMatrix<-data.frame(samples=colnames(exon_cm), condition=factor(c(rep("Normal",4), rep("Tumor",4)), levels=c("Normal", "Tumor")))
rownames(designMatrix)<-designMatrix$samples
BPPARAM=MulticoreParam(20)
geneExon2<-geneExon
colnames(geneExon2)[2]<-"isoform_id"
rownames(geneExon2)<-rownames(exon_cm)
totalCounts<-totalGeneCounts(exon_cm, geneExon2, BPPARAM)

ratioThres<-0.05
countThres<-1

idxLowRat<-filterLowRatio(iso_cm=cbind(exon_cm, totalCounts), designMatrix, ratioThres, countThres, BPPARAM)

exon_cm_express<-exon_cm[-idxLowRat,]

geneExon2<-geneExon2[rownames(exon_cm_express),]
##Limma
library(limma);library(edgeR)
y.all <- DGEList(counts=exon_cm_express, genes=geneExon2)
y.all<- calcNormFactors(y.all)
condition<-factor(c("Normal","Normal","Normal","Normal","Tumor","Tumor","Tumor","Tumor"), levels=c("Normal", "Tumor"))
design <- model.matrix(~ condition)
v <- voom(y.all,design,plot=FALSE)
fit <- lmFit(v, design)
ex <- diffSplice(fit[,"conditionTumor"], geneid = "gene_id", exonid = "exon")
# Total number of exons:  35363 
# Total number of genes:  12861 
# Number of genes with 1 exon:  4547 
# Mean number of exons in a gene:  3 
# Max number of exons in a gene:  13 

DSRes<-topSplice(ex, test="simes", n=length(ex))
ASGenes<-DSRes[DSRes$FDR < 0.05,]
table(ASGenes$gene_id %in% unique(iso_info$gene_id[iso_info$DIEDS | iso_info$DS]))
# FALSE  TRUE 
#      18   142 

LimmaGenes<-ASGenes$gene_id

##EDGER

y <- estimateDisp(y.all, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)
qlf <- glmQLFTest(fit, coef=2)
qlf<- diffSpliceDGE(fit, coef=2, geneid="gene_id", exonid="exon")
# Total number of exons:  15376 
# Total number of genes:  6880 
# Number of genes with 1 exon:  3642 
# Mean number of exons in a gene:  2 
# Max number of exons in a gene:  13 
DSResEdgeR<-topSpliceDGE(qlf,  n=length(ex))
ASGenesEdgeR<-DSResEdgeR[DSResEdgeR$FDR < 0.05,]
table(unique(ASGenesEdgeR$gene_id) %in% unique(iso_info$gene_id[ iso_info$DIEDS | iso_info$DS]))
# FALSE  TRUE 
#     5    54 

edgeRGenes<-ASGenesEdgeR$gene_id

##SpliceVariants
#spliceVariants
y.all <- DGEList(counts=exon_cm_express, genes= geneExon2,group=condition )
# Max number of exons in a gene:  10
SV<-spliceVariants(y.all, geneExon2$gene_id)
SVResults<-SV$table
SVResults$FDR<-p.adjust(SVResults$PValue, "fdr")
table(SVResults$FDR < 0.05)
# FALSE  TRUE 
#12666   195 

table(rownames(SVResults[SVResults$FDR < 0.05,]) %in% unique(iso_info$gene_id[ iso_info$DIEDS | iso_info$DS]))
# 
# FALSE  TRUE 
#    29   166 

edgeRSVGenes<-rownames(SVResults[SVResults$FDR < 0.05,])

####
exon_cpm<-round(t(t(exon_cm)/colSums(exon_cm)*1000000))
geneExon2<-geneExon
colnames(geneExon2)[2]<-"isoform_id"
rownames(geneExon2)<-rownames(exon_cpm)

totalCounts2<-totalGeneCounts(exon_cpm, geneExon2, BPPARAM)


ratioThres<-0.05
countThres<-1

idxLowRat2<-filterLowRatio(iso_cm=cbind(exon_cpm, totalCounts2), designMatrix, ratioThres, countThres, BPPARAM)


myResults<-NBRTest(exon_cpm, totalCounts2,idxLowRat2, geneExon2, designMatrix,family="negativebinomial",test="F", BPPARAM=BPPARAM)

DEISO<-as.character(myResults$iso[which(myResults$FDR < 0.05 )])
NBGenes<-unique(myResults[myResults$geneFDR< 0.05,"gene"])

table(NBGenes %in%  unique(iso_info$gene_id[iso_info$DS | iso_info$DIEDS]))
# FALSE  TRUE 
#    5   125 

myResultsQB<-NBRTest(exon_cpm, totalCounts2,idxLowRat2, geneExon2, designMatrix,family="quasibinomial",test="F", BPPARAM=BPPARAM)

DEISOQB<-as.character(myResultsQB$iso[which(myResultsQB$FDR < 0.05 )])
QBNGenes<-unique(myResultsQB[myResultsQB$geneFDR< 0.05,"gene"])

table(QBNGenes %in%  unique(iso_info$gene_id[iso_info$DS | iso_info$DIEDS]))
# FALSE  TRUE 
#         275   699 
save.image("sim1IsoCounts.RData", compress="xz")

