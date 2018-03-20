# profilesSimulation.R
# This script was designed to perform ten simulations of one experimental 
# scenario were 4 replicates per condition and 10% of genes with differential
# splicing were considered.  
library("ggplot2");library("gridExtra");library("BiocParallel")
#setting path to files
setwd("/path_to_simulation_scenario")
# read file containing gene/isoforms definition
g2t<-read.delim("/path2geneIsoFile/gene2transc.txt", header=FALSE)
head(g2t,3)
#             gene_id     transcript_id
# 92274 ENSG00000000003 ENST00000373020
# 92275 ENSG00000000003 ENST00000496771
# 92276 ENSG00000000003 ENST00000494424
# gene
## load RSEM quantification for all samples generated in step 3
#gene counts
gene_files<-paste("SRR0576",c(49:52,31,43,45,48), ".genes.results",  sep="")
gene_res<-lapply(gene_files, function(gene_f){
  return(read.delim(gene_f))
})
names(gene_res)<-sapply(gene_files,function(gene_f){
  return(strsplit(gene_f, split="[.]")[[1]][1])
})
gene_cm<-as.data.frame(do.call(cbind, lapply(1:length(gene_res), function(y){
  gene_sample<-gene_res[[y]]
  return(gene_sample$expected_count)
})))
rownames(gene_cm)<-gene_res[[1]]$gene_id
colnames(gene_cm)<-names(gene_res)
gene_tpm<-as.data.frame(do.call(cbind, lapply(1:length(gene_res), function(y){
  gene_sample<-gene_res[[y]]
  return(gene_sample$TPM )
})))
rownames(gene_tpm)<-gene_res[[1]]$gene_id
colnames(gene_tpm)<-names(gene_res)
dim(gene_cm)
# [1] 44464    8
# isoform counts
iso_files<-paste("SRR0576",c(49:52,31,43,45,48), ".isoforms.results", sep="")
iso_res<-lapply(iso_files, function(iso_f){
  return(read.delim(iso_f))
})
names(iso_res)<-sapply(iso_files,function(iso_f){
  return(strsplit(iso_f, split="[.]")[[1]][1])
})
iso_cm<-as.data.frame(do.call(cbind, lapply(1:length(iso_res), function(x){
  iso_sample<-iso_res[[x]]
  return(iso_sample$expected_count)
})))
rownames(iso_cm)<-iso_res[[1]]$transcript_id
colnames(iso_cm)<-c(paste("C1R", c(1:4), sep=""), paste("C2R", c(1:4), sep=""))
iso_tpm<-as.data.frame(do.call(cbind, lapply(1:length(iso_res), function(x){
  iso_sample<-iso_res[[x]]
  return(iso_sample$TPM)
})))
rownames(iso_tpm)<-iso_res[[1]]$transcript_id
colnames(iso_tpm)<-c(paste("C1R", c(1:4), sep=""), paste("C2R", c(1:4), sep=""))
dim(iso_cm)
# [1] 169481     8

##Filter not-expressed genes
index_g<-unique(do.call(c,sapply(1:nrow(gene_tpm), function(x){
   if(any(gene_tpm[x,1:4] ==0) | any(gene_tpm[x,5:8] ==0) ) {	
    return(x)
    }else return(NULL)
})))
length(index_g)

# [1] 28280
# # filter
gene_cm<-gene_cm[-index_g,]
dim(gene_cm)
# [1] 16184    8
gene_tpm<-gene_tpm[-index_g,]
dim(gene_tpm)
# [1] 16184    8

g2t_ok<-g2t[(g2t$gene_id%in% rownames(gene_cm)),]
iso_cm<-iso_cm[rownames(iso_cm) %in% g2t_ok$transcript_id,]
dim(iso_cm)
# [1] 119965     8
iso_tpm<-iso_tpm[rownames(iso_tpm) %in% rownames(iso_cm),]

pr_comp_y<-prcomp(iso_cm)
resumen <- summary(pr_comp_y)
labX <- signif((resumen$importance[2,1])*100, 3)
labY <- signif((resumen$importance[2,2])*100, 3)
PCAdata<-cbind(as.data.frame(pr_comp_y$x), condition=factor(c(rep("Normal",4), 
    rep("Tumor",4)), levels=c("Normal", "Tumor")))

g<-ggplot(PCAdata, aes(x=PC1, y=PC2, fill=condition, color=condition))+labs(
    title="PCA of isoform counts", x=paste("PC1 (",labX,"%)", sep=""), y=paste(
    "PC2 (",labY,"%)", sep=""))+geom_point(aes(fill = condition), size = 10)
biplot(pr_comp_y,cex=c(1,0.001),xlabs=colnames(iso_cm),ylabs=rep("",
    nrow(iso_cm)),var.axes=FALSE)
# x11(type="cairo")
#1 ggplot(log2(iso_cm+1), aes(as.factor(colnames(iso_cm))))+geom_violin()


#Characterization of the number of isoforms per gene
genes<-as.character(unique(g2t_ok$gene_id))
genes<-sort(genes)
numb_iso<-do.call(c, bplapply(genes, function(gen){
  return(length(g2t[g2t$gene_id == gen,"transcript_id"]))
  }, BPPARAM=MulticoreParam(20)))
#construct a table with gene information
  
gene_info<-data.frame(gene_id=genes, numb_iso=numb_iso)
freq<-data.frame(table(gene_info$numb_iso[!gene_info$numb_iso ==1]), accum=
    cumsum(table(gene_info$numb_iso[!gene_info$numb_iso ==1]))/sum(table(
    gene_info$numb_iso[!gene_info$numb_iso ==1])))
colnames(freq)<-c("Iso_num", "freq", "cum_rel_freq")
freq[,"score"]<-cut(freq[,"cum_rel_freq"], breaks=c(0,0.33,0.67,1),
    include.lowest=FALSE, right=TRUE)
# x11(type="cairo")

g1<-ggplot(freq,aes(x=Iso_num, y=freq, fill=as.factor(score)))+geom_bar(
    stat="identity")+guides(fill=guide_legend(title="intervals"))+theme(
    axis.text=element_text(size=12), axis.title=element_text(size=14), 
    legend.text=element_text(size=14))

g2<-ggplot(freq,aes(x=Iso_num, y=cum_rel_freq, fill=as.factor(score)))+
    geom_bar(stat="identity")+scale_fill_hue()+guides(fill=guide_legend(title=
    "intervals"))+theme(axis.text=element_text(size=12), axis.title=
    element_text(size=14), legend.text=element_text(size=14))
p<-arrangeGrob(g1,g2, nrow=1, ncol=2)
plot(p)
# genes were clustered in four groups according to the number of annotated transcripts #that they have: 
# "0” = 1 annotated transcript
# "I” = 2-4 annotated transcripts
# "II”= 5-9 annotated transcripts
# "III”= >9 annotated transcripts
# add the gene groups to the gene_info data.frame
gene_info$group<-cut(gene_info[,"numb_iso"], breaks=c(1,2,5,10,Inf), 
    include.lowest=TRUE, right=FALSE)
levels(gene_info$group)<-c("0","I", "II", "III")
# add quantification along replicates
gene_info$meanC1<-rowMeans(gene_cm[,1:4])
gene_info$meanC2<-rowMeans(gene_cm[,5:8])
gene_info$mean<-rowMeans(gene_cm)
gene_info$varC1<-apply(gene_cm[,1:4],1,var)
gene_info$varC2<-apply(gene_cm[,5:8],1,var)
gene_info$var<-apply(gene_cm,1,var)

# create a table with isoform information
iso_info<-merge(x=g2t_ok ,y=gene_info, by="gene_id", sort=FALSE)
iso_info<-iso_info[order(iso_info$gene_id, iso_info$transcript_id),
    c(2,1,3:ncol(iso_info))]
iso_info$iso_meanC1<-rowMeans(iso_cm[,1:4])
iso_info$iso_meanC2<-rowMeans(iso_cm[,5:8])
iso_info$iso_mean<-rowMeans(iso_cm)
iso_info$iso_varC1<-apply(iso_cm[,1:4],1,var)
iso_info$iso_varC2<-apply(iso_cm[,5:8],1,var)
iso_info$iso_var<-apply(iso_cm,1,var)
iso_info$ratiosC1<-iso_info$iso_meanC1/iso_info$meanC1
iso_info$ratiosC2<-iso_info$iso_meanC2/iso_info$meanC2
# determining the major isoform for each gene
id_mayor<-do.call(c, bplapply(1:length(unique(g2t_ok$gene_id)), function(x){
    maximo<-which.max(iso_info$iso_mean[iso_info$gene_id == as.character(unique(
        g2t_ok$gene_id)[x])]) 
    return(as.character(iso_info$transcript_id[iso_info$gene_id ==as.character(
        unique(g2t_ok$gene_id)[x])])[maximo])
}, BPPARAM=MulticoreParam(22)))
g2t_ok$mayor<-g2t_ok$transcript_id %in% id_mayor
ggplot(iso_info[iso_info$transcript_id %in% id_mayor,], aes(x=group, y=
    iso_mean/mean,fill=as.factor(group)))+geom_boxplot()+scale_fill_hue()+
    guides(fill=guide_legend(title="group"))+labs(title="")+theme(axis.text=
    element_text(size=12), axis.title=element_text(size=14), legend.text=
    element_text(size=14))


##  Selecting target genes simulated as differentially spliced
# We only consider those genes having at least a mean counts of 50 in at
# least one condition
index_exp<-unique(do.call(c,sapply(1:nrow(gene_tpm), function(x){
    if(mean(gene_tpm[x,1:4] >=50) | mean(gene_tpm[x,5:8] >=50) ) {	
    return(x)
    }else return(NULL)
})))
names_NSAgenes<-as.character(gene_info$gene_id[index_exp][gene_info[
    index_exp,]$group == "0"])
names_SAgenesgI<-as.character(gene_info$gene_id[index_exp][gene_info[
    index_exp,]$group == "I"])
names_SAgenesgII<-as.character(gene_info$gene_id[index_exp][gene_info[
    index_exp,]$group == "II"])
names_SAgenesgIII<-as.character(gene_info$gene_id[index_exp][ gene_info[

length(names_SAgenesgI)+length(names_SAgenesgII)+length(names_SAgenesgIII)
# [1] 4129
#But the total of simulated genes are
dim(gene_info)
#  16184     9
## We take the 10% of genes, 1680, to be DS. 
# We have three groups of genes to simulate differential splicing, we take 560
# from each group.
DSI <- sample( 1:length(names_SAgenesgI),560)
DSgI<-names_SAgenesgI[DSI]
DSII <- sample( 1:length(names_SAgenesgII),560)
DSgII<-names_SAgenesgII[DSII]
DSIII <- sample( 1:length(names_SAgenesgIII),560)
DSgIII<-names_SAgenesgIII[DSIII]
DS_index<-data.frame(DSgI=DSgI, DSgII=DSgII, DSgIII=DSgIII)
DSgenes<-c(DSgI, DSgII, DSgIII)

##  SIMULATION PROFILES

# Only Differential Splicing, we sample 480 of the 1680 DS genes 
DSg_index<-sample(1:length(DSgenes), 480)
DSg_genes<-DSgenes[DSg_index]
# Differential Splicing combined with differentia expression, 
# we take the remaining 1200 DS genes 
DIEDSg_genes<-DSgenes[-DSg_index]

# Add the group information at the gene info matrix 
gene_info$DS<-gene_info$gene_id %in% DSg_genes
gene_info$DIEDS<-gene_info$gene_id %in% DIEDSg_genes


## 1) Differential splicing
# simulated changes
ratios<-c("0.5-0.7","0.7-0.5","0.5-0.9","0.9-0.5","0.3-0.5", "0.5-0.3",
    "0.3-0.7","0.7-0.3")

# 120 genes for each ratio group
a<-sample(DSg_genes, 120)
b<-sample(DSg_genes[!DSg_genes %in% a],120)
d<-sample(DSg_genes[!DSg_genes %in% a & !DSg_genes %in% b],120)
e<-DSg_genes[!DSg_genes %in% a & !DSg_genes %in% b & !DSg_genes %in% d]
DS_genesdf<-data.frame(a=a, b=b,d=d, e=e)

gene_info$DSgroup<-0
gene_info$DSgroup[gene_info$gene_id %in% DS_genesdf$a]<-c(rep(ratios[1],60),
    rep(ratios[2], 60))
gene_info$DSgroup[gene_info$gene_id %in% DS_genesdf$b]<-(c(rep(ratios[3],60), 
    rep(ratios[4], 60)))
gene_info$DSgroup[gene_info$gene_id %in% DS_genesdf$d]<-(c(rep(ratios[5],60), 
    rep(ratios[6], 60)))
gene_info$DSgroup[gene_info$gene_id %in% DS_genesdf$e]<-(c(rep(ratios[7],60), 
    rep(ratios[8], 60)))

## 2) Differential isoforms expression and differential splicing
#folds and ratios simulated 
folds<-c("2-0.5-0.7","0.5-0.5-0.7","2-0.5-0.9", "0.5-0.5-0.9", "2-0.3-0.5",
    "0.5-0.3-0.5","2-0.3-0.7", "0.5-0.3-0.7", "2-0.5-0.3", "0.5-0.5-0.3", 
    "2-0.7-0.3", "0.5-0.7-0.3", "4-0.5-0.7", "0.25-0.5-0.7", "4-0.5-0.9", 
    "0.25-0.5-0.9", "4-0.7-0.5", "0.25-0.7-0.5", "4-0.9-0.5", "0.25-0.9-0.5")
# 240 genes for each fold-ratio comb 
a<-sample(DIEDSg_genes, 240)
b<-sample(DIEDSg_genes[!DIEDSg_genes %in% a],240)
d<-sample(DIEDSg_genes[!DIEDSg_genes %in% a & !DIEDSg_genes %in% b],240)
e<-sample(DIEDSg_genes[!DIEDSg_genes %in% a & !DIEDSg_genes %in% b & !DIEDSg_genes %in% d],240)
f<-DIEDSg_genes[!DIEDSg_genes %in% a & !DIEDSg_genes %in% b & !DIEDSg_genes %in% 
    d & !DIEDSg_genes %in% e]
DIEDSg_genesdf<-data.frame(a=a, b=b,d=d,e=e, f=f)
# 
gene_info$DIEDSgroup<-factor(rep(0, nrow(gene_info)), levels=c(0,folds))
gene_info$DIEDSgroup[gene_info$gene_id %in% DIEDSg_genesdf$a]<-c(rep(folds[1],60),
    rep(folds[2], 60), rep(folds[3], 60),rep(folds[4], 60))
gene_info$DIEDSgroup[gene_info$gene_id %in% DIEDSg_genesdf$b]<-c(rep(folds[5],60),
    rep(folds[6], 60), rep(folds[7], 60),rep(folds[8], 60))
gene_info$DIEDSgroup[gene_info$gene_id %in% DIEDSg_genesdf$d]<-c(rep(folds[9],60),
    rep(folds[10], 60), rep(folds[11], 60),rep(folds[12], 60))
gene_info$DIEDSgroup[gene_info$gene_id %in% DIEDSg_genesdf$e]<-c(rep(folds[13],60),
    rep(folds[14], 60), rep(folds[15], 60),rep(folds[16], 60))
gene_info$DIEDSgroup[gene_info$gene_id %in% DIEDSg_genesdf$f]<-c(rep(folds[17],60),
    rep(folds[18], 60), rep(folds[19], 60),rep(folds[20], 60))

##Creating an isoform level data matrix
iso_info<-merge(x=gene_info[,c(1,10:13)], y=iso_info, by="gene_id")
head(iso_info)
iso_info<-iso_info[,c(6,1,7:ncol(iso_info), 2:5)]

## Modifying distributional parameters values
iso_info$meanC1Sim<-iso_info$iso_meanC2
iso_info$meanC2Sim<-iso_info$iso_meanC2
iso_info$varC1Sim<-iso_info$iso_varC2
iso_info$varC2Sim<-iso_info$iso_varC2
iso_info$mayor<-(iso_info$transcript_id %in% g2t_ok$transcript_id[g2t_ok$mayor])
# we add an attribute flag indicating if the gene was up or down regulated
iso_info$updown<-0

# the genes that no are differential expressed should have the same distribution in both conditions

## 1) Differential splicing
#0.5-0.7
downgene<-as.character(iso_info$gene_id[iso_info$DSgroup == "0.5-0.7" & 
    iso_info$mayor])
# in this case the mean and variance that are scaled were computed along two conditions and for each GENE

iso_info$meanC1Sim[iso_info$DSgroup == "0.5-0.7" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DSgroup == "0.5-0.7" & iso_info$mayor]*0.5
iso_info$varC1Sim[iso_info$DSgroup == "0.5-0.7" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DSgroup == "0.5-0.7" & iso_info$mayor]*0.5^2

iso_info$meanC2Sim[iso_info$DSgroup == "0.5-0.7" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DSgroup == "0.5-0.7" & iso_info$mayor]*0.7
iso_info$varC2Sim[iso_info$DSgroup == "0.5-0.7" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DSgroup == "0.5-0.7" & iso_info$mayor]*0.7^2
iso_info$updown[iso_info$DSgroup == "0.5-0.7" & iso_info$mayor]<-1

# now should set the expression values to the rest!!
downinfo<-do.call(rbind,lapply(as.character(downgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DSgroup == "0.5-0.7" & 
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DSgroup == "0.5-0.7" & !iso_info$mayor 
    & iso_info$gene_id == gen]*(1-0.5)/(numb_iso-1)
  varC1_sim<-iso_info$varC2[iso_info$DSgroup == "0.5-0.7" & !iso_info$mayor & 
    iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1))^2

  meanC2_sim<-iso_info$meanC2[iso_info$DSgroup == "0.5-0.7" & !iso_info$mayor & 
    iso_info$gene_id == gen]*(1-0.7)/(numb_iso-1)
  varC2_sim<-iso_info$varC2[iso_info$DSgroup == "0.5-0.7" & !iso_info$mayor &
    iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1))^2
  
  index<-which(iso_info$DSgroup == "0.5-0.7" & !iso_info$mayor & 
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info[downinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-downinfo[,-1]
iso_info$updown[downinfo$index]<--1

#0.7-0.5

upgene<-as.character(iso_info$gene_id[iso_info$DSgroup == "0.7-0.5" & 
    iso_info$mayor]  )
iso_info$meanC1Sim[iso_info$DSgroup == "0.7-0.5" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DSgroup == "0.7-0.5" & iso_info$mayor]  *0.7
iso_info$varC1Sim[iso_info$DSgroup == "0.7-0.5" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DSgroup == "0.7-0.5" & iso_info$mayor]  *0.7^2
iso_info$meanC2Sim[iso_info$DSgroup == "0.7-0.5" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DSgroup == "0.7-0.5" & iso_info$mayor]  *0.5
iso_info$varC2Sim[iso_info$DSgroup == "0.7-0.5" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DSgroup == "0.7-0.5" & iso_info$mayor]  *0.5^2
iso_info$updown[iso_info$DSgroup == "0.7-0.5" & iso_info$mayor]<- -1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DSgroup == "0.7-0.5" & 
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DSgroup == "0.7-0.5" & !iso_info$mayor & 
    iso_info$gene_id == gen]*(1-0.7)/(numb_iso-1)
  varC1_sim<-iso_info$varC2[iso_info$DSgroup == "0.7-0.5" & !iso_info$mayor & 
    iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1))^2

  meanC2_sim<-iso_info$meanC2[iso_info$DSgroup == "0.7-0.5" & !iso_info$mayor &
    iso_info$gene_id == gen]*(1-0.5)/(numb_iso-1)
  varC2_sim<-iso_info$varC2[iso_info$DSgroup == "0.7-0.5" & !iso_info$mayor & 
    iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1))^2
  
  index<-which(iso_info$DSgroup == "0.7-0.5" & !iso_info$mayor &
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim, 
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info$updown[upinfo$index]<-1

iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]
  
#0.5-0.9
downgene<-as.character(iso_info$gene_id[iso_info$DSgroup == "0.5-0.9" & iso_info$mayor])

iso_info$meanC1Sim[iso_info$DSgroup == "0.5-0.9" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DSgroup == "0.5-0.9" & iso_info$mayor]*0.5
iso_info$varC1Sim[iso_info$DSgroup == "0.5-0.9" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DSgroup == "0.5-0.9" & iso_info$mayor]*0.5^2

iso_info$meanC2Sim[iso_info$DSgroup == "0.5-0.9" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DSgroup == "0.5-0.9" & iso_info$mayor]*0.9
iso_info$varC2Sim[iso_info$DSgroup == "0.5-0.9" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DSgroup == "0.5-0.9" & iso_info$mayor]*0.9^2
iso_info$updown[iso_info$DSgroup == "0.5-0.9" & iso_info$mayor]<-1

# now should set the expression values to the rest!!
downinfo<-do.call(rbind,lapply(as.character(downgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DSgroup == "0.5-0.9" & 
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DSgroup == "0.5-0.9" & !iso_info$mayor &
    iso_info$gene_id == gen]*(1-0.5)/(numb_iso-1)
  varC1_sim<-iso_info$varC2[iso_info$DSgroup == "0.5-0.9" & !iso_info$mayor &
    iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1))^2

  meanC2_sim<-iso_info$meanC2[iso_info$DSgroup == "0.5-0.9" & !iso_info$mayor &
    iso_info$gene_id == gen]*(1-0.9)/(numb_iso-1)
  varC2_sim<-iso_info$varC2[iso_info$DSgroup == "0.5-0.9" & !iso_info$mayor &
    iso_info$gene_id == gen]*((1-0.9)/(numb_iso-1))^2
  
  index<-which(iso_info$DSgroup == "0.5-0.9" & !iso_info$mayor &
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info[downinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim",
    "varC2Sim")]<-downinfo[,-1]
iso_info$updown[downinfo$index]<- -1

#0.9-0.5

upgene<-as.character(iso_info$gene_id[iso_info$DSgroup == "0.9-0.5" &
    iso_info$mayor]  )
iso_info$meanC1Sim[iso_info$DSgroup == "0.9-0.5" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DSgroup == "0.9-0.5" & iso_info$mayor]  *0.9
iso_info$varC1Sim[iso_info$DSgroup == "0.9-0.5" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DSgroup == "0.9-0.5" & iso_info$mayor]  *0.8^2
iso_info$meanC2Sim[iso_info$DSgroup == "0.9-0.5" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DSgroup == "0.9-0.5" & iso_info$mayor]  *0.5
iso_info$varC2Sim[iso_info$DSgroup == "0.9-0.5" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DSgroup == "0.9-0.5" & iso_info$mayor]  *0.5^2
iso_info$updown[iso_info$DSgroup == "0.9-0.5" & iso_info$mayor]<- -1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DSgroup == "0.9-0.5" & 
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DSgroup == "0.9-0.5" & !iso_info$mayor &
    iso_info$gene_id == gen]*(1-0.9)/(numb_iso-1)
  varC1_sim<-iso_info$varC2[iso_info$DSgroup == "0.9-0.5" & !iso_info$mayor & 
    iso_info$gene_id == gen]*((1-0.9)/(numb_iso-1))^2

  meanC2_sim<-iso_info$meanC2[iso_info$DSgroup == "0.9-0.5" & !iso_info$mayor &
    iso_info$gene_id == gen]*(1-0.5)/(numb_iso-1)
  varC2_sim<-iso_info$varC2[iso_info$DSgroup == "0.9-0.5" & !iso_info$mayor &
    iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1))^2
  
  index<-which(iso_info$DSgroup == "0.9-0.5" & !iso_info$mayor & 
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim, 
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info$updown[upinfo$index]<-1

iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-
    upinfo[,-1]
  
#0.3-0.5
table(gene_info$DSgroup=="0.3-0.5")
# FALSE  TRUE 
#16124    60 

downgene<-as.character(iso_info$gene_id[iso_info$DSgroup == "0.3-0.5" & 
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DSgroup == "0.3-0.5" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DSgroup == "0.3-0.5" & iso_info$mayor]*0.3
iso_info$varC1Sim[iso_info$DSgroup == "0.3-0.5" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DSgroup == "0.3-0.5" & iso_info$mayor]*0.3^2

iso_info$meanC2Sim[iso_info$DSgroup == "0.3-0.5" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DSgroup == "0.3-0.5" & iso_info$mayor]*0.5
iso_info$varC2Sim[iso_info$DSgroup == "0.3-0.5" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DSgroup == "0.3-0.5" & iso_info$mayor]*0.5^2
iso_info$updown[iso_info$DSgroup == "0.3-0.5" & iso_info$mayor]<-1

# now should set the expression values to the rest!!
downinfo<-do.call(rbind,lapply(as.character(downgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DSgroup == "0.3-0.5" &    
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DSgroup == "0.3-0.5" & !iso_info$mayor &
    iso_info$gene_id == gen]*(1-0.3)/(numb_iso-1)
  varC1_sim<-iso_info$varC2[iso_info$DSgroup == "0.3-0.5" & !iso_info$mayor &
    iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1))^2

  meanC2_sim<-iso_info$meanC2[iso_info$DSgroup == "0.3-0.5" & !iso_info$mayor & 
    iso_info$gene_id == gen]*(1-0.5)/(numb_iso-1)
  varC2_sim<-iso_info$varC2[iso_info$DSgroup == "0.3-0.5" & !iso_info$mayor &
    iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1))^2
  
  index<-which(iso_info$DSgroup == "0.3-0.5" & !iso_info$mayor &
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info[downinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-downinfo[,-1]
iso_info$updown[downinfo$index]<- -1

# 0.5-0.3

upgene<-as.character(iso_info$gene_id[iso_info$DSgroup == "0.5-0.3" & 
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DSgroup == "0.5-0.3" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DSgroup == "0.5-0.3" & iso_info$mayor]*0.5
iso_info$varC1Sim[iso_info$DSgroup == "0.5-0.3" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DSgroup == "0.5-0.3" & iso_info$mayor]*0.5^2
iso_info$meanC2Sim[iso_info$DSgroup == "0.5-0.3" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DSgroup == "0.5-0.3" & iso_info$mayor]*0.3
iso_info$varC2Sim[iso_info$DSgroup == "0.5-0.3" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DSgroup == "0.5-0.3" & iso_info$mayor]*0.3^2
iso_info$updown[iso_info$DSgroup == "0.5-0.3" & iso_info$mayor]<- -1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DSgroup == "0.5-0.3" &
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DSgroup == "0.5-0.3" & !iso_info$mayor &
    iso_info$gene_id == gen]*(1-0.5)/(numb_iso-1)
  varC1_sim<-iso_info$varC2[iso_info$DSgroup == "0.5-0.3" & !iso_info$mayor & 
    iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1))^2

  meanC2_sim<-iso_info$meanC2[iso_info$DSgroup == "0.5-0.3" & !iso_info$mayor &
    iso_info$gene_id == gen]*(1-0.3)/(numb_iso-1)
  varC2_sim<-iso_info$varC2[iso_info$DSgroup == "0.5-0.3" & !iso_info$mayor &
    iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1))^2
  
  index<-which(iso_info$DSgroup == "0.5-0.3" & !iso_info$mayor & 
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info$updown[upinfo$index]<-1

iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]


#0.3-0.7
table(gene_info$DSgroup=="0.3-0.7")
# FALSE  TRUE 
# 16124    60 

downgene<-as.character(iso_info$gene_id[iso_info$DSgroup == "0.3-0.7" & 
    iso_info$mayor])

iso_info$meanC1Sim[iso_info$DSgroup == "0.3-0.7" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DSgroup == "0.3-0.7" & iso_info$mayor]*0.3
iso_info$varC1Sim[iso_info$DSgroup == "0.3-0.7" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DSgroup == "0.3-0.7" & iso_info$mayor]*0.3^2
iso_info$meanC2Sim[iso_info$DSgroup == "0.3-0.7" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DSgroup == "0.3-0.7" & iso_info$mayor]*0.7
iso_info$varC2Sim[iso_info$DSgroup == "0.3-0.7" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DSgroup == "0.3-0.7" & iso_info$mayor]*0.7^2

iso_info$updown[iso_info$DSgroup == "0.3-0.7" & iso_info$mayor]<-1

# now should set the expression values to the rest!!
downinfo<-do.call(rbind,lapply(as.character(downgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DSgroup == "0.3-0.7" & 
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DSgroup == "0.3-0.7" & !iso_info$mayor &
    iso_info$gene_id == gen]*(1-0.3)/(numb_iso-1)
  varC1_sim<-iso_info$varC2[iso_info$DSgroup == "0.3-0.7" & !iso_info$mayor &
    iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1))^2

  meanC2_sim<-iso_info$meanC2[iso_info$DSgroup == "0.3-0.7" & !iso_info$mayor &
    iso_info$gene_id == gen]*(1-0.7)/(numb_iso-1)
  varC2_sim<-iso_info$varC2[iso_info$DSgroup == "0.3-0.7" & !iso_info$mayor &
    iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1))^2
  
  index<-which(iso_info$DSgroup == "0.3-0.7" & !iso_info$mayor & 
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info[downinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-downinfo[,-1]
iso_info$updown[downinfo$index]<-1

# 0.7-0.3

upgene<-as.character(iso_info$gene_id[iso_info$DSgroup == "0.7-0.3" & 
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DSgroup == "0.7-0.3" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DSgroup == "0.7-0.3" & iso_info$mayor]*0.7
iso_info$varC1Sim[iso_info$DSgroup == "0.7-0.3" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DSgroup == "0.7-0.3" & iso_info$mayor]*0.7^2

iso_info$meanC2Sim[iso_info$DSgroup == "0.7-0.3" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DSgroup == "0.7-0.3" & iso_info$mayor]*0.3
iso_info$varC2Sim[iso_info$DSgroup == "0.7-0.3" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DSgroup == "0.7-0.3" & iso_info$mayor]*0.3^2
iso_info$updown[iso_info$DSgroup == "0.7-0.3" & iso_info$mayor]<--1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DSgroup == "0.7-0.3" & 
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DSgroup == "0.7-0.3" & !iso_info$mayor &
    iso_info$gene_id == gen]*(1-0.7)/(numb_iso-1)
  varC1_sim<-iso_info$varC2[iso_info$DSgroup == "0.7-0.3" & !iso_info$mayor &
    iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1))^2

  meanC2_sim<-iso_info$meanC2[iso_info$DSgroup == "0.7-0.3" & !iso_info$mayor &
    iso_info$gene_id == gen]*(1-0.3)/(numb_iso-1)
  varC2_sim<-iso_info$varC2[iso_info$DSgroup == "0.7-0.3" & !iso_info$mayor &
    iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1))^2
  
  index<-which(iso_info$DSgroup == "0.7-0.3" & !iso_info$mayor & 
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info$updown[upinfo$index]<-1

iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]

## 2) Differential isoform expression and differential splicing

# 2-0.5-0.7
upgene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "2-0.5-0.7" & 
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "2-0.5-0.7" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "2-0.5-0.7" & iso_info$mayor]*0.5
iso_info$varC1Sim[iso_info$DIEDSgroup == "2-0.5-0.7" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "2-0.5-0.7" & iso_info$mayor]*0.5^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "2-0.5-0.7" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "2-0.5-0.7" & iso_info$mayor]*(2*0.7)
iso_info$varC2Sim[iso_info$DIEDSgroup == "2-0.5-0.7" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "2-0.5-0.7" & iso_info$mayor]*(2*0.7)^2
iso_info$updown[iso_info$DIEDSgroup == "2-0.5-0.7" & iso_info$mayor]<-1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "2-0.5-0.7" &
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "2-0.5-0.7" & 
    !iso_info$mayor & iso_info$gene_id == gen]*(1-0.5)/(numb_iso-1)
  varC1_sim<-iso_info$varC2[iso_info$DIEDSgroup == "2-0.5-0.7" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1))^2

  meanC2_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "2-0.5-0.7" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1)*2)
  varC2_sim<-iso_info$varC2[iso_info$DIEDSgroup == "2-0.5-0.7" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1)*2)^2
  
  index<-which(iso_info$DIEDSgroup == "2-0.5-0.7" & !iso_info$mayor & 
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))

iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]
iso_info$updown[upinfo$index]<-1

# 0.5-0.5-0.7
upgene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "0.5-0.5-0.7" &
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "0.5-0.5-0.7" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.5-0.7" & iso_info$mayor]*2*0.5
iso_info$varC1Sim[iso_info$DIEDSgroup == "0.5-0.5-0.7" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.5-0.7" & iso_info$mayor]*(2*0.5)^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "0.5-0.5-0.7" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.5-0.7" & iso_info$mayor]*(0.7)
iso_info$varC2Sim[iso_info$DIEDSgroup == "0.5-0.5-0.7" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.5-0.7" & iso_info$mayor]*(0.7)^2

iso_info$updown[iso_info$DIEDSgroup == "0.5-0.5-0.7" & iso_info$mayor]<--1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "0.5-0.5-0.7" & 
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.5-0.7" &
    !iso_info$mayor & iso_info$gene_id == gen]*(1-0.5)/(numb_iso-1)*2
  varC1_sim<-iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.5-0.7" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1)*2)^2
    meanC2_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.5-0.7" &
        !iso_info$mayor & iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1))
  varC2_sim<-iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.5-0.7" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1))^2
  
  index<-which(iso_info$DIEDSgroup == "0.5-0.5-0.7" & !iso_info$mayor &
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info$updown[upinfo$index]<--1
iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]

# 2-0.5-0.9
upgene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "2-0.5-0.9" &
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "2-0.5-0.9" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "2-0.5-0.9" & iso_info$mayor]*0.5
iso_info$varC1Sim[iso_info$DIEDSgroup == "2-0.5-0.9" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "2-0.5-0.9" & iso_info$mayor]*0.5^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "2-0.5-0.9" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "2-0.5-0.9" & iso_info$mayor]*(2*0.9)
iso_info$varC2Sim[iso_info$DIEDSgroup == "2-0.5-0.9" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "2-0.5-0.9" & iso_info$mayor]*(2*0.9)^2
iso_info$updown[iso_info$DIEDSgroup == "2-0.5-0.9" & iso_info$mayor]<-1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "2-0.5-0.9" & 
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "2-0.5-0.9" &
    !iso_info$mayor & iso_info$gene_id == gen]*(1-0.5)/(numb_iso-1)
  varC1_sim<-iso_info$varC2[iso_info$DIEDSgroup == "2-0.5-0.9" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1))^2

  meanC2_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "2-0.5-0.9" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.9)/(numb_iso-1)*2)
  varC2_sim<-iso_info$varC2[iso_info$DIEDSgroup == "2-0.5-0.9" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.9)/(numb_iso-1)*2)^2
  
  index<-which(iso_info$DIEDSgroup == "2-0.5-0.9" & !iso_info$mayor & 
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))

iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]
iso_info$updown[upinfo$index]<-1

# 0.5-0.5-0.9
upgene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "0.5-0.5-0.9" & 
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "0.5-0.5-0.9" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.5-0.9" & iso_info$mayor]*2*0.5
iso_info$varC1Sim[iso_info$DIEDSgroup == "0.5-0.5-0.9" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.5-0.9" & iso_info$mayor]*(2*0.5)^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "0.5-0.5-0.9" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.5-0.9" & iso_info$mayor]*(0.9)
iso_info$varC2Sim[iso_info$DIEDSgroup == "0.5-0.5-0.9" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.5-0.9" & iso_info$mayor]*(0.9)^2

iso_info$updown[iso_info$DIEDSgroup == "0.5-0.5-0.9" & iso_info$mayor]<--1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "0.5-0.5-0.9" &
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.5-0.9" & 
    !iso_info$mayor & iso_info$gene_id == gen]*(1-0.5)/(numb_iso-1)*2
  varC1_sim<-iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.5-0.9" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1)*2)^2
    meanC2_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.5-0.9" & 
        !iso_info$mayor & iso_info$gene_id == gen]*((1-0.9)/(numb_iso-1))
  varC2_sim<-iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.5-0.9" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.9)/(numb_iso-1))^2
  
  index<-which(iso_info$DIEDSgroup == "0.5-0.5-0.9" & !iso_info$mayor & 
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info$updown[upinfo$index]<--1
iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]

# 2-0.3-0.5
upgene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "2-0.3-0.5" &
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "2-0.3-0.5" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "2-0.3-0.5" & iso_info$mayor]*0.3
iso_info$varC1Sim[iso_info$DIEDSgroup == "2-0.3-0.5" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "2-0.3-0.5" & iso_info$mayor]*0.3^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "2-0.3-0.5" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "2-0.3-0.5" & iso_info$mayor]*(2*0.5)
iso_info$varC2Sim[iso_info$DIEDSgroup == "2-0.3-0.5" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "2-0.3-0.5" & iso_info$mayor]*(2*0.5)^2
iso_info$updown[iso_info$DIEDSgroup == "2-0.3-0.5" & iso_info$mayor]<-1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "2-0.3-0.5" &
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "2-0.3-0.5" & 
    !iso_info$mayor & iso_info$gene_id == gen]*(1-0.3)/(numb_iso-1)
  varC1_sim<-iso_info$varC2[iso_info$DIEDSgroup == "2-0.3-0.5" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1))^2

  meanC2_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "2-0.3-0.5" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1)*2)
  varC2_sim<-iso_info$varC2[iso_info$DIEDSgroup == "2-0.3-0.5" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1)*2)^2
  
  index<-which(iso_info$DIEDSgroup == "2-0.3-0.5" & !iso_info$mayor &
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))

iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]
iso_info$updown[upinfo$index]<-1

# 0.5-0.3-0.5
upgene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "0.5-0.3-0.5" &
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "0.5-0.3-0.5" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.3-0.5" & iso_info$mayor]*2*0.3
iso_info$varC1Sim[iso_info$DIEDSgroup == "0.5-0.3-0.5" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.3-0.5" & iso_info$mayor]*(2*0.3)^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "0.5-0.3-0.5" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.3-0.5" & iso_info$mayor]*(0.5)
iso_info$varC2Sim[iso_info$DIEDSgroup == "0.5-0.3-0.5" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.3-0.5" & iso_info$mayor]*(0.5)^2

iso_info$updown[iso_info$DIEDSgroup == "0.5-0.3-0.5" & iso_info$mayor]<--1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "0.5-0.3-0.5" &
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.3-0.5" &
    !iso_info$mayor & iso_info$gene_id == gen]*(1-0.3)/(numb_iso-1)*2
  varC1_sim<-iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.3-0.5" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1)*2)^2
    meanC2_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.3-0.5" &
        !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1))
  varC2_sim<-iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.3-0.5" &
  !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1))^2
  
  index<-which(iso_info$DIEDSgroup == "0.5-0.3-0.5" & !iso_info$mayor &
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info$updown[upinfo$index]<--1
iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]

# 2-0.3-0.7
upgene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "2-0.3-0.7" &
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "2-0.3-0.7" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "2-0.3-0.7" & iso_info$mayor]*0.3
iso_info$varC1Sim[iso_info$DIEDSgroup == "2-0.3-0.7" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "2-0.3-0.7" & iso_info$mayor]*0.3^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "2-0.3-0.7" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "2-0.3-0.7" & iso_info$mayor]*(2*0.7)
iso_info$varC2Sim[iso_info$DIEDSgroup == "2-0.3-0.7" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "2-0.3-0.7" & iso_info$mayor]*(2*0.7)^2
iso_info$updown[iso_info$DIEDSgroup == "2-0.3-0.7" & iso_info$mayor]<-1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "2-0.3-0.7" &
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "2-0.3-0.7" &
    !iso_info$mayor & iso_info$gene_id == gen]*(1-0.3)/(numb_iso-1)
  varC1_sim<-iso_info$varC2[iso_info$DIEDSgroup == "2-0.3-0.7" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1))^2

  meanC2_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "2-0.3-0.7" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1)*2)
  varC2_sim<-iso_info$varC2[iso_info$DIEDSgroup == "2-0.3-0.7" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1)*2)^2
  
  index<-which(iso_info$DIEDSgroup == "2-0.3-0.7" & !iso_info$mayor & 
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
  varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))

iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]
iso_info$updown[upinfo$index]<-1

# 0.5-0.3-0.7
upgene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "0.5-0.3-0.7" &
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "0.5-0.3-0.7" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.3-0.7" & iso_info$mayor]*2*0.3
iso_info$varC1Sim[iso_info$DIEDSgroup == "0.5-0.3-0.7" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.3-0.7" & iso_info$mayor]*(2*0.3)^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "0.5-0.3-0.7" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.3-0.7" & iso_info$mayor]*(0.7)
iso_info$varC2Sim[iso_info$DIEDSgroup == "0.5-0.3-0.7" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.3-0.7" & iso_info$mayor]*(0.7)^2

iso_info$updown[iso_info$DIEDSgroup == "0.5-0.3-0.7" & iso_info$mayor]<--1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "0.5-0.3-0.7" & 
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.3-0.7" & 
    !iso_info$mayor & iso_info$gene_id == gen]*(1-0.3)/(numb_iso-1)*2
  varC1_sim<-iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.3-0.7" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1)*2)^2
    meanC2_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.3-0.7" & 
        !iso_info$mayor & iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1))
  varC2_sim<-iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.3-0.7" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1))^2
  
  index<-which(iso_info$DIEDSgroup == "0.5-0.3-0.7" & !iso_info$mayor & 
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info$updown[upinfo$index]<--1
iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]

# 2-0.5-0.3
upgene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "2-0.5-0.3" & 
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "2-0.5-0.3" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "2-0.5-0.3" & iso_info$mayor]*0.5
iso_info$varC1Sim[iso_info$DIEDSgroup == "2-0.5-0.3" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "2-0.5-0.3" & iso_info$mayor]*0.5^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "2-0.5-0.3" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "2-0.5-0.3" & iso_info$mayor]*(2*0.3)
iso_info$varC2Sim[iso_info$DIEDSgroup == "2-0.5-0.3" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "2-0.5-0.3" & iso_info$mayor]*(2*0.3)^2
iso_info$updown[iso_info$DIEDSgroup == "2-0.5-0.3" & iso_info$mayor]<-1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "2-0.5-0.3" & 
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "2-0.5-0.3" &
    !iso_info$mayor & iso_info$gene_id == gen]*(1-0.5)/(numb_iso-1)
  varC1_sim<-iso_info$varC2[iso_info$DIEDSgroup == "2-0.5-0.3" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1))^2

  meanC2_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "2-0.5-0.3" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1)*2)
  varC2_sim<-iso_info$varC2[iso_info$DIEDSgroup == "2-0.5-0.3" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1)*2)^2
  
  index<-which(iso_info$DIEDSgroup == "2-0.5-0.3" & !iso_info$mayor & iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))

iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]
iso_info$updown[upinfo$index]<-1

# 0.5-0.5-0.3
upgene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "0.5-0.5-0.3" &
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "0.5-0.5-0.3" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.5-0.3" & iso_info$mayor]*2*0.5
iso_info$varC1Sim[iso_info$DIEDSgroup == "0.5-0.5-0.3" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.5-0.3" & iso_info$mayor]*(2*0.5)^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "0.5-0.5-0.3" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.5-0.3" & iso_info$mayor]*(0.3)
iso_info$varC2Sim[iso_info$DIEDSgroup == "0.5-0.5-0.3" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.5-0.3" & iso_info$mayor]*(0.3)^2

iso_info$updown[iso_info$DIEDSgroup == "0.5-0.5-0.3" & iso_info$mayor]<--1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "0.5-0.5-0.3" & 
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.5-0.3" & 
    !iso_info$mayor & iso_info$gene_id == gen]*(1-0.5)/(numb_iso-1)*2
  varC1_sim<-iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.5-0.3" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1)*2)^2
    meanC2_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.5-0.3" & 
        !iso_info$mayor & iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1))
  varC2_sim<-iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.5-0.3" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1))^2
  
  index<-which(iso_info$DIEDSgroup == "0.5-0.5-0.3" & !iso_info$mayor &
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info$updown[upinfo$index]<--1
iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]

# 2-0.7-0.3
upgene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "2-0.7-0.3" & 
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "2-0.7-0.3" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "2-0.7-0.3" & iso_info$mayor]*0.7
iso_info$varC1Sim[iso_info$DIEDSgroup == "2-0.7-0.3" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "2-0.7-0.3" & iso_info$mayor]*0.7^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "2-0.7-0.3" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "2-0.7-0.3" & iso_info$mayor]*(2*0.3)
iso_info$varC2Sim[iso_info$DIEDSgroup == "2-0.7-0.3" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "2-0.7-0.3" & iso_info$mayor]*(2*0.3)^2
iso_info$updown[iso_info$DIEDSgroup == "2-0.7-0.3" & iso_info$mayor]<-1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "2-0.7-0.3" &
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "2-0.7-0.3" &
    !iso_info$mayor & iso_info$gene_id == gen]*(1-0.7)/(numb_iso-1)
  varC1_sim<-iso_info$varC2[iso_info$DIEDSgroup == "2-0.7-0.3" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1))^2

  meanC2_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "2-0.7-0.3" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1)*2)
  varC2_sim<-iso_info$varC2[iso_info$DIEDSgroup == "2-0.7-0.3" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1)*2)^2
  
  index<-which(iso_info$DIEDSgroup == "2-0.7-0.3" & !iso_info$mayor &
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))

iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]
iso_info$updown[upinfo$index]<-1

# 0.5-0.7-0.3
upgene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "0.5-0.7-0.3" &
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "0.5-0.7-0.3" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.7-0.3" & iso_info$mayor]*2*0.7
iso_info$varC1Sim[iso_info$DIEDSgroup == "0.5-0.7-0.3" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.7-0.3" & iso_info$mayor]*(2*0.7)^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "0.5-0.7-0.3" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.7-0.3" & iso_info$mayor]*(0.3)
iso_info$varC2Sim[iso_info$DIEDSgroup == "0.5-0.7-0.3" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.7-0.3" & iso_info$mayor]*(0.3)^2

iso_info$updown[iso_info$DIEDSgroup == "0.5-0.7-0.3" & iso_info$mayor]<--1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "0.5-0.7-0.3" &
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.7-0.3" &
    !iso_info$mayor & iso_info$gene_id == gen]*(1-0.7)/(numb_iso-1)*2
  varC1_sim<-iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.7-0.3" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1)*2)^2
    meanC2_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "0.5-0.7-0.3" & 
        !iso_info$mayor & iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1))
  varC2_sim<-iso_info$varC2[iso_info$DIEDSgroup == "0.5-0.7-0.3" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1))^2
  
  index<-which(iso_info$DIEDSgroup == "0.5-0.7-0.3" & !iso_info$mayor &
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info$updown[upinfo$index]<--1
iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]

# 4-0.5-0.7
upgene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "4-0.5-0.7" &
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "4-0.5-0.7" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "4-0.5-0.7" & iso_info$mayor]*0.5
iso_info$varC1Sim[iso_info$DIEDSgroup == "4-0.5-0.7" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "4-0.5-0.7" & iso_info$mayor]*0.5^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "4-0.5-0.7" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "4-0.5-0.7" & iso_info$mayor]*(4*0.7)
iso_info$varC2Sim[iso_info$DIEDSgroup == "4-0.5-0.7" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "4-0.5-0.7" & iso_info$mayor]*(4*0.7)^2
iso_info$updown[iso_info$DIEDSgroup == "4-0.5-0.7" & iso_info$mayor]<-1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "4-0.5-0.7" &
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "4-0.5-0.7" &
    !iso_info$mayor & iso_info$gene_id == gen]*(1-0.5)/(numb_iso-1)
  varC1_sim<-iso_info$varC2[iso_info$DIEDSgroup == "4-0.5-0.7" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1))^2

  meanC2_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "4-0.5-0.7" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1)*4)
  varC2_sim<-iso_info$varC2[iso_info$DIEDSgroup == "4-0.5-0.7" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1)*4)^2
  
  index<-which(iso_info$DIEDSgroup == "4-0.5-0.7" & !iso_info$mayor &
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))

iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]
iso_info$updown[upinfo$index]<-1

# 0.25-0.5-0.7
upgene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "0.25-0.5-0.7" &
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "0.25-0.5-0.7" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "0.25-0.5-0.7" & iso_info$mayor]*4*0.5
iso_info$varC1Sim[iso_info$DIEDSgroup == "0.25-0.5-0.7" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "0.25-0.5-0.7" & iso_info$mayor]*(4*0.5)^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "0.25-0.5-0.7" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "0.25-0.5-0.7" & iso_info$mayor]*(0.7)
iso_info$varC2Sim[iso_info$DIEDSgroup == "0.25-0.5-0.7" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "0.25-0.5-0.7" & iso_info$mayor]*(0.7)^2

iso_info$updown[iso_info$DIEDSgroup == "0.25-0.5-0.7" & iso_info$mayor]<--1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "0.25-0.5-0.7" &
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "0.25-0.5-0.7" & 
    !iso_info$mayor & iso_info$gene_id == gen]*(1-0.5)/(numb_iso-1)*4
  varC1_sim<-iso_info$varC2[iso_info$DIEDSgroup == "0.25-0.5-0.7" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1)*4)^2
    meanC2_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "0.25-0.5-0.7" &
        !iso_info$mayor & iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1))
  varC2_sim<-iso_info$varC2[iso_info$DIEDSgroup == "0.25-0.5-0.7" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1))^2
  
  index<-which(iso_info$DIEDSgroup == "0.25-0.5-0.7" & !iso_info$mayor & 
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info$updown[upinfo$index]<--1
iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]

# 4-0.5-0.9
upgene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "4-0.5-0.9" &
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "4-0.5-0.9" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "4-0.5-0.9" & iso_info$mayor]*0.5
iso_info$varC1Sim[iso_info$DIEDSgroup == "4-0.5-0.9" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "4-0.5-0.9" & iso_info$mayor]*0.5^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "4-0.5-0.9" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "4-0.5-0.9" & iso_info$mayor]*(4*0.9)
iso_info$varC2Sim[iso_info$DIEDSgroup == "4-0.5-0.9" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "4-0.5-0.9" & iso_info$mayor]*(4*0.9)^2
iso_info$updown[iso_info$DIEDSgroup == "4-0.5-0.9" & iso_info$mayor]<-1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "4-0.5-0.9" &
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "4-0.5-0.9" &
    !iso_info$mayor & iso_info$gene_id == gen]*(1-0.5)/(numb_iso-1)
  varC1_sim<-iso_info$varC2[iso_info$DIEDSgroup == "4-0.5-0.9" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1))^2

  meanC2_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "4-0.5-0.9" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.9)/(numb_iso-1)*4)
  varC2_sim<-iso_info$varC2[iso_info$DIEDSgroup == "4-0.5-0.9" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.9)/(numb_iso-1)*4)^2
  
  index<-which(iso_info$DIEDSgroup == "4-0.5-0.9" & !iso_info$mayor &
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))

iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]
iso_info$updown[upinfo$index]<-1

# 0.25-0.5-0.9
upgene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "0.25-0.5-0.9" &
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "0.25-0.5-0.9" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "0.25-0.5-0.9" & iso_info$mayor]*4*0.5
iso_info$varC1Sim[iso_info$DIEDSgroup == "0.25-0.5-0.9" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "0.25-0.5-0.9" & iso_info$mayor]*(4*0.5)^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "0.25-0.5-0.9" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "0.25-0.5-0.9" & iso_info$mayor]*(0.9)
iso_info$varC2Sim[iso_info$DIEDSgroup == "0.25-0.5-0.9" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "0.25-0.5-0.9" & iso_info$mayor]*(0.9)^2

iso_info$updown[iso_info$DIEDSgroup == "0.25-0.5-0.9" & iso_info$mayor]<--1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "0.25-0.5-0.9" &
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "0.25-0.5-0.9" &
    !iso_info$mayor & iso_info$gene_id == gen]*(1-0.5)/(numb_iso-1)*4
  varC1_sim<-iso_info$varC2[iso_info$DIEDSgroup == "0.25-0.5-0.9" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1)*4)^2
    meanC2_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "0.25-0.5-0.9" &
        !iso_info$mayor & iso_info$gene_id == gen]*((1-0.9)/(numb_iso-1))
  varC2_sim<-iso_info$varC2[iso_info$DIEDSgroup == "0.25-0.5-0.9" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.9)/(numb_iso-1))^2
  
  index<-which(iso_info$DIEDSgroup == "0.25-0.5-0.9" & !iso_info$mayor &
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info$updown[upinfo$index]<--1
iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]

# 4-0.7-0.5
upgene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "4-0.7-0.5" &
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "4-0.7-0.5" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "4-0.7-0.5" & iso_info$mayor]*0.7
iso_info$varC1Sim[iso_info$DIEDSgroup == "4-0.7-0.5" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "4-0.7-0.5" & iso_info$mayor]*0.7^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "4-0.7-0.5" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "4-0.7-0.5" & iso_info$mayor]*(4*0.5)
iso_info$varC2Sim[iso_info$DIEDSgroup == "4-0.7-0.5" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "4-0.7-0.5" & iso_info$mayor]*(4*0.5)^2
iso_info$updown[iso_info$DIEDSgroup == "4-0.7-0.5" & iso_info$mayor]<-1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "4-0.7-0.5" &
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "4-0.7-0.5" &
    !iso_info$mayor & iso_info$gene_id == gen]*(1-0.7)/(numb_iso-1)
  varC1_sim<-iso_info$varC2[iso_info$DIEDSgroup == "4-0.7-0.5" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1))^2

  meanC2_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "4-0.7-0.5" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1)*4)
  varC2_sim<-iso_info$varC2[iso_info$DIEDSgroup == "4-0.7-0.5" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1)*4)^2
  
  index<-which(iso_info$DIEDSgroup == "4-0.7-0.5" & !iso_info$mayor &
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))

iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]
iso_info$updown[upinfo$index]<-1

# 0.25-0.7-0.5
upgene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "0.25-0.7-0.5" &
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "0.25-0.7-0.5" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "0.25-0.7-0.5" & iso_info$mayor]*4*0.7
iso_info$varC1Sim[iso_info$DIEDSgroup == "0.25-0.7-0.5" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "0.25-0.7-0.5" & iso_info$mayor]*(4*0.7)^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "0.25-0.7-0.5" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "0.25-0.7-0.5" & iso_info$mayor]*(0.5)
iso_info$varC2Sim[iso_info$DIEDSgroup == "0.25-0.7-0.5" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "0.25-0.7-0.5" & iso_info$mayor]*(0.5)^2

iso_info$updown[iso_info$DIEDSgroup == "0.25-0.7-0.5" & iso_info$mayor]<--1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "0.25-0.7-0.5" &
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "0.25-0.7-0.5" & 
    !iso_info$mayor & iso_info$gene_id == gen]*(1-0.7)/(numb_iso-1)*4
  varC1_sim<-iso_info$varC2[iso_info$DIEDSgroup == "0.25-0.7-0.5" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1)*4)^2
    meanC2_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "0.25-0.7-0.5" &
        !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1))
  varC2_sim<-iso_info$varC2[iso_info$DIEDSgroup == "0.25-0.7-0.5" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1))^2
  
  index<-which(iso_info$DIEDSgroup == "0.25-0.7-0.5" & !iso_info$mayor &
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info$updown[upinfo$index]<--1
iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]

# 4-0.9-0.5
upgene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "4-0.9-0.5" &
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "4-0.9-0.5" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "4-0.9-0.5" & iso_info$mayor]*0.9
iso_info$varC1Sim[iso_info$DIEDSgroup == "4-0.9-0.5" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "4-0.9-0.5" & iso_info$mayor]*0.9^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "4-0.9-0.5" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "4-0.9-0.5" & iso_info$mayor]*(4*0.5)
iso_info$varC2Sim[iso_info$DIEDSgroup == "4-0.9-0.5" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "4-0.9-0.5" & iso_info$mayor]*(4*0.5)^2
iso_info$updown[iso_info$DIEDSgroup == "4-0.9-0.5" & iso_info$mayor]<-1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "4-0.9-0.5" & 
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "4-0.9-0.5" & 
    !iso_info$mayor & iso_info$gene_id == gen]*(1-0.9)/(numb_iso-1)
  varC1_sim<-iso_info$varC2[iso_info$DIEDSgroup == "4-0.9-0.5" & 
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.9)/(numb_iso-1))^2

  meanC2_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "4-0.9-0.5" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1)*4)
  varC2_sim<-iso_info$varC2[iso_info$DIEDSgroup == "4-0.9-0.5" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1)*4)^2
  
  index<-which(iso_info$DIEDSgroup == "4-0.9-0.5" & !iso_info$mayor & 
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
    varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))

iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]
iso_info$updown[upinfo$index]<-1

# 0.25-0.9-0.5
upgene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "0.25-0.9-0.5" &
    iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "0.25-0.9-0.5" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "0.25-0.9-0.5" & iso_info$mayor]*4*0.9
iso_info$varC1Sim[iso_info$DIEDSgroup == "0.25-0.9-0.5" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "0.25-0.9-0.5" & iso_info$mayor]*(4*0.9)^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "0.25-0.9-0.5" & iso_info$mayor]<-
    iso_info$meanC2[iso_info$DIEDSgroup == "0.25-0.9-0.5" & iso_info$mayor]*(0.5)
iso_info$varC2Sim[iso_info$DIEDSgroup == "0.25-0.9-0.5" & iso_info$mayor]<-
    iso_info$varC2[iso_info$DIEDSgroup == "0.25-0.9-0.5" & iso_info$mayor]*(0.5)^2

iso_info$updown[iso_info$DIEDSgroup == "0.25-0.9-0.5" & iso_info$mayor]<--1

upinfo<-do.call(rbind,lapply(as.character(upgene), function(gen){
  numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "0.25-0.9-0.5" &
    !iso_info$mayor & iso_info$gene_id == gen])
  meanC1_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "0.25-0.9-0.5" &
    !iso_info$mayor & iso_info$gene_id == gen]*(1-0.9)/(numb_iso-1)*4
  varC1_sim<-iso_info$varC2[iso_info$DIEDSgroup == "0.25-0.9-0.5" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.9)/(numb_iso-1)*4)^2
    meanC2_sim<-iso_info$meanC2[iso_info$DIEDSgroup == "0.25-0.9-0.5" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1))
  varC2_sim<-iso_info$varC2[iso_info$DIEDSgroup == "0.25-0.9-0.5" &
    !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1))^2
  
  index<-which(iso_info$DIEDSgroup == "0.25-0.9-0.5" & !iso_info$mayor &
    iso_info$gene_id == gen)
  return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim,
  varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info$updown[upinfo$index]<--1
iso_info[upinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-upinfo[,-1]

## Generating replicates values
# NB parameters
meansC1<-iso_info$meanC1Sim
meansC2<-iso_info$meanC2Sim
f.loess.C1<- loess(varC1 ~ mediasC1, data.frame(mediasC1=meansC1, varC1=iso_info$varC1Sim), degree=2)
f.loess.C2<- loess(varC2 ~ mediasC2, data.frame(mediasC2=meansC2, varC2=iso_info$varC2Sim), degree=2)
# variance prediction
var.predict1 <- predict(f.loess.C1, data.frame(mediasC1=meansC1))
var.predict2 <- predict(f.loess.C2, data.frame(mediasC2=meansC2))
size<-meansC1^2/(var.predict1-meansC1)
prob<-size/(size+meansC1)

numberExperimentRealizations<-10
for(nsim in 1:numberExperimentRealizations){
setwd(paste("/path_to_sim", nsim))
prob[size<0]<-size[size<0]<-1
C1repeticiones10<-t(sapply(1:nrow(iso_info), function(iso){
    iso_cts_C1<-t(sapply(1:1, function(x){
      cond1<-rnbinom(n=4, size=size[iso], prob=prob[iso])
      return(c(cond1))

    }))
  
    cero_index_cts<-which(is.na(iso_cts_C1[,1]) | is.na(iso_cts_C1[,2]) |   
        is.na(iso_cts_C1[,3]) | is.na(iso_cts_C1[,4]))
    iso_cts_C1[cero_index_cts,]<-0
    return(colMeans(iso_cts_C1))
  }))
iso_cts_C1<-as.data.frame(C1repeticiones10)
rownames(iso_cts_C1)<-iso_info$transcript_id
names(iso_cts_C1)<-c("C1R1", "C1R2", "C1R3", "C1R4")

# Condition2
size<-meansC2^2/(var.predict2-meansC2)
prob<-size/(size+meansC2)
prob[size<0]<-size[size<0]<-1
C2repeticiones10<-t(sapply(1:nrow(iso_info), function(iso){
    iso_cts_C2<-t(sapply(1:1, function(x){
      cond2<-rnbinom(n=4, size=size[iso], prob=prob[iso])
      return(c(cond2))

    }))
  
    cero_index_cts<-which(is.na(iso_cts_C2[,1]) | is.na(iso_cts_C2[,2]) | 
        is.na(iso_cts_C2[,3]) | is.na(iso_cts_C2[,4]))
    iso_cts_C2[cero_index_cts,]<-0
    return(colMeans(iso_cts_C2))
  }))
iso_cts_C2<-as.data.frame(C2repeticiones10)
rownames(iso_cts_C2)<-iso_info$transcript_id
names(iso_cts_C2)<-c("C2R1", "C2R2", "C2R3", "C2R4")


iso_cm_sim<-data.frame(iso_cts_C1, iso_cts_C2)
names(iso_cm_sim)<-c("C1R1", "C1R2", "C1R3", "C1R4", "C2R1", "C2R2", "C2R3", "C2R4")
rownames(iso_cm_sim)<-rownames(iso_cm)

## Generate files containing simulation profiles 
transcript_ids<-as.character(iso_res[[1]]$transcript_id)
index<-match(transcript_ids, rownames(iso_cm_sim), nomatch=0)
iso_cm_sim_complete<-data.frame(matrix(0,nrow=length(transcript_ids), 
    ncol=ncol(iso_cm_sim)))
iso_cm_sim_complete[index!= 0,]<-iso_cm_sim[index[index !=0],]
rownames(iso_cm_sim_complete)<-transcript_ids
names(iso_cm_sim_complete)<-names(iso_cm_sim)
genes<-as.list(unique(as.character(iso_info$gene_id)))
sampleIndex<-1:8
effect_length<-do.call(cbind,bplapply(sampleIndex, function(sample){
  
  return(iso_res[[sample]][,"effective_length"])
  
  }, BPPARAM=MulticoreParam(16)))
  
effect_length<-as.data.frame(effect_length)
rownames(effect_length)<-rownames(iso_cm_sim_complete)
t_length<-do.call(cbind,bplapply(sampleIndex, function(sample){
  
  return(iso_res[[sample]][,"length"])
  
  }, BPPARAM=MulticoreParam(18)))
  
t_length<-as.data.frame(t_length)
rownames(t_length)<-rownames(iso_cm_sim_complete)
scale_factors<-sapply(1:ncol(iso_cm_sim_complete), function(i){
    
  allratio<-sum(iso_cm_sim_complete[effect_length[,i] > 0,i]/effect_length[
        effect_length[,i] > 0,i])

  return(allratio)
})


TPM<-do.call(rbind,bplapply(as.list(transcript_ids), function(iso){
    
    tpm<-iso_cm_sim_complete[rownames(iso_cm_sim_complete) == iso,]/
        effect_length[rownames(effect_length) == iso,]
    if(any(is.na(tpm)) | any(tpm == Inf)){
      tpm[is.na(tpm) | tpm == Inf]<-0
      }
    tpm<-tpm*10^6/scale_factors
    return(tpm)

}, BPPARAM=MulticoreParam(18)))
sample_sizes<-colSums(iso_cm_sim)
FPKM<-sapply(1:ncol(iso_cm_sim_complete), function(sample){
  FPKM<-(iso_cm_sim_complete[,sample]) / (effect_length[,sample])*10^9/
    sample_sizes[sample]
  FPKM[is.na(FPKM) | FPKM == Inf]<-0
  return(FPKM)
  })
FPKM<-as.data.frame(FPKM)
names(FPKM)<-names(TPM)
rownames(FPKM)<-rownames(TPM)
index<-match(transcript_ids,g2t$transcript_id, nomatch=0)
g2t_ord<-g2t[index,]
# writing isoform files
sim_data<-bplapply(1:ncol(iso_cm_sim_complete), function(sample){
  sample_data<-data.frame(transcript_id=transcript_ids, gene_id=g2t_ord$gene_id,
  length=as.character(round(t_length[,sample],2)),effective_length= 
  as.character(round(effect_length[,sample],2)), expected_count=as.character(
  round(iso_cm_sim_complete[,sample],2)), TPM=as.character(round(TPM[,sample],
  2)),FPKM=as.character(round(FPKM[,sample],2)),IsoPct=as.character(round(
  ratios_complete[,sample]*100,2)))
  return(sample_data)
}, BPPARAM=MulticoreParam(10))

names(sim_data)<-names(iso_res)[sampleIndex]

lapply(1:length(sim_data), function(sim){
  
    write.table(sim_data[[sim]], file=paste(names(sim_data)[sim], 
        "_sim_iso.results", sep=""), row.names=FALSE, col.names=TRUE,dec=".", 
        quote=FALSE, sep="\t")
  })
}  
write.table(gene_info, "gene_info.tab", sep="\t", row.names=FALSE, 
    col.names=TRUE, quote=FALSE)  
write.table(iso_info, "iso_info.tab", sep="\t", row.names=FALSE,
    col.names=TRUE, quote=FALSE)  
save.image("simulation.RData", compress="xz")
