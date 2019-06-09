###############################################
## NBSplice performance at the isoform level ##
###############################################
rm(list=ls())
library("ggplot2")
library("gridExtra")
library("BiocParallel")
library("tximport")
library("NBSplice")
library("iCOBRA")
library("parallel")
# setting the amount of simulated replications
n<-10
setwd("path_to_simulated_samples")


mclapply(1:n, function(n.sim){
    setwd(paste("sim", n.sim,"/", sep=""))
    # loading simulated data information and matrices
    load(paste("iso_info_final_sim", n.sim,".RData", sep=""))
    load(paste("iso_cm_sim", n.sim,".RData", sep=""))
    load( file=paste("gene_iso_sim", n.sim,".RData", sep=""))
    DSgenes<-unique(iso_info[ iso_info$DS  | iso_info$DIEDS , "gene_id"])
    # identifiyng the genes simulated as differentially spliced
    DSgenes<-unique(iso_info[ iso_info$DS  | iso_info$DIEDS , "gene_id"])
    id_mayorDS<-as.character(g2t_ok$transcript_id[g2t_ok$gene_id %in% DSgenes & iso_info$mayor])
    # building the desing Matrix
    sampleData<-data.frame(samples=colnames(iso_cm), condition=factor(c(rep("Normal",4), rep("Tumor",4)), levels=c("Normal", "Tumor")))
    geneIso<-iso_info[,c("gene_id", "transcript_id")]
    colnames(geneIso)[2]<-"isoform_id"
    rownames(sampleData)<-sampleData$samples
  
  
    ## Loading NBSplice results
    load(file=paste("NBSpliceResSim", n.sim, ".RData", sep=""))
    DSGeneRes<-unique(results(myDSRes, filter = FALSE)[,c("iso", "FDR")])
    
    load(file=paste("NBSpliceResWFSim", n.sim, ".RData", sep=""))
    DSGeneResWF<-unique(results(myDSResWF, filter = FALSE)[,c("iso", "FDR")])

    load(file=paste("NBSpliceResWFRatSim", n.sim, ".RData", sep=""))
    DSGeneResWFR<-unique(results(myDSResWFR, filter = FALSE)[,c("iso", "FDR")])
    
    load(file=paste("NBSpliceResWFCountSim", n.sim, ".RData", sep=""))
    DSGeneResWFC<-unique(results(myDSResWFC, filter = FALSE)[,c("iso", "FDR")])
    
    load(file=paste("NBSpliceRes2Sim", n.sim, ".RData", sep=""))
    DSGeneRes2<-unique(results(myDSRes2, filter = FALSE)[,c("iso", "FDR")])
    
    load(file=paste("NBSpliceRes3Sim", n.sim, ".RData", sep=""))
    DSGeneRes3<-unique(results(myDSRes3, filter = FALSE)[,c("iso", "FDR")])
    
    load(file=paste("NBSpliceRes4Sim", n.sim, ".RData", sep=""))
    DSGeneRes4<-unique(results(myDSRes4, filter = FALSE)[,c("iso", "FDR")])

    load(file=paste("NBSpliceRes5Sim", n.sim, ".RData", sep=""))
    DSGeneRes5<-unique(results(myDSRes5, filter = FALSE)[,c("iso", "FDR")])

    load(file=paste("NBSpliceRes6Sim", n.sim, ".RData", sep=""))
    DSGeneRes6<-unique(results(myDSRes6, filter = FALSE)[,c("iso", "FDR")])
    
    # obtaining the adjusted p-values at the isoform level indicating occurrence of differential transcript usage
    padj.gene <- data.frame(row.names=unique(iso_info$transcript_id[which(!is.na(iso_info$transcript_id))]))
    padj.gene$NBSplice <- DSGeneRes$FDR[match(rownames(padj.gene), DSGeneRes$iso)]
    padj.gene$NBSpliceWF <- DSGeneResWF$FDR[match(rownames(padj.gene), DSGeneResWF$iso)]
    padj.gene$NBSpliceWFR <- DSGeneResWFR$FDR[match(rownames(padj.gene), DSGeneResWFR$iso)]
    padj.gene$NBSpliceWFC <- DSGeneResWFC$FDR[match(rownames(padj.gene), DSGeneResWFC$iso)]
    padj.gene$NBSplice2 <- DSGeneRes2$FDR[match(rownames(padj.gene), DSGeneRes2$iso)]
    padj.gene$NBSplice3 <- DSGeneRes3$FDR[match(rownames(padj.gene), DSGeneRes3$iso)]
    padj.gene$NBSplice4 <- DSGeneRes4$FDR[match(rownames(padj.gene), DSGeneRes4$iso)]
    padj.gene$NBSplice5 <- DSGeneRes5$FDR[match(rownames(padj.gene), DSGeneRes5$iso)]
    padj.gene$NBSplice6 <- DSGeneRes6$FDR[match(rownames(padj.gene), DSGeneRes6$iso)]
    
    # obtaining the true state of each transcript
    truth<-data.frame(status=as.numeric(rownames(padj.gene) %in% DSIso),
                      row.names=rownames(padj.gene), geneGroup=iso_info$group[match(rownames(padj.gene), iso_info$transcript_id)])
        
    nas <- apply(padj.gene, 1, function(x) all(is.na(x)))
    padj.gene <- padj.gene[!nas,]
    truth <- truth[!nas,,drop=FALSE]
    # saving differential transcript usage adjusted p-values
    DSResults<-cbind(padj.gene, truth$status)
    save(DSResults, file=paste("Qvals_transcript_sim", n.sim, ".RData",sep=""), compress = "xz")
    cd.gene <- COBRAData(padj=padj.gene, truth=truth)
    cp <- calculate_performance(cd.gene,
                                binary_truth="status",
                                aspect=c("fdrtpr","fdrtprcurve"),
                                thrs=c(.01,.05,.1))
    dfIso<-cp@fdrtpr
    save(dfIso, file=paste("CobraPerfIsoSim", n.sim, ".RData", sep=""), compress = "xz")
}, mc.cores=16)

