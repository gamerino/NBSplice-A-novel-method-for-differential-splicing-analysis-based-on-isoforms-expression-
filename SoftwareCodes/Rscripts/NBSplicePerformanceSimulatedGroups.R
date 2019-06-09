######################################################################
## Computing performance measures on simualted groups and subgroups ##
######################################################################

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
    # identifiying simulated groups and subgroups
    iso_info$simGroup<-""
    iso_info$simGroup[iso_info$DS]<-"DS"
    iso_info$simGroup[iso_info$DIEDS]<-"DIEDS"
    iso_info$simSubGroup<-""
    iso_info$simSubGroup[iso_info$DIEDS]<-iso_info$DIEDSgroup[iso_info$DIEDS]
    iso_info$simSubGroup[iso_info$DS]<-iso_info$DSgroup[iso_info$DS]
    
    # building the desing Matrix
    sampleData<-data.frame(samples=colnames(iso_cm), condition=factor(c(rep("Normal",4), rep("Tumor",4)), levels=c("Normal", "Tumor")))
    geneIso<-iso_info[,c("gene_id", "transcript_id")]
    colnames(geneIso)[2]<-"isoform_id"
    rownames(sampleData)<-sampleData$samples

    ## loading the NBSplice results with its nine configurations
    load(file=paste("NBSpliceResSim", n.sim, ".RData", sep=""))
    DSGeneRes<-unique(results(myDSRes, filter = FALSE)[,c("gene", "geneFDR")])
    
    load(file=paste("NBSpliceResWFSim", n.sim, ".RData", sep=""))
    DSGeneResWF<-unique(results(myDSResWF, filter = FALSE)[,c("gene", "geneFDR")])

    load(file=paste("NBSpliceResWFRatSim", n.sim, ".RData", sep=""))
    DSGeneResWFR<-unique(results(myDSResWFR, filter = FALSE)[,c("gene", "geneFDR")])
    
    load(file=paste("NBSpliceResWFCountSim", n.sim, ".RData", sep=""))
    DSGeneResWFC<-unique(results(myDSResWFC, filter = FALSE)[,c("gene", "geneFDR")])
    
    load(file=paste("NBSpliceRes2Sim", n.sim, ".RData", sep=""))
    DSGeneRes2<-unique(results(myDSRes2, filter = FALSE)[,c("gene", "geneFDR")])
    
    load(file=paste("NBSpliceRes3Sim", n.sim, ".RData", sep=""))
    DSGeneRes3<-unique(results(myDSRes3, filter = FALSE)[,c("gene", "geneFDR")])
    
    load(file=paste("NBSpliceRes4Sim", n.sim, ".RData", sep=""))
    DSGeneRes4<-unique(results(myDSRes4, filter = FALSE)[,c("gene", "geneFDR")])

    load(file=paste("NBSpliceRes5Sim", n.sim, ".RData", sep=""))
    DSGeneRes5<-unique(results(myDSRes5, filter = FALSE)[,c("gene", "geneFDR")])

    load(file=paste("NBSpliceRes6Sim", n.sim, ".RData", sep=""))
    DSGeneRes6<-unique(results(myDSRes6, filter = FALSE)[,c("gene", "geneFDR")])
    
    # Extracting Differential splicing analysis results at the gene level
    padj.gene <- data.frame(row.names=unique(iso_info$gene_id[which(!is.na(iso_info$gene_id))]))
    padj.gene$NBSplice <- DSGeneRes$geneFDR[match(rownames(padj.gene), DSGeneRes$gene)]
    padj.gene$NBSpliceWF <- DSGeneResWF$geneFDR[match(rownames(padj.gene), DSGeneResWF$gene)]
    padj.gene$NBSpliceWFR <- DSGeneResWFR$geneFDR[match(rownames(padj.gene), DSGeneResWFR$gene)]
    padj.gene$NBSpliceWFC <- DSGeneResWFC$geneFDR[match(rownames(padj.gene), DSGeneResWFC$gene)]
    padj.gene$NBSplice2 <- DSGeneRes2$geneFDR[match(rownames(padj.gene), DSGeneRes2$gene)]
    padj.gene$NBSplice3 <- DSGeneRes3$geneFDR[match(rownames(padj.gene), DSGeneRes3$gene)]
    padj.gene$NBSplice4 <- DSGeneRes4$geneFDR[match(rownames(padj.gene), DSGeneRes4$gene)]
    padj.gene$NBSplice5 <- DSGeneRes5$geneFDR[match(rownames(padj.gene), DSGeneRes5$gene)]
    padj.gene$NBSplice6 <- DSGeneRes6$geneFDR[match(rownames(padj.gene), DSGeneRes6$gene)]
    
    # Computing performance results at the gene level on the simulated gene groups DS and DIEDS
    truth<-data.frame(status=as.numeric(rownames(padj.gene) %in% DSgenes),
                      row.names=rownames(padj.gene), 
                      geneGroup=geneInfo$simGroup[match(rownames(padj.gene), geneInfo$gene_id)],
                      geneSubGroup=geneInfo$simSubGroup[match(rownames(padj.gene), geneInfo$gene_id)])
        
    nas <- apply(padj.gene, 1, function(x) all(is.na(x)))
    padj.gene <- padj.gene[!nas,]
    truth <- truth[!nas,,drop=FALSE]
    DSResults<-cbind(padj.gene, truth$status)
    save(DSResults, file=paste("Qvals_geneGroup_sim", n.sim, ".RData",sep=""), compress = "xz")
    cd.gene <- COBRAData(padj=padj.gene, truth=truth)
    cp <- calculate_performance(cd.gene,thrs=c(0.01,0.05,0.1),
                                binary_truth="status",
                                aspect=c("tpr"), splv ="geneGroup")
    dfGeneG<-cp@tpr
    save(dfGeneG, file=paste("CobraPerfGeneGSim", n.sim, ".RData", sep=""), compress = "xz")
    cp <- calculate_performance(cd.gene,thrs=c(0.01,0.05,0.1),
                                binary_truth="status",
                                aspect=c("tpr"), splv ="geneSubGroup",maxsplit=length(levels(truth$geneSubGroup)))
    dfGeneSG<-cp@tpr
    save(dfGeneSG, file=paste("CobraPerfGeneSGSim", n.sim, ".RData", sep=""), compress = "xz")
    
}, mc.cores=4)

