##################################################
## Differential splicing analysis with NBSplice ##
##################################################

# run one time for each simulation
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
    # loading the expression, simulation information, and gene-isoform relation matrices.
    load(paste("iso_info_final_sim", n.sim,".RData", sep=""))
    load(paste("iso_cm_sim", n.sim,".RData", sep=""))
    load(paste("gene_iso_sim", n.sim,".RData", sep=""))
    # identifiyng the genes simulated as differentially spliced
    DSgenes<-unique(iso_info[ iso_info$DS  | iso_info$DIEDS , "gene_id"])
    id_mayorDS<-as.character(g2t_ok$transcript_id[g2t_ok$gene_id %in% DSgenes & iso_info$mayor])
    # building the desing Matrix
    sampleData<-data.frame(samples=colnames(iso_cm), condition=factor(c(rep("Normal",4), rep("Tumor",4)), levels=c("Normal", "Tumor")))
    geneIso<-iso_info[,c("gene_id", "transcript_id")]
    colnames(geneIso)[2]<-"isoform_id"
    rownames(sampleData)<-sampleData$samples
    
    #### RUNNING NBSplice with different parameters configuration. ####
    BPPARAM<-SerialParam()
    colName<-"condition"
    test<-"F"

    ## NBSplice without filtering, i.e. cT=0, rT=0 or do not calling the `buildLowExpIdx()` function
    
    myIsoDSWF<-IsoDataSet(isoCounts = iso_cm, colName =colName, geneIso = geneIso, experimentData = sampleData,BPPARAM = BPPARAM)
    # ratioThres<-0
    # countThres<-0
    # myIsoDSWF<-buildLowExpIdx(myIsoDSWF, colName=colName, ratioThres = ratioThres, countThres = countThres, BPPARAM = BPPARAM)
    # 
    # lowExpID<-lowExpIndex(myIsoDSWF)
    myDSResWF<-NBTest(object = myIsoDSWF, colName = colName, test=test, BPPARAM=BPPARAM)
    NBGenesWF<-GetDSGenes(myDSResWF)
    DSGeneResWF<-unique(results(myDSResWF)[,c("gene", "geneFDR")])
    save(myDSResWF, file=paste("NBSpliceResWFSim", n.sim, ".RData", sep=""), compress="xz")
    
    ## NBSplice cT=1
    ratioThres<-0
    countThres<-1
    myIsoDSWFR<-buildLowExpIdx(myIsoDSWF, ratioThres = ratioThres, countThres = countThres, colName=colName, BPPARAM = BPPARAM)
    myDSResWFR<-NBTest(object = myIsoDSWFR, colName = colName, test=test, BPPARAM=BPPARAM)
    NBGenesWFR<-GetDSGenes(myDSResWFR)
    DSGeneResWFR<-unique(results(myDSResWFR)[,c("gene", "geneFDR")])
    save(myDSResWFR, file=paste("NBSpliceResWFRatSim", n.sim, ".RData", sep=""), compress="xz")
    
    # NBSplice rT=0.01
    ratioThres<-0.01
    countThres<-0
    myIsoDSWFC<-buildLowExpIdx(myIsoDSWF, ratioThres = ratioThres, countThres = countThres, colName=colName, BPPARAM = BPPARAM)
    myDSResWFC<-NBTest(object = myIsoDSWFC, colName = colName, test=test, BPPARAM=BPPARAM)
    NBGenesWFC<-GetDSGenes(myDSResWFC)
    DSGeneResWFC<-unique(results(myDSResWFC)[,c("gene", "geneFDR")])
    save(myDSResWFC, file=paste("NBSpliceResWFCountSim", n.sim, ".RData", sep=""), compress="xz")
    
    # NBSplice rT=0.01;cT=1
    myIsoDS<-IsoDataSet(isoCounts = iso_cm, colName =colName, geneIso = geneIso, experimentData = sampleData,BPPARAM = BPPARAM)
    myIsoDS<-buildLowExpIdx(myIsoDS, colName=colName, BPPARAM = BPPARAM)
    myDSRes<-NBTest(object = myIsoDS, colName = colName, test=test, BPPARAM=BPPARAM)
    save(myDSRes, file=paste("NBSpliceResSim", n.sim, ".RData", sep=""), compress="xz") 
    NBGenes<-GetDSGenes(myDSRes)
    DSGeneRes<-unique(results(myDSRes)[,c("gene", "geneFDR")])
    
    # NBSplice rT=0.05;cT=1
    ratioThres<-0.05
    countThres<-1
    myIsoDS2<-buildLowExpIdx(myIsoDSWF, ratioThres = ratioThres, countThres = countThres, colName=colName, BPPARAM = BPPARAM)
    myDSRes2<-NBTest(object = myIsoDS2, colName = colName, test=test, BPPARAM=BPPARAM)
    NBGenes2<-GetDSGenes(myDSRes2)
    DSGeneRes2<-unique(results(myDSRes2)[,c("gene", "geneFDR")])
    save(myDSRes2, file=paste("NBSpliceRes2Sim", n.sim, ".RData", sep=""), compress="xz")
    
    # NBSplice rT=0.01;cT=2
    ratioThres<-0.01
    countThres<-2
    myIsoDS3<-buildLowExpIdx(myIsoDSWF, ratioThres = ratioThres, countThres = countThres, colName=colName, BPPARAM = BPPARAM)
    myDSRes3<-NBTest(object = myIsoDS3, colName = colName, test=test, BPPARAM=BPPARAM)
    NBGenes3<-GetDSGenes(myDSRes3)
    DSGeneRes3<-unique(results(myDSRes3)[,c("gene", "geneFDR")])
    save(myDSRes3, file=paste("NBSpliceRes3Sim", n.sim, ".RData", sep=""), compress="xz")
    
    # NBSplice rT=0.05;cT=2
    ratioThres<-0.05
    countThres<-2
    myIsoDS4<-buildLowExpIdx(myIsoDSWF, ratioThres = ratioThres, countThres = countThres, colName=colName, BPPARAM = BPPARAM)
    myDSRes4<-NBTest(object = myIsoDS4, colName = colName, test=test, BPPARAM=BPPARAM)
    NBGenes4<-GetDSGenes(myDSRes4)
    DSGeneRes4<-unique(results(myDSRes4)[,c("gene", "geneFDR")])
    save(myDSRes4, file=paste("NBSpliceRes4Sim", n.sim, ".RData", sep=""), compress="xz")
    
    # NBSplice rT=0.1;cT=1
    ratioThres<-0.1
    countThres<-1
    myIsoDS5<-buildLowExpIdx(myIsoDSWF, ratioThres = ratioThres, countThres = countThres, colName=colName, BPPARAM = BPPARAM)
    myDSRes5<-NBTest(object = myIsoDS5, colName = colName, test=test, BPPARAM=BPPARAM)
    NBGenes5<-GetDSGenes(myDSRes5)
    DSGeneRes5<-unique(results(myDSRes5)[,c("gene", "geneFDR")])
    save(myDSRes5, file=paste("NBSpliceRes5Sim", n.sim, ".RData", sep=""), compress="xz")

    # NBSplice rT=0.1;cT=2
    ratioThres<-0.1
    countThres<-2
    myIsoDS6<-buildLowExpIdx(myIsoDSWF, ratioThres = ratioThres, countThres = countThres, colName=colName, BPPARAM = BPPARAM)
    myDSRes6<-NBTest(object = myIsoDS6, colName = colName, test=test, BPPARAM=BPPARAM)
    NBGenes6<-GetDSGenes(myDSRes6)
    DSGeneRes6<-unique(results(myDSRes6)[,c("gene", "geneFDR")])
    save(myDSRes6, file=paste("NBSpliceRes6Sim", n.sim, ".RData", sep=""), compress="xz")
    # obtaining the adjusted p-values at the gene level indicating occurrence of differential splicing
    padj.gene <- data.frame(row.names=unique(iso_info$gene_id[which(!is.na(iso_info$transcript_id))]))
    padj.gene$NBSplice <- DSGeneRes$geneFDR[match(rownames(padj.gene), DSGeneRes$gene)]
    padj.gene$NBSpliceWF <- DSGeneResWF$geneFDR[match(rownames(padj.gene), DSGeneResWF$gene)]
    padj.gene$NBSpliceWFR <- DSGeneResWFR$geneFDR[match(rownames(padj.gene), DSGeneResWFR$gene)]
    padj.gene$NBSpliceWFC <- DSGeneResWFC$geneFDR[match(rownames(padj.gene), DSGeneResWFC$gene)]
    padj.gene$NBSplice2 <- DSGeneRes2$geneFDR[match(rownames(padj.gene), DSGeneRes2$gene)]
    padj.gene$NBSplice3 <- DSGeneRes3$geneFDR[match(rownames(padj.gene), DSGeneRes3$gene)]
    padj.gene$NBSplice4 <- DSGeneRes4$geneFDR[match(rownames(padj.gene), DSGeneRes4$gene)]
    padj.gene$NBSplice5 <- DSGeneRes5$geneFDR[match(rownames(padj.gene), DSGeneRes5$gene)]
    padj.gene$NBSplice6 <- DSGeneRes6$geneFDR[match(rownames(padj.gene), DSGeneRes6$gene)]
    # obtaining the true state of each gene
    truth<-data.frame(status=as.numeric(rownames(padj.gene) %in% DSgenes),
                      row.names=rownames(padj.gene))    
    nas <- apply(padj.gene, 1, function(x) all(is.na(x)))
    padj.gene <- padj.gene[!nas,]
    truth <- truth[!nas,,drop=FALSE]
    DSResults<-cbind(padj.gene, truth$status)
    # saving differential-splicing adjusted p-values
    save(DSResults, file=paste("Qvals_sim", n.sim, ".RData",sep=""), compress = "xz")
    # obtaining performance measures for all NBSplice configurations
    cd.gene <- COBRAData(padj=padj.gene, truth=truth)
    cp <- calculate_performance(cd.gene,
                                binary_truth="status",
                                aspect=c("fdrtpr","fdrtprcurve"),
                                thrs=c(.01,.05,.1))
    df<-cp@fdrtpr
    #saving results
    save(df, file=paste("CobraPerfSim", n.sim, ".RData", sep=""), compress = "xz")
    save.image(paste("sim", n.sim, "IsoCounts.RData", sep=""), compress="xz")
    save(iso_cm, file=paste("expressionMatrixSim", n.sim, ".RData", sep=""), compress="xz")
    save(iso_info, file=paste("isoInfoSim", n.sim, ".RData", sep=""), compress="xz")
    
}, mc.cores=16)

