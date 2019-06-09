###############################################
## Simulated counts analysis and exploration ##
###############################################
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
    
    ## NBSplice
    BPPARAM<-SerialParam()
    colName<-"condition"
    # loading NBSplice results using the optimal configuration
    load(paste("NBSpliceResSim", n.sim, ".RData", sep=""))
    DSRes <- NBSplice::results(myDSRes, filter=FALSE)
    genesEval<-unique(DSRes[,c("gene", "geneFDR")])
    # extracting the index of low-expressed isoforms
    lowExpIdx<-lowExpIndex(myDSRes)
       

    ##DRIMSeq with the optimal NBSplice filtering strategy
    
    counts <- data.frame(gene_id=geneIso$gene_id[-lowExpIdx],feature_id=geneIso$isoform_id[-lowExpIdx])
    counts<-cbind(counts,iso_cm[-lowExpIdx,])
    samps<-sampleData
    colnames(samps)[1]<-"sample_id"
    d <- dmDSdata(counts=counts, samples=samps)
    design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))
    d <- dmPrecision(d, design=design_full)    
    d <- dmFit(d, design=design_full)
    d <- dmTest(d, coef="conditionTumor")
    save(d, file=paste("DRIMSeq_Sim",n.sim,"_", ".RData", sep=""), compress="xz")

    res <- DRIMSeq::results(d)
    res.txp <- DRIMSeq::results(d, level="feature")
    DRIMGenes<-unique(res$gene_id[!is.na(res$adj_pvalue) & res$adj_pvalue<0.05])
    
    ##DRIMSeq with the SwimmF filtering strategy   
    counts2<- data.frame(gene_id=geneIso$gene_id,feature_id=geneIso$isoform_id)
    counts2<-cbind(counts2,iso_cm)
    samps<-sampleData
    colnames(samps)[1]<-"sample_id"
    d2 <- dmDSdata(counts=counts2, samples=samps)
    n <- 8
    n.small <- 4
    # Filtering low-reliable transcripts with the SwimmF strategy
    d2 <- dmFilter(d2,
                   min_samps_feature_expr=n.small, min_feature_expr=10,
                   min_samps_feature_prop=n.small, min_feature_prop=0.1,
                   min_samps_gene_expr=n, min_gene_expr=10)
    design_full <- model.matrix(~condition, data=DRIMSeq::samples(d2))
    
    d2 <- dmPrecision(d2, design=design_full)
    d2 <- dmFit(d2, design=design_full)
    d2 <- dmTest(d2, coef="conditionTumor")
    save(d2, file=paste("DRIMSeq_FiltSwim_Sim",n.sim,"_", ".RData", sep=""), compress="xz")
    
    res2 <- DRIMSeq::results(d2)
    res.txp2 <- DRIMSeq::results(d2, level="feature")
    DRIMGenes2<-unique(res2$gene_id[!is.na(res2$adj_pvalue) & res2$adj_pvalue<0.05])

    ##DEXSeq with the NBSplice filtering strategy
    
    count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
    sample.data <- DRIMSeq::samples(d)
    dxd <- DEXSeqDataSet(countData=count.data,
                         sampleData=sample.data,
                         design=~sample + exon + condition:exon,
                         featureID=counts(d)$feature_id,
                         groupID=counts(d)$gene_id)
    dxd <- estimateSizeFactors(dxd)
    dxd <- estimateDispersions(dxd, quiet=TRUE)
    dxd <- testForDEU(dxd, reducedModel=~sample + exon)  
    save(dxd, file=paste("DEXSeqSim", n.sim, ".RData", sep=""), compress="xz")

    dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
    qval <- perGeneQValue(dxr)
    dxr.g <- data.frame(gene=names(qval),qval)    
    columns <- c("featureID","groupID","pvalue", "padj")
    dxr <- as.data.frame(dxr[,columns])
    
    ## DEXSeq with the SwimmF filtering strategy    
    count.data2 <- round(as.matrix(counts(d2)[,-c(1:2)]))
    sample.data <- DRIMSeq::samples(d2)
    dxd2 <- DEXSeqDataSet(countData=count.data2,
                         sampleData=sample.data,
                         design=~sample + exon + condition:exon,
                         featureID=counts(d2)$feature_id,
                         groupID=counts(d2)$gene_id)
    dxd2 <- estimateSizeFactors(dxd2)
    dxd2 <- estimateDispersions(dxd2, quiet=TRUE)
    dxd2 <- testForDEU(dxd2, reducedModel=~sample + exon)
    save(dxd2, file=paste("DEXSeqFiltSwim_Sim", n.sim, ".RData", sep=""),	compress="xz")
    dxr2 <- DEXSeqResults(dxd2, independentFiltering=FALSE)
    qval2 <- perGeneQValue(dxr2)
    dxr.g2 <- data.frame(gene=names(qval2),qval2)
    dxr2 <- as.data.frame(dxr2[,columns])
    
    ## NBSplice with the SwimmF filtering strategy
    geneIso<-data.frame(isoform_id=counts(d2)$feature_id,gene_id=counts(d2)$gene_id)
    rownames(geneIso)<-geneIso$isoform_id
    colnames(count.data2)<-rownames(sample.data)<-sample.data$sample_id
    rownames(count.data2)<-geneIso$isoform_id
    myIsoDS2<-IsoDataSet(isoCounts = count.data2, experimentData = sample.data, colName = colName, 
                         geneIso = geneIso, BPPARAM = BPPARAM)
    test<-"F"
    myDSRes2<-NBTest(object = myIsoDS2, colName = colName,test =test, BPPARAM = BPPARAM)
    DSRes2 <- NBSplice::results(myDSRes2, filter=FALSE)
    genesEval2<-unique(DSRes2[,c("gene", "geneFDR")])
    
    padj.gene <- data.frame(row.names=unique(iso_info$gene_id[which(!is.na(iso_info$transcript_id))]))
    # NBSplice filtering
    padj.gene$DRIMSeq <- res$adj_pvalue[match(rownames(padj.gene), res$gene_id)]
    padj.gene$DEXSeq <- dxr.g$qval[match(rownames(padj.gene), dxr.g$gene)]
    padj.gene$NBSplice <- genesEval$geneFDR[match(rownames(padj.gene), genesEval$gene)]
    # SwimmF filtering
    padj.gene$DRIMSeq2 <- res2$adj_pvalue[match(rownames(padj.gene), res2$gene_id)]
    padj.gene$DEXSeq2 <- dxr.g2$qval[match(rownames(padj.gene), dxr.g2$gene)]
    padj.gene$NBSplice2 <- genesEval2$geneFDR[match(rownames(padj.gene), genesEval2$gene)]

    # defining the true state of each gene
    truth<-data.frame(status=as.numeric(rownames(padj.gene) %in% DSgenes),
                      row.names=rownames(padj.gene))
        
    nas <- apply(padj.gene, 1, function(x) all(is.na(x)))
    padj.gene <- padj.gene[!nas,]
    truth <- truth[!nas,,drop=FALSE]
    DSResults<-cbind(padj.gene, truth$status)
    save(DSResults, file=paste("Qvals_FiltSwimm_sim", n.sim, ".RData",sep=""), compress = "xz")
    # obtaining performance measures for all NBSplice configurations
    cd.gene <- COBRAData(padj=padj.gene, truth=truth)
    cp <- calculate_performance(cd.gene,
                                binary_truth="status",
                                aspect=c("fdrtpr","fdrtprcurve"),
                                thrs=c(.01,.05,.1))
    df<-cp@fdrtpr
    save(df, file=paste("CobraPerfDTUMethodsFiltSwimmSim", n.sim, "RData", sep=""), compress = "xz")
    
}, mc.cores=16)

