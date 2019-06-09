library("tidyverse")
library("ggplot2")
library("cowplot")

setwd("path_DS_analysis_results")
stop()
rm(list=ls())
# Number of experiment realizations
n<-10

# Loading NBSplice results over its nine configurations at the gene level
load("sim1/CobraPerfSim1.RData")
dfAll<-df
dfAll$sim<-1
for(i in 2:n){
    load(paste("sim", i, "/CobraPerfSim",i,".RData", sep=""))
    df$sim<-i
    dfAll<-rbind(dfAll, df)
}
dfAll$sim<-factor(dfAll$sim)
levels(dfAll$sim)<-paste("sim",1:10, sep="")
# Correcting the negative counts becacuse iCOBRA do not consider NA's
dfAll$FNAll<-dfAll$TP/dfAll$TPR-dfAll$TP
dfAll$Fscore<-2*dfAll$TP/(2*dfAll$TP+dfAll$FP+dfAll$FNAll)
dfAll$TNAll<-12602-dfAll$TP-dfAll$FNAll-dfAll$FP
dfAll$Acc<-(dfAll$TP+dfAll$TNAll)/12602

# Figure 1
dfAll$method<-as.factor(dfAll$method)
levels(dfAll$method)<-c("rT=0.01;cT=1","rT=0.05;cT=1","rT=0.01;cT=2","rT=0.05;cT=2","rT=0.1;cT=1","rT=0.1;cT=2","WF","rT=0.01","cT=1")
dfAll$method<-factor(dfAll$method, levels = c("WF","rT=0.01","cT=1","rT=0.01;cT=1","rT=0.01;cT=2","rT=0.05;cT=1","rT=0.05;cT=2","rT=0.1;cT=1","rT=0.1;cT=2"))
thresholds<-c(0.01,0.05,0.1)
levels(dfAll$thr)<-thresholds
p0<-ggplot(dfAll[dfAll$thr=="0.05", ], aes(x=method, y=Acc, color=method, fill=method))+geom_boxplot(alpha=0.5)+labs(x="", y="Accuracy")+
    scale_color_manual(breaks=levels(dfAll$method), values=plotcolors, name="NBSplice configuration")+
    theme_bw()+theme(legend.position = "none", axis.text.x=element_blank() )+guides(col=guide_legend(ncol=3))+scale_fill_manual(breaks=levels(
        dfAll$method), values=plotcolors, name="NBSplice configuration")+scale_y_continuous(breaks = seq(0.8:1, by=0.05), limits = c(0.8,1))

p1<-ggplot(dfAll[dfAll$thr=="0.05", ], aes(x=method, y=TPR, color=method, fill=method))+geom_boxplot(alpha=0.5)+scale_color_manual(breaks=levels(dfAll$method), values=plotcolors, name="NBSplice configuration")+
    theme_bw()+scale_fill_manual(breaks=levels(
        dfAll$method), values=plotcolors, name="NBSplice configuration")+theme(
            legend.position = "bottom",axis.text.x=element_blank() )+guides(col=guide_legend(ncol=9,title.hjust=0.5, title.position="top"), fill=guide_legend(ncol=9,title.hjust=0.5, title.position="top"))+ labs(x="",y="Sensitivity")+scale_y_continuous(breaks = seq(0.4:1, by=0.1), limits = c(0.4,1))

p2<-ggplot(dfAll[dfAll$thr=="0.05", ], aes(x=method, y=1-FDR, color=method, fill=method))+geom_boxplot(alpha=0.5)+scale_color_manual(breaks=levels(dfAll$method), values=plotcolors, name="NBSplice configuration")+
    theme_bw()+labs(x="", y="Precision")+theme(legend.position = "none", axis.text.x=element_blank() 
)+guides(col=guide_legend(ncol=9))+scale_fill_manual(breaks=levels(
    dfAll$method), values=plotcolors, name="NBSplice configuration")+scale_y_continuous(
        breaks = seq(0.4:1, by=0.1), limits = c(0.4,1))+geom_hline(yintercept= 0.95,linetype ="dashed", color="#666666")

p3<-ggplot(dfAll[dfAll$thr=="0.05", ], aes(x=method, y=Fscore, color=method, fill=method))+geom_boxplot(alpha=0.5)+ scale_color_manual(breaks=levels(dfAll$method), values=plotcolors, name="NBSplice configuration")+
    theme_bw()+labs(x="", y="F-score")+theme(legend.position = "none", axis.text.x=element_blank() )+guides(
        col=guide_legend(ncol=3))+scale_fill_manual(breaks=levels(
            dfAll$method), values=plotcolors, name="NBSplice configuration")+scale_y_continuous(breaks = seq(0.4:1, by=0.1), limits = c(0.4,1))

legend <- get_legend(p1+theme(legend.position="bottom", 
                              legend.text = element_text(size=12),
                              legend.title = element_text(size=12),
                              legend.spacing.x = unit(0.25, 'cm'), 
                              legend.spacing.y = unit(0.25, 'cm'),
                              legend.key.height = unit(0.6, "cm"),
                              legend.key.width = unit(2.5, "cm"),
                              )+guides(fill=guide_legend("NBSplice\nconfigurations",ncol=3),
                                       color=guide_legend("NBSplice\nconfigurations",ncol=3)))
                            

pF<-plot_grid(p0+scale_y_continuous(limits=c(0.85,0.95), expand = c(0.1, 0))+theme(axis.title = element_text(size=12)),
               p1+scale_y_continuous(limits=c(0.4,1), expand = c(0, 0))+theme(axis.title = element_text(size=12)),
               p2+scale_y_continuous(limits=c(0.4,1), expand = c(0, 0))+theme(axis.title = element_text(size=12)),
               p3+scale_y_continuous(limits=c(0.4,1), expand = c(0, 0))+theme(axis.title = element_text(size=12)),
               nrow=1, labels = c("A", "B", "C", "D"))
pF1<-plot_grid(pF,legend,ncol=1, rel_heights = c(1,0.15))

ggplot2::ggsave(pF1, file="overallPerformanceProstThr005V3.pdf", height=8, width=9, dpi=450)

# Summarizing performance measures results
dfMean<-summarise(group_by(dfAll, thr, method), FDR=mean(FDR), TPR=mean(TPR), Fscore=mean(Fscore), Acc=mean(Acc))
dfMean$satis<-dfMean$FDR<=as.numeric(do.call(c, lapply(as.character(dfMean$thr), function(x){strsplit(x, "thr")[[1]][2]})))
dfMean$satis[dfMean$satis]<-"yes"
dfMean$satis[dfMean$satis=="FALSE"]<-"no"
dfMean$method<-as.factor(dfMean$method)
levels(dfMean$method)<-c("rT=0.01;cT=1","rT=0.05;cT=1","rT=0.01;cT=2","rT=0.05;cT=2","rT=0.1;cT=1","rT=0.1;cT=2","WF","rT=0.01","cT=1")
dfMean$method<-factor(as.character(dfMean$method), levels=c("WF","rT=0.01","cT=1","rT=0.01;cT=1","rT=0.01;cT=2","rT=0.05;cT=1","rT=0.05;cT=2","rT=0.1;cT=1","rT=0.1;cT=2"))
dfMean$method.satis<-paste0(dfMean$satis)
dfMean$method.satis[dfMean$method.satis=="yes"]<-as.character(dfMean$method[dfMean$method.satis=="yes"])
dfMean$method.satis<-as.factor(dfMean$method.satis)

## Plotting summarized results
nthr<-length(levels(dfMean$thr))
nlevs <- length(unique(dfMean$method))
plotcolors<-c("WF"="black","rT=0.01"="mediumpurple3","cT=1"="royalblue4","rT=0.01;cT=1"="firebrick3","rT=0.01;cT=2"="deeppink2","rT=0.05;cT=1"="darkorange2","rT=0.05;cT=2"="gold3","rT=0.1;cT=1"="steelblue2","rT=0.1;cT=2"="palegreen4")
plotFillColors<-c(plotcolors, "no"="white")
thresholds<-c(0.01,0.05,0.1)
levels(dfMean$thr)<-thresholds
dfMean$satis<-as.factor(dfMean$satis)

# Figure 2
pF2<-ggplot(dfMean, aes(x=FDR, y=TPR, group=method))+geom_point(size=5,
aes(fill=method.satis,colour=method, shape=thr, ), alpha=0.7, stroke=1
)+geom_line(aes(group=method,colour=method), size=1,linetype ="dashed")+ geom_vline(xintercept= thresholds,linetype ="dotted"
)+ylim(c(0,0.75))+xlim(0,0.15)+scale_shape_manual(values = rep(21, nthr),guide = FALSE
) +scale_fill_manual(values=plotFillColors,guide=FALSE, name="")+scale_colour_manual(
    values=plotcolors, name="NBSplice configuration")+theme_bw()+theme(legend.position = "bottom")+guides(col=guide_legend(ncol=3))+labs(y="Sensitivity")
ggsave(pF2, file="ROCNBSpliceConfigs.pdf", height = 5, width=6.2, dpi=500)
# Wilcoxon tests
dfAux<-dfAll[dfAll$thr==0.05,]
meth<-levels(dfAux$method)
AccComp<-as.data.frame(do.call(rbind, lapply(1:(length(meth)-1),function(i){
    comp<-NULL
    for(j in (i+1):length(meth)){
        pval<-wilcox.test(dfAux$Acc[dfAux$method==meth[i]], dfAux$Acc[dfAux$method==meth[j]], exact=F, paired=T)$p.value
        names(pval)<-meth[j]
        comp<-c(comp, pval)
    }
   return(cbind(meth[i], names(comp), comp))
})))

wilcox.test(dfAux$Acc[dfAux$method=="rT=0.05;cT=1"], dfAux$Acc[dfAux$method=="rT=0.05;cT=2"], exact=F, paired=T)
wilcox.test(dfAux$FDR[dfAux$method=="rT=0.01"], dfAux$FDR[dfAux$method=="rT=0.01;cT=1"], exact=F, paired=T)


# Loading NBSplice results over its nine configurations at the isoform level ##

load("sim1/CobraPerfIsoSim1.RData")
dfIso$method<-factor(dfIso$method)
dfAllIso<-dfIso
dfAllIso$sim<-1
for(i in 2:n){
    load(paste("sim", i, "/CobraPerfIsoSim",i,".RData", sep=""))
    dfIso$sim<-i
    dfIso$method<-factor(dfIso$method)
dfAllIso<-rbind(dfAllIso, dfIso)
}
dfAllIso$sim<-factor(dfAllIso$sim)
levels(dfAllIso$sim)<-paste("sim",1:10, sep="")
# correcting iCOBRA results to consider transcript identified as NAs
dfAllIso$FNAll<-dfAllIso$TP/dfAllIso$TPR-dfAllIso$TP
dfAllIso$Fscore<-2*dfAllIso$TP/(2*dfAllIso$TP+dfAllIso$FP+dfAllIso$FNAll)
dfAllIso$TNAll<-104829-dfAllIso$TP-dfAllIso$FNAll-dfAllIso$FP
dfAllIso$Acc<-(dfAllIso$TP+dfAllIso$TNAll)/104829
dfAllIso$P<-dfAllIso$TP+dfAllIso$FP

dfAllIso$method<-as.factor(dfAllIso$method)
levels(dfAllIso$method)<-c("rT=0.01;cT=1","rT=0.05;cT=1","rT=0.01;cT=2","rT=0.05;cT=2","rT=0.1;cT=1","rT=0.1;cT=2","WF","rT=0.01","cT=1")
dfAllIso$method<-factor(dfAllIso$method, levels = c("WF","rT=0.01","cT=1","rT=0.01;cT=1","rT=0.01;cT=2","rT=0.05;cT=1","rT=0.05;cT=2","rT=0.1;cT=1","rT=0.1;cT=2"))
thresholds<-c(0.01,0.05,0.1)
levels(dfAllIso$thr)<-thresholds

dfMeanIso<-as.data.frame(summarise(group_by(dfAllIso, thr, method), 
                            TP=mean(TP), FP=mean(FP), 
                            TOT_CALLED=mean(TOT_CALLED), 
                            P=mean(TP+FP), Acc=mean(Acc), TPR=mean(TPR), FDR=mean(FDR), Fscore=mean(Fscore)))

dfSDIso<-as.data.frame(summarise(group_by(dfAllIso, thr, method), 
                              TP=sd(TP), FP=sd(FP), 
                              TOT_CALLED=sd(TOT_CALLED), 
                              P=sd(P), Acc=sd(Acc), 
                              TPR=sd(TPR), FDR=sd(FDR), Fscore=sd(Fscore)))

for(i in c(7:10)){
    dfMeanIso[,i]<-round(dfMeanIso[,i],3)
    dfSDIso[,i]<-round(dfSDIso[,i],3)
    
}
for(i in c(3:6)){
    dfSDIso[,i]<-round(dfSDIso[,i])
    dfMeanIso[,i]<-round(dfMeanIso[,i])
}

dfMeanIso[dfMeanIso$thr=="0.05", c("method", "TOT_CALLED", "P", "TP")]
dfSDIso[dfSDIso$thr=="0.05", c("method", "TOT_CALLED", "P", "TP")]

df[df$thr==0.05, c("method", "TOT_CALLED", "P", "TP")]
dfSD[dfSD$thr==0.05, c("method", "TOT_CALLED", "P", "TP")]

dfAux<-dfAll[dfAll$thr==0.05,c("method", "TOT_CALLED", "P", "TP")]
dfAuxIso<-dfAllIso[dfAllIso$thr==0.05,c("method", "TOT_CALLED", "P", "TP")]
# Comparing results
for(i in meth){
  for(j in meth[meth !=i]){
    print(paste("Total:", i, "vs", j,wilcox.test(dfAux$TOT_CALLED[dfAux$method==i], dfAux$TOT_CALLED[dfAux$method==j], exact=F, paired=T)$p.value, sep=" "))
    print(paste("Total:", i, "vs", j,wilcox.test(dfAux$P[dfAux$method==i], dfAux$P[dfAux$method==j], exact=F, paired=T), exact=F, paired=T)$p.value, sep=" "))
    print(paste("Total:", i, "vs", j,wilcox.test(dfAux$TP[dfAux$method==i], dfAux$TP[dfAux$method==j], exact=F, paired=T), exact=F, paired=T)$p.value, sep=" "))
  }
}

# Loading NBSplice results over its nine configurations and for the simulated groups and soubgroups of genes

rm(list=ls())
n<-10
# loading the simulated groups' results
load("sim1/CobraPerfGeneGSim1.RData")
dfAll<-dfGeneG
dfAll$sim<-1
for(i in 2:n){
    load(paste("sim", i, "/CobraPerfGeneGSim",i,".RData", sep=""))
    dfGeneG$sim<-i
    dfAll<-rbind(dfAll, dfGeneG)
    
}
dfAll$sim<-factor(dfAll$sim)
levels(dfAll$sim)<-paste("sim",1:10, sep="")
dfAll$FNAll<-dfAll$TP/dfAll$TPR-dfAll$TP
dfAll$TNAll<-12602-dfAll$TP-dfAll$FNAll-dfAll$FP

dfMean<-summarise(group_by(dfAll, thr, method, splitval), TPR=mean(TPR))
dfMean$method<-as.factor(dfMean$method)
levels(dfMean$method)<-c("rT=0.01;cT=1","rT=0.05;cT=1","rT=0.01;cT=2","rT=0.05;cT=2","rT=0.1;cT=1","rT=0.1;cT=2","WF","rT=0.01","cT=1")
dfMean$method<-factor(as.character(dfMean$method), levels=c("WF","rT=0.01","cT=1","rT=0.01;cT=1","rT=0.01;cT=2","rT=0.05;cT=1","rT=0.05;cT=2","rT=0.1;cT=1","rT=0.1;cT=2"))

nlevs <- length(unique(dfMean$method))
plotcolors<-c("WF"="black","rT=0.01"="mediumpurple3","cT=1"="royalblue4","rT=0.01;cT=1"="firebrick3","rT=0.01;cT=2"="deeppink2","rT=0.05;cT=1"="darkorange2","rT=0.05;cT=2"="gold3","rT=0.1;cT=1"="steelblue2","rT=0.1;cT=2"="palegreen4")
plotFillColors<-c(plotcolors, "no"="white")
thresholds<-c(0.01,0.05,0.1)
levels(dfMean$thr)<-thresholds

dfAll$method<-as.factor(dfAll$method)
levels(dfAll$method)<-c("rT=0.01;cT=1","rT=0.05;cT=1","rT=0.01;cT=2","rT=0.05;cT=2","rT=0.1;cT=1","rT=0.1;cT=2","WF","rT=0.01","cT=1")
dfAll$method<-factor(dfAll$method, levels = c("WF","rT=0.01","cT=1","rT=0.01;cT=1","rT=0.01;cT=2","rT=0.05;cT=1","rT=0.05;cT=2","rT=0.1;cT=1","rT=0.1;cT=2"))
thresholds<-c(0.01,0.05,0.1)
levels(dfAll$thr)<-thresholds
levels(dfAll$splitval)[2:3]<-c("DIEDS", "DS")
dfAll$splitval<-factor(as.character(dfAll$splitval), levels=c("DS", "DIEDS", "overall"))

# keeping results for 0.05 significance threshold
dfAux<-dfAll[dfAll$splitval != "overall" &  dfAll$thr=="0.05",]
meth<-levels(dfAux$method)

# loading the simulated subgroups' results
load("sim1/CobraPerfGeneSGSim1.RData")
dfAllSG<-dfGeneSG
dfAllSG$sim<-1
for(i in 2:n){
    load(paste("sim", i, "/CobraPerfGeneSGSim",i,".RData", sep=""))
    dfGeneSG$sim<-i
    dfAllSG<-rbind(dfAllSG, dfGeneSG)
}
dfAllSG$sim<-factor(dfAllSG$sim)
levels(dfAllSG$sim)<-paste("sim",1:10, sep="")

dfAllSG$FNAll<-dfAllSG$TP/dfAllSG$TPR-dfAllSG$TP
dfMean<-summarise(group_by(dfAllSG, thr, method, splitval), TPR=mean(TPR))
dfMean$method<-as.factor(dfMean$method)
levels(dfMean$method)<-c("rT=0.01;cT=1","rT=0.05;cT=1","rT=0.01;cT=2","rT=0.05;cT=2","rT=0.1;cT=1","rT=0.1;cT=2","WF","rT=0.01","cT=1")
dfMean$method<-factor(as.character(dfMean$method), levels=c("WF","rT=0.01","cT=1","rT=0.01;cT=1","rT=0.01;cT=2","rT=0.05;cT=1","rT=0.05;cT=2","rT=0.1;cT=1","rT=0.1;cT=2"))

nlevs <- length(unique(dfMean$method))
plotcolors<-c("WF"="black","rT=0.01"="mediumpurple3","cT=1"="royalblue4","rT=0.01;cT=1"="firebrick3","rT=0.01;cT=2"="deeppink2","rT=0.05;cT=1"="darkorange2","rT=0.05;cT=2"="gold3","rT=0.1;cT=1"="steelblue2","rT=0.1;cT=2"="palegreen4")
plotFillColors<-c(plotcolors, "no"="white")
thresholds<-c(0.01,0.05,0.1)
levels(dfMean$thr)<-thresholds

dfMean[dfMean$splitval != "overall" & dfMean$thr == "0.05",]
dfAllSG$method<-as.factor(dfAllSG$method)
levels(dfAllSG$method)<-c("rT=0.01;cT=1","rT=0.05;cT=1","rT=0.01;cT=2","rT=0.05;cT=2","rT=0.1;cT=1","rT=0.1;cT=2","WF","rT=0.01","cT=1")
dfAllSG$method<-factor(dfAllSG$method, levels = c("WF","rT=0.01","cT=1","rT=0.01;cT=1","rT=0.01;cT=2","rT=0.05;cT=1","rT=0.05;cT=2","rT=0.1;cT=1","rT=0.1;cT=2"))
thresholds<-c(0.01,0.05,0.1)
levels(dfAllSG$thr)<-thresholds
levels(dfAllSG$splitval)[2:29]<-do.call(rbind, strsplit(levels(dfAllSG$splitval)[2:29], split=":"))[,2]
dfAllSG$SimSubGroup<-as.character(dfAllSG$splitval)
dfAllSG$SimSubGroup[dfAllSG$SimSubGroup =="0.5-0.3"]<-"0.3-0.5"
dfAllSG$SimSubGroup[dfAllSG$SimSubGroup=="0.7-0.3"]<-"0.3-0.7"
dfAllSG$SimSubGroup[dfAllSG$SimSubGroup =="0.7-0.5"]<-"0.5-0.7" 
dfAllSG$SimSubGroup[dfAllSG$SimSubGroup =="0.9-0.5"]<- "0.5-0.9"
dfAllSG$SimSubGroup[dfAllSG$SimSubGroup =="0.5-0.5-0.3"]<-"2-0.3-0.5" 
dfAllSG$SimSubGroup[dfAllSG$SimSubGroup =="0.5-0.7-0.3"]<-"2-0.3-0.7"
dfAllSG$SimSubGroup[dfAllSG$SimSubGroup =="0.5-0.3-0.5"]<-"2-0.5-0.3"
dfAllSG$SimSubGroup[dfAllSG$SimSubGroup =="0.5-0.3-0.7" ]<-"2-0.7-0.3"  

dfAllSG$SimSubGroup[dfAllSG$SimSubGroup =="0.25-0.5-0.7" ]<-"4-0.7-0.5"  
dfAllSG$SimSubGroup[dfAllSG$SimSubGroup =="0.25-0.7-0.5" ]<-"4-0.5-0.7"  
dfAllSG$SimSubGroup[dfAllSG$SimSubGroup =="0.25-0.9-0.5" ]<-"4-0.5-0.9"   
dfAllSG$SimSubGroup[dfAllSG$SimSubGroup =="0.25-0.5-0.9" ]<-"4-0.9-0.5" 

dfAllSG$SimSubGroup<-factor(as.character(dfAllSG$SimSubGroup), levels=c("overall", 
"0.3-0.5", "0.3-0.7", "0.5-0.7", "0.5-0.9","0.5-0.5-0.7","0.5-0.5-0.9"
,"2-0.5-0.3","2-0.3-0.5","2-0.3-0.7","2-0.7-0.3","2-0.5-0.7","2-0.5-0.9", "4-0.5-0.7",
"4-0.7-0.5", "4-0.5-0.9","4-0.9-0.5" ))

groupSubgroup<-data.frame(group=c(rep("DS",4), rep("DIEDS", 12)), subgroup=levels(dfAllSG$SimSubGroup)[2:17])

dfAllSG$SimGroup<-factor(as.character(groupSubgroup[match(dfAllSG$SimSubGroup,groupSubgroup$subgroup),"group"]), levels=c("DS", "DIEDS"))
    
# Exploring results of groups and subgroups for the optimal NBSplice configuration at a significance threshold of 0.05
# Figure 3
p1<-ggplot(dfAll[dfAll$splitval != "overall" & dfAll$thr=="0.05" & dfAll$method=="rT=0.01;cT=1",], 
           aes(x=splitval, y=TPR, color=splitval, fill=splitval))+geom_boxplot(alpha=0.5)+
    ylim(c(0.45,0.6))+scale_fill_manual(values=c("DS"="red3", "DIEDS"="royalblue4"))+scale_color_manual(values=c("DS"="red3", "DIEDS"="royalblue4"))+
    theme_bw()+theme(legend.position = "none",axis.text.x = element_text(angle=90),strip.text=element_text(size=12), 
                     panel.background=element_rect(fill="white"), 
                     panel.grid.minor=element_line(colour="black", linetype="dashed"),
                     panel.border = element_rect( colour = "black", fill=NA))+guides(
                         col=guide_legend(ncol=3))+ labs(x="Simulation Group", y="Sensitivity", color="Group",fill="Group")
levels(dfAllSG$SimSubGroup)
cols<-c("0.3-0.5" ="#FF6633",  "0.3-0.7" ="#FF0000",  "0.5-0.7"="#CC0000", "0.5-0.9"="#990033",
        "0.5-0.5-0.7"="#66CCCC","0.5-0.5-0.9" ="#00CCCC",
        "2-0.5-0.3"="#33CCFF", "2-0.3-0.5"="#0099CC", "2-0.3-0.7"="#0099FF" ,
        "2-0.7-0.3"="#0066FF", "2-0.5-0.7"="#0033FF",    "2-0.5-0.9"="#0033CC",
        "4-0.5-0.7" ="#000099",  "4-0.7-0.5"="#000066" ,"4-0.5-0.9"="#003366" , 
        "4-0.9-0.5"="#000033" 
        )
p2<-ggplot(dfAllSG[dfAllSG$SimSubGroup != "overall" & (dfAllSG$SimGroup != "") & dfAllSG$thr=="0.05" & dfAllSG$method=="rT=0.01;cT=1",], aes(
    x=SimSubGroup, y=TPR, color=SimSubGroup, fill=SimSubGroup))+geom_boxplot(alpha=0.5)+theme_bw()+theme(
                legend.position = "none",axis.text.x = element_text(angle = 90),strip.text=element_text(size=12), 
                panel.background=element_rect(fill="white"), 
                panel.grid.minor=element_line(colour="black", linetype="dashed"),
                panel.border = element_rect( colour = "black", fill=NA))+guides(col=guide_legend(ncol=9))+ labs(
                    x="Simulation Subgroup", y="Sensitivity")+facet_grid(~SimGroup,
            scales="free_x", space="free_x")+ylim(0.1,0.95)+scale_fill_manual(values=cols)+scale_color_manual(values=cols)
p2<-p2+theme(legend.position="none")
p1<-p1+ylim(0.1,0.95)

TPRSimGroupsBEST<-plot_grid(plot_grid(ggplot(),p1, ggplot(), nrow =3, rel_heights = c(0.04,0.91,0.04))
                        ,p2,nrow=1, rel_widths = c(0.2,0.8), labels=c("A", "B"))
ggplot2::ggsave(TPRSimGroupsBEST, file="TPRSimGroupsOptimu.pdf", height = 7, width=11, dpi=400)

subGro<-as.character(unique(dfAux$SimSubGroup))

for (i in 1:length(subGro)){
    for(j in 1:length(subGro)){
        gi<-unique(dfAux$SimGroup[dfAux$method=="rT=0.01;cT=1" & dfAux$splitval == subGro[i]])
        gj<-unique(dfAux$SimGroup[dfAux$method=="rT=0.01;cT=1" & dfAux$splitval == subGro[j]])
        
        if(i !=j & gi ==gj ){
    print(paste(subGro[i], "-",subGro[j], ": ", round(wilcox.test(dfAux$TPR[dfAux$method=="rT=0.01;cT=1" & dfAux$splitval == subGro[i]],
                      dfAux$TPR[dfAux$method=="rT=0.01;cT=1" & dfAux$splitval == subGro[j]], 
                      exact=F, paired=T)$p.value,3), sep=""))
    }}}
    
#######################################################
## NBSplice, DEXSeq, and DRIMSeq results' comparison ##
#######################################################

n<-10
load("sim1/CobraPerfDTUMethodsFiltSwimmSim1RData")
df$method<-factor(df$method)
dfAll<-df
dfAll$sim<-1

for(i in 2:n){
    load(paste("sim", i, "/CobraPerfDTUMethodsFiltSwimmSim",i,"RData", sep=""))
    df$sim<-i
    df$method<-factor(df$method)
    dfAll<-rbind(dfAll, df)
}
dfAll$sim<-factor(dfAll$sim)
levels(dfAll$sim)<-paste("sim",1:10, sep="")
dfAll$FNAll<-dfAll$TP/dfAll$TPR-dfAll$TP
dfAll$Fscore<-2*dfAll$TP/(2*dfAll$TP+dfAll$FP+dfAll$FNAll)
dfAll$TNAll<-12602-dfAll$TP-dfAll$FNAll-dfAll$FP

dfAll$Acc<-(dfAll$TP+dfAll$TNAll)/12602
levels(dfAll$method)<-c("DEXSeq:rT=0.01;cT=1","DEXSeq:SwimmF","DRIMSeq:rT=0.01;cT=1","DRIMSeq:SwimmF","NBSplice:rT=0.01;cT=1","NBSplice:SwimmF")

thresholds<-c(0.01,0.05,0.1)
levels(dfAll$thr)<-thresholds
# Figure 4
p4<-ggplot(dfAll[dfAll$thr==0.05,], aes(x=method, y=Acc, color=method, fill=method))+geom_boxplot(alpha=0.5)+labs(x="", y="Accuracy")+scale_color_manual(breaks=levels(dfAll$method), values=plotcolors, name="Methods")+scale_fill_manual(breaks=levels(dfAll$method), values=plotcolors, name="Methods")+
    theme_bw()+theme(legend.position = "none", axis.text.x=element_blank())+guides(col=guide_legend(ncol=3))

p5<-ggplot(dfAll[dfAll$thr==0.05,], aes(x=method, y=TPR, color=method, fill=method))+geom_boxplot(alpha=0.5)+
scale_color_manual(breaks=levels(dfAll$method), values=plotcolors, name="Methods")+scale_fill_manual(breaks=levels(dfAll$method), values=plotcolors, name="Methods")+
    theme_bw()+theme(legend.position = "none", axis.text.x=element_blank())+guides(col=guide_legend(ncol=3))+ labs(x="", y="Sensitivity")

p6<-ggplot(dfAll[dfAll$thr==0.05,], aes(x=method, y=1-FDR, color=method, fill=method))+geom_boxplot(alpha=0.5)+
scale_color_manual(breaks=levels(dfAll$method), values=plotcolors, name="Methods")+scale_fill_manual(breaks=levels(dfAll$method), values=plotcolors, name="Methods")+geom_hline(yintercept= 0.95,linetype ="dashed", color="#666666")+
    theme_bw()+labs(x="", y="Precision")+theme(legend.position = "bottom", axis.text.x=element_blank())+guides(col=guide_legend(ncol=3))

p7<-ggplot(dfAll[dfAll$thr==0.05,], aes(x=method, y=Fscore, color=method, fill=method))+geom_boxplot(alpha=0.5)+
scale_color_manual(breaks=levels(dfAll$method), values=plotcolors, name="Methods")+scale_fill_manual(breaks=levels(dfAll$method), values=plotcolors, name="Methods")+
    theme_bw()+labs(x="", y="F-score")+theme(legend.position = "none", axis.text.x=element_blank())+guides(col=guide_legend(ncol=3))
p4<-p4+scale_y_continuous(limits = c(0.85,0.95), expand = c(0, 0))
p5<-p5+scale_y_continuous(limits = c(0.25,1),expand = c(0, 0))
p6<-p6+scale_y_continuous(limits = c(0.25,1),expand = c(0, 0))
p7<-p7+scale_y_continuous(limits = c(0.25,1),expand = c(0, 0))
legend <- get_legend(p6)
p6<-p6+theme(legend.position = "none")
p2<-plot_grid(p4,p5,p6,p7, nrow=1, labels = c("A", "B", "C", "D"))
pF2<-plot_grid(p2,legend,ncol=1, rel_heights = c(1,0.1))

ggplot2::ggsave(pF2, file="overallPerformanceDTUMethodSwimmFProst.pdf", height=4.5, width=8, dpi=450)
# Figure 5
dfMean<-as.data.frame(summarise(group_by(dfAll, thr, method), FDR=mean(FDR), TPR=mean(TPR), Fscore=mean(Fscore), Acc=mean(Acc)))
dfMean$satis<-dfMean$FDR<=(as.numeric(do.call(c, lapply(as.character(dfMean$thr), function(x){strsplit(x, "thr")[[1]][2]}))))
dfMean$satis[dfMean$satis]<-"yes"
dfMean$satis[dfMean$satis=="FALSE"]<-"no"
dfMean$method<-as.factor(dfMean$method)
levels(dfMean$method)<-c("DEXSeq:rT=0.01;cT=1","DEXSeq:SwimmF","DRIMSeq:rT=0.01;cT=1","DRIMSeq:SwimmF","NBSplice:rT=0.01;cT=1","NBSplice:SwimmF")
dfMean$method.satis<-paste0(dfMean$satis)
dfMean$method.satis[dfMean$method.satis=="yes"]<-as.character(dfMean$method[dfMean$method.satis=="yes"])
dfMean$method.satis<-as.factor(dfMean$method.satis)


nthr<-length(levels(dfMean$thr))
nlevs <- length(unique(dfMean$method))
plotcolors<-c("DEXSeq:rT=0.01;cT=1"="blue","DEXSeq:SwimmF"="deepskyblue3","DRIMSeq:rT=0.01;cT=1"="green4","DRIMSeq:SwimmF"="olivedrab3","NBSplice:rT=0.01;cT=1"="firebrick3","NBSplice:SwimmF"="hotpink1")
plotFillColors<-c(plotcolors, "no"="white")
thresholds<-c(0.01,0.05,0.1)
levels(dfMean$thr)<-thresholds
dfMean$satis<-as.factor(dfMean$satis)

g9<-ggplot(dfMean, aes(x=FDR, y=TPR, group=method))+geom_point(size=5,
                                                               aes(fill=method.satis,colour=method, shape=thr, ), alpha=0.7, stroke=1
)+geom_line(aes(group=method,colour=method), size=1,linetype ="dashed")+ geom_vline(xintercept= thresholds,linetype ="dotted"
)+ylim(c(0,0.75))+xlim(0,0.35)+scale_shape_manual(values = rep(21, nthr),guide = FALSE
) +scale_fill_manual(values=plotFillColors,guide=FALSE, name="")+scale_colour_manual(
    values=plotcolors, name="Methods")+theme_bw()+theme(legend.position = "bottom")+guides(col=guide_legend(ncol=3))
ggsave(g9, file="ROCDTUMethSwimmFiltProst.pdf", height = 5, width=6.2, dpi=500)


dfAux<-dfAll[dfAll$thr==0.05, ]
for(i in meth){
    for(j in meth[meth !=i]){
        print(paste(i, "vs", j, wilcox.test(dfAux$Acc[dfAux$method==i], dfAux$Acc[dfAux$method==j], exact=F, paired=T)$p.value,sep=" "))
        print(paste(i, "vs", j,wilcox.test(dfAux$TPR[dfAux$method==i], dfAux$TPR[dfAux$method==j], exact=F, paired=T)$p.value,sep=" "))
        print(paste(i, "vs", j,wilcox.test(dfAux$FDR[dfAux$method==i], dfAux$FDR[dfAux$method==j], exact=F, paired=T)$p.value,sep=" "))
        print(paste(i, "vs", j,wilcox.test(dfAux$Fscore[dfAux$method==i], dfAux$Fscore[dfAux$method==j], exact=F, paired=T)$p.value,sep=" "))
        
        
    }
}

for(i in c("DEXSeq:rT=0.01;cT=1", "DRIMSeq:rT=0.01;cT=1" )){
        print(paste("Acc: NBSplice:rT=0.01;cT=1 vs", i, wilcox.test(dfAux$Acc[dfAux$method=="NBSplice:rT=0.01;cT=1"], dfAux$Acc[dfAux$method==i], exact=F, paired=T, alternative="greater")$p.value,sep=" "))
        print(paste("Sen: NBSplice:rT=0.01;cT=1 vs", i,wilcox.test(dfAux$TPR[dfAux$method=="NBSplice:rT=0.01;cT=1"], dfAux$TPR[dfAux$method==i], exact=F, paired=T, alternative="greater")$p.value,sep=" "))
        print(paste("Pre: NBSplice:rT=0.01;cT=1 vs", i,wilcox.test(1-dfAux$FDR[dfAux$method=="NBSplice:rT=0.01;cT=1"], 1-dfAux$FDR[dfAux$method==i], exact=F, paired=T, alternative="greater")$p.value,sep=" "))
        print(paste("F: NBSplice:rT=0.01;cT=1 vs", i,wilcox.test(dfAux$Fscore[dfAux$method=="NBSplice:rT=0.01;cT=1"], dfAux$Fscore[dfAux$method==i], exact=F, paired=T, alternative="greater")$p.value,sep=" "))
}

