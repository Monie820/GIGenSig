###################################################################################################
############################ Genome instability related genes identification ######################
###################################################################################################
load("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_clin2_Binary_matrix2_Metabric_exp2.RData")
Binary_sum2=apply(Binary_matrix2,2,function(x) sum(x)) #299
min(Binary_sum2) # 0
max(Binary_sum2) # 2056
GS_pos=order(Binary_sum2)[1:75]  ##genomic stable-like group (GS)
GU_pos=order(Binary_sum2)[225:299]  ##genomic unstable-like group (GU)
GS_sam_raw=names(Binary_sum2[GS_pos])
GU_sam_raw=names(Binary_sum2[GU_pos])
pValues=c()
FCs=c()
for(i in 1:dim(Metabric_exp2)[1]){
  pValue=t.test(Metabric_exp2[i,GU_pos],Metabric_exp2[i,GS_pos])$p.value
  FC=mean(as.numeric(Metabric_exp2[i,GU_pos]))-mean(as.numeric(Metabric_exp2[i,GS_pos]))
  pValues=c(pValues,pValue)
  FCs=c(FCs,FC)
}
stats_all=cbind(rownames(Metabric_exp2),pValues,FCs)
stats_all=as.data.frame(stats_all,stringsAsFactors = F)
colnames(stats_all)=c("gene","pValue","logFC")
DEGs_2 <- stats_all[as.numeric(stats_all$pValue) < 0.05 & (as.numeric(stats_all$logFC) < -1 |as.numeric(stats_all$logFC) > 1),]
write.table(DEGs_2,"/Users/username/input/cBioPortal/brca_metabric/Analysis/DEGs_2.txt",sep="\t",col.names=T,row.names=F,quote=F)
dim(DEGs_2)
# 111 DGEs


#### heatmap analysis(add PAM50 subtype)
clin3=read.table("/Users/username/input/cBioPortal/brca_metabric/data_clinical_patient.txt",skip=4,header = TRUE, sep="\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
pos=match(Metabric_clin2$PATIENT_ID,clin3$PATIENT_ID)
Metabric_clin3=cbind(Metabric_clin2,clin3[pos,15])
library(gplots)
deg2=read.table("/Users/username/input/cBioPortal/brca_metabric/Analysis/DEGs_2_sort.txt", header=TRUE, sep="\t")
deg2_pos=match(deg2$gene,rownames(Metabric_exp2))
deg2_expValue=Metabric_exp2[deg2_pos,]
deg2_expValue=scale(deg2_expValue) #z-score normal
heat.exp.value=as.matrix(deg2_expValue)
mycol <- colorpanel(100000,"blue", "white", "red")
group<- as.character(Metabric_clin3[,8])
group[group=="Basal"]="#7FFFD4"
group[group=="claudin-low"]="#DEB887"
group[group=="Her2"]="#6495ED"
group[group=="LumA"]="#CAFF70"
group[group=="Normal"]="#FF8C00"
patientcolors=group
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/DEGs_2_heatmap.pdf",width=7, height=5)
heatmap.2(heat.exp.value,
          col=mycol, ColSideColors=patientcolors, key=FALSE, symkey=FALSE, density.info="none", 
          trace="none", labCol=NA, labRow=NA, Colv=TRUE, Rowv=TRUE, dendrogram="both")
dev.off()

#heatmap analysis(add raw top 25% GU/GS sample bar)
deg2=read.table("/Users/username/input/cBioPortal/brca_metabric/Analysis/DEGs_2_sort.txt", header=TRUE, sep="\t")
deg2_pos=match(deg2$gene,rownames(Metabric_exp2))
deg2_expValue=Metabric_exp2[deg2_pos,]
deg2_expValue=scale(deg2_expValue) #z-score normal
heat.exp.value=as.matrix(deg2_expValue)
mycol <- colorpanel(100000,"blue", "white", "red")
group<- as.character(Metabric_clin3[,1])
GS_sam_raw=names(Binary_sum2[GS_pos])
GU_sam_raw=names(Binary_sum2[GU_pos])
group[match(GS_sam_raw,rownames(Metabric_clin3))]="GS"
group[match(GU_sam_raw,rownames(Metabric_clin3))]="GU"
group[(group!="GS") & (group!="GU")]="Other"
group[group=="GU"]="red"
group[group=="GS"]="blue"
group[group=="Other"]="grey"
patientcolors=group
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/DEGs_2_heatmap_GUGS.pdf",width=7, height=5)
heatmap.2(heat.exp.value,
          col=mycol, ColSideColors=patientcolors, key=FALSE, symkey=FALSE, density.info="none", 
          trace="none", labCol=NA, labRow=NA, Colv=TRUE, Rowv=TRUE, dendrogram="both")
dev.off()




########### boxplot for mutation number and driver expression between GU and GS groups ################
library(lattice)
hc <- hclust(dist(t(heat.exp.value)))
memb <- cutree(hc, k = 2)
GU_sam=names(which(memb==1))
GS_sam=names(which(memb==2))

pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/MutNO_FOXM1_exp.pdf",width=7, height=7)
par(mfrow = c(1, 2))
GU_count=Binary_sum2[match(GU_sam,names(Binary_sum2))]
GS_count=Binary_sum2[match(GS_sam,names(Binary_sum2))]
ddf_count = data.frame(count = c(GU_count,GS_count), class = c(rep("GU-like group",length(GU_sam)),rep("GS-like group",length(GS_sam))))
boxplot(count ~ class, data = ddf_count, lwd = 0.5, ylab = 'Mutation number')
stripchart(count ~ class, vertical = TRUE, data = ddf_count, 
           method = "jitter", add = TRUE, pch = 20, col = c("blue","red"),cex=0.6)
wilcox.test(GU_count,GS_count)$p.value   ## P = 1.16e-19
text(1.5,2000,"P < 0.001")

GU_FOXM1_exp=as.numeric(Metabric_exp2[which(rownames(Metabric_exp2)=="FOXM1"),match(GU_sam,colnames(Metabric_exp2))])
GS_FOXM1_exp=as.numeric(Metabric_exp2[which(rownames(Metabric_exp2)=="FOXM1"),match(GS_sam,colnames(Metabric_exp2))])
ddf_FOXM1_exp = data.frame(exp = c(GU_FOXM1_exp,GS_FOXM1_exp), class = c(rep("GU-like group",length(GU_sam)),rep("GS-like group",length(GS_sam))))
boxplot(exp ~ class, data = ddf_FOXM1_exp, lwd = 0.5, ylab = 'FOXM1 epression')
stripchart(exp ~ class, vertical = TRUE, data = ddf_FOXM1_exp, 
           method = "jitter", add = TRUE, pch = 20, col = c("blue","red"),cex=0.6)
wilcox.test(GU_FOXM1_exp,GS_FOXM1_exp)$p.value  ## P = 2.82e-25
text(1.5,10,"P < 0.001")
dev.off()

################## function enrichment analysis by DAVID #########################
bp <- read.table("/Users/username/input/cBioPortal/brca_metabric/Analysis/DAVID/GO_BP.txt",sep="\t",header=TRUE)
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/DAVID/GO_BP.pdf")
bp <- as.matrix(bp[1:28,c(2,5)])
GO <- c()
for(i in nrow(bp):1){
  str <- strsplit(bp[i,1],"~")
  GO <- rbind(GO,c(str[[1]][2],-log10(as.numeric(bp[i,2]))))
}

plot.new()
par(new=TRUE)
par(mar=c(5,15,4,2)) 
par(mgp=c(2,0.5,0))
barplot(as.numeric(GO[,2]),width=1,space=0.4,xaxt="n",xlab="-log(P value)",names.arg=GO[,1],cex.names=0.8,
        cex.axis=0.6,horiz=TRUE,col="lightblue",offset=0,las=2,xlim=c(0,max(as.numeric(GO[,2]))+0.2))
axis(1,at=seq(0,max(as.numeric(GO[,2])),1),las=1)
#box()
dev.off()


cc <- read.table("/Users/username/input/cBioPortal/brca_metabric/Analysis/DAVID/GO_CC.txt",sep="\t",header=TRUE)
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/DAVID/GO_CC.pdf",width=7, height=9)
par(mfrow = c(2,1))
cc <- as.matrix(cc[1:16,c(2,5)])
GO <- c()
for(i in nrow(cc):1){
  str <- strsplit(cc[i,1],"~")
  GO <- rbind(GO,c(str[[1]][2],-log10(as.numeric(cc[i,2]))))
}
plot.new()
par(new=TRUE)
par(mar=c(5,17,4,2)) 
par(mgp=c(2,0.5,0))
barplot(as.numeric(GO[,2]),width=1,space=0.4,xaxt="n",xlab="-log(P value)",names.arg=GO[,1],cex.names=0.8,
        cex.axis=0.6,horiz=TRUE,col="lightblue",offset=0,las=2,xlim=c(0,max(as.numeric(GO[,2]))+0.2))
axis(1,at=seq(0,max(as.numeric(GO[,2])),1),las=1)
#box()
dev.off()


mf <- read.table("/Users/username/input/cBioPortal/brca_metabric/Analysis/DAVID/GO_MF.txt",sep="\t",header=TRUE)
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/DAVID/GO_MF.pdf",width=7, height=8.5)
par(mfrow = c(3,1))
mf <- as.matrix(mf[1:10,c(2,5)])
GO <- c()
for(i in nrow(mf):1){
  str <- strsplit(mf[i,1],"~")
  GO <- rbind(GO,c(str[[1]][2],-log10(as.numeric(mf[i,2]))))
}
plot.new()
par(new=TRUE)
par(mar=c(5,19,4,2)) 
par(mgp=c(2,0.5,0))
barplot(as.numeric(GO[,2]),width=1,space=0.4,xaxt="n",xlab="-log(P value)",names.arg=GO[,1],cex.names=0.8,
        cex.axis=0.6,horiz=TRUE,col="lightblue",offset=0,las=2,xlim=c(0,max(as.numeric(GO[,2]))+0.2))
axis(1,at=seq(0,max(as.numeric(GO[,2])),1),las=1)
#box()
dev.off()


################## function enrichment analysis by clusterProfiler #########################
####ID transfer
library(stringr)
library(org.Hs.eg.db)
x <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(x)
# Convert to a list
ID2symbol <- as.list(x[mapped_genes])
# For the reverse map:
x <- org.Hs.egSYMBOL2EG
# Get the entrez gene identifiers that are mapped to a gene symbol
mapped_genes <- mappedkeys(x)
# Convert to a list
symbol2ID <- as.list(x[mapped_genes])

dotplot_ylab <- function(x,width,top){
  x.df <- as.data.frame(x)
  if(nrow(x.df)>top){
    x.df <- x.df[1:top,]
  }
  x.df$GeneRatio <-  unlist(lapply(as.list(x.df$GeneRatio), function(x) eval(parse(text=x)))) 
  x.df$Description <- factor(x.df$Description,levels = x.df$Description[order(x.df$Count,decreasing = F)])
  #x.df <- x.df[order(x.df$Count,decreasing = F),]
  p<-ggplot(x.df, aes(x=GeneRatio, y=Description,color=p.adjust)) + 
    geom_point(aes(size = Count))+scale_color_gradient(low="red", high="blue")+scale_y_discrete(labels=function(x) str_wrap(x,width=width))
}

library(gplots)
library(ggplot2)
library(clusterProfiler)
outpath="/Users/username/input/cBioPortal/brca_metabric/Analysis"
deg2=read.table("/Users/username/input/cBioPortal/brca_metabric/Analysis/DEGs_2_sort.txt", header=TRUE, sep="\t")
tmp.geneset <- as.vector(deg2$gene)
tmp.geneset2id <- as.vector(unlist((sapply(tmp.geneset,function(x){symbol2ID[[x]]}))))

for(goCate in c("BP","CC","MF")){
  tmp.enrichGOTerms <- enrichGO(tmp.geneset2id,OrgDb="org.Hs.eg.db",ont=goCate)
  tmp.enrichGOTerms.df <- as.data.frame(tmp.enrichGOTerms)
  if(nrow(tmp.enrichGOTerms.df)>0){
    gsea.folder=paste(outpath,"/","gsea/",sep="")
    gsea.fig.folder=paste(outpath,"/","gsea/figure/",sep="")
    if(!dir.exists(gsea.folder)){dir.create(gsea.folder,recursive=T)}
    if(!dir.exists(gsea.fig.folder)){dir.create(gsea.fig.folder,recursive=T)}
    write.table(tmp.enrichGOTerms.df,paste(gsea.folder,"_GO_",goCate,".txt",sep=""),col.names = T,sep="\t",row.names = F,quote = F)
    pdf(paste(gsea.fig.folder,"_GO_",goCate,".pdf",sep=""))
    #print(DOSE::dotplot(tmp.enrichGOTerms,font.size=6))
    print(dotplot_ylab(tmp.enrichGOTerms,width = 50,top = 10))
    dev.off()
  }
}
tmp.enrichKEGG <- enrichKEGG(tmp.geneset2id,organism = "hsa")
tmp.enrichKEGG.df <- as.data.frame(tmp.enrichKEGG)
if(nrow(tmp.enrichKEGG.df) > 0){
  gsea.folder=paste(outpath,"/","gsea/",sep="")
  gsea.fig.folder=paste(outpath,"/","gsea/figure/",sep="")
  if(!dir.exists(gsea.folder)){dir.create(gsea.folder,recursive=T)}
  if(!dir.exists(gsea.fig.folder)){dir.create(gsea.fig.folder,recursive=T)}
  write.table(tmp.enrichKEGG.df,paste(gsea.folder,"_KEGG.txt",sep=""),col.names = T,sep="\t",row.names = F,quote = F)
  pdf(paste(gsea.fig.folder,"_KEGG.pdf",sep=""))
  #print(DOSE::dotplot(tmp.enrichKEGG,font.size=6))
  print(dotplot_ylab(tmp.enrichKEGG,width = 50,top = 10))
  dev.off()
}
