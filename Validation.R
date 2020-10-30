############################# #####################################################
############################# Shanghai TNBC validation analysis  ##################
############################# #####################################################
# clinical
clin=read.table("/Users/username/input/Shanghai_TNBC/CNV/GSE118527_series_matrix.txt",skip = 25)
cnv_sam=as.character(t(clin[11,-1])) #424
tumor_sam=as.character(t(clin[11,which(clin[12,]=="tissue: Primary triple negative breast cancer")])) #401
blood_sam=as.character(t(clin[11,which(clin[12,]=="tissue: Peripheral white blood cell")])) #23
germline_sam=cnv_sam[match(blood_sam,cnv_sam)+1]

sra=read.table("/Users/username/input/Shanghai_TNBC/SRRlist/SraRunTable.txt",header=TRUE,sep="\t")
wes_sam=sra[which(sra[,1]=="WXS"),c(12,16)] #279
rna_sam=sra[which(sra[,1]=="RNA-Seq"),c(12,16)] #448
inter_sam=intersect(intersect(wes_sam[,2],rna_sam[,2]),cnv_sam) #235

wes_interpos=match(inter_sam,wes_sam[,2])
wes_inter=wes_sam[wes_interpos,]#235
wes_inter=wes_inter[order(wes_inter[,1]),]
write.table(wes_inter,"/Users/username/input/Shanghai_TNBC/SRRlist/WES_inter_srrFUS.txt",sep="\t",col.names=F,row.names=F,quote=F)

rna_interpos=match(inter_sam,rna_sam[,2])
rna_inter=rna_sam[rna_interpos,]#235
rna_inter=rna_inter[order(rna_inter[,1]),]
write.table(rna_inter,"/Users/username/input/Shanghai_TNBC/SRRlist/RNA_inter_srrFUS.txt",sep="\t",col.names=F,row.names=F,quote=F)

clin=read.table("/Users/username/input/Shanghai_TNBC/Clinical/Sample_information.txt",header = TRUE, sep="\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
Shanghai_clin=clin[match(rna_inter[,2],clin[,1]),c(1,38,37,7,31)]
write.table(Shanghai_clin,"/Users/username/input/Shanghai_TNBC/Shanghai_clin.txt",sep="\t", row.names=F, col.names=T, quote=F)


# expression
library(edgeR)
rnaseqv2_file=read.delim("/Users/username/input/Shanghai_TNBC/RNASeq/gene_count_matrix_Shanghai.txt",row.names=1)
PCG_gene_Anno=All_gene_Anno[All_gene_Anno$V11=="protein_coding",]
gene_pos=na.omit(match(rownames(rnaseqv2_file),PCG_gene_Anno[,10]))
gene_pos2=match(PCG_gene_Anno[gene_pos,10],rownames(rnaseqv2_file))
rnaseqv2_matrix=rnaseqv2_file[gene_pos2,]  # 20327   235
gene_length=PCG_gene_Anno[gene_pos,5]-PCG_gene_Anno[gene_pos,4]
rnaseqv2_matrix_rpkm=rpkm(y = rnaseqv2_matrix, gene.length = gene_length)
colnames(rnaseqv2_matrix_rpkm)=Shanghai_clin[,1] ## 20327   235
rownames(rnaseqv2_matrix_rpkm)=PCG_gene_Anno[gene_pos,9]
write.table(rnaseqv2_matrix_rpkm,"/Users/username/input/Shanghai_TNBC/Shanghai_expression.txt",sep="\t", row.names=T, col.names=T, quote=F)


################# Shanghai risk score difference for different stages ###############
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/Shanghai_scores_forGrade.pdf")
Shanghai_surv_score_2=Shanghai_surv_score[-which(is.na(Shanghai_surv_score$Grade)==TRUE),]
Shanghai_surv_score_2[Shanghai_surv_score_2$Grade==2,5]="G2"
Shanghai_surv_score_2[Shanghai_surv_score_2$Grade=="2 to 3",5]="G2 to G3"
Shanghai_surv_score_2[Shanghai_surv_score_2$Grade==3,5]="G3"
del_pos1= which(Shanghai_surv_score_2$scores > 10)
del_pos2= which(Shanghai_surv_score_2$scores < -8)
Shanghai_surv_score_2=Shanghai_surv_score_2[-c(del_pos1,del_pos2),]
ddf_count = data.frame(exp = -(Shanghai_surv_score_2$scores), class = Shanghai_surv_score_2$Grade)
boxplot(exp ~ class, data = ddf_count, lwd = 0.5, ylab = 'Overall risk score')
stripchart(exp ~ class, vertical = TRUE, data = ddf_count, 
           method = "jitter", add = TRUE, pch = 20, col = c("blue","orange","red"),cex=0.7)
text(1.5,12,"P = 0.001")
wilcox.test(exp ~ class, data = ddf_count)
wilcox.test(exp ~ class, data = ddf_count,subset = class %in% c("G2", "G2 to G3"))
wilcox.test(exp ~ class, data = ddf_count,subset = class %in% c("G2", "G3"))
wilcox.test(exp ~ class, data = ddf_count,subset = class %in% c("G2 to G3", "G3"))
dev.off()

############################# #####################################################
############################# GSE21653 TNBC validation analysis  ##################
############################# #####################################################
clin = read.table("/Users/username/input/cBioPortal/brca_metabric/Analysis/GEO_validation/GSE21653_clin.txt",header = TRUE, sep="\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
clinFeatureMatrix = clin[,c(1,13,12,6,7,8,9,10)]
colnames(clinFeatureMatrix) = c("sample","dfs_mouth","dfs_status","grade","er","pr","erbb2","p53")
featureMatrix = clinFeatureMatrix[,5:7]
tnbcTag = apply(featureMatrix,1,function(x) all(x==0))
tnbc_sam = clinFeatureMatrix$sample[which(tnbcTag =="TRUE")]
length(tnbc_sam) #87
clin_tnbc=clinFeatureMatrix[match(tnbc_sam,clinFeatureMatrix$sample),]

exp=read.table("/Users/username/input/cBioPortal/brca_metabric/Analysis/GEO_validation/GSE21653_series_matrix.txt",skip = 90,header=TRUE,sep="\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
probe_anno=read.table("/Users/username/input/cBioPortal/brca_metabric/Analysis/GEO_validation/GPL570-55999.txt",skip = 16,header=TRUE,sep="\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
probe_gene=probe_anno[,c(1,11)]
gene_names=apply(probe_gene,1,function(x) unlist(strsplit(x[2]," /// "))[1])
probe_gene2=data.frame(probe_gene[,1],gene_names)
colnames(probe_gene2)=c("probe","gene")
exp[,1]=probe_gene2[match(exp$ID_REF,probe_gene2$probe),2]
exp2=aggregate(x = exp[,-1], by = list(exp$ID_REF), FUN = "mean")
colnames(exp2)[1]=c("Gene")
HG_U133_Plus_2_exp=exp2[,-1]
rownames(HG_U133_Plus_2_exp)=exp2$Gene

############################# ##############################################
############################# GSE25066 TNBC validation analysis  ###########
############################# ##############################################
clin = read.table("/Users/username/input/cBioPortal/brca_metabric/Analysis/GEO_validation/GSE25066_clin.txt",header = TRUE, sep="\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
clinFeatureMatrix = clin[,c(1,15,14,10,11,5,6,7)]
colnames(clinFeatureMatrix) = c("sample","dfs_days","dfs_status","stage","grade","er","pr","erbb2")
featureMatrix = clinFeatureMatrix[,6:8]
tnbcTag = apply(featureMatrix,1,function(x) all(x == "N"))
tnbc_sam = clinFeatureMatrix$sample[which(tnbcTag =="TRUE")]
length(tnbc_sam) #178
clin_tnbc=clinFeatureMatrix[match(tnbc_sam,clinFeatureMatrix$sample),]

exp=read.table("/Users/username/input/cBioPortal/brca_metabric/Analysis/GEO_validation/GSE25066_series_matrix.txt",skip = 87,header=TRUE,sep="\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
probe_anno=read.table("/Users/username/input/cBioPortal/brca_metabric/Analysis/GEO_validation/GPL96-57554.txt",skip = 16,header=TRUE,sep="\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
probe_gene=probe_anno[,c(1,11)]
gene_names=apply(probe_gene,1,function(x) unlist(strsplit(x[2]," /// "))[1])
probe_gene2=data.frame(probe_gene[,1],gene_names)
colnames(probe_gene2)=c("probe","gene")
exp[,1]=probe_gene2[match(exp$ID_REF,probe_gene2$probe),2]
exp2=aggregate(x = exp[,-1], by = list(exp$ID_REF), FUN = "mean")
colnames(exp2)[1]=c("Gene")
HG_U133A_exp=exp2[,-1]
rownames(HG_U133A_exp)=exp2$Gene

###############################################################################
####### boxplot for TFF3 expression in different grade group (G1/2, G3) #######
###############################################################################
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/GSE21653andGSE31448andGSE25066_TFF3_exp_forGrade.pdf",width=9, height=4)
par(mfrow = c(1, 3))
GSE21653_surv_2=GSE21653_surv[-which(is.na(GSE21653_surv$grade)==TRUE),]
GSE21653_surv_2[GSE21653_surv_2$grade==1,4]="G1/2"
GSE21653_surv_2[GSE21653_surv_2$grade==2,4]="G1/2"
GSE21653_surv_2[GSE21653_surv_2$grade==3,4]="G3"
ddf_count = data.frame(exp = GSE21653_surv_2$TFF3, class = GSE21653_surv_2$grade)
boxplot(exp ~ class, data = ddf_count, lwd = 0.5, ylab = 'TFF3 expression')
stripchart(exp ~ class, vertical = TRUE, data = ddf_count, 
           method = "jitter", add = TRUE, pch = 20, col = c("red","blue"),cex=0.7)
text(1.5,12,"P = 0.003")
wilcox.test(exp ~ class, data = ddf_count)

GSE31448_surv_2=GSE31448_surv[-which(is.na(GSE31448_surv$grade)==TRUE),]
GSE31448_surv_2[GSE31448_surv_2$grade==1,4]="G1/2"
GSE31448_surv_2[GSE31448_surv_2$grade==2,4]="G1/2"
GSE31448_surv_2[GSE31448_surv_2$grade==3,4]="G3"
ddf_count = data.frame(exp = GSE31448_surv_2$TFF3, class = GSE31448_surv_2$grade)
boxplot(exp ~ class, data = ddf_count, lwd = 0.5, ylab = 'TFF3 expression')
stripchart(exp ~ class, vertical = TRUE, data = ddf_count, 
           method = "jitter", add = TRUE, pch = 20, col = c("red","blue"),cex=0.7)
text(1.5,12,"P = 0.004")
wilcox.test(exp ~ class, data = ddf_count)

GSE25066_surv_2=GSE25066_surv[-c(which(is.na(GSE25066_surv$grade)==TRUE),which(GSE25066_surv$grade=="Indeterminate")),]
GSE25066_surv_2[GSE25066_surv_2$grade==1,5]="G1/2"
GSE25066_surv_2[GSE25066_surv_2$grade==2,5]="G1/2"
GSE25066_surv_2[GSE25066_surv_2$grade==3,5]="G3"
ddf_count = data.frame(exp = GSE25066_surv_2$TFF3, class = GSE25066_surv_2$grade)
boxplot(exp ~ class, data = ddf_count, lwd = 0.5, ylab = 'TFF3 expression')
stripchart(exp ~ class, vertical = TRUE, data = ddf_count, 
           method = "jitter", add = TRUE, pch = 20, col = c("red","blue"),cex=0.7)
text(1.5,12,"P = 0.019")
wilcox.test(exp ~ class, data = ddf_count)
dev.off()

#################################################################################
##### boxplot for FOXM1 expression according to high and low  TFF3 epxression ####
#################################################################################
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/GSE21653andGSE31448andGSE25066_FOXM1_expByTFF3.pdf",width=9, height=4)
par(mfrow = c(1, 3))
GSE21653_surv_2=GSE21653_surv
FOXM1_class=rep(0,length(GSE21653_surv_2$TFF3))
highExp_pos=which(GSE21653_surv_2$TFF3 >= median(GSE21653_surv_2$TFF3))
FOXM1_class[highExp_pos]=2
FOXM1_class[-highExp_pos]=1
ddf_count = data.frame(exp = GSE21653_surv_2$FOXM1, class = FOXM1_class)
boxplot(exp ~ class, data = ddf_count, lwd = 0.5, ylab = 'FOXM1 expression')
stripchart(exp ~ class, vertical = TRUE, data = ddf_count, 
           method = "jitter", add = TRUE, pch = 20, col = c("red","blue"),cex=0.7)
text(1.5,11,"P < 0.001")
wilcox.test(exp ~ class, data = ddf_count)

GSE31448_surv_2=GSE31448_surv
FOXM1_class=rep(0,length(GSE31448_surv_2$TFF3))
highExp_pos=which(GSE31448_surv_2$TFF3 >= median(GSE31448_surv_2$TFF3))
FOXM1_class[highExp_pos]=2
FOXM1_class[-highExp_pos]=1
ddf_count = data.frame(exp = GSE31448_surv_2$FOXM1, class = FOXM1_class)
boxplot(exp ~ class, data = ddf_count, lwd = 0.5, ylab = 'FOXM1 expression')
stripchart(exp ~ class, vertical = TRUE, data = ddf_count, 
           method = "jitter", add = TRUE, pch = 20, col = c("red","blue"),cex=0.7)
text(1.5,11,"P < 0.001")
wilcox.test(exp ~ class, data = ddf_count)

GSE25066_surv_2=GSE25066_surv
FOXM1_class=rep(0,length(GSE25066_surv_2$TFF3))
highExp_pos=which(GSE25066_surv_2$TFF3 >= median(GSE25066_surv_2$TFF3))
FOXM1_class[highExp_pos]=2
FOXM1_class[-highExp_pos]=1
ddf_count = data.frame(exp = GSE25066_surv_2$FOXM1, class = FOXM1_class)
boxplot(exp ~ class, data = ddf_count, lwd = 0.5, ylab = 'FOXM1 expression')
stripchart(exp ~ class, vertical = TRUE, data = ddf_count, 
           method = "jitter", add = TRUE, pch = 20, col = c("red","blue"),cex=0.7)
text(1.5,11,"P = 0.003")
wilcox.test(exp ~ class, data = ddf_count)
dev.off()

