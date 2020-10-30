
###########################################################################################
################################### Starting data processing analysis #####################
###########################################################################################
library(readr)

geneAnnoAll=read_delim("/Users/username/input/GENCODE/gencode.v19.annotation.gtf",delim="\t",skip=5,col_names=FALSE)
geneAnnoAll=as.data.frame(geneAnnoAll)

## geneAnno
geneAnno=geneAnnoAll[geneAnnoAll[,3]=="gene",]
geneList=apply(geneAnno,1,function(x) unlist(strsplit(x[9],"\""))[10])
geneType=apply(geneAnno,1,function(x) unlist(strsplit(x[9],"\""))[6])
PCG_Anno=data.frame(geneAnno[which(geneType=="protein_coding"),1:8],geneList[which(geneType=="protein_coding")],geneType[which(geneType=="protein_coding")])
write.table(PCG_Anno,"/Users/username/input/GENCODE/PCG_Anno.txt",sep="\t",col.names=F,row.names=F,quote=F)

geneList_name=apply(geneAnno,1,function(x) unlist(strsplit(x[9],"\""))[10])
geneList_ensg=apply(geneAnno,1,function(x) unlist(strsplit(x[9],"\""))[2])
geneList_ensg2=apply(as.matrix(geneList_ensg),1,function(x) unlist(strsplit(x[1],"[.]"))[1])
geneType_all=apply(geneAnno,1,function(x) unlist(strsplit(x[9],"\""))[6])
All_gene_Anno=data.frame(geneAnno[,1:8],geneList_name,geneList_ensg2,geneType_all)
write.table(All_gene_Anno,"/Users/username/input/GENCODE/All_gene_Anno.txt",sep="\t",col.names=F,row.names=F,quote=F)

## transcriptAnno
transcriptAnno=geneAnnoAll[geneAnnoAll[,3]=="transcript",]
transcriptType=apply(geneAnno,1,function(x) unlist(strsplit(x[9],"\""))[6])
geneList=apply(transcriptAnno,1,function(x) unlist(strsplit(x[9],"\""))[10])
transcriptList=apply(transcriptAnno,1,function(x) unlist(strsplit(x[9],"\""))[4])
transcript_Anno=data.frame(transcriptAnno[which(transcriptType=="protein_coding"),1:8],geneList[which(transcriptType=="protein_coding")],transcriptType[which(transcriptType=="protein_coding")],transcriptList[which(transcriptType=="protein_coding")])
write.table(transcript_Anno,"/Users/username/input/GENCODE/transcript_Anno.txt",sep="\t",col.names=F,row.names=F,quote=F)


###########################################################################################
################################### cBioPortal METABRIC TNBC samples ######################
###########################################################################################
# clinical  (N=2509)
library(readr)
clin=read.table("/Users/username/input/cBioPortal/brca_metabric/data_clinical_sample.txt",skip=4,header = TRUE, sep="\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
featureRow=c("ER_STATUS","HER2_STATUS","PR_STATUS")
featureMatrix=clin[,match(featureRow,colnames(clin))]
clin2=read.table("/Users/username/input/cBioPortal/brca_metabric/data_clinical_sample.txt",skip=4,header = TRUE, sep="\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
featureRow2=c("TUMOR_STAGE","GRADE")
featureMatrix2=clin2[,match(featureRow2,colnames(clin2))]
clin3=read.table("/Users/username/input/cBioPortal/brca_metabric/data_clinical_patient.txt",skip=4,header = TRUE, sep="\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
featureRow3=c("PATIENT_ID","OS_MONTHS","OS_STATUS","AGE_AT_DIAGNOSIS","INFERRED_MENOPAUSAL_STATE")
featureMatrix3=clin3[,match(featureRow3,colnames(clin3))]

# tnbc samples (N=320(299 have expr))
tnbcTag=apply(featureMatrix,1,function(x) all(x=="Negative"))
tnbc_pos=which(tnbcTag=="TRUE")
tnbc_sam=clin[tnbc_pos,1] ##320
Metabric_clin=cbind(featureMatrix3[tnbc_pos,],featureMatrix2[tnbc_pos,])

# expression
rnaseqv2_file=read.delim("/Users/username/input/cBioPortal/brca_metabric/Analysis/data_expression_median.txt",row.names=1) #24368  1905
sample_pos=na.omit(match(tnbc_sam,colnames(rnaseqv2_file))) 
PCG_pos=na.omit(match(PCG_Anno[,9],rownames(rnaseqv2_file)))
length(sample_pos) ##299
Metabric_exp=rnaseqv2_file[PCG_pos,sample_pos] ##  16331   299
write.table(Metabric_exp,"/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_exp.txt",sep="\t",col.names=T,row.names=T,quote=F)

tnbc_pos2=match(colnames(Metabric_exp),tnbc_sam)
Metabric_clin2=Metabric_clin[tnbc_pos2,]
write.table(Metabric_clin2,"/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_clin.txt",sep="\t",col.names=T,row.names=F,quote=F)



# CNV
cnv_file=read.delim("/Users/username/input/cBioPortal/brca_metabric/Analysis/data_CNA.txt",row.names=1)
sample_pos2=na.omit(match(colnames(Metabric_exp),colnames(cnv_file)))
length(sample_pos2) ##299
Metabric_cnv=cnv_file[,sample_pos2]
Metabric_cnv2=Metabric_cnv
Metabric_cnv2[Metabric_cnv==1]=0
Metabric_cnv2[Metabric_cnv==-1]=0
Metabric_cnv2[Metabric_cnv==2]=1
Metabric_cnv2[Metabric_cnv==-2]=1
Metabric_cnv2[is.na(Metabric_cnv)]=0  ### 22544   299
write.table(Metabric_cnv2,"/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_cnv.txt",sep="\t",col.names=T,row.names=T,quote=F)

# mutation
mut_file=read.delim("/Users/username/input/cBioPortal/brca_metabric/Analysis/data_mutations_mskcc.txt",header=TRUE,sep="\t",fill=TRUE)
sample_pos3=unlist(apply(as.matrix(colnames(Metabric_exp)),1,function(x) grep(x[1],mut_file[,17])))
Metabric_mut=mut_file[sample_pos3,] ## 2284   45
write.table(Metabric_mut,"/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_mut.txt",sep="\t",col.names=T,row.names=F,quote=F)
sample_pos4=unique(na.omit(match(Metabric_mut[,1],rownames(Metabric_cnv2))))
sample_pos5=unlist(apply(as.matrix(rownames(Metabric_cnv2)[sample_pos4]),1,function(x) which(Metabric_mut[,1]==x[1])))
Metabric_mut2=Metabric_mut[sample_pos5,] ## 2234   45

# Binary matrix of somatic genetic alteration
for(i in 1:dim(Metabric_mut2)[1]){
  row_pos=match(Metabric_mut2[i,1],rownames(Metabric_cnv2))
  col_pos=match(Metabric_mut2[i,17],colnames(Metabric_cnv2))
  Metabric_cnv2[row_pos,col_pos]=1
}
Binary_matrix=Metabric_cnv2  ##22544   299
write.table(Binary_matrix,"/Users/username/input/cBioPortal/brca_metabric/Analysis/Binary_matrix.txt",sep="\t",col.names=T,row.names=T,quote=F)

gene_inter=intersect(rownames(Binary_matrix),rownames(Metabric_exp))
Binary_matrix2=Binary_matrix[match(gene_inter,rownames(Binary_matrix)),] # 15872   299
Metabric_exp2=Metabric_exp[match(gene_inter,rownames(Metabric_exp)),] # 15872   299

# training and test analysis 
Metabric_clin2=read.table("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_clin.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
Metabric_clin2[Metabric_clin2[,3]=="LIVING",3]=0
Metabric_clin2[Metabric_clin2[,3]=="DECEASED",3]=1
rownames(Metabric_clin2)=colnames(Metabric_exp2)

Metabric_clin2_0=Metabric_clin2[Metabric_clin2$OS_STATUS==0,]
Metabric_clin2_1=Metabric_clin2[Metabric_clin2$OS_STATUS==1,]
train_sam0=sample(rownames(Metabric_clin2_0),69)
train_sam1=sample(rownames(Metabric_clin2_1),81)
train_sam=c(train_sam0,train_sam1)
Metabric_sam_train_pos=match(train_sam,rownames(Metabric_clin2))

Metabric_clin2_train=Metabric_clin2[Metabric_sam_train_pos,]
Metabric_clin2_test=Metabric_clin2[-Metabric_sam_train_pos,]

Metabric_exp2_train=Metabric_exp2[,Metabric_sam_train_pos]
Metabric_exp2_test=Metabric_exp2[,-Metabric_sam_train_pos]

Binary_matrix2_train=Binary_matrix2[,Metabric_sam_train_pos]
Binary_matrix2_test=Binary_matrix2[,-Metabric_sam_train_pos]

save(Metabric_clin2,Binary_matrix2,Metabric_exp2,file="/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_clin2_Binary_matrix2_Metabric_exp2.RData")
load(file="/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_clin2_Binary_matrix2_Metabric_exp2.RData")

save(Metabric_exp2_train,Metabric_exp2_test,Binary_matrix2_train,Binary_matrix2_test,Metabric_clin2_train,Metabric_clin2_test,file="/Users/username/input/cBioPortal/brca_metabric/Analysis/Expression_mut_clin2.RData")
load(file="/Users/username/input/cBioPortal/brca_metabric/Analysis/Expression_mut_clin2.RData")
