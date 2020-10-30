###################################################################################################
################ Genome instability-derived gene signature (GIGenSig) identification ##############
###################################################################################################
#Univariate cox regression analysis: identify survival associated DEG from 112 DEGs
load(file="/Users/maoni/Documents/TNBC2/input/cBioPortal/brca_metabric/Analysis/Expression_mut_clin2.RData")
deg2=read.table("/Users/username/input/cBioPortal/brca_metabric/Analysis/DEGs_2.txt", header=TRUE, sep="\t")
Metabric_exp_111_train=Metabric_exp2_train[match(DEGs_2$gene,rownames(Metabric_exp2_train)),]
Metabric_exp_111_2_train=t(Metabric_exp_111_train)
Metabric_surv2_train=cbind(Metabric_clin2_train,Metabric_exp_111_2_train)
library(splines)
library(survival)
Metabric_uni_stats=c()
for(i in 8:dim(Metabric_surv2_train)[2]){
  surv.stat=summary(coxph(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ Metabric_surv2_train[,i], data = Metabric_surv2_train))
  Coefficient=surv.stat[[7]][1]
  HR=surv.stat[[7]][2]
  pValue=surv.stat[[7]][5]
  lower_CI=surv.stat[[8]][3]
  upper_CI=surv.stat[[8]][4]
  Metabric_stat=c(Coefficient,HR,lower_CI,upper_CI,pValue)
  Metabric_uni_stats=rbind(Metabric_uni_stats,Metabric_stat)
}
rownames(Metabric_uni_stats)=colnames(Metabric_exp_111_2_train)
sig_pos2=which(Metabric_uni_stats[,5]<0.05)
DEG_surv=Metabric_uni_stats[sig_pos2,]
DEG_surv_Anno=All_gene_Anno[match(rownames(DEG_surv),All_gene_Anno[,9]),]
DEG_surv2_111_train=cbind(DEG_surv_Anno[,c(9,1,4,5)],DEG_surv)
colnames(DEG_surv2_111_train)=c("gene","chr","start","end","Coefficient","HR","lower_CI","upper_CI","pValue")
write.table(DEG_surv2_111_train,"/Users/username/input/cBioPortal/brca_metabric/Analysis/DEG_surv2_uni_stats_train.txt",sep="\t",row.names=F,col.names=T,quote=F)
dim(DEG_surv2_111_train) # 11


#Multivariate cox regression analysis
DEG_surv2_pos=match(DEG_surv2_111_train$gene,colnames(Metabric_surv2_train))
Metabric_surv3_train=Metabric_surv2_train[,c(1:7,DEG_surv2_pos)]
library(splines)
library(survival)
Metabric_mul2_stats2_train=c()
surv.stat2=summary(coxph(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ 
                           ART3 + CD52 + CD79A + FZD9 + GABRP + IRF8 + ITM2A + PRKCB + SOX10 + TFF3 + VGLL1,
                         data = Metabric_surv3_train))
mul_stats2=surv.stat2[[7]]
Metabric_mul2_stats2_train=mul_stats2
write.table(Metabric_mul2_stats2_train,"/Users/username/input/cBioPortal/brca_metabric/Analysis/DEG_surv2_mul_stats2_train.txt",sep="\t",row.names=T,col.names=T,quote=F)


#survival analysis (train plot mul2 Coeff surv curve)
Metabric_exp_111_train=Metabric_exp2_train[match(DEGs_2$gene,rownames(Metabric_exp2_train)),]
Metabric_exp_111_2_train=t(Metabric_exp_111_train) #150 111
scores=apply(Metabric_exp_111_2_train,1,function(x) sum(Metabric_mul2_stats2_train[,1]*x[sig_pos2]))
classify=rep(0,length(scores))
pos1=which(scores < median(scores))
pos2=which(scores >= median(scores))
median_scores = median(scores)  ## median_scores = -2.321193
classify[pos1]=1
classify[pos2]=2
Metabric_surv_score_train=cbind(Metabric_clin2_train,scores,classify)
surv <- survfit(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ classify, data = Metabric_surv_score_train)
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_surv2_mul2_train.pdf")
plot(surv,mark.time=TRUE,mark=3,lty=1,lwd=3,col=3:2,xlab="Years",ylab="Survival",xscale=12)
pp<-survdiff(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ classify, data = Metabric_surv_score_train)
p.val <- 1 - pchisq(pp$chisq, length(pp$n) - 1) 
p.val
text(240,0.7,"Log-rank p = 2.66e-4")
dev.off()

#survival analysis (test plot mul2 surv curve)
Metabric_exp_111_test=Metabric_exp2_test[match(DEGs_2$gene,rownames(Metabric_exp2_test)),]
Metabric_exp_111_2_test=t(Metabric_exp_111_test)
scores=apply(Metabric_exp_111_2_test,1,function(x) sum(Metabric_mul2_stats2_train[,1]*x[sig_pos2]))
classify=rep(0,length(scores))
pos1=which(scores < median_scores)
pos2=which(scores >= median_scores)
classify[pos1]=1
classify[pos2]=2
Metabric_surv_score_test=cbind(Metabric_clin2_test,scores,classify)
surv <- survfit(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ classify, data = Metabric_surv_score_test)
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_surv2_mul2_test.pdf")
plot(surv,mark.time=TRUE,mark=3,lty=1,lwd=3,col=3:2,xlab="Years",ylab="Survival",xscale=12)
pp<-survdiff(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ classify, data = Metabric_surv_score_test)
p.val <- 1 - pchisq(pp$chisq, length(pp$n) - 1) 
p.val
text(240,0.7,"Log-rank p = 2.45e-2")
dev.off()


#survival analysis (all cBioPortal plot mul2 surv curve)
Metabric_exp_111=Metabric_exp2[match(DEGs_2$gene,rownames(Metabric_exp2)),]
Metabric_exp_111_2=t(Metabric_exp_111)
scores=apply(Metabric_exp_111_2,1,function(x) sum(Metabric_mul2_stats2_train[,1]*x[sig_pos2]))
classify=rep(0,length(scores))
pos1=which(scores < median_scores)
pos2=which(scores >= median_scores)
classify[pos1]=1
classify[pos2]=2
Metabric_surv_score=cbind(Metabric_clin2,scores,classify)
surv <- survfit(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ classify, data = Metabric_surv_score)
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_surv2_mul2_total.pdf")
plot(surv,mark.time=TRUE,mark=3,lty=1,lwd=3,col=3:2,xlab="Years",ylab="Survival",xscale=12)
pp<-survdiff(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ classify, data = Metabric_surv_score)
p.val <- 1 - pchisq(pp$chisq, length(pp$n) - 1) 
p.val
text(240,0.7,"Log-rank p = 2.57e-5") #significant
dev.off()

save(Metabric_exp_111_2_train,Metabric_exp_111_2_test,Metabric_exp_111_2,file="/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_exp_111_2_train_test_total.RData")
load("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_exp_111_2_train_test_total.RData")

save(Metabric_surv_score_train,Metabric_surv_score_test,Metabric_surv_score,file="/Users/username/input/cBioPortal/brca_metabric/Analysis/Train_test_total_scores.RData")
load("/Users/username/input/cBioPortal/brca_metabric/Analysis/Train_test_total_scores.RData")
