#################### ROC comparing analysis with another three signatures ####################
load("/Users/username/input/cBioPortal/brca_metabric/Analysis/Train_test_total_scores.RData")
Mansour_signature = c("ACSM4","SPDYC")
Wang_signature = c("EOMES", "FA2H", "GSPT1", "RASGRP1", "SOD2")
Kim_signature = c("CEBPD", "MMP20", "WLS", "ASF1A", "ASPSCR1", "CHAF1B", "DNMT1", "GINS2") 

library(splines)
library(survival)
Metabric_exp_112=Metabric_exp2[match(Mansour_signature,rownames(Metabric_exp2)),]
Metabric_exp_112_2=t(Metabric_exp_112)
Metabric_surv2=cbind(Metabric_clin2,Metabric_exp_112_2)
Metabric_mul2_stats2=c()
surv.stat2=summary(coxph(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ CEBPD + ASPSCR1 + CHAF1B, data = Metabric_surv2))
mul_stats2=surv.stat2[[7]]
Metabric_mul2_stats2=mul_stats2
scores=apply(Metabric_exp_112_2,1,function(x) sum(Metabric_mul2_stats2[,1]*x))
classify=rep(0,length(scores))
pos1=which(scores < median(scores))
pos2=which(scores >= median(scores))
classify[pos1]=1
classify[pos2]=2
Metabric_surv_score_Mansour=cbind(Metabric_clin2,scores,classify)
surv <- survfit(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ classify, data = Metabric_surv_score_Mansour)
#pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_surv2_mul2_total.pdf")
plot(surv,mark.time=TRUE,mark=3,lty=1,lwd=3,col=3:2,xlab="Years",ylab="Survival",xscale=12)
pp<-survdiff(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ classify, data = Metabric_surv_score_Mansour)
p.val <- 1 - pchisq(pp$chisq, length(pp$n) - 1) 
p.val
text(240,0.7,"Log-rank p = 0.009") #significant
#dev.off()

# ROC combine plot
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/ROC_compare_at5years.pdf")

total=Metabric_surv_score
nobs <- NROW(total)
cutoff <- 12*5
totalscore= survivalROC(Stime=total$OS_MONTHS, status=total$OS_STATUS, marker = total$scores, predict.time = cutoff,span = 0.25*nobs^(-0.20))
plot(totalscore$FP, totalscore$TP, type="l",col="red",lwd=2, xlim=c(0,1), ylim=c(0,1),xlab="1-Specificity",ylab="Sensitivity")
text(0.6,0.4,paste("GIGenSig ROC-AUC at 5 years = ",round(totalscore$AUC,3)))
abline(0,1)

par(new=TRUE)
Wang=Metabric_surv_score_Wang
nobs <- NROW(Wang)
cutoff <- 12*5
Wangscore= survivalROC(Stime=Wang$OS_MONTHS, status=Wang$OS_STATUS, marker = Wang$scores, predict.time = cutoff,span = 0.25*nobs^(-0.20))
plot(Wangscore$FP, Wangscore$TP, type="l",col="blue",lwd=2, xlim=c(0,1), ylim=c(0,1),xlab="1-Specificity",ylab="Sensitivity")
text(0.6,0.35,paste("WangGenSig ROC-AUC at 5 years = ",round(Wangscore$AUC,3)))
#abline(0,1)

par(new=TRUE)
Kim=Metabric_surv_score_Kim
nobs <- NROW(Kim)
cutoff <- 12*5
Kimscore= survivalROC(Stime=Kim$OS_MONTHS, status=Kim$OS_STATUS, marker = Kim$scores, predict.time = cutoff,span = 0.25*nobs^(-0.20))
plot(Kimscore$FP, Kimscore$TP, type="l",col="orange",lwd=2, xlim=c(0,1), ylim=c(0,1),xlab="1-Specificity",ylab="Sensitivity")
text(0.6,0.30,paste("KimGenSig ROC-AUC at 5 years = ",round(Kimscore$AUC,3)))
#abline(0,1)

par(new=TRUE)
Mansour=Metabric_surv_score_Mansour
nobs <- NROW(Mansour)
cutoff <- 12*5
Mansourscore= survivalROC(Stime=Mansour$OS_MONTHS, status=Mansour$OS_STATUS, marker = Mansour$scores, predict.time = cutoff,span = 0.25*nobs^(-0.20))
plot(Mansourscore$FP, Mansourscore$TP, type="l",col="grey",lwd=2, xlim=c(0,1), ylim=c(0,1),xlab="1-Specificity",ylab="Sensitivity")
text(0.6,0.25,paste("MansourGenSig ROC-AUC at 5 years = ",round(Mansourscore$AUC,3)))
#abline(0,1)

dev.off()
