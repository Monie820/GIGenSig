############################################################################################################
##################################### Stratification analysis ##############################################
############################################################################################################
########################## Clinical information for train ,test, and total Metabric datasets ###############
load("/Users/username/input/cBioPortal/brca_metabric/Analysis/Train_test_total_scores.RData")
Metabric_surv_score[which(Metabric_surv_score$TUMOR_STAGE==1),6]="I/II"
Metabric_surv_score[which(Metabric_surv_score$TUMOR_STAGE==2),6]="I/II"
Metabric_surv_score[which(Metabric_surv_score$TUMOR_STAGE==3),6]="III/IV"
Metabric_surv_score[which(Metabric_surv_score$GRADE==1),7]="G1/2"
Metabric_surv_score[which(Metabric_surv_score$GRADE==2),7]="G1/2"
Metabric_surv_score[which(Metabric_surv_score$GRADE==3),7]="G3"

age_len=c(length(which(Metabric_surv_score$AGE_AT_DIAGNOSIS < 55)),length(which(Metabric_surv_score$AGE_AT_DIAGNOSIS >= 55)))
menopausal_len=c(length(which(Metabric_surv_score$INFERRED_MENOPAUSAL_STATE == "Pre")),length(which(Metabric_surv_score$INFERRED_MENOPAUSAL_STATE == "Post")))
stage_len=c(length(which(Metabric_surv_score$TUMOR_STAGE == "I/II")),length(which(Metabric_surv_score$TUMOR_STAGE == "III/IV")))
grade_len=c(length(which(Metabric_surv_score$GRADE == "G1/2")),length(which(Metabric_surv_score$GRADE == "G3")))
status_len=c(length(which(Metabric_surv_score$OS_STATUS == 0)),length(which(Metabric_surv_score$OS_STATUS == 1)))

total_infor=c(age_len,menopausal_len,stage_len,grade_len,status_len)

clinical_infor=data.frame(train_infor,train_infor/150*100,test_infor,test_infor/149*100,total_infor,total_infor/299*100)
write.table(clinical_infor,"/Users/username/input/cBioPortal/brca_metabric/Analysis/Clinical_information.txt",sep="\t",row.names=F)
clinical_infor2=data.frame(train_infor,test_infor,total_infor)

Age <- chisq.test(clinical_infor2[1:2,])$p.value
Menopausal <- chisq.test(clinical_infor2[3:4,])$p.value
Stage <- chisq.test(clinical_infor2[5:6,])$p.value
Grade <- chisq.test(clinical_infor2[7:8,])$p.value
Status <- chisq.test(clinical_infor2[9:10,])$p.value
P_value <- c(Age,Menopausal,Stage,Grade,Status)

############# Univariate and multivariate analysis for train ,test, and total Metabric datasets ##########
load("/Users/username/input/cBioPortal/brca_metabric/Analysis/Train_test_total_scores.RData")
### uni train
Metabric_surv_score_train[which(Metabric_surv_score_train$TUMOR_STAGE==1),6]="I/II"
Metabric_surv_score_train[which(Metabric_surv_score_train$TUMOR_STAGE==2),6]="I/II"
Metabric_surv_score_train[which(Metabric_surv_score_train$TUMOR_STAGE==3),6]="III/IV"
Metabric_surv_score_train[which(Metabric_surv_score_train$GRADE==1),7]="G1/2"
Metabric_surv_score_train[which(Metabric_surv_score_train$GRADE==2),7]="G1/2"
Metabric_surv_score_train[which(Metabric_surv_score_train$GRADE==3),7]="G3"
surv.stats=c()
for(i in c(8,4,5,6,7)){
  surv.stat=summary(coxph(Surv(OS_MONTHS,as.numeric(OS_STATUS))~Metabric_surv_score_train[,i],Metabric_surv_score_train))
  hr=surv.stat[[7]][2] ## HR
  pvalue=surv.stat[[7]][5] ## p
  ci_lower=surv.stat[[8]][3]
  ci_upper=surv.stat[[8]][4] ## 95%CI
  stats=cbind(hr,ci_lower,ci_upper,pvalue)
  surv.stats=rbind(surv.stats,stats)
}
uni_surv_stats_train=surv.stats

### uni test
Metabric_surv_score_test[which(Metabric_surv_score_test$TUMOR_STAGE==1),6]="I/II"
Metabric_surv_score_test[which(Metabric_surv_score_test$TUMOR_STAGE==2),6]="I/II"
Metabric_surv_score_test[which(Metabric_surv_score_test$TUMOR_STAGE==3),6]="III/IV"
Metabric_surv_score_test[which(Metabric_surv_score_test$GRADE==1),7]="G1/2"
Metabric_surv_score_test[which(Metabric_surv_score_test$GRADE==2),7]="G1/2"
Metabric_surv_score_test[which(Metabric_surv_score_test$GRADE==3),7]="G3"
surv.stats=c()
for(i in c(8,4,5,6,7)){
  surv.stat=summary(coxph(Surv(OS_MONTHS,as.numeric(OS_STATUS))~Metabric_surv_score_test[,i],Metabric_surv_score_test))
  hr=surv.stat[[7]][2] ## HR
  pvalue=surv.stat[[7]][5] ## p
  ci_lower=surv.stat[[8]][3]
  ci_upper=surv.stat[[8]][4] ## 95%CI
  stats=cbind(hr,ci_lower,ci_upper,pvalue)
  surv.stats=rbind(surv.stats,stats)
}
uni_surv_stats_test=surv.stats

### uni total
Metabric_surv_score[which(Metabric_surv_score$TUMOR_STAGE==1),6]="I/II"
Metabric_surv_score[which(Metabric_surv_score$TUMOR_STAGE==2),6]="I/II"
Metabric_surv_score[which(Metabric_surv_score$TUMOR_STAGE==3),6]="III/IV"
Metabric_surv_score[which(Metabric_surv_score$GRADE==1),7]="G1/2"
Metabric_surv_score[which(Metabric_surv_score$GRADE==2),7]="G1/2"
Metabric_surv_score[which(Metabric_surv_score$GRADE==3),7]="G3"
surv.stats=c()
for(i in c(8,4,5,6,7)){
  surv.stat=summary(coxph(Surv(OS_MONTHS,as.numeric(OS_STATUS))~Metabric_surv_score[,i],Metabric_surv_score))
  hr=surv.stat[[7]][2] ## HR
  pvalue=surv.stat[[7]][5] ## p
  ci_lower=surv.stat[[8]][3]
  ci_upper=surv.stat[[8]][4] ## 95%CI
  stats=cbind(hr,ci_lower,ci_upper,pvalue)
  surv.stats=rbind(surv.stats,stats)
}
uni_surv_stats_total=surv.stats

### mul train
surv.stat=summary(coxph(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ scores + AGE_AT_DIAGNOSIS + 
                          INFERRED_MENOPAUSAL_STATE + TUMOR_STAGE + GRADE, Metabric_surv_score_train))
hr=surv.stat[[7]][,2]  ## HR
pvalue=surv.stat[[7]][,5] ## p
ci_lower=surv.stat[[8]][,3]
ci_upper=surv.stat[[8]][,4] ## 95%CI
mul_surv_stats_train=cbind(hr,ci_lower,ci_upper,pvalue)

### mul test
surv.stat=summary(coxph(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ scores + AGE_AT_DIAGNOSIS + 
                          INFERRED_MENOPAUSAL_STATE + TUMOR_STAGE + GRADE, Metabric_surv_score_test))
hr=surv.stat[[7]][,2]  ## HR
pvalue=surv.stat[[7]][,5] ## p
ci_lower=surv.stat[[8]][,3]
ci_upper=surv.stat[[8]][,4] ## 95%CI
mul_surv_stats_test=cbind(hr,ci_lower,ci_upper,pvalue)

### mul total
surv.stat=summary(coxph(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ scores + AGE_AT_DIAGNOSIS + 
                          INFERRED_MENOPAUSAL_STATE + TUMOR_STAGE + GRADE, Metabric_surv_score))
hr=surv.stat[[7]][,2]  ## HR
pvalue=surv.stat[[7]][,5] ## p
ci_lower=surv.stat[[8]][,3]
ci_upper=surv.stat[[8]][,4] ## 95%CI
mul_surv_stats_total=cbind(hr,ci_lower,ci_upper,pvalue)

uni=rbind(uni_surv_stats_train,uni_surv_stats_test,uni_surv_stats_total)
mul=rbind(mul_surv_stats_train,mul_surv_stats_test,mul_surv_stats_total)
uni_mul=cbind(uni,mul)
write.table(uni_mul,"/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_train_test_total_stats.txt",sep="\t",quote=FALSE)



############################# Stratification survival plot analysis ##################################
load("/Users/username/input/cBioPortal/brca_metabric/Analysis/Train_test_total_scores.RData")
#survival analysis for age
age_pos1=which(Metabric_surv_score$AGE_AT_DIAGNOSIS < 55)
age_pos2=which(Metabric_surv_score$AGE_AT_DIAGNOSIS >= 55)
Metabric_surv_score_less55=Metabric_surv_score[age_pos1,]
Metabric_surv_score_more55=Metabric_surv_score[age_pos2,]

surv <- survfit(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ classify, data = Metabric_surv_score_less55)
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_surv_score_less55.pdf")
plot(surv,mark.time=TRUE,mark=3,lty=1,lwd=3,col=3:2,xlab="Years",ylab="Survival",xscale=12)
pp<-survdiff(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ classify, data = Metabric_surv_score_less55)
p.val <- 1 - pchisq(pp$chisq, length(pp$n) - 1) 
p.val
text(240,0.7,"Log-rank p = 0.019") #significant
dev.off()

surv <- survfit(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ classify, data = Metabric_surv_score_more55)
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_surv_score_more55.pdf")
plot(surv,mark.time=TRUE,mark=3,lty=1,lwd=3,col=3:2,xlab="Years",ylab="Survival",xscale=12)
pp<-survdiff(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ classify, data = Metabric_surv_score_more55)
p.val <- 1 - pchisq(pp$chisq, length(pp$n) - 1) 
p.val
text(240,0.7,"Log-rank p = 0.002") #significant
dev.off()

#survival analysis for manopausal
menopausal_pos1=which(Metabric_surv_score$INFERRED_MENOPAUSAL_STATE == "Pre")
menopausal_pos2=which(Metabric_surv_score$INFERRED_MENOPAUSAL_STATE == "Post")
Metabric_surv_score_Pre=Metabric_surv_score[menopausal_pos1,]
Metabric_surv_score_Post=Metabric_surv_score[menopausal_pos2,]

surv <- survfit(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ classify, data = Metabric_surv_score_Pre)
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_surv_score_Pre.pdf")
plot(surv,mark.time=TRUE,mark=3,lty=1,lwd=3,col=3:2,xlab="Years",ylab="Survival",xscale=12)
pp<-survdiff(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ classify, data = Metabric_surv_score_Pre)
p.val <- 1 - pchisq(pp$chisq, length(pp$n) - 1) 
p.val
text(240,0.7,"Log-rank p = 0.074") #significant
dev.off()

surv <- survfit(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ classify, data = Metabric_surv_score_Post)
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_surv_score_Post.pdf")
plot(surv,mark.time=TRUE,mark=3,lty=1,lwd=3,col=3:2,xlab="Years",ylab="Survival",xscale=12)
pp<-survdiff(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ classify, data = Metabric_surv_score_Post)
p.val <- 1 - pchisq(pp$chisq, length(pp$n) - 1) 
p.val
text(240,0.7,"Log-rank p < 0.001") #significant
dev.off()


#survival analysis for stage
stage_pos1=which(Metabric_surv_score$TUMOR_STAGE == "I/II")
stage_pos2=which(Metabric_surv_score$TUMOR_STAGE == "III/IV")
Metabric_surv_score_lowStage=Metabric_surv_score[stage_pos1,]
Metabric_surv_score_highStage=Metabric_surv_score[stage_pos2,]

surv <- survfit(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ classify, data = Metabric_surv_score_lowStage)
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_surv_score_lowStage.pdf")
plot(surv,mark.time=TRUE,mark=3,lty=1,lwd=3,col=3:2,xlab="Years",ylab="Survival",xscale=12)
pp<-survdiff(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ classify, data = Metabric_surv_score_lowStage)
p.val <- 1 - pchisq(pp$chisq, length(pp$n) - 1) 
p.val
text(240,0.7,"Log-rank p < 0.001") #significant
dev.off()

surv <- survfit(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ classify, data = Metabric_surv_score_highStage)
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_surv_score_highStage.pdf")
plot(surv,mark.time=TRUE,mark=3,lty=1,lwd=3,col=3:2,xlab="Years",ylab="Survival",xscale=12)
pp<-survdiff(Surv(OS_MONTHS,as.numeric(OS_STATUS)) ~ classify, data = Metabric_surv_score_highStage)
p.val <- 1 - pchisq(pp$chisq, length(pp$n) - 1) 
p.val
text(240,0.7,"Log-rank p = 0.387") #significant
dev.off()

# Number at risk
Low_risk=Metabric_surv_score_highStage[which(Metabric_surv_score_highStage$classify==1),]
Low_risk_surv_0 = length(Low_risk[which(Low_risk$OS_MONTHS >= 0),3])
Low_risk_surv_5 = length(Low_risk[which(Low_risk$OS_MONTHS >= 5*12),3])
Low_risk_surv_10 = length(Low_risk[which(Low_risk$OS_MONTHS >= 10*12),3])
Low_risk_surv_15 = length(Low_risk[which(Low_risk$OS_MONTHS >= 15*12),3])
Low_risk_surv_20 = length(Low_risk[which(Low_risk$OS_MONTHS >= 20*12),3])
Low_risk_surv_25 = length(Low_risk[which(Low_risk$OS_MONTHS >= 25*12),3])

High_risk=Metabric_surv_score_highStage[which(Metabric_surv_score_highStage$classify==2),]
High_risk_surv_0 = length(High_risk[which(High_risk$OS_MONTHS >= 0),3])
High_risk_surv_5 = length(High_risk[which(High_risk$OS_MONTHS >= 5*12),3])
High_risk_surv_10 = length(High_risk[which(High_risk$OS_MONTHS >= 10*12),3])
High_risk_surv_15 = length(High_risk[which(High_risk$OS_MONTHS >= 15*12),3])
High_risk_surv_20 = length(High_risk[which(High_risk$OS_MONTHS >= 20*12),3])
High_risk_surv_25 = length(High_risk[which(High_risk$OS_MONTHS >= 25*12),3])

Number_at_risk_highStage=rbind(c(Low_risk_surv_0,Low_risk_surv_5,Low_risk_surv_10,Low_risk_surv_15,Low_risk_surv_20,Low_risk_surv_25),
                               c(High_risk_surv_0,High_risk_surv_5,High_risk_surv_10,High_risk_surv_15,High_risk_surv_20,High_risk_surv_25))

Number_at_risk_Stratification=rbind(Number_at_risk_less55, Number_at_risk_more55, Number_at_risk_Pre, Number_at_risk_Post,
                                    Number_at_risk_lowStage, Number_at_risk_highStage)
write.table(Number_at_risk_Stratification,"/Users/username/input/cBioPortal/brca_metabric/Analysis/Number_at_risk_Stratification.txt",sep="\t",col.names=F,row.names=F,quote=F)



##########################  Number at risk ##################################
# for train test total
load("/Users/username/input/cBioPortal/brca_metabric/Analysis/Train_test_total_scores.RData")

Low_risk=Metabric_surv_score[which(Metabric_surv_score$classify==1),]
Low_risk_surv_0 = length(Low_risk[which(Low_risk$OS_MONTHS >= 0),3])
Low_risk_surv_5 = length(Low_risk[which(Low_risk$OS_MONTHS >= 5*12),3])
Low_risk_surv_10 = length(Low_risk[which(Low_risk$OS_MONTHS >= 10*12),3])
Low_risk_surv_15 = length(Low_risk[which(Low_risk$OS_MONTHS >= 15*12),3])
Low_risk_surv_20 = length(Low_risk[which(Low_risk$OS_MONTHS >= 20*12),3])
Low_risk_surv_25 = length(Low_risk[which(Low_risk$OS_MONTHS >= 25*12),3])

High_risk=Metabric_surv_score[which(Metabric_surv_score$classify==2),]
High_risk_surv_0 = length(High_risk[which(High_risk$OS_MONTHS >= 0),3])
High_risk_surv_5 = length(High_risk[which(High_risk$OS_MONTHS >= 5*12),3])
High_risk_surv_10 = length(High_risk[which(High_risk$OS_MONTHS >= 10*12),3])
High_risk_surv_15 = length(High_risk[which(High_risk$OS_MONTHS >= 15*12),3])
High_risk_surv_20 = length(High_risk[which(High_risk$OS_MONTHS >= 20*12),3])
High_risk_surv_25 = length(High_risk[which(High_risk$OS_MONTHS >= 25*12),3])

Number_at_risk_total=rbind(c(Low_risk_surv_0,Low_risk_surv_5,Low_risk_surv_10,Low_risk_surv_15,Low_risk_surv_20,Low_risk_surv_25),
                           c(High_risk_surv_0,High_risk_surv_5,High_risk_surv_10,High_risk_surv_15,High_risk_surv_20,High_risk_surv_25))

Number_at_risk=rbind(Number_at_risk_train,Number_at_risk_test,Number_at_risk_total)
write.table(Number_at_risk,"/Users/username/input/cBioPortal/brca_metabric/Analysis/Number_at_risk.txt",sep="\t",col.names=F,row.names=F,quote=F)


##########################  Survvial ROC for train, test and total dataset ##################################
load("/Users/username/input/cBioPortal/brca_metabric/Analysis/Train_test_total_scores.RData")
library(survivalROC)
## trainSCORE, METHOD = NNE
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/Train_survival_at5years.pdf")
train=Metabric_surv_score_train
nobs <- NROW(train)
cutoff <- 12*5
trainscore= survivalROC(Stime=train$OS_MONTHS, status=train$OS_STATUS, marker = train$scores,
                        predict.time = cutoff,span = 0.25*nobs^(-0.20))
plot(trainscore$FP, trainscore$TP, type="l",col="red",lwd=2, xlim=c(0,1), ylim=c(0,1),xlab="1-Specificity",ylab="Sensitivity")
text(0.6,0.4,paste("ROC-AUC at 5 years = ",round(trainscore$AUC,3)))
abline(0,1)
dev.off()


pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/TestandTotal_survival_at5years.pdf",width=14,height=7)
par(mfrow = c(1,2))
test=Metabric_surv_score_test
nobs <- NROW(test)
cutoff <- 12*5
testscore= survivalROC(Stime=test$OS_MONTHS, status=test$OS_STATUS, marker = test$scores,
                       predict.time = cutoff,span = 0.25*nobs^(-0.20))
plot(testscore$FP, testscore$TP, type="l",col="red",lwd=2, xlim=c(0,1), ylim=c(0,1),xlab="1-Specificity",ylab="Sensitivity")
text(0.6,0.4,paste("ROC-AUC at 5 years = ",round(testscore$AUC,3)))
abline(0,1)

total=Metabric_surv_score
nobs <- NROW(total)
cutoff <- 12*5
totalscore= survivalROC(Stime=total$OS_MONTHS, status=total$OS_STATUS, marker = total$scores,
                        predict.time = cutoff,span = 0.25*nobs^(-0.20))
plot(totalscore$FP, totalscore$TP, type="l",col="red",lwd=2, xlim=c(0,1), ylim=c(0,1),xlab="1-Specificity",ylab="Sensitivity")
text(0.6,0.4,paste("ROC-AUC at 5 years = ",round(totalscore$AUC,3)))
abline(0,1)

dev.off()

