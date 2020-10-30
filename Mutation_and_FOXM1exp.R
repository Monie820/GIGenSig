##########################################################################################
############################ Mutation and FOXM1 expression pattern #######################
##########################################################################################
load("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_exp_111_2_train_test_total.RData")
load("/Users/username/input/cBioPortal/brca_metabric/Analysis/Train_test_total_scores.RData")
################# Train: expression patterns for 11 genes in GIGenSig ####################
### expression patterns
Metabric_exp_111_train=Metabric_exp2_train[match(DEGs_2$gene,rownames(Metabric_exp2_train)),]
Metabric_exp_111_2_train=t(Metabric_exp_111_train)
scores_pos=order(Metabric_surv_score_train$scores)
mul2_pos=match(rownames(Metabric_mul2_stats2_train),rownames(Metabric_exp_111_train))
Metabric_mul2_exp=Metabric_exp_111_train[mul2_pos,scores_pos]
library(gplots)
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_mul2_exp_train.pdf",width=7, height=3)
Metabric_mul2_exp=scale(Metabric_mul2_exp) #z-score normal
heat.exp.value=as.matrix(Metabric_mul2_exp)
mycol <- colorpanel(100000,"blue", "white", "red")
heatmap.2(heat.exp.value,
          col=mycol, key=FALSE, symkey=FALSE, density.info="none", trace="none", 
          labCol=NA, labRow=rownames(heat.exp.value),cexRow=0.7, Colv=FALSE, Rowv=TRUE, dendrogram="none")
dev.off()

### the distribution of mutation/FOXM1 expxpression
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_mul2_mutNOExp_train.pdf",width=7, height=7)
par(mfrow = c(2,1))
scores_pos=order(Metabric_surv_score_train$scores)
Metabric_surv_score_sort=Metabric_surv_score_train[scores_pos,]
Metabric_surv_score_sort[,1]=rownames(Metabric_exp_111_2_train)[scores_pos]

Low_count=Binary_sum2[match(Metabric_surv_score_sort[Metabric_surv_score_sort$classify==1,1],names(Binary_sum2))]
High_count=Binary_sum2[match(Metabric_surv_score_sort[Metabric_surv_score_sort$classify==2,1],names(Binary_sum2))]
x=c(Low_count,High_count)
col=c(rep("blue",length(Low_count)),rep("red",length(High_count)))
plot(x,pch=16,cex=0.5,col=col,xlab="",ylab="Mutation count") # scatter plot

Low_FOXM1_exp=as.numeric(Metabric_exp2[which(rownames(Metabric_exp2)=="FOXM1"),match(Metabric_surv_score_sort[Metabric_surv_score_sort$classify==1,1],colnames(Metabric_exp2))])
High_FOXM1_exp=as.numeric(Metabric_exp2[which(rownames(Metabric_exp2)=="FOXM1"),match(Metabric_surv_score_sort[Metabric_surv_score_sort$classify==2,1],colnames(Metabric_exp2))])
x = c(Low_FOXM1_exp,High_FOXM1_exp)
col=c(rep("blue",length(Low_FOXM1_exp)),rep("red",length(High_FOXM1_exp)))
barplot(x,width=0.5,space=0.8,border=NA,col=col, xlab="",ylab="FOXM1 expression")  # barplot
box()
dev.off()

### boxplot of mutation count/FOXM1 expxpression #######
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_mul2_mutNOExpBoxplot_train.pdf",width=7, height=7)
par(mfrow = c(1, 2))
ddf_count = data.frame(count = c(Low_count,High_count), class = c(rep("Low-risk",length(Low_count)),rep("High-risk",length(High_count))))
boxplot(count ~ class, data = ddf_count, lwd = 0.5, ylab = 'Mutation number')
stripchart(count ~ class, vertical = TRUE, data = ddf_count, 
           method = "jitter", add = TRUE, pch = 20, col = c("red","blue"),cex=0.7)
wilcox.test(Low_count,High_count)$p.value
text(1.5,2000,"P < 0.001")

ddf_FOXM1_exp = data.frame(exp = c(Low_FOXM1_exp,High_FOXM1_exp), class = c(rep("Low-risk",length(Low_FOXM1_exp)),rep("High-risk",length(High_FOXM1_exp))))
boxplot(exp ~ class, data = ddf_FOXM1_exp, lwd = 0.5, ylab = 'FOXM1 epression')
stripchart(exp ~ class, vertical = TRUE, data = ddf_FOXM1_exp, 
           method = "jitter", add = TRUE, pch = 20, col = c("red","blue"),cex=0.7)
wilcox.test(Low_FOXM1_exp,High_FOXM1_exp)$p.value
text(1.5,10,"P = 0.01")
dev.off()


############### Test: expression patterns for 12 Metabric_mul2_stats2 ####################
###expression patterns
Metabric_exp_111_test=Metabric_exp2_test[match(DEGs_2$gene,rownames(Metabric_exp2_test)),]
Metabric_exp_111_2_test=t(Metabric_exp_111_test)
scores_pos=order(Metabric_surv_score_test$scores)
mul2_pos=match(rownames(Metabric_mul2_stats2_train),rownames(Metabric_exp_111_test))
Metabric_mul2_exp=Metabric_exp_111_test[mul2_pos,scores_pos]
library(gplots)
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_mul2_exp_test.pdf",width=7, height=3)
Metabric_mul2_exp=scale(Metabric_mul2_exp) #z-score normal
heat.exp.value=as.matrix(Metabric_mul2_exp)
mycol <- colorpanel(100000,"blue", "white", "red")
heatmap.2(heat.exp.value,
          col=mycol, key=FALSE, symkey=FALSE, density.info="none", trace="none", 
          labCol=NA, labRow=rownames(heat.exp.value),cexRow=0.7, Colv=FALSE, Rowv=TRUE, dendrogram="none")
dev.off()

###the distribution of mutation/FOXM1 expression
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_mul2_mutNOExp_test.pdf",width=7, height=7)
par(mfrow = c(2,1))
scores_pos=order(Metabric_surv_score_test$scores)
Metabric_surv_score_sort=Metabric_surv_score_test[scores_pos,]
Metabric_surv_score_sort[,1]=rownames(Metabric_exp_111_2_test)[scores_pos]

Low_count=Binary_sum2[match(Metabric_surv_score_sort[Metabric_surv_score_sort$classify==1,1],names(Binary_sum2))]
High_count=Binary_sum2[match(Metabric_surv_score_sort[Metabric_surv_score_sort$classify==2,1],names(Binary_sum2))]
x=c(Low_count,High_count)
col=c(rep("blue",length(Low_count)),rep("red",length(High_count)))
plot(x,pch=16,cex=0.5,col=col,xlab="",ylab="Mutation count") # scatter plot

Low_FOXM1_exp=as.numeric(Metabric_exp2[which(rownames(Metabric_exp2)=="FOXM1"),match(Metabric_surv_score_sort[Metabric_surv_score_sort$classify==1,1],colnames(Metabric_exp2))])
High_FOXM1_exp=as.numeric(Metabric_exp2[which(rownames(Metabric_exp2)=="FOXM1"),match(Metabric_surv_score_sort[Metabric_surv_score_sort$classify==2,1],colnames(Metabric_exp2))])
x = c(Low_FOXM1_exp,High_FOXM1_exp)
col=c(rep("blue",length(Low_FOXM1_exp)),rep("red",length(High_FOXM1_exp)))
barplot(x,width=0.5,space=0.8,border=NA,col=col, xlab="",ylab="FOXM1 expression")  # barplot
box()
dev.off()

### boxplot of mutation count/FOXM1 expxpression 
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_mul2_mutNOExpBoxplot_test.pdf",width=7, height=7)
par(mfrow = c(1, 2))
ddf_count = data.frame(count = c(Low_count,High_count), class = c(rep("Low-risk",length(Low_count)),rep("High-risk",length(High_count))))
boxplot(count ~ class, data = ddf_count, lwd = 0.5, ylab = 'Mutation number')
stripchart(count ~ class, vertical = TRUE, data = ddf_count, 
           method = "jitter", add = TRUE, pch = 20, col = c("red","blue"),cex=0.7)
wilcox.test(Low_count,High_count)$p.value
text(1.5,2000,"P < 0.001")

ddf_FOXM1_exp = data.frame(exp = c(Low_FOXM1_exp,High_FOXM1_exp), class = c(rep("Low-risk",length(Low_FOXM1_exp)),rep("High-risk",length(High_FOXM1_exp))))
boxplot(exp ~ class, data = ddf_FOXM1_exp, lwd = 0.5, ylab = 'FOXM1 epression')
stripchart(exp ~ class, vertical = TRUE, data = ddf_FOXM1_exp, 
           method = "jitter", add = TRUE, pch = 20, col = c("red","blue"),cex=0.7)
wilcox.test(Low_FOXM1_exp,High_FOXM1_exp)$p.value
text(1.5,10,"P = 0.001")
dev.off()


############### Total: expression patterns for 12 Metabric_mul2_stats2 ####################
###expression patterns
Metabric_exp_111=Metabric_exp2[match(DEGs_2$gene,rownames(Metabric_exp2)),]
Metabric_exp_111_2=t(Metabric_exp_111)
scores_pos=order(Metabric_surv_score$scores)
mul2_pos=match(rownames(Metabric_mul2_stats2_train),rownames(Metabric_exp_111))
Metabric_mul2_exp=Metabric_exp_111[mul2_pos,scores_pos]
library(gplots)
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_mul2_exp_total.pdf",width=7, height=3)
Metabric_mul2_exp=scale(Metabric_mul2_exp) #z-score normal
heat.exp.value=as.matrix(Metabric_mul2_exp)
mycol <- colorpanel(100000,"blue", "white", "red")
heatmap.2(heat.exp.value,
          col=mycol, key=FALSE, symkey=FALSE, density.info="none", trace="none", 
          labCol=NA, labRow=rownames(heat.exp.value),cexRow=0.7, Colv=FALSE, Rowv=TRUE, dendrogram="none")
dev.off()

###the distribution of mutation/FOXM1 expxpression
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_mul2_mutNOExp_total.pdf",width=7, height=7)
par(mfrow = c(2,1))
scores_pos=order(Metabric_surv_score$scores)
Metabric_surv_score_sort=Metabric_surv_score[scores_pos,]
Metabric_surv_score_sort[,1]=rownames(Metabric_exp_111_2)[scores_pos]

Low_count=Binary_sum2[match(Metabric_surv_score_sort[Metabric_surv_score_sort$classify==1,1],names(Binary_sum2))]
High_count=Binary_sum2[match(Metabric_surv_score_sort[Metabric_surv_score_sort$classify==2,1],names(Binary_sum2))]
x=c(Low_count,High_count)
col=c(rep("blue",length(Low_count)),rep("red",length(High_count)))
plot(x,pch=16,cex=0.5,col=col,xlab="",ylab="Mutation count") # scatter plot

Low_FOXM1_exp=as.numeric(Metabric_exp2[which(rownames(Metabric_exp2)=="FOXM1"),match(Metabric_surv_score_sort[Metabric_surv_score_sort$classify==1,1],colnames(Metabric_exp2))])
High_FOXM1_exp=as.numeric(Metabric_exp2[which(rownames(Metabric_exp2)=="FOXM1"),match(Metabric_surv_score_sort[Metabric_surv_score_sort$classify==2,1],colnames(Metabric_exp2))])
x = c(Low_FOXM1_exp,High_FOXM1_exp)
col=c(rep("blue",length(Low_FOXM1_exp)),rep("red",length(High_FOXM1_exp)))
barplot(x,width=0.5,space=0.8,border=NA,col=col, xlab="",ylab="FOXM1 expression")  # barplot
box()
dev.off()

### boxplot of mutation count/FOXM1 expxpression 
pdf("/Users/username/input/cBioPortal/brca_metabric/Analysis/Metabric_mul2_mutNOExpBoxplot_total.pdf",width=7, height=7)
par(mfrow = c(1, 2))
ddf_count = data.frame(count = c(Low_count,High_count), class = c(rep("Low-risk",length(Low_count)),rep("High-risk",length(High_count))))
boxplot(count ~ class, data = ddf_count, lwd = 0.5, ylab = 'Mutation number')
stripchart(count ~ class, vertical = TRUE, data = ddf_count, 
           method = "jitter", add = TRUE, pch = 20, col = c("red","blue"),cex=0.7)
wilcox.test(Low_count,High_count)$p.value
text(1.5,2000,"P < 0.001")

ddf_FOXM1_exp = data.frame(exp = c(Low_FOXM1_exp,High_FOXM1_exp), class = c(rep("Low-risk",length(Low_FOXM1_exp)),rep("High-risk",length(High_FOXM1_exp))))
boxplot(exp ~ class, data = ddf_FOXM1_exp, lwd = 0.5, ylab = 'FOXM1 epression')
stripchart(exp ~ class, vertical = TRUE, data = ddf_FOXM1_exp, 
           method = "jitter", add = TRUE, pch = 20, col = c("red","blue"),cex=0.7)
wilcox.test(Low_FOXM1_exp,High_FOXM1_exp)$p.value
text(1.5,10,"P < 0.001")
dev.off()
