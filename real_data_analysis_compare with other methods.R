
#### In real data analysis, use "rda", "e1071" to compare the performance ####

# 20231103: change chosen pairs which is greater than 29.
#           combine two R codes into one code

set.seed(34)

#### read data and clean ####

library(JGL)
library(NetDA)
library(XICOR)
library(dplyr)
library(glasso)
library(GGally)
library(network)
library(sna)
library(scales)
library(ggplot2)
library(sqldf)

read.table("GCM_Training.csv",header=F, sep = ",")->DATA

data_original = NULL


for(i in 1:144)  {
  data_original = cbind(data_original, DATA[,(2*i-1)])
}

dim(data_original)[1] -> p
dim(data_original)[2] -> n


gene = NULL

for(j in 1:p) {
  
  vt = data_original[j,]
  
  vt = (vt - mean(vt))/sqrt(var(vt))
  
  gene = rbind(gene,vt)
  
}

X = t(gene)

y = as.factor(c(rep(1,8),rep(2,8),rep(3,8),rep(4,8),rep(5,16),rep(6,8),rep(7,8),rep(8,8),rep(9,24),rep(10,8),
                rep(11,8),rep(12,8),rep(13,8),rep(14,16)))

levels(y) = c("BR","PR","LU","CO","LY","BL","ML","UT","LE","RE","PA","OV","ME","CNS")

level_y = c("BR","PR","LU","CO","LY","BL","ML","UT","LE","RE","PA","OV","ME","CNS")

########################################################################################

test_data <- read.table("GCM_Test.csv",header=F, sep = ",")

test_original = NULL
for(i in 1:54){
  test_original = cbind(test_original, test_data[,(2*i-1)])
}
dim(test_original)[2] -> n_test
gene_test = NULL

for(j in 1:p) {
  
  vt_test = test_original[j,]
  
  vt_test = (vt_test - mean(vt_test))/sqrt(var(vt_test))
  
  gene_test = rbind(gene_test,vt_test)
}

X_test = t(gene_test)

y_test = as.factor(c(rep(1,3),rep(2,4),rep(3,2),rep(4,4),rep(5,4),rep(6,4),rep(7,2),rep(8,3),rep(9,6),rep(10,4),
                     rep(11,3),rep(12,6),rep(13,3),rep(14,6)))
levels(y_test) = c("BR","PR","LU","CO","LY","BL","ML","UT","LE","RE","PA","OV","ME","CNS")

level_y_test = c("BR","PR","LU","CO","LY","BL","ML","UT","LE","RE","PA","OV","ME","CNS")

########################################################################################

#### precision matrix ####

y = as.factor(c(rep(1,8),rep(2,8),rep(3,8),rep(4,8),rep(5,16),rep(6,8),rep(7,8),rep(8,8),rep(9,24),rep(10,8),
                rep(11,8),rep(12,8),rep(13,8),rep(14,16)))

classes <- length(unique(y))
load("del_ind.RData")

library("rda")
library("e1071")
library("glmnet")
library("pemultinom")

#### rda ####   
data1 <- as.data.frame(cbind(y,X))

fit <- rda(t(X),y,alpha=0.5,prior=rep(1/14,14),delta=1.7)
predict_rda <- predict.rda(object=fit,x=t(X),y=y,xnew=t(X_test))

comb <- cbind(y_test,predict_rda)

# rda performance 
TP <- NULL; FP <- NULL; FN <- NULL
PRE <- NULL; REC <- NULL
for(i in 1:classes){
  TP[i] <- length(which(comb[,1]==i & comb[,2]==i))
  FP[i] <- length(which(comb[,1]!=i & comb[,2]==i))
  FN[i] <- length(which(comb[,1]==i & comb[,2]!=i))
  PRE[i] <- TP[i]/(TP[i]+FP[i])
  REC[i] <- TP[i]/(TP[i]+FN[i])
}
PRE_micro <- sum(TP)/sum(TP+FP)
REC_micro <- sum(TP)/sum(TP+FN)
F_micro <- 2 * (PRE_micro * REC_micro) / (PRE_micro + REC_micro)
PRE_macro <- mean(PRE)
REC_macro <- mean(REC)
F_macro <- 2 * (PRE_macro * REC_macro) / (PRE_macro + REC_macro)
PRE_micro; REC_micro; F_micro
PRE_macro; REC_macro; F_macro
print("rda")

#### e1071 #### 

#### svm ####
fit <- svm(X,y)
predict_svm <- predict(fit,X_test)

comb <- cbind(y_test,predict_svm)

# svm performance 
TP <- NULL; FP <- NULL; FN <- NULL
PRE <- NULL; REC <- NULL
for(i in 1:classes){
  TP[i] <- length(which(comb[,1]==i & comb[,2]==i))
  FP[i] <- length(which(comb[,1]!=i & comb[,2]==i))
  FN[i] <- length(which(comb[,1]==i & comb[,2]!=i))
  PRE[i] <- TP[i]/(TP[i]+FP[i])
  REC[i] <- TP[i]/(TP[i]+FN[i])
}
PRE_micro <- sum(TP)/sum(TP+FP)
REC_micro <- sum(TP)/sum(TP+FN)
F_micro <- 2 * (PRE_micro * REC_micro) / (PRE_micro + REC_micro)
PRE_macro <- mean(PRE)
REC_macro <- mean(REC)
F_macro <- 2 * (PRE_macro * REC_macro) / (PRE_macro + REC_macro)
PRE_micro; REC_micro; F_micro
PRE_macro; REC_macro; F_macro
print("svm")


#### naiveBayes ####
fit <- naiveBayes(X, y, data=data1)
predict_Bayes <- predict(fit,X_test)

comb <- cbind(y_test,predict_Bayes)

# naiveBayes performance 
TP <- NULL; FP <- NULL; FN <- NULL
PRE <- NULL; REC <- NULL
for(i in 1:classes){
  TP[i] <- length(which(comb[,1]==i & comb[,2]==i))
  FP[i] <- length(which(comb[,1]!=i & comb[,2]==i))
  FN[i] <- length(which(comb[,1]==i & comb[,2]!=i))
  PRE[i] <- TP[i]/(TP[i]+FP[i])
  REC[i] <- TP[i]/(TP[i]+FN[i])
}
PRE_micro <- sum(TP)/sum(TP+FP)
REC_micro <- sum(TP)/sum(TP+FN)
F_micro <- 2 * (PRE_micro * REC_micro) / (PRE_micro + REC_micro)
PRE_macro <- mean(PRE)
REC_macro <- mean(REC)
F_macro <- 2 * (PRE_macro * REC_macro) / (PRE_macro + REC_macro)
PRE_micro; REC_micro; F_micro
PRE_macro; REC_macro; F_macro
print("naiveBayes")


#### glmnet multinomial logistic regression ####

cvfit <- cv.glmnet(X, y, family = "multinomial", alpha = 1)
fit <- glmnet(X, y, family = "multinomial", alpha = 1)
predict_logistic <- predict(fit, X_test, s = cvfit$lambda.min, type = "class")
comb <- cbind(y_test, as.numeric(predict_logistic))

# logistic performance 
TP <- NULL; FP <- NULL; FN <- NULL
PRE <- NULL; REC <- NULL
for(i in 1:classes){
  TP[i] <- length(which(comb[,1]==i & comb[,2]==i))
  FP[i] <- length(which(comb[,1]!=i & comb[,2]==i))
  FN[i] <- length(which(comb[,1]==i & comb[,2]!=i))
  PRE[i] <- TP[i]/(TP[i]+FP[i])
  REC[i] <- TP[i]/(TP[i]+FN[i])
}
PRE_micro <- sum(TP)/sum(TP+FP)
REC_micro <- sum(TP)/sum(TP+FN)
F_micro <- 2 * (PRE_micro * REC_micro) / (PRE_micro + REC_micro)
PRE_macro <- mean(PRE)
REC_macro <- mean(REC)
F_macro <- 2 * (PRE_macro * REC_macro) / (PRE_macro + REC_macro)
PRE_micro; REC_micro; F_micro
PRE_macro; REC_macro; F_macro
print("glmnet_logistic")


#### before do sparse logistic, do feature selection ####
X <- X[,-del_ind]
pairs = NULL
thre = 0.4
XI = NULL; pair = NULL
p = dim(X)[2]

for(i in 1:p){
  x = max(xicor(X[,i],as.numeric(y)),xicor(as.numeric(y),X[,i]))*1
  if(x > thre) {
    pair = rbind(pair,i)
    XI = c(XI,x)
  }
}
pairs <- cbind(pair,XI)
pairs ## xypairs_linear.RData

pairs <- as.data.frame(pairs)

# label the homogeneous variables
pairs_ord <- sqldf("
                   SELECT * FROM pairs
                   ORDER BY XI DESC
                  ")
### change the dimension which will be greater than 29.
# pairs_chosen <- pairs_ord[1:ceiling(n/log(n)),]
pairs_chosen <- pairs_ord[1:35,]


#### sparse multinomial logistic regression model ####
# eg. Abramovich et al. 2021; Tian et al. 2023 
#colnames(X) <- paste("X",1:15015,sep="")
#fit <- cv.pemultinom(X, y, ncores=2)
fit <- cv.pemultinom(X[,pairs_chosen[,1]], y, ncores=2)
beta <- fit$beta.min
predict_Tian_logistic <- predict_pemultinom(fit$beta.min, ref=3, xnew=X_test[,pairs_chosen[,1]], type="class") 
comb <- cbind(y_test, as.numeric(predict_Tian_logistic))

# Tian logistic performance 
TP <- NULL; FP <- NULL; FN <- NULL
PRE <- NULL; REC <- NULL
for(i in 1:classes){
  TP[i] <- length(which(comb[,1]==i & comb[,2]==i))
  FP[i] <- length(which(comb[,1]!=i & comb[,2]==i))
  FN[i] <- length(which(comb[,1]==i & comb[,2]!=i))
  PRE[i] <- TP[i]/(TP[i]+FP[i])
  REC[i] <- TP[i]/(TP[i]+FN[i])
}
PRE_micro <- sum(TP)/sum(TP+FP)
REC_micro <- sum(TP)/sum(TP+FN)
F_micro <- 2 * (PRE_micro * REC_micro) / (PRE_micro + REC_micro)
PRE_macro <- mean(PRE)
REC_macro <- mean(REC)
F_macro <- 2 * (PRE_macro * REC_macro) / (PRE_macro + REC_macro)
PRE_micro; REC_micro; F_micro
PRE_macro; REC_macro; F_macro
print("Tian_logistic")

