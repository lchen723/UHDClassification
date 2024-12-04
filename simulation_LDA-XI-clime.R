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
library(clime)

n=500
p=3000
PRE_ave = NULL; REC_ave = NULL; F_ave = NULL

strue=diag(1,p,p)
strue[2,1]=strue[1,2]=1
strue[3,1]=strue[1,3]=1;strue[3,2]=strue[2,3]=1
strue[4,3]=strue[3,4]=1;strue[4,6]=strue[6,4]=1
strue[5,2]=strue[2,5]=1;strue[5,6]=strue[6,5]=1
strue[7,3]=strue[3,7]=1;strue[7,5]=strue[5,7]=1

#repeat 100 times
#Theta_ng=array(dim = c(p,p,100))
for (q in 1:1){
 # cat("\r",round(q/100*100,2), '%     ')
  
  #e3?a¡V?a?¡Ân*pc??c?ce?¢G
  x1=runif(n,-1,1)
  
  x2=6*cos(x1)+runif(n,-1,1) #x2a¡¦?x1a??e¡X?
  x3=5*sin(x1)+x2+rnorm(n,0,1) #x3a¡¦?x1a£á?x2a??e¡X?
  
  x6=runif(n,-1,1)
  
  x4=5*cos(x3*x6)+3*x3+3*x6+rnorm(n,0,1) #x4a¡¦?x3a£á?x6a??e¡X?
  x5=0.05*(x2+x6)^3+rnorm(n,0,1) #x5a¡¦?x2a£á?x6a??e¡X?
  x7=6*cos(0.2*(x3+log(abs(5*x5)+1)))+runif(n,-1,1) #x7a¡¦?x3a£á?x5a??e¡X?
  
  x=data.frame()
  for (i in 1:p){
    if (i==1)
      x[1:n,i]=x1
    else if (i==2)
      x[1:n,i]=x2
    else if (i==3)
      x[1:n,i]=x3
    else if (i==4)
      x[1:n,i]=x4
    else if (i==5)
      x[1:n,i]=x5
    else if (i==6)
      x[1:n,i]=x6
    else if (i==7)
      x[1:n,i]=x7
    else
      x[1:n,i]=rnorm(n,0,1)
  }
  x=matrix(unlist(x),ncol=p)

prop = exp(rowSums(x)) / (1 + exp(rowSums(x)))
X = x
y = rbinom(n,1,prop)


########################################################

pairs = NULL
thre = 0.1
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
pairs_chosen <- pairs_ord[1:35,] #pairs_chosen <- pairs_ord[1:50,]
colnames(pairs_chosen) <- c("X", "XI")
#pairs_chosen_end <- sqldf("
#                   SELECT * FROM pairs_chosen
#                   ORDER BY X
#                  ")

distinct_var <- NULL; seq_vars <- NULL
distinct_var <- pairs_chosen$X
seq_vars <- cbind((1:length(distinct_var)), distinct_var)

# 29 new variables forms X_new and p_new
#X_new = X
X_new <- X[,distinct_var] # (selected variables)
p_new <- length(distinct_var)

# compute corr(pairs) with sparsity
pair = NULL; XI = NULL; thre = 0.05
for(i in 1:p_new) {
  for(j in 1:p_new) {
    x = max(xicor(X_new[,i],X_new[,j]),xicor(X_new[,j],X_new[,i]))*1
    if(i<j && x > thre) {
      pair = rbind(pair,c(i,j))
      XI = c(XI,x)
    }
    
  }
}
pairs <- cbind(pair,XI)
dim(pairs)[1] ## xxpairs_linear.RData


# find out the overall XI 
pairs_new <- pairs
pairs_new2 <- cbind(pairs_new[,2],pairs_new[,1],pairs_new[,3])
XI <- diag(dim(seq_vars)[1])
pairs_end <- rbind(pairs_new,pairs_new2)
for(i in 1:(dim(pairs_end)[1])){	
	XI[pairs_end[i,1],pairs_end[i,2]] <- pairs_end[i,3]
}
XI  ## XI_linear.RData
XI[which(abs(XI)<0.25)]=0

#### step 2 ####
clime=clime(XI,lambda = 0.001,sigma = TRUE)
theta_hat <- matrix(unlist(clime$Omegalist),nrow = p_new)
theta_hat[which(abs(theta_hat)<0.45)]=0
   ## theta_hat_linear.RData


pi <- NULL
mu <- NULL
classes <- length(unique(y))
delta <- matrix(rep(0,classes*length(y)),ncol=length(y))

for(i in 1:classes){
	mu[[i]] <- colMeans(X_new[which(y==(i-1)),1:35])
}

for(i in 1:classes){
	pi <-  c(length(which(y == 0))/length(y),length(which(y == 1))/length(y))
	for(j in 1:length(y)){
		delta[i,j] <- log(pi[i]) - 0.5 * (matrix(mu[[i]],nrow=1) %*% theta_hat %*% matrix(mu[[i]],ncol=1)) + as.numeric(X_new[j, 1:35]) %*% theta_hat %*% matrix(mu[[i]],ncol=1)
	# x is the obs. ready to estimate
	}
}
dim(delta)


class_est <- NULL
for(j in 1:length(y)){
	class_est[j] <- which(delta[,j]==max(delta[,j]))
}
class_est <- class_est-1


comb <- cbind(y, class_est) # real classes vs. estimate classes
  ## comb_linear_train.RData
which(class_est!=y) # misclassification obs.
length(which(class_est!=y))

#### true positive, false positive, false negative ####
head(comb)
colnames(comb)

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

PRE_ave = c(PRE_ave,PRE_micro)
REC_ave = c(REC_ave,REC_micro)
F_ave = c(F_ave,F_micro)

SPE = 0 #e??0
SEN = 0 #0

for(i in 1:dim(theta_hat)[2]) {
  for(j in 1:dim(theta_hat)[2])   {
    if(strue[i,j]!=0 && theta_hat[i,j] !=0) {
      SPE = SPE + 1  }
    if(strue[i,j]==0 && theta_hat[i,j] ==0) {
      SEN = SEN + 1  }
    
  }
}




}

PRE_ave; REC_ave; F_ave



#n0
SPE / sum((strue[1:dim(theta_hat)[2],1:dim(theta_hat)[2]]!=0)*1)
#0
SEN / (dim(theta_hat)[2]*dim(theta_hat)[2] - sum((strue[1:dim(theta_hat)[2],1:dim(theta_hat)[2]]!=0)*1))
