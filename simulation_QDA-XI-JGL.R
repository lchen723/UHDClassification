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

n=200
p=10000
PRE_ave = NULL; REC_ave = NULL; F_ave = NULL

strue=diag(1,p,p)
strue[2,1]=strue[1,2]=1
strue[3,1]=strue[1,3]=1;strue[3,2]=strue[2,3]=1
strue[4,3]=strue[3,4]=1;strue[4,6]=strue[6,4]=1
strue[5,2]=strue[2,5]=1;strue[5,6]=strue[6,5]=1
strue[7,3]=strue[3,7]=1;strue[7,5]=strue[5,7]=1

#repeat 100 times
#Theta_ng=array(dim = c(p,p,100))
#for (q in 1:100){
 # cat("\r",round(q/100*100,2), '%     ')
  
  #e3?a!V?a?!An*pc??c?ce?¢FG
  x1=runif(n,-1,1)
  
  x2=6*cos(x1)+runif(n,-1,1) #x2a!|?x1a??e!X?
  x3=5*sin(x1)+x2+rnorm(n,0,1) #x3a!|?x1a¢Ga?x2a??e!X?
  
  x6=runif(n,-1,1)
  
  x4=5*cos(x3*x6)+3*x3+3*x6+rnorm(n,0,1) #x4a!|?x3a¢Ga?x6a??e!X?
  x5=0.05*(x2+x6)^3+rnorm(n,0,1) #x5a!|?x2a¢Ga?x6a??e!X?
  x7=6*cos(0.2*(x3+log(abs(5*x5)+1)))+runif(n,-1,1) #x7a!|?x3a¢Ga?x5a??e!X?
  
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
y = y + 1


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
pairs ## xypairs_quad.RData

pairs <- as.data.frame(pairs)

pairs_ord <- sqldf("
                   SELECT * FROM pairs
                   ORDER BY XI DESC
                  ")
### change the dimension which will be greater than 29.
# pairs_chosen <- pairs_ord[1:ceiling(n/log(n)),]
pairs_chosen <- pairs_ord[1:35,]
colnames(pairs_chosen) <- c("X", "XI")
pairs_chosen_end <- sqldf("
                   SELECT * FROM pairs_chosen
                   ORDER BY X
                  ")
distinct_var <- NULL; seq_vars <- NULL
distinct_var <- pairs_chosen_end$X
seq_vars <- cbind((1:length(distinct_var)), distinct_var)

X_new <- X[,distinct_var]
p_new <- length(distinct_var)
classes <- length(unique(y))
pairs_i <- NULL
for(k in 1:classes){
  pair = NULL; XI = NULL
  for(i in 1:p_new) {
    for(j in 1:p_new) {
      x = max(xicor(X_new[which(y==k),i],X_new[which(y==k),j]),xicor(X_new[which(y==k),j],X_new[which(y==k),i]))*1
      if(i>j && x > 0) {
        pair = rbind(pair,c(i,j))
        XI = c(XI,x)
      }
    }
  }
  pairs_i[[k]] <- cbind(pair,XI)
} 
## xxpairs_i_quad_update.RData

pairs_new <- NULL; pairs_new2 <- NULL
XI <- NULL; pairs_end <- NULL



for(k in 1:classes){
  pairs_new[[k]] <- pairs_i[[k]][(which(pairs_i[[k]][,3] > 0.1)),] 
  pairs_new2[[k]] <- cbind(pairs_new[[k]][,2],pairs_new[[k]][,1],pairs_new[[k]][,3])
  XI[[k]] <- diag(dim(seq_vars)[1])
  pairs_end[[k]] <- rbind(pairs_new[[k]],pairs_new2[[k]])
  for(i in 1:(dim(pairs_end[[k]])[1])){	
    XI[[k]][pairs_end[[k]][i,1],pairs_end[[k]][i,2]] <- pairs_end[[k]][i,3]
  }
}
## XI_quad_update.RData


X_end <- NULL
for(k in 1:classes){
  X_end[[k]] <- X_new[which(y==k),]
}

###############


quad_model <- JGL_new(Y=X_end,XI=XI, lambda1=0.001, lambda2=0.002)

pp <- quad_model$theta ## theta_hat_quad.RData
quad_model$diag.theta.unconnected
quad_model$connected
conn <- which(quad_model$connected==TRUE)


mu <- NULL
classes <- length(unique(y))
phi <- matrix(rep(0,classes*length(y)),ncol=length(y))

for(i in 1:classes){
  mu[[i]] <- colMeans(X_new[which(y==unique(y)[i]),])
}

for(i in 1:classes){
  pi[i] <- summary(y)[i]/length(y)
  for(j in 1:length(y)){
#    phi[i,j] <- log(pi[i]) + 0.5 * log(det(XI[[i]])) - 0.5 * (X_new[j,conn] - (matrix(mu[[i]][conn],nrow=1))) %*% XI[[i]] %*% (t(t(X_new[j,conn])) - matrix(mu[[i]][conn],ncol=1))
    phi[i,j] <- log(pi[i]) + 0.5 * log(det(quad_model$theta[[i]])) - 0.5 * (X_new[j,conn] - (matrix(mu[[i]][conn],nrow=1))) %*% quad_model$theta[[i]] %*% (t(t(X_new[j,conn])) - matrix(mu[[i]][conn],ncol=1))
    # x is the obs. ready to estimate
  }
}


class_est <- NULL
for(j in 1:length(y)){
  class_est[j] <- which(phi[,j]==max(phi[,j]))
}
class_est <- as.factor(class_est)


comb <- cbind(y, class_est) # real classes vs. estimate classes
  ## comb_quad_train.RData
which(class_est!=y) # misclassification obs.
length(which(class_est!=y))

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
PRE_micro; REC_micro; F_micro

SPE = 0 #e??0
SEN = 0 #0
theta_hat = quad_model$theta[[1]]
for(i in 1:dim(theta_hat)[2]) {
  for(j in 1:dim(theta_hat)[2])   {
    if(strue[i,j]!=0 && theta_hat[i,j] !=0) {
      SPE = SPE + 1  }
    if(strue[i,j]==0 && theta_hat[i,j] ==0) {
      SEN = SEN + 1  }
    }
      }


#n0
SPE / sum((strue[1:dim(theta_hat)[2],1:dim(theta_hat)[2]]!=0)*1)
#0
SEN / (dim(theta_hat)[2]*dim(theta_hat)[2] - sum((strue[1:dim(theta_hat)[2],1:dim(theta_hat)[2]]!=0)*1))



