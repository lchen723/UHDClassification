
#### linear discriminent function ####
#### difference: 0. remove repeated variables (the vars that corr>=0.97)
####             1. compute corr(Y, each X) first
####             2. choose those ind. variables
####             3. compute pairs
set.seed(34)

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

X = t(gene)#

y = as.factor(c(rep(1,8),rep(2,8),rep(3,8),rep(4,8),rep(5,16),rep(6,8),rep(7,8),rep(8,8),rep(9,24),rep(10,8),
rep(11,8),rep(12,8),rep(13,8),rep(14,16)))

levels(y) = c("BR","PR","LU","CO","LY","BL","ML","UT","LE","RE","PA","OV","ME","CNS")

level_y = c("BR","PR","LU","CO","LY","BL","ML","UT","LE","RE","PA","OV","ME","CNS")

# del_ind <- pairs[which(pairs[,3]>=0.97),2]
# write.csv(del_ind,"del_ind.csv")
del_ind <- read.csv("del_ind.csv"); del_ind <- del_ind[,2]
# X <- X[,-del_ind]

###############################


############
## Step 1 ##
############

# 
# write.csv(pairs,"pairs.csv")
# pairs <- read.csv("pairs.csv"); pairs <- pairs[,-1]; pairs <- pairs[-which(pairs[,3]>=0.97),]
# pairs <- pairs[which(pairs[,3]>0.78),] 
# Here, "pairs" is the xicor between X_i and X_j. 


# find out the homogeneous variables
#### corr(y,each x) ####
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
pairs_chosen <- pairs_ord[1:35,] #pairs_chosen <- pairs_ord[1:50,]
colnames(pairs_chosen) <- c("X", "XI")
pairs_chosen_end <- sqldf("
                   SELECT * FROM pairs_chosen
                   ORDER BY X
                  ")
distinct_var <- NULL; seq_vars <- NULL
distinct_var <- pairs_chosen_end$X
seq_vars <- cbind((1:length(distinct_var)), distinct_var)

# 29 new variables forms X_new and p_new
X_new <- X[,distinct_var] # (selected variables)
p_new <- length(distinct_var)

# compute corr(pairs) with sparsity
pair = NULL; XI = NULL; thre = 0.1
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
pairs_new <- pairs # 594*3
pairs_new2 <- cbind(pairs_new[,2],pairs_new[,1],pairs_new[,3]) # 594*3
XI <- diag(dim(seq_vars)[1]) # 35*35
pairs_end <- rbind(pairs_new,pairs_new2) # 1188*3
for(i in 1:(dim(pairs_end)[1])){	
	XI[pairs_end[i,1],pairs_end[i,2]] <- pairs_end[i,3]
}
XI  ## XI_linear.RData

#### step 2 ####
theta_hat <- glasso(XI,rho=.01)$wi # 35*35
   ## theta_hat_linear.RData

#### step 3 ####
pi <- NULL
mu <- NULL
classes <- length(unique(y))
delta <- matrix(rep(0,classes*length(y)),ncol=length(y)) # 14*144

for(i in 1:classes){
	mu[[i]] <- colMeans(X_new[which(y==unique(y)[i]),])
}  # compute the mu of each class

for(i in 1:classes){
	pi[i] <- summary(y)[i]/length(y)  # number of each class/ total sample size
	for(j in 1:length(y)){
		delta[i,j] <- log(pi[i]) - 0.5 * (matrix(mu[[i]],nrow=1) %*% theta_hat %*% matrix(mu[[i]],ncol=1)) + X_new[j,] %*% theta_hat %*% matrix(mu[[i]],ncol=1)
	# x is the obs. ready to estimate
	}
}
dim(delta)


class_est <- NULL
for(j in 1:length(y)){
	class_est[j] <- which(delta[,j]==max(delta[,j])) 
}  # choose the max delta among classes for all samples (the class with highest probability)

class_est <- as.factor(class_est)
levels(class_est) = c("BR","PR","LU","CO","LY","BL","ML","UT","LE","RE","PA","OV","ME","CNS")

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
} # values above for each class

PRE_micro <- sum(TP)/sum(TP+FP) # values above for overall samples
REC_micro <- sum(TP)/sum(TP+FN)
F_micro <- 2 * (PRE_micro * REC_micro) / (PRE_micro + REC_micro)
PRE_macro <- mean(PRE)
REC_macro <- mean(REC)
F_macro <- 2 * (PRE_macro * REC_macro) / (PRE_macro + REC_macro)
PRE_micro; REC_micro; F_micro
PRE_macro; REC_macro; F_macro

# plot precision matrix
net = theta_hat
net = network(net, directed = FALSE)
network.vertex.names(net) = paste0("X",distinct_var)
ggnet2(net, model="target",mode="circle",size=3,label.size=3,color="color",label=TRUE,node.color = "lightgray")


# heat map
Y <- c("BR","PR","LE","CO","LU","BL","CNS","UT","LY","RE","PA","OV","ME","ML")
Fitted <- paste0(Y)
data <- cbind(expand.grid(Fitted=Fitted, Y=Y),expand.grid(Fitted2=1:14, Y2=1:14))

comb <- as.data.frame(comb)
tmp <- sqldf("
              SELECT class_est, count(y) as total_count
              FROM comb
              GROUP BY class_est
            ")
tmp2 <- sqldf("
              SELECT class_est, y, count() as group_count
              FROM comb
              GROUP BY class_est, y
             ")
tmp3 <- sqldf("
              SELECT tmp2.class_est as class_est, tmp2.y as y, 
                     tmp2.group_count as group_count, tmp.total_count 
              FROM tmp2
              LEFT JOIN tmp
              ON tmp2.class_est = tmp.class_est
              ")
tmp3$ratio <- tmp3$group_count/tmp3$total_count 

data_new <- sqldf("
                  SELECT data.Fitted as Fitted, data.Y as Y, tmp3.ratio as ratio 
                  FROM data
                  LEFT JOIN tmp3
                  ON data.Fitted2 = tmp3.class_est and data.Y2 = tmp3.y
                  ")
data_new[which(is.na(data_new$ratio)),3] <- 0

head(data_new)

ggplot(data_new, aes(Y, Fitted, fill= ratio)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="blue")
