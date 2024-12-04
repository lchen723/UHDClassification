
#### LOOK FOR XI_i for quadratic discrimination ####

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

#### precision matrix ####

y = as.factor(c(rep(1,8),rep(2,8),rep(3,8),rep(4,8),rep(5,16),rep(6,8),rep(7,8),rep(8,8),rep(9,24),rep(10,8),
                rep(11,8),rep(12,8),rep(13,8),rep(14,16)))

classes <- length(unique(y))
load("del_ind.RData")

# find out the homogeneous variables
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
pairs ## xypairs_quad.RData

pairs <- as.data.frame(pairs)


# label the homogeneous variables
pairs_ord <- sqldf("
                   SELECT * FROM pairs
                   ORDER BY XI DESC
                  ")
### change the dimension which will be greater than 29.
pairs_chosen <- pairs_ord[1:35,]
colnames(pairs_chosen) <- c("X", "XI")
pairs_chosen_end <- sqldf("
                   SELECT * FROM pairs_chosen
                   ORDER BY X
                  ")
distinct_var <- NULL; seq_vars <- NULL
distinct_var <- pairs_chosen_end$X
seq_vars <- cbind((1:length(distinct_var)), distinct_var)

# 29 new variables forms X_new and p_new
X_new <- X[,distinct_var]
p_new <- length(distinct_var)
# find out the XI of each class
y_ind = as.factor(c(rep(1,8),rep(2,8),rep(3,8),rep(4,8),rep(5,16),rep(6,8),rep(7,8),rep(8,8),rep(9,24),rep(10,8),
                rep(11,8),rep(12,8),rep(13,8),rep(14,16)))

pairs_i <- NULL
for(k in 1:classes){
  pair = NULL; XI = NULL
  for(i in 1:p_new) {
    for(j in 1:p_new) {
      x = max(xicor(X_new[which(y_ind==k),i],X_new[which(y_ind==k),j]),xicor(X_new[which(y_ind==k),j],X_new[which(y_ind==k),i]))*1
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
  X_end[[k]] <- X_new[which(y_ind==k),]
}

# We now get matrix: "XI_i" and "X_end" 


#### ICDA_RDA quadratic part ####

# We should still take linear part and lots of JGL functions into consideration.


#### Functions we changed ####
# 1. JGL

#### RUNNING STEPS ####
# 0. above code (to find "XI_i")
# 1. screen.fgl 
# 2. screen.ggl 
# 3. admm.iters.unconnected
# 4. penalty.as.matrix 
# 5. admm.iters
# 6. JGL (changed)  no.1~5 are prepared for no.6

quad_model <- JGL_new(Y=X_end,XI=XI, lambda1=0.1, lambda2=0.2)


pp <- quad_model$theta ## theta_hat_quad_update.RData
quad_model$diag.theta.unconnected
quad_model$connected
conn <- which(quad_model$connected==TRUE)

str(quad_model)
print.jgl(quad_model)

mu <- NULL
classes <- length(unique(y))
phi <- matrix(rep(0,classes*length(y)),ncol=length(y))

for(i in 1:classes){
  mu[[i]] <- colMeans(X_new[which(y==unique(y)[i]),])
}

for(i in 1:classes){
  pi[i] <- summary(y)[i]/length(y)
  for(j in 1:length(y)){
    phi[i,j] <- log(pi[i]) + 0.5 * log(det(pp[[i]])) - 0.5 * (X_new[j,conn] - (matrix(mu[[i]][conn],nrow=1))) %*% pp[[i]] %*% (t(t(X_new[j,conn])) - matrix(mu[[i]][conn],ncol=1))
    # x is the obs. ready to estimate
  }
}

class_est <- NULL
for(j in 1:length(y)){
  class_est[j] <- which(phi[,j]==max(phi[,j]))
}
class_est <- as.factor(class_est)
levels(class_est) = c("BR","PR","LU","CO","LY","BL","ML","UT","LE","RE","PA","OV","ME","CNS")

levels(y) = c("BR","PR","LU","CO","LY","BL","ML","UT","LE","RE","PA","OV","ME","CNS")
level_y = c("BR","PR","LU","CO","LY","BL","ML","UT","LE","RE","PA","OV","ME","CNS")

comb <- cbind(y, class_est) # real classes vs. estimate classes
  ## comb_quad_train.RData
which(class_est!=y) # misclassification obs.
length(which(class_est!=y))

#### statistics ####
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
PRE_macro <- mean(PRE)
REC_macro <- mean(REC)
F_macro <- 2 * (PRE_macro * REC_macro) / (PRE_macro + REC_macro)
PRE_micro; REC_micro; F_micro
PRE_macro; REC_macro; F_macro


# add at 20241005
pai_quad = NULL; tmp_pai_quad = NULL
for(k in 1:classes){
  for(i in 2:35){
    for(j in 1:(i-1)){
      tmp <- ifelse(abs(pp[[k]][i,j])==0, 0, 1)
      tmp_pai_quad <- rbind(tmp_pai_quad,c(j,i,tmp))
    }
  }
  pai_quad[[k]] <- tmp_pai_quad
  tmp_pai_quad = NULL
}
pai_quad[[1]][,3]-pai_quad[[5]][,3]

pai <- NULL
for(k in 1:classes){
  pai <- cbind(pai, pai_quad[[k]][,3])
}
same_ind <- which(rowSums(pai)==14) ; length(same_ind)
for(i in 1:595){
  cc <- ifelse(rowSums(pai)==14,1,0)
}
cc2 <- cbind(1:595, cc)

for(i in 1:classes){
  net = matrix(rep(1,35*35),ncol=35)
  #net = pp[[i]]
  net = network(net, directed = FALSE)
  network.vertex.names(net) = paste0("X",distinct_var)
  edges <- which(pai_quad[[i]][,3]==1) # all edges in quad (theta_hat) (ith class)
  net %e% "weight" <- ifelse((pai_quad[[i]][,3]==1 & cc==1),"lightgray",ifelse((pai_quad[[i]][,3]==1 & rowSums(pai)==1),"red","white"))
  del_ind <- which(net %e% "weight"=="white")
  delete.edges(net, del_ind)
  pl <- ggnet2(net,model="target",mode="circle",size=3,label.size=3,color="color",label=TRUE,node.color = "lightgray",edge.color="weight")
  plot(pl)
}
for(i in 1:classes){
  print(which(rowSums(pai)<14 & pai_quad[[i]][,3]==1))
}
which(rowSums(pai)<14 & pai_quad[[1]][,3]==1)

tmp5<- cbind(pai_quad[[1]][,3],pai_quad[[2]][,3],)

####
net = theta_hat
net = network(net, directed = FALSE)
network.vertex.names(net) = paste0("X",distinct_var)

linear_edges <- which(pai_linear[,3]==1) # all edges in linear (theta_hat)
clime_edges <- which(pai_clime[,3]==1) # all edge in clime 
# linear plot
net %e% "weight" <- ifelse((pai_linear[linear_edges,3]==1 & tmp2[linear_edges]==1), "red",ifelse(pai_clime[linear_edges,3]==1,"lightgray","white"))
ggnet2(net,model="target",mode="circle",size=3,label.size=3,color="color",label=TRUE,node.color = "lightgray",edge.color="weight")
# clime plot
del_ind <- which(ifelse((pai_linear[linear_edges,3]==1 & tmp2[linear_edges]==1), "red",ifelse(pai_clime[linear_edges,3]==1,"lightgray","white"))=="red")
delete.edges(net, del_ind)
ggnet2(net,model="target",mode="circle",size=3,label.size=3,color="color",label=TRUE,node.color = "lightgray",edge.color="weight")
####


pl <- ggnet2(net, model="target", mode="circle",size=3,label.size=3,color="color",label=TRUE,node.color = "lightgray")
plot(pl)
which((pai_quad[[1]][,3]-pai_quad[[12]][,3])!=0)


# plot precision matrix
for(k in 1:14){
  net = pp[[k]]
  net = network(net, directed = FALSE)
  network.vertex.names(net) = paste0("X",distinct_var)
  pl <- ggnet2(net, model="target", mode="circle",size=3,label.size=3,color="color",label=TRUE,node.color = "lightgray")
  plot(pl)
}

pl <- ggnet2(net,model="target",mode="circle",size=3,label.size=3,color="color",label=TRUE,node.color = "lightgray")

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

remove(comb)
