###################
## Yeast dataset ##
##  BR benchmark ##
###################

packages=c("mldr","stats","ROSE","gridExtra","Metrics")
lapply(packages,FUN=library,character.only=TRUE)

#####Obtain clusters#####
yeast=read.csv('yeast.csv',header = TRUE)
label.yeast=yeast[,104:117] # ignore the last class

time0.cluster=Sys.time()
D=1-cor(label.yeast)

#hierarchical clustering
hc.complete <-hclust(as.dist(D),method="complete")
plot(hc.complete)
result.clusters <- cutree(hc.complete,k=8);result.clusters



####create label subsets####
n=dim(label.yeast)[1]
x.matrix <- as.data.frame(yeast[,1:103])
x.matrix <- apply(x.matrix, 2, unlist) #1 indicates rows, unlist makes rows into single vectors
x.matrix=apply(x.matrix,2,as.numeric)
y <- as.data.frame(label.yeast)
label.subset <- matrix(0,nrow=n,ncol=8)
for (i in 1:n){
  if(label.yeast[i,1]==1 && label.yeast[i,2]==1){
    label.subset[i,1] = 1
  } 
  if(label.yeast[i,3]==1&&label.yeast[i,4]==1){
    label.subset[i,2] = 1
  } 
  if(label.yeast[i,5]==1 &&label.yeast[i,6]==1){
    label.subset[i,3] = 1
  } 
  if(label.yeast[i,7]==1 && label.yeast[i,8]==1){
    label.subset[i,4] = 1
  } 
  if(label.yeast[i,9]==1 ){ 
    label.subset[i,5] = 1
  } 
  if(label.yeast[i,10]==1 && label.yeast[i,11]==1 ){
    label.subset[i,6] = 1
  } 
  if (label.yeast[i,12]==1 && label.yeast[i,13]==1){
    label.subset[i,7] = 1
  } 
  if(label.yeast[i,14]==1){ 
    label.subset[i,8] = 1
  }
}
freq=apply(label.subset,2,sum)
freq=freq/n
names(label.subset)<-c("labelset1","labelset2","labelset3","labelset4","labelset5","labelset6",
                       "labelset7","labelset8")
label.index <- list(c(1,2),c(3,4),c(5,6),c(7,8),c(9),c(10,11),c(12,13),c(14))
time.total.cluster=Sys.time()-time0.cluster

####cross-validation####


#######################################
##create index subsets for 10-fold CV##
#######################################

##Note: this functoin, fn.CV, is created s.t. each fold contains 1/K obs from the true 1's from the given label subset
##use this function to realize 10-fold CV for predicting each label subset
fn.CV <- function(label.subset)
{
  K <- 10 #10-fold CV as an example
  n.fold <- floor(n/K) #size of each fold, rounded down to the nearest integer (so the last fold might be larger in size)
  index.fold <- list()
  
  id.labelsubset <- sample(which(label.subset==1),length(which(label.subset==1)),replace=FALSE) #id's for 1's of the given label subset, already shuffled
  n.labelsubset.fold <- floor(length(id.labelsubset)/K) #number of indexes in each fold that need to be taken from the true 1's of the given label subset
  id.remain <- sample(setdiff(1:n, id.labelsubset),n-length(id.labelsubset),replace=FALSE) #shuffle the remaining id's that correspond to 0's for the label subset
  n.remain <- n.fold-n.labelsubset.fold
  
  for(i in 1:K)
  { 
    if(i<K)
    {
      index.fold[[i]] <- c(id.labelsubset[((i-1)*n.labelsubset.fold+1):(i*n.labelsubset.fold)],id.remain[((i-1)*n.remain+1):(i*n.remain)])
    }else
    {
      index.fold[[i]] <- c(id.labelsubset[((K-1)*n.labelsubset.fold+1):length(id.labelsubset)],id.remain[((K-1)*n.remain+1):length(id.remain)])
    }
  }
  
  return(index.fold)
}





####benchmark####
pred.label <- matrix(0,nrow=n,ncol=dim(y)[2])

accuracy =  recall = fmeasure = precision = 0
K = 10 #use 10-fold CV
th.fix = 0.5 #threshold

time0 <- Sys.time()

for(j in 1:dim(label.yeast)[2])
{
  set.seed(j) #control randomness
  index.fold = fn.CV(label.yeast[,j]) #create the 10-fold index sets for label subset j
  
  for(k in 1:K)
  {
    mydata_j <- as.data.frame(cbind(label.yeast[,j],x.matrix))
    names(mydata_j)[1] <- "V1"
    
    train.id=c()
    for (t in setdiff(1:10,k)){
      train.id=c(train.id,index.fold[[t]])
    }
    
    fit <- glm(V1~.,data=mydata_j[train.id,],family=binomial)
    pred <- predict(fit,newdata = mydata_j[index.fold[[k]],],type="response")
    Y.hat <- ifelse(pred>th.fix,1,0)
    id1.set <- which(Y.hat==1)
    pred.label[index.fold[[k]][id1.set],j] <- 1 
    
  }
  
  
  print(k) #print the curent fold to track the computation 
}

pred.label.meta=matrix(0,nrow=n,ncol=dim(y)[2])

for(j in 1:dim(label.yeast)[2])
{
  set.seed(j) #control randomness
  index.fold = fn.CV(label.yeast[,j]) #create the 10-fold index sets for label subset j
  
  for(k in 1:K)
  {
    mydata_j <- as.data.frame(cbind(label.yeast[,j],pred.label))
    names(mydata_j)[1] <- "V1"
    
    train.id=c()
    for (t in setdiff(1:10,k)){
      train.id=c(train.id,index.fold[[t]])
    }
    
    fit <- glm(V1~.,data=mydata_j[train.id,],family=binomial)
    pred <- predict(fit,newdata = mydata_j[index.fold[[k]],],type="response")
    Y.hat <- ifelse(pred>th.fix,1,0)
    id1.set <- which(Y.hat==1)
    pred.label.meta[index.fold[[k]][id1.set],j] <- 1 
    
  }
  
  
  print(k) #print the curent fold to track progress 
}

time.total <- Sys.time() - time0

pred.label=pred.label.meta

accuracy = mldr::accuracy(label.yeast,pred.label,undefined_value="diagnose")
precision = mldr::precision(label.yeast,pred.label,undefined_value="diagnose")
recall = mldr::recall(label.yeast,pred.label,undefined_value="diagnose")
fmeasure = fmeasure(label.yeast,pred.label,undefined_value="diagnose")
HL = hamming_loss(label.yeast,pred.label)

accuracy;precision;recall;fmeasure;HL

mac_accuracy = mean ( apply(label.yeast == pred.label, 2, sum)/dim(label.yeast)[1] )
mac_precision = macro_precision(label.yeast,pred.label)
mac_recall = macro_recall(label.yeast,pred.label)
mac_fmeasure = macro_fmeasure(label.yeast,pred.label)

mic_accuracy = mean ( apply(label.yeast == pred.label, 1, sum)/dim(label.yeast)[2] )
mic_precision = micro_precision(label.yeast,pred.label)
mic_recall = micro_recall(label.yeast,pred.label)
mic_fmeasure = micro_fmeasure(label.yeast,pred.label)

mac_accuracy;mac_precision;mac_recall;mac_fmeasure
mic_accuracy;mic_precision;mic_recall;mic_fmeasure

