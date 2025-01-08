#############################
##      Yeast dataset      ##
## resampled training sets ##
#############################

packages=c("mldr","mldr.datasets","stats","ROSE","gridExtra","Metrics")
lapply(packages,FUN=library,character.only=TRUE)

#####Obtain clusters#####
yeast=yeast()
label.yeast=yeast$dataset[,104:117] 
D=1-cor(label.yeast)

#hierarchical clustering
hc.complete <-hclust(as.dist(D),method="complete")
plot(hc.complete)
result.clusters <- cutree(hc.complete,k=8);result.clusters



####create label subsets####
n=dim(label.yeast)[1]
x.matrix <- as.data.frame(yeast$dataset[,1:103])
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
freq=apply(label.subset,2,sum);freq
freq=freq/n
names(label.subset)<-c("labelset1","labelset2","labelset3","labelset4","labelset5","labelset6",
                       "labelset7","labelset8")
label.index <- list(c(1,2),c(3,4),c(5,6),c(7,8),c(9),c(10,11),c(12,13),c(14))


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



####resampling threshold####
pred.label.under <- matrix(0,nrow=n,ncol=dim(label.yeast)[2])
pred.label.over <- matrix(0,nrow=n,ncol=dim(label.yeast)[2])
pred.label.smote <- matrix(0,nrow=n,ncol=dim(label.yeast)[2])


accuracy.under =  recall.under = fmeasure.under = precision.under =
accuracy.over =  recall.over = fmeasure.over = precision.over =
accuracy.smote =  recall.smote = fmeasure.smote = precision.smote = 0


K = 10 #use 10-fold CV
th.fix = 0.5 #threshold

time0 <- Sys.time()
for(k in 1:K)
{
  
  for(j in 1:dim(label.subset)[2])
  {
    set.seed(j) #control randomness
    index.fold = fn.CV(label.subset[,j])
    
    if(length(label.index[[j]])>1) #if the labelset contains more than one label
    {
      mydata_j <- as.data.frame(cbind(label.subset[,j],x.matrix))
      names(mydata_j)[1] <- "V1"
      #undersample: 
      mydata_j.under <- ovun.sample(as.factor(V1)~., data=mydata_j[-index.fold[[k]],],method="under")$data
      #oversample: 
      mydata_j.over <- ovun.sample(as.factor(V1)~., data=mydata_j[-index.fold[[k]],],method="over")$data
      #smote:
      mydata_j.smote <- ROSE(as.factor(V1)~.,data=mydata_j[-index.fold[[k]],],seed=123)$data
      
      fit.under <- glm(V1~.,data=mydata_j.under,family=binomial)
      fit.over <- glm(V1~.,data=mydata_j.over,family=binomial)
      fit.smote <- glm(V1~.,data=mydata_j.smote,family=binomial)
      
      
      pred.under <- predict(fit.under,newdata = mydata_j[index.fold[[k]],],type="response") #type=response returns the predicted probabilities of having this label subset
      pred.over <- predict(fit.over,newdata = mydata_j[index.fold[[k]],],type="response")
      pred.smote <- predict(fit.smote,newdata = mydata_j[index.fold[[k]],],type="response")
      
      Y.hat.under <- ifelse(pred.under>th.fix,1,0)
      Y.hat.over <- ifelse(pred.over>th.fix,1,0)
      Y.hat.smote <- ifelse(pred.smote>th.fix,1,0)
      
      id1.set.under <- which(Y.hat.under==1)
      id1.set.over <- which(Y.hat.over==1)
      id1.set.smote <- which(Y.hat.smote==1)
      
      
      pred.label.under[index.fold[[k]][id1.set.under],label.index[[j]]] <- 1 
      pred.label.over[index.fold[[k]][id1.set.over],label.index[[j]]] <- 1
      pred.label.smote[index.fold[[k]][id1.set.smote],label.index[[j]]] <- 1
      
      id0.set.under <- which(Y.hat.under !=1)
      id0.set.over <- which(Y.hat.over !=1)
      id0.set.smote <- which(Y.hat.smote !=1)
      
      for(i in 1:length(label.index[[j]]))
      {
        mydata.temp <- as.data.frame(cbind(y[,label.index[[j]][i]],x.matrix)) 
        names(mydata.temp)[1] <- "V1"
        temp <- glm(V1~.,data=mydata.temp[-index.fold[[k]],],family=binomial)
        
        pred.temp.under <- predict(temp,newdata=mydata.temp[index.fold[[k]][id0.set.under],],type="response")
        pred.temp.over <- predict(temp,newdata=mydata.temp[index.fold[[k]][id0.set.over],],type="response")
        pred.temp.smote <- predict(temp,newdata=mydata.temp[index.fold[[k]][id0.set.smote],],type="response")
        
        Y.hat.temp.under <- ifelse(pred.temp.under>th.fix,1,0)   
        Y.hat.temp.over <- ifelse(pred.temp.over>th.fix,1,0)
        Y.hat.temp.smote <- ifelse(pred.temp.smote>th.fix,1,0)
        
        pred.label.under[index.fold[[k]][id0.set.under[which(Y.hat.temp.under==1)]],label.index[[j]][i]] <- 1 
        pred.label.over[index.fold[[k]][id0.set.over[which(Y.hat.temp.over==1)]],label.index[[j]][i]] <- 1 
        pred.label.smote[index.fold[[k]][id0.set.smote[which(Y.hat.temp.smote==1)]],label.index[[j]][i]] <- 1 
        
      }
      #}
    }else #if the labelset is a singleton
    {
      mydata_j <- as.data.frame(cbind(label.subset[,j],x.matrix))
      names(mydata_j)[1] <- "V1"
      fit <- glm(V1~.,data=mydata_j[-index.fold[[k]],],family=binomial)
      pred <- predict(fit,newdata = mydata_j[index.fold[[k]],],type="response")
      Y.hat <- ifelse(pred>th.fix,1,0)
      id1.set <- which(Y.hat==1)
      pred.label.under[index.fold[[k]][id1.set],label.index[[j]]] = 1
        pred.label.over[index.fold[[k]][id1.set],label.index[[j]]]= 1
        pred.label.smote[index.fold[[k]][id1.set],label.index[[j]]] <- 1 
    }
    
  }
  
  
  print(k) #print the curent fold to track progress 
}

time.total <- Sys.time() - time0

#undersample metrics
accuracy.under = mldr::accuracy(label.yeast,pred.label.under,undefined_value="diagnose")
precision.under = mldr::precision(label.yeast,pred.label.under,undefined_value="diagnose")
recall.under = mldr::recall(label.yeast,pred.label.under,undefined_value="diagnose")
fmeasure.under = fmeasure(label.yeast,pred.label.under,undefined_value="diagnose")
HL.under = hamming_loss(label.yeast,pred.label.under)
accuracy.under;precision.under;recall.under;fmeasure.under;HL.under

mac_accuracy.under = mean ( apply(label.yeast == pred.label.under, 2, sum)/dim(label.yeast)[1] )
mac_precision.under = macro_precision(label.yeast,pred.label.under)
mac_recall.under = macro_recall(label.yeast,pred.label.under)
mac_fmeasure.under = macro_fmeasure(label.yeast,pred.label.under)

mic_accuracy.under = mean ( apply(label.yeast == pred.label.under, 1, sum)/dim(label.yeast)[2] )
mic_precision.under = micro_precision(label.yeast,pred.label.under)
mic_recall.under = micro_recall(label.yeast,pred.label.under)
mic_fmeasure.under = micro_fmeasure(label.yeast,pred.label.under)

mac_accuracy.under;mac_precision.under;mac_recall.under;mac_fmeasure.under
mic_accuracy.under;mic_precision.under;mic_recall.under;mic_fmeasure.under

#oversample metrics
accuracy.over = mldr::accuracy(label.yeast,pred.label.over,undefined_value="diagnose")
precision.over = mldr::precision(label.yeast,pred.label.over,undefined_value="diagnose")
recall.over = mldr::recall(label.yeast,pred.label.over,undefined_value="diagnose")
fmeasure.over = fmeasure(label.yeast,pred.label.over,undefined_value="diagnose")
HL.over = hamming_loss(label.yeast,pred.label.over)
accuracy.over;precision.over;recall.over;fmeasure.over;HL.over

mac_accuracy.over = mean ( apply(label.yeast == pred.label.over, 2, sum)/dim(label.yeast)[1] )
mac_precision.over = macro_precision(label.yeast,pred.label.over)
mac_recall.over = macro_recall(label.yeast,pred.label.over)
mac_fmeasure.over = macro_fmeasure(label.yeast,pred.label.over)

mic_accuracy.over = mean ( apply(label.yeast == pred.label.over, 1, sum)/dim(label.yeast)[2] )
mic_precision.over = micro_precision(label.yeast,pred.label.over)
mic_recall.over = micro_recall(label.yeast,pred.label.over)
mic_fmeasure.over = micro_fmeasure(label.yeast,pred.label.over)

mac_accuracy.over;mac_precision.over;mac_recall.over;mac_fmeasure.over
mic_accuracy.over;mic_precision.over;mic_recall.over;mic_fmeasure.over

#smote-sample metrics
accuracy.smote = mldr::accuracy(label.yeast,pred.label.smote,undefined_value="diagnose")
precision.smote = mldr::precision(label.yeast,pred.label.smote,undefined_value="diagnose")
recall.smote = mldr::recall(label.yeast,pred.label.smote,undefined_value="diagnose")
fmeasure.smote = fmeasure(label.yeast,pred.label.smote,undefined_value="diagnose")
HL.smote = hamming_loss(label.yeast,pred.label.smote)
accuracy.smote;precision.smote;recall.smote;fmeasure.smote;HL.smote

mac_accuracy.smote = mean ( apply(label.yeast == pred.label.smote, 2, sum)/dim(label.yeast)[1] )
mac_precision.smote = macro_precision(label.yeast,pred.label.smote)
mac_recall.smote = macro_recall(label.yeast,pred.label.smote)
mac_fmeasure.smote = macro_fmeasure(label.yeast,pred.label.smote)

mic_accuracy.smote = mean ( apply(label.yeast == pred.label.smote, 1, sum)/dim(label.yeast)[2] )
mic_precision.smote = micro_precision(label.yeast,pred.label.smote)
mic_recall.smote = micro_recall(label.yeast,pred.label.smote)
mic_fmeasure.smote = micro_fmeasure(label.yeast,pred.label.smote)

mac_accuracy.smote;mac_precision.smote;mac_recall.smote;mac_fmeasure.smote
mic_accuracy.smote;mic_precision.smote;mic_recall.smote;mic_fmeasure.smote

