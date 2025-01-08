##################
## Yeast dataset##
## CV threshold ##
##################

packages=c("mldr","mldr.datasets","stats","ROSE","gridExtra","Metrics")
lapply(packages,FUN=library,character.only=TRUE)

#####Obtain clusters#####
yeast=yeast()
label.yeast=yeast$dataset[,104:117] # ignore the last class
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
freq=apply(label.subset,2,sum)
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

###################################
##create a function to compute HL## 
## for a given threshold via CV  ##
###################################
K <- 10 #10-fold CV

HL_labelset = function(data_labelset, th)
{
  HL = 0 #initialize value of Hamming Loss
  index.fold <- fn.CV(data_labelset[,1]) #data_labelset[,1] is the 0/1 column for the given label subset
  names(data_labelset[,1]) <- "V1"
  pred.Y <- rep(0,dim(data_labelset)[1])
  
  for(i in 1:K) #10 is the num of folds in CV
  {
    fit = glm(V1~., data=data_labelset[-index.fold[[i]],]) #fit a logistic regression on delete-i-fold training set for the given label set
    pred.prob = predict(fit, newdata = data_labelset[index.fold[[i]],], type="response") #predict on the removed ith fold
    pred.Y[index.fold[[i]]] = ifelse(pred.prob>th, 1, 0)
  }
  HL <- hamming_loss(data_labelset[,1],pred.Y)
  return(HL)
}

F_labelset = function(data_labelset,th)
{
  F_measure = 0 #initialize value of Hamming Loss
  index.fold <- fn.CV(data_labelset[,1]) #data_labelset[,1] is the 0/1 column for the given label subset
  names(data_labelset[,1]) <- "V1"
  pred.Y <- rep(0,dim(data_labelset)[1])
  for(i in 1:K) #10 is the num of folds in CV
  {
    fit = glm(V1~., data=data_labelset[-index.fold[[i]],]) #fit a logistic regression on delete-i-fold training set for the given label set
    pred.prob = predict(fit, newdata = data_labelset[index.fold[[i]],], type="response") #predict on the removed ith fold
    pred.Y[index.fold[[i]]] = ifelse(pred.prob>th, 1, 0)
    
  }
  prec <- length(which((data_labelset[,1]+pred.Y)==2))/sum(pred.Y)
  rec <- length(which((data_labelset[,1]+pred.Y)==2))/sum(data_labelset[,1])
  F_measure <- 2*prec*rec/(prec+rec)
  
  return(F_measure)
}

#obtain the best threshold for each subset
threshold = seq(0.05,0.75,0.05) #sequence of thresholds to choose from
cv_threshold.F = cv_threshold.HL=rep(0,dim(label.subset)[2])

for(j in 1:dim(label.subset)[2])
{ # iterate over all label subsets
  
  if(length(label.index[[j]])>1)
  { 
    # check if the subset is a singleton
    mydata_j <- as.data.frame(cbind(label.subset[,j],x.matrix))
    
    vec_F=vec_HL = rep(0,length(threshold))
    for (i in 1:length(threshold))
    { 
      vec_HL[i]=HL_labelset(mydata_j,threshold[i])
      vec_F[i]=F_labelset(mydata_j,threshold[i])
    }
    thres_index.HL=which.min(vec_HL) #minimize HL
    thres_index.F=which.max(vec_F) #maximize F-measure
    cv_threshold.HL[j] = threshold[thres_index.HL]
    cv_threshold.F[j] = threshold[thres_index.F]
  }                                #end: check if the subset is a singleton
  else
  { 
    cv_threshold.F[j] = cv_threshold.HL[j]=0.5
  }
  
} # end: iterate over all label subset
cv_threshold.F;cv_threshold.HL
cv.threshold.F = cv_threshold.F
cv.threshold.HL = cv_threshold.HL


####cross-validated threshold####
pred.label.F =  pred.label.HL <- matrix(0,nrow=n,ncol=14)
#oversample <-

accuracy.F =  recall.F = fmeasure.F = precision.F = 0
accuracy.HL =  recall.HL = fmeasure.HL = precision.HL = 0

th.fix = 0.5 #threshold

time0 <- Sys.time()
for(k in 1:K)
{
  
  for(j in 1:dim(label.subset)[2])
  {
    set.seed(j) #control randomness
    index.fold = fn.CV(label.subset[,j]) #create the 10-fold index sets for label subset j
    th.F.j <- cv.threshold.F[j] #CV threshold for label subset j based on opt F-measure
    th.HL.j <- cv.threshold.HL[j] #CV threshold for label subset j based on opt HL
    
    if(length(label.index[[j]])>1) #if the labelset contains more than one label
    {
      mydata_j <- as.data.frame(cbind(label.subset[,j],x.matrix))
      names(mydata_j)[1] <- "V1"
      fit <- glm(V1~.,data=mydata_j[-index.fold[[k]],],family=binomial)
      pred <- predict(fit,newdata = mydata_j[index.fold[[k]],],type="response") #type=response returns the predicted probabilities of having this label subset
      
      
      Y.hat.F <- ifelse(pred>th.F.j,1,0)
      Y.hat.HL <- ifelse(pred>th.HL.j,1,0)
      
      id1.set.F <- which(Y.hat.F==1)
      pred.label.F[index.fold[[k]][id1.set.F],label.index[[j]]] <- 1 
      id1.set.HL <- which(Y.hat.HL==1)
      pred.label.HL[index.fold[[k]][id1.set.HL],label.index[[j]]] <- 1
      
      
      id0.set.F <- which(Y.hat.F!=1)
      id0.set.HL <- which(Y.hat.HL!=1)
      
      for(i in 1:length(label.index[[j]]))
      {
        mydata.temp <- as.data.frame(cbind(y[,label.index[[j]][i]],x.matrix)) 
        names(mydata.temp)[1] <- "V1"
        temp <- glm(V1~.,data=mydata.temp[-index.fold[[k]],],family=binomial)
        
        if (length(id0.set.F) != 0){
        pred.temp.F <- predict(temp,newdata=mydata.temp[index.fold[[k]][id0.set.F],],type="response")
        Y.hat.temp.F <- ifelse(pred.temp.F>th.fix,1,0) 
        pred.label.F[index.fold[[k]][id0.set.F[which(Y.hat.temp.F==1)]],label.index[[j]][i]] <- 1 
        }
        
        if (length(id0.set.HL) != 0){
        pred.temp.HL <- predict(temp,newdata=mydata.temp[index.fold[[k]][id0.set.HL],],type="response")
        Y.hat.temp.HL <- ifelse(pred.temp.HL>th.fix,1,0) 
        pred.label.HL[index.fold[[k]][id0.set.HL[which(Y.hat.temp.HL==1)]],label.index[[j]][i]] <- 1
        }
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
      pred.label.F[index.fold[[k]][id1.set],label.index[[j]]] <- 1 
      pred.label.HL[index.fold[[k]][id1.set],label.index[[j]]] <- 1 
    }
    
  }
  
  
  print(k) #print the curent fold to track progress
}

time.total <- Sys.time() - time0

accuracy.F = mldr::accuracy(label.yeast,pred.label.F,undefined_value="diagnose")
precision.F = mldr::precision(label.yeast,pred.label.F,undefined_value="diagnose")
recall.F = mldr::recall(label.yeast,pred.label.F,undefined_value="diagnose")
fmeasure.F = fmeasure(label.yeast,pred.label.F,undefined_value="diagnose");fmeasure.F
HL.F = hamming_loss(label.yeast,pred.label.F)

#CV by optimizing F-measure
accuracy.F;precision.F;recall.F;fmeasure.F;HL.F

mac_accuracy.F = mean ( apply(label.yeast == pred.label.F, 2, sum)/dim(label.yeast)[1] )
mac_precision.F = macro_precision(label.yeast,pred.label.F)
mac_recall.F = macro_recall(label.yeast,pred.label.F)
mac_fmeasure.F = macro_fmeasure(label.yeast,pred.label.F)

mic_accuracy.F = mean ( apply(label.yeast == pred.label.F, 1, sum)/dim(label.yeast)[2] )
mic_precision.F = micro_precision(label.yeast,pred.label.F)
mic_recall.F = micro_recall(label.yeast,pred.label.F)
mic_fmeasure.F = micro_fmeasure(label.yeast,pred.label.F)

mac_accuracy.F;mac_precision.F;mac_recall.F;mac_fmeasure.F
mic_accuracy.F;mic_precision.F;mic_recall.F;mic_fmeasure.F


accuracy.HL = mldr::accuracy(label.yeast,pred.label.HL,undefined_value="diagnose")
#do the same for the rest of the measures
precision.HL = mldr::precision(label.yeast,pred.label.HL,undefined_value="diagnose")
recall.HL = mldr::recall(label.yeast,pred.label.HL,undefined_value="diagnose")
fmeasure.HL = fmeasure(label.yeast,pred.label.HL,undefined_value="diagnose");fmeasure.HL
HL.HL = hamming_loss(label.yeast,pred.label.HL)

#CV by optimizing HL
accuracy.HL;precision.HL;recall.HL;fmeasure.HL;HL.HL

mac_accuracy.HL = mean ( apply(label.yeast == pred.label.HL, 2, sum)/dim(label.yeast)[1] )
mac_precision.HL = macro_precision(label.yeast,pred.label.HL)
mac_recall.HL = macro_recall(label.yeast,pred.label.HL)
mac_fmeasure.HL = macro_fmeasure(label.yeast,pred.label.HL)

mic_accuracy.HL = mean ( apply(label.yeast == pred.label.HL, 1, sum)/dim(label.yeast)[2] )
mic_precision.HL = micro_precision(label.yeast,pred.label.HL)
mic_recall.HL = micro_recall(label.yeast,pred.label.HL)
mic_fmeasure.HL = micro_fmeasure(label.yeast,pred.label.HL)

mac_accuracy.HL;mac_precision.HL;mac_recall.HL;mac_fmeasure.HL
mic_accuracy.HL;mic_precision.HL;mic_recall.HL;mic_fmeasure.HL

