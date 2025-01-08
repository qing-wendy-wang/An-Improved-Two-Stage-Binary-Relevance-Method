##Birds Dataset##
##  CC method  ##

packages=c("mldr","mldr.datasets","stats","ROSE","gridExtra","Metrics","utiml","gdata")
lapply(packages,FUN=library,character.only=TRUE)

#####function for cross-validation#####
fn.CV <- function(label.subset)
{
  K <- 10 #10-fold CV as an example
  index.fold <- list()
  
  if (sum(label.subset==1) <10){
    
    id.subset = which(label.subset==1)  #the positive cases in the dataset
    
    n.1.fold = ceiling(n/K*0.02)  #the number of positive cases needed in each fold
    n.1.total = ceiling(n/K*0.02)*K   #the number of positive cases in all folds
    id.1 = gdata::resample(id.subset,n.1.total,replace=TRUE) #I would recommend that we sample with replacement for all the needed positive cases
    # id.1 = sample(id.1,length(id.1),replace=FALSE) #this line is not necessary, as we already did random sampling in the previous line; the order is already arbitrary
    
    n.remain.total = n-length(id.subset)
    id.remain <- sample(setdiff(1:n,id.subset),n.remain.total,replace=FALSE)
    n.remain.fold = floor(n.remain.total/K)
    
    for(i in 1:K)
    {
      if(i<K)
      {index.fold[[i]] <- c(id.1[ ((i-1)*n.1.fold+1) :(i*n.1.fold)],
                            id.remain[ ((i-1)*n.remain.fold+1) :(i*n.remain.fold)])    }
      else
      {index.fold[[i]] <- c(id.1[ ((i-1)*n.1.fold+1) :length(id.1)]
                            ,id.remain[((K-1)*n.remain.fold+1):length(id.remain)])    }
    }
  }
  else
  { 
    n.fold = floor(n/K)
    id.1 <- sample(which(label.subset==1),length(which(label.subset==1)),replace=FALSE)
    
    n.labelsubset.fold <- floor(length(id.1)/K)
    #number of indexes in each fold that need to be taken from the true 1's of the given label subset
    id.remain <- sample(setdiff(1:n, id.1),n-length(id.1),replace=FALSE)
    #shuffle the remaining id's that correspond to 0's for the label subset
    n.remain <- n.fold-n.labelsubset.fold
    
    for(i in 1:K)
    {
      if(i<K)
      {
        index.fold[[i]] <- c(id.1[((i-1)* n.labelsubset.fold+1):(i*n.labelsubset.fold)],id.remain[((i-1)*n.remain+1):(i*n.remain)])
      }else
      {
        index.fold[[i]] <- c(id.1[((K-1)* n.labelsubset.fold+1):length(id.1)],id.remain[((K-1)*n.remain+1):length(id.remain)])
      }
    }
    
  }
  
  return(index.fold)
} 





n=dim(birds$dataset)[1]
x=birds$dataset[,1:260]
label.birds = birds$dataset[,261:279]

accuracy =  recall = fmeasure = precision = 0
K = 10 #use 10-fold CV
th.fix = 0.5 #thresholdv

pred.label <- matrix(0,nrow=n,ncol=dim(label.birds)[2])
num_columns <- dim(birds$dataset)[2]
num_labels <- dim(label.birds)[2]

time0 <- Sys.time()


x_incre = birds$dataset[,1:260]
for (j in 1:num_labels) {
  
  set.seed(j) #control randomness
  index.fold = fn.CV(label.birds[,j]) #create the 10-fold index sets for label subset j
  
  current_label <- label.birds[, j]
  mydata_j <- as.data.frame(cbind(label.birds[,j],x))
  names(mydata_j)[1] <- "V1"
  
  if (j==1){
    mydata_j_test <- as.data.frame(cbind(label.birds[,j],x_incre))
    names(mydata_j_test)[1] <- "V1"
  }else{
    mydata_j_test <- as.data.frame(cbind(label.birds[,j],x_incre))
    names(mydata_j_test)[1] <- "V1"
  }
  
  for (k in 1: K){
    train.id=c()
    for (t in setdiff(1:10,k)){
      train.id=c(train.id,index.fold[[t]])
    }
    
    # Fit logistic regression model
    fit <- glm(V1~.,data=mydata_j[train.id,],family=binomial)
    pred <- predict(fit,newdata = mydata_j_test[index.fold[[k]],],type="response")
    
    Y.hat <- ifelse(pred>th.fix,1,0)
    id1.set <- which(Y.hat==1)
    pred.label[index.fold[[k]][id1.set],j] <- 1 
  }
  

  
  x[,paste0("pred_","label",j)]=label.birds[,j]
  x_incre[,paste0("pred_","label",j)]=pred.label[,j]
  print(j)
}



time.total <- Sys.time() - time0

accuracy = mldr::accuracy(label.birds,pred.label,undefined_value="diagnose")
precision = mldr::precision(label.birds,pred.label,undefined_value="diagnose")
recall = mldr::recall(label.birds,pred.label,undefined_value="diagnose")
fmeasure = fmeasure(label.birds,pred.label,undefined_value="diagnose")
HL = hamming_loss(label.birds,pred.label)

accuracy;precision;recall;fmeasure;HL

mac_accuracy = mean ( apply(label.birds == pred.label, 2, sum)/dim(label.birds)[1] )
mac_precision = macro_precision(label.birds,pred.label)
mac_recall = macro_recall(label.birds,pred.label)
mac_fmeasure = macro_fmeasure(label.birds,pred.label)

mic_accuracy = mean ( apply(label.birds == pred.label, 1, sum)/dim(label.birds)[2] )
mic_precision = micro_precision(label.birds,pred.label)
mic_recall = micro_recall(label.birds,pred.label)
mic_fmeasure = micro_fmeasure(label.birds,pred.label)

mac_accuracy;mac_precision;mac_recall;mac_fmeasure
mic_accuracy;mic_precision;mic_recall;mic_fmeasure


