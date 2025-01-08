###################
## Yeast dataset ##
##  CC method    ##
###################
packages=c("mldr","mldr.datasets","stats","ROSE","gridExtra","Metrics","utiml")
lapply(packages,FUN=library,character.only=TRUE)

#####function for cross-validation#####
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





yeast=yeast()
n=dim(yeast$dataset)[1]
x=yeast$dataset[,1:103]
label.yeast = yeast$dataset[,104:117]

accuracy =  recall = fmeasure = precision = 0
K = 10 #use 10-fold CV
th.fix = 0.5 #thresholdv

pred.label <- matrix(0,nrow=n,ncol=dim(label.yeast)[2])
num_columns <- dim(yeast$dataset)[2]
num_labels <- dim(label.yeast)[2]

time0 <- Sys.time()

x_incre = yeast$dataset[,1:103]
for (j in 1:num_labels) {
  
  set.seed(j) #control randomness
  index.fold = fn.CV(label.yeast[,j]) #create the 10-fold index sets for label subset j
  
  current_label <- label.yeast[, j]
  mydata_j <- as.data.frame(cbind(label.yeast[,j],x))
  names(mydata_j)[1] <- "V1"
  
  mydata_j_test <- as.data.frame(cbind(label.yeast[,j],x_incre))
  names(mydata_j_test)[1] <- "V1"
  
  for (k in 1: K){
    # Fit logistic regression model
    fit <- glm(V1~.,data=mydata_j[-index.fold[[k]],],family=binomial)
    pred <- predict(fit,newdata = mydata_j_test[index.fold[[k]],],type="response")
    Y.hat <- ifelse(pred>th.fix,1,0)
    id1.set <- which(Y.hat==1)
    pred.label[index.fold[[k]][id1.set],j] <- 1 
  }
  
  #update x (feature)
  x[,paste0("pred_","label",j)]=label.yeast[,j]
  x_incre[,paste0("pred_","label",j)]=pred.label[,j]
  print(j)
  
}



time.total <- Sys.time() - time0

accuracy = mldr::accuracy(label.yeast,pred.label,undefined_value="diagnose")
precision = mldr::precision(label.yeast,pred.label,undefined_value="diagnose")
recall = mldr::recall(label.yeast,pred.label,undefined_value="diagnose")
fmeasure = fmeasure(label.yeast,pred.label,undefined_value="diagnose");fmeasure
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

