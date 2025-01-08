##################
#Birds dataset#
##################
#default proposal

packages=c("mldr","mldr.datasets","stats","ROSE","gridExtra","Metrics","gdata")
lapply(packages,FUN=library,character.only=TRUE)
birds=read.csv('birds_copy.csv',header=TRUE)

#####Obtain clusters#####
label.birds <- birds[,261:279]  #645-by-19 matrix of observed labels
n=dim(label.birds)[1]


D <- 1-(cor(label.birds)) #correlation | cor() returns pairwise correlation between columns

hc.complete <-hclust(as.dist(D),method="complete")
plot(hc.complete)
result.clusters <- cutree(hc.complete,k=12);result.clusters



####create label subsets####
y= birds[,261:279]
x.matrix <- birds[,1:260]
x.matrix=apply(x.matrix,2,as.numeric)

label.subset <- matrix(0,nrow=n,ncol=12) 
for (i in 1:n){
  if (y[i,1]==1 && y[i,8]==1)
  {
    label.subset[i,1] = 1
  } 
  if (y[i,2]==1 && y[i,12]==1) 
  {
    label.subset[i,2] = 1
  } 
  if (y[i,3]==1)
  {
    label.subset[i,3]=1
  } 
  if (y[i,4]==1)
  {
    label.subset[i,4]=1
  } 
  if (y[i,5]==1 && y[i,7]==1)
  {
    label.subset[i,5]=1
  }
  if ( y[i,6]==1 && y[i,13]==1)
  {
    label.subset[i,6]=1
  }
  if (y[i,9]==1 && y[i,11]==1 && y[i,15]==1)
  {
    label.subset[i,7]=1
  }
  if (y[i,10]==1)
  {
    label.subset[i,8]=1
  }
  if (y[i,14]==1 && y[i,17]==1)
  {
    label.subset[i,9]=1
  }
  if (y[i,16]==1)
  {
    label.subset[i,10]=1
  }
  if (y[i,18]==1)
  {
    label.subset[i,11]=1
  }
  if (y[i,19]==1 )
  {
    label.subset[i,12]=1
  }
}

freq=apply(label.subset,2,sum)
freq/n;freq
label.subset=as.data.frame(label.subset)
names(label.subset) <- c("labelset1","labelset2","labelset3","labelset4",
                         "labelset5","labelset6","labelset7","labelset8", 
                         "labelset9","labelset10","labelset11","labelset12")
label.index <- list( c(1,8), c(2,12),c(3),c(4),c(5,7),c(6,13),
                     c(9,11,15),c(10),c(14,17),c(16),c(18),c(19))


####cross-validation####


#######################################
##create index subsets for 10-fold CV##
#######################################

##Note: this functoin, fn.CV, is created s.t. each fold contains 1/K obs from the true 1's from the given label subset
##use this function to realize 10-fold CV for predicting each label subset
fn.CV <- function(label.subset)
{
  K <- 10 #10-fold CV as an example
  index.fold <- list()
  
  if (sum(label.subset==1) <10){
    #case 1: 1's less than the number of folds. We make sure each fold contains AT LEAST n/K*0.02 positive cases
    
    id.subset = which(label.subset==1)  #the positive cases in the dataset
    
    n.1.fold = ceiling(n/K*0.02)  #the number of positive cases needed in each fold
    n.1.total = ceiling(n/K*0.02)*K   #the number of positive cases in all folds
    id.1 = gdata::resample(id.subset,n.1.total,replace=TRUE) #I would recommend that we sample with replacement for all the needed positive cases
    
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

####fixed threshold####
pred.label <- matrix(0,nrow=n,ncol=dim(label.birds)[2])
accuracy =  recall = fmeasure = precision = 0
K = 10 #use 10-fold CV
th.fix = 0.5 #threshold

time0 <- Sys.time()
for(k in 1:K)
{
  for(j in 1:dim(label.subset)[2])
  {
    set.seed(j) #control randomness
    index.fold = fn.CV(label.subset[,j]) #create the 10-fold index sets for label subset j
    train.id=c()
      for (t in setdiff(1:10,k)){
        train.id=c(train.id,index.fold[[t]])
      }
    
    if(length(label.index[[j]])>1) #the labelset contains more than one label
    {
      mydata_j <- as.data.frame(cbind(label.subset[,j],x.matrix))
      names(mydata_j)[1] <- "V1"
      
      fit <- glm(V1~.,data=mydata_j[train.id,],family=binomial)
      
      pred <- predict(fit,newdata = mydata_j[index.fold[[k]],],type="response") #type=response returns the predicted probabilities of having this label subset
      
      Y.hat <- ifelse(pred>th.fix,1,0)
      
      id1.set <- which(Y.hat==1)
      pred.label[index.fold[[k]][id1.set],label.index[[j]]] <- 1 
     
      
      id0.set <- which(Y.hat!=1)
      
      for(i in 1:length(label.index[[j]]))
      {
        mydata.temp <- as.data.frame(cbind(y[,label.index[[j]][i]],x.matrix)) 
        names(mydata.temp)[1] <- "V1"
        temp <- glm(V1~.,data=mydata.temp[train.id,],family=binomial)
        pred.temp <- predict(temp,newdata=mydata.temp[index.fold[[k]][id0.set],],type="response")
        Y.hat.temp <- ifelse(pred.temp>th.fix,1,0)    
        pred.label[index.fold[[k]][id0.set[which(Y.hat.temp==1)]],label.index[[j]][i]] <- 1 
      }
      #}
    }else #if the labelset is a singleton
    {
      mydata_j <- as.data.frame(cbind(label.subset[,j],x.matrix))
      names(mydata_j)[1] <- "V1"
      fit <- glm(V1~.,data=mydata_j[train.id,],family=binomial)
      pred <- predict(fit,newdata = mydata_j[index.fold[[k]],],type="response")
      Y.hat <- ifelse(pred>th.fix,1,0)
      id1.set <- which(Y.hat==1)
      pred.label[index.fold[[k]][id1.set],label.index[[j]]] <- 1 
    }
    
  }
  
  print(k) #print the current fold to track progress
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









