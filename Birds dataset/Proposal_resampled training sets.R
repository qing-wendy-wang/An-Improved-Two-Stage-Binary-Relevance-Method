#########################################
##            Birds dataset            ##
##Porposal with resampled training sets##
#########################################

packages=c("mldr","mldr.datasets","stats","ROSE","gridExtra","Metrics","gdata")
lapply(packages,FUN=library,character.only=TRUE)

birds=read.csv('birds_copy.csv',header=TRUE)

#####Obtain clusters#####
birds=read.csv('birds_copy.csv',header=TRUE)
label.birds <- birds[,261:279]  #645-by-19 matrix of observed labels
n=dim(label.birds)[1]

D <- 1-(cor(label.birds)) #correlation | cor() returns pairwise correlation between columns

hc.complete <-hclust(as.dist(D),method="complete")
plot(hc.complete) #dendrogram
rect.hclust(hc.complete,k=11)
plot(hc.complete$height,type='b',pch=19,lwd=1.5, xlab="Number of clusters", ylab="The optimal between-cluster dissimilarity at merging",xaxt='n')
axis(1, at=1:18, labels=18:1)

result.clusters <- cutree(hc.complete,k=11)
result.clusters



####create label subsets####
y= birds[,261:279]
x.matrix <- birds[,1:260]
x.matrix=apply(x.matrix,2,as.numeric)

label.subset <- matrix(0,nrow=n,ncol=11) 
for (i in 1:n)
{
  if (y[i,1]==1 && y[i,8]==1)
  {
    label.subset[i,1] = 1
  } 
  if (y[i,2]==1 && y[i,12]==1) 
  {
    label.subset[i,2] = 1
  } 
  if (y[i,3]==1 && y[i,4]==1)
  {
    label.subset[i,3]=1
  } 
  if (y[i,5]==1 && y[i,7]==1)
  {
    label.subset[i,4]=1
  }
  if ( y[i,6]==1 && y[i,13]==1)
  {
    label.subset[i,5]=1
  }
  if (y[i,9]==1 && y[i,11]==1 && y[i,15]==1)
  {
    label.subset[i,6]=1
  }
  if (y[i,10]==1)
  {
    label.subset[i,7]=1
  }
  if (y[i,14]==1 && y[i,17]==1)
  {
    label.subset[i,8]=1
  }
  if (y[i,16]==1)
  {
    label.subset[i,9]=1
  }
  if (y[i,18]==1)
  {
    label.subset[i,10]=1
  }
  if (y[i,19]==1 )
  {
    label.subset[i,11]=1
  }
}

freq=apply(label.subset,2,sum)
freq/n;freq
label.subset=as.data.frame(label.subset)
names(label.subset) <- c("labelset1","labelset2","labelset3","labelset4",
                         "labelset5","labelset6","labelset7","labelset8", 
                         "labelset9","labelset10","labelset11")
label.index <- list( c(1,8), c(2,12),c(3,4),c(5,7),c(6,13),
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
      
      if (sum(label.subset==1) <10)
      {
        
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

    

####resample####
pred.label.under <- matrix(0,nrow=n,ncol=dim(label.birds)[2])
pred.label.over <- matrix(0,nrow=n,ncol=dim(label.birds)[2])
pred.label.smote <- matrix(0,nrow=n,ncol=dim(label.birds)[2])


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
    train.id=c()
    for (t in setdiff(1:10,k)){
      train.id=c(train.id,index.fold[[t]])
    }
    
    if(length(label.index[[j]])>1 ) #if the labelset contains more than one label
    {
      mydata_j <- as.data.frame(cbind(label.subset[,j],x.matrix))
      names(mydata_j)[1] <- "V1"
      
      #undersample: 
      mydata_j.under <- ovun.sample(as.factor(V1)~., data=mydata_j[train.id,],method="under")$data
      #oversample: 
      mydata_j.over <- ovun.sample(as.factor(V1)~., data=mydata_j[train.id,],method="over")$data
      #smote:
      mydata_j.smote <- ROSE(as.factor(V1)~.,data=mydata_j[train.id,],seed=123)$data
      
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
        temp <- glm(V1~.,data=mydata.temp[train.id,],family=binomial)
        
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
      fit <- glm(V1~.,data=mydata_j[train.id,],family=binomial)
      pred <- predict(fit,newdata = mydata_j[index.fold[[k]],],type="response")
      Y.hat <- ifelse(pred>th.fix,1,0)
      id1.set <- which(Y.hat==1)
      pred.label.under[index.fold[[k]][id1.set],label.index[[j]]] =
      pred.label.over[index.fold[[k]][id1.set],label.index[[j]]]=
        pred.label.smote[index.fold[[k]][id1.set],label.index[[j]]] <- 1 
    }
    
  }
  
  
  print(k) #print the curent fold to track progress 
}

time.total <- Sys.time() - time0

#undersample metrics
accuracy.under = mldr::accuracy(label.birds,pred.label.under,undefined_value="diagnose")
precision.under = mldr::precision(label.birds,pred.label.under,undefined_value="diagnose")
recall.under = mldr::recall(label.birds,pred.label.under,undefined_value="diagnose")
fmeasure.under = fmeasure(label.birds,pred.label.under,undefined_value="diagnose")
HL.under = hamming_loss(label.birds,pred.label.under)
accuracy.under;precision.under;recall.under;fmeasure.under;HL.under

mac_accuracy.under = mean ( apply(label.birds == pred.label.under, 2, sum)/dim(label.birds)[1] )
mac_precision.under = macro_precision(label.birds,pred.label.under)
mac_recall.under = macro_recall(label.birds,pred.label.under)
mac_fmeasure.under = macro_fmeasure(label.birds,pred.label.under)

mic_accuracy.under = mean ( apply(label.birds == pred.label.under, 1, sum)/dim(label.birds)[2] )
mic_precision.under = micro_precision(label.birds,pred.label.under)
mic_recall.under = micro_recall(label.birds,pred.label.under)
mic_fmeasure.under = micro_fmeasure(label.birds,pred.label.under)

mac_accuracy.under;mac_precision.under;mac_recall.under;mac_fmeasure.under
mic_accuracy.under;mic_precision.under;mic_recall.under;mic_fmeasure.under


#oversample metrics
accuracy.over = mldr::accuracy(label.birds,pred.label.over,undefined_value="diagnose")
precision.over = mldr::precision(label.birds,pred.label.over,undefined_value="diagnose")
recall.over = mldr::recall(label.birds,pred.label.over,undefined_value="diagnose")
fmeasure.over = fmeasure(label.birds,pred.label.over,undefined_value="diagnose")
HL.over = hamming_loss(label.birds,pred.label.over)
accuracy.over;precision.over;recall.over;fmeasure.over;HL.over

mac_accuracy.over = mean ( apply(label.birds == pred.label.over, 2, sum)/dim(label.birds)[1] )
mac_precision.over = macro_precision(label.birds,pred.label.over)
mac_recall.over = macro_recall(label.birds,pred.label.over)
mac_fmeasure.over = macro_fmeasure(label.birds,pred.label.over)

mic_accuracy.over = mean ( apply(label.birds == pred.label.over, 1, sum)/dim(label.birds)[2] )
mic_precision.over = micro_precision(label.birds,pred.label.over)
mic_recall.over = micro_recall(label.birds,pred.label.over)
mic_fmeasure.over = micro_fmeasure(label.birds,pred.label.over)

mac_accuracy.over;mac_precision.over;mac_recall.over;mac_fmeasure.over
mic_accuracy.over;mic_precision.over;mic_recall.over;mic_fmeasure.over

#smote-sample metrics
accuracy.smote = mldr::accuracy(label.birds,pred.label.smote,undefined_value="diagnose")
precision.smote = mldr::precision(label.birds,pred.label.smote,undefined_value="diagnose")
recall.smote = mldr::recall(label.birds,pred.label.smote,undefined_value="diagnose")
fmeasure.smote = fmeasure(label.birds,pred.label.smote,undefined_value="diagnose")
HL.smote = hamming_loss(label.birds,pred.label.smote)
accuracy.smote;precision.smote;recall.smote;fmeasure.smote;HL.smote

mac_accuracy.smote = mean ( apply(label.birds == pred.label.smote, 2, sum)/dim(label.birds)[1] )
mac_precision.smote = macro_precision(label.birds,pred.label.smote)
mac_recall.smote = macro_recall(label.birds,pred.label.smote)
mac_fmeasure.smote = macro_fmeasure(label.birds,pred.label.smote)

mic_accuracy.smote = mean ( apply(label.birds == pred.label.smote, 1, sum)/dim(label.birds)[2] )
mic_precision.smote = micro_precision(label.birds,pred.label.smote)
mic_recall.smote = micro_recall(label.birds,pred.label.smote)
mic_fmeasure.smote = micro_fmeasure(label.birds,pred.label.smote)

mac_accuracy.smote;mac_precision.smote;mac_recall.smote;mac_fmeasure.smote
mic_accuracy.smote;mic_precision.smote;mic_recall.smote;mic_fmeasure.smote

