# Taken from CRAN and modified as required for this Project

################################################################################### 
## singleIter function that is called by adaSample()
## modified to include functionality for RF and ntree, mtry parameter support
###################################################################################
singleIter <- function(Ps, Ns, dat, test=NULL, pos.probs=NULL, una.probs=NULL, classifier="svm", sampleFactor, seed, mtry=dim(dat)[2], ntree=100) {
  set.seed(seed);
  
  positive.train <- c()
  positive.cls <- c()
  
  # determine the proper sample size for creating a balanced dataset
  sampleN <- ifelse(length(Ps) < length(Ns), length(Ps), length(Ns))
  
  # bootstrap sampling to build the positive training set (labeled as 'P')
  idx.pl <- unique(sample(x=Ps, size=sampleFactor*sampleN, replace=TRUE, prob=pos.probs[Ps]))
  positive.train <- dat[idx.pl,]
  positive.cls <- rep("P", nrow(positive.train))
  
  # bootstrap sampling to build the "unannotate" or "negative" training set (labeled as 'N')
  idx.dl <- unique(sample(x=Ns, size=sampleFactor*sampleN, replace=TRUE, prob=una.probs[Ns]))
  unannotate.train <- dat[idx.dl,]
  unannotate.cls <- rep("N", nrow(unannotate.train))
  
  # combine data
  train.sample <- rbind(positive.train, unannotate.train)
  rownames(train.sample) <- NULL;
  cls <- as.factor(c(positive.cls, unannotate.cls))
  
  # training svm classifier
  if (classifier == "svm") {
    model.svm <- svm(cls~.,train.sample, probability=TRUE, scale=TRUE);
    svm.pred <- c();
    if (is.null(test)) {
      svm.pred <- predict(model.svm, dat, decision.values=TRUE, probability=TRUE);
    } else {
      svm.pred <- predict(model.svm, test, decision.values=TRUE, probability=TRUE);
    }
    return(attr(svm.pred,"probabilities"));
    
  } else if (classifier == "knn") {
    # training knn classifier
    if (is.null(test)) {
      knn.fit <- knn(train.sample, dat, cl=cls, k=5, prob=TRUE)
      
      p <- attr(knn.fit, "prob")
      idx <- which(knn.fit == "N")
      p[idx] <- 1- p[idx]
      knn.pred <- cbind(p, 1 - p)
      colnames(knn.pred) <- c("P", "N")
      rownames(knn.pred) <- rownames(dat)
      return(knn.pred)
    } else {
      test.mat <- test
      rownames(test.mat) <- NULL
      knn.fit <- knn(train.sample, test.mat, cl=cls, k=5, prob=TRUE)
      
      p <- attr(knn.fit, "prob")
      idx <- which(knn.fit == "N")
      p[idx] <- 1- p[idx]
      knn.pred <- cbind(p, 1 - p)
      colnames(knn.pred) <- c("P", "N")
      rownames(knn.pred) <- rownames(test)
      return(knn.pred)
    }
  } else if (classifier == "logit") {
    logit.model <- glm(cls~., family=binomial(link='logit'), data=data.frame(train.sample, cls))
    if (is.null(test)) {
      p <- predict(logit.model, newdata=data.frame(dat), type='response')
      logit.pred <- cbind(p, 1-p)
      colnames(logit.pred) <- c("P", "N")
      rownames(logit.pred) <- rownames(dat)
      return(logit.pred)
    } else {
      test.mat <- data.frame(test)
      rownames(test.mat) <- NULL
      colnames(test.mat) <- colnames(dat)
      p <- predict(logit.model, newdata=test.mat, type='response')
      logit.pred <- cbind(p, 1-p)
      colnames(logit.pred) <- c("P", "N")
      rownames(logit.pred) <- rownames(test)
      return(logit.pred)
    }
  } else if (classifier == "lda") {
    lda.model <- MASS::lda(cls~., data=data.frame(train.sample, cls))
    if (is.null(test)) {
      lda.pred <- predict(lda.model, data.frame(dat))$posterior
      colnames(lda.pred) <- c("N", "P")
      rownames(lda.pred) <- rownames(dat)
      return(lda.pred)
    } else {
      test.mat <- data.frame(test)
      rownames(test.mat) <- NULL
      colnames(test.mat) <- colnames(dat)
      lda.pred <- predict(lda.model, test.mat)$posterior
      colnames(lda.pred) <- c("N", "P")
      rownames(lda.pred) <- rownames(test)
      return(lda.pred)
    }
  } else if (classifier == "rf") { ## Add in random forest model
    rf.model <- randomForest::randomForest(cls~., data=data.frame(train.sample, cls), mtry=mtry, ntree=ntree)
    if (is.null(test)){
      rf.pred <- predict(rf.model, newdata=data.frame(dat),type="prob")
      rf.pred <- data.matrix(data.frame(rf.pred))
      colnames(rf.pred) <- c("N","P")
      rownames(rf.pred) <- rownames(dat)
      
      return(rf.pred)
    } else {
      test.mat <- data.frame(test)
      rownames(test.mat) <- NULL
      colnames(test.mat) <- colnames(dat)
      rf.pred <- predict(rf.model, newdata=test.mat, type='prob')
      rf.pred <- data.matrix(data.frame(rf.pred))
      colnames(rf.pred) <- c("N", "P")
      rownames(rf.pred) <- rownames(test)
      return(rf.pred)
    }
  }
  
}


########################################
# adaSample function (modified)
########################################
adaSample <- function(Ps, Ns, train.mat, test.mat, classifier="svm", s=1, C=1, sampleFactor=1, eps=0.01, mtry=dim(train.mat)[2], ntree=100) {
  
  # checking the input
  if(ncol(train.mat) != ncol(test.mat)) {stop("train.mat and test.mat do not have the same number of columns")}
  
  # initialize sampling probablity
  pos.probs <- rep(1, length(Ps))
  una.probs <- rep(1, length(Ns))
  names(pos.probs) <- Ps
  names(una.probs) <- Ns
  
  iterdiff <- eps+0.0001
  
  i = 0
  
  
  while (iterdiff >= eps & i < 5) {
    #while (i < 5) {
    # update count
    # i <- i + 1
    # training the predictive model
    
    model <- singleIter(Ps=Ps, Ns=Ns, dat=train.mat, pos.probs=pos.probs,
                        una.probs=una.probs, seed=i, classifier=classifier, sampleFactor=sampleFactor, mtry=mtry, ntree=ntree)
    
    
    
    # Calculate average difference in probabilities before and after loop before updating.
    iterdiff <- (sum(abs(pos.probs - model[Ps,"P"])) + sum(abs(una.probs - model[Ns,"N"]))) / (length(pos.probs) + length(una.probs))
    
    
    # update probability arrays
    pos.probs <- model[Ps, "P"]
    una.probs <- model[Ns, "N"]
    
    
    
    i <- i+1
  }
  
  pred <- singleIter(Ps=Ps, Ns=Ns, dat=train.mat, test=test.mat,
                     pos.probs=pos.probs, una.probs=una.probs, seed=s, classifier=classifier, sampleFactor=sampleFactor, mtry=mtry, ntree=ntree)
  
  # if C is greater than 1, create an ensemble
  if (C > 1){
    for (j in 2:C){
      pred <- pred + singleIter(Ps=Ps, Ns=Ns, dat=train.mat, test=test.mat,
                                pos.probs=pos.probs, una.probs=una.probs, seed=j, classifier=classifier, sampleFactor=sampleFactor, mtry=mtry, ntree=ntree)
    }
    pred <- pred/C
  }
  
  return(pred)
}


#################################### 
## Adasample Wrapper Function
####################################

adaSampleWrapperRF <- function(dat, cls, cls.known.neg, folds, mtry=dim(dat)[2], ntree=100){
  
  pred.known <- cls.known <- c()
  pred.all <- cls.all <-  c() # results for all labels, including unknown labels (assumed to be negative)
  pred.known.all <- cls.known.all <-  c() # results for known labels
  
  # For each fold:
  n=1
  for (fold in folds){
    
    if (n %% 5 == 0 | n == 1){
      print(paste("Running fold: ", n, sep=''))
    }
    
    #length(fold[1])
    #length(dat.extra.ada.akt_cls)
    #length(dat.extra.ada.akt_cls[c(-1,-4,-7)])
    
    trn.cls <- cls[-fold]
    tst.cls <- cls[fold] # need to use to estimate sensitivity
    
    tst.cls.known.neg <- cls.known.neg[fold] # Need to use to estimate lower bound for specificity
    
    train.mat <- dat[-fold,]
    test.mat <- dat[fold,]
    # cls <- dat.extra.ada.mtor_cls
    
    # Index positive and negative instances, this assumes all unlabeled samples are of the negative class
    Ps <- rownames(train.mat)[which(trn.cls == 1)]
    Ns <- rownames(train.mat)[which(trn.cls == 0)]
    
    ## Perform AdaSampling on the Training Set and provide prediction on the test set
    # require(AdaSampling)
    pred.prob <- adaSample(Ps, Ns, train.mat, test=test.mat, classifier="rf", C=50, mtry=mtry, ntree=ntree)
    
    ## Find which pred.probs were known positive class samples and which were known negative class samples
    # Known positive class examples are tst.cls
    pred.pos <- pred.prob[which(tst.cls == 1), "P"]
    # Known negative class examples are from tst.cls.mTOR
    pred.neg <- pred.prob[which(tst.cls.known.neg == 1), "P"]
    # Known actual classes stored in pred.prob.known
    pred.prob.known <- c(pred.pos, pred.neg)
    
    # Store the predictions made by AdaSampling and the Known actual classes into an array
    pred.known <- c(pred.known, ifelse(pred.prob.known > 0.5, 1, 0))
    
    
    cls.known <- c(cls.known,c(rep(1,length(which(tst.cls==1))), rep(0,length(which(tst.cls.known.neg==1)))   ))
    
    # Store all prediction probabilities of the test sets into an extending array
    pred.known.all <- c(pred.known.all, pred.known)
    cls.known.all <- c(cls.known.all, cls.known)
    
    ### Store everything including unknown labels
    pred.all <- c(pred.all, ifelse(pred.prob[, "P"] > 0.5, 1, 0))
    cls.all <- c(cls.all, as.numeric(tst.cls)-1)
    
    
    n <- n+1
  }
  
  
  # Sensitivity estimate using known samples
  known.scores <- scorer(pred.known.all,cls.known.all)
  
  # Specificity, using all samples
  allSample.scores <- scorer(pred.all,cls.all)
  
  print("Finished, returning results")
  
  return(c(known.scores,allSample.scores))
}


################################################################################################
## Scoring Function for ada sampling
## For sensitivity, specificity, precision, F1 score, Geometric Mean, accuracy
## Valid inputs for scoretypes is a vector containing
################################################################################################

scorer <- function(pred, actual, scoretypes = c('sen','spe','pre','F1','GM','acc')){
  # Note this function ASSUMES classification is binary and positive class is 1 while negative class is 0.
  # Therefore cannot use factors such as 'pos' and 'neg', they must be converted to 1 and 0.
  # pred and actual should be arrays of predictions and actual classes respectively
  TP <-  TN <-  FP <-  FN <-  0
  TP = sum(pred == actual & pred == 1)
  TN = sum(pred == actual & pred == 0)
  FP = sum(pred != actual & pred == 1)
  FN = sum(pred != actual & pred == 0)
  
  #resultDF <- data.frame(c("TP","TN","FP","FN"),c(TP,TN,FP,FN))
  #names(resultDF) <- c("Metric","Score")
  resultDF <- data.frame()
  
  
  
  for (type in scoretypes){
    if (type == 'sen'){
      res <- TP / (TP+FN)
      tempDF <- data.frame(type,res)
      names(tempDF) <- c("Metric","Score")
      resultDF <- rbind(resultDF,tempDF)
    }
    else if (type == 'spe')
    {
      res <- TN/(TN+FP)
      tempDF <- data.frame(type,res)
      names(tempDF) <- c("Metric","Score")
      resultDF <- rbind(resultDF,tempDF)
    }
    else if (type == 'pre')
    {
      res <- TP/(TP+FP)
      tempDF <- data.frame(type,res)
      names(tempDF) <- c("Metric","Score")
      resultDF <- rbind(resultDF,tempDF)
    }
    else if (type == 'F1')
    {
      res <- 2*TP/(2*TP+FP+FN)
      tempDF <- data.frame(type,res)
      names(tempDF) <- c("Metric","Score")
      resultDF <- rbind(resultDF,tempDF)
    }
    else if (type == 'GM')
    {
      res <- sqrt((TP/ (TP+FN))*(TP/ (TP+FP)))
      tempDF <- data.frame(type,res)
      names(tempDF) <- c("Metric","Score")
      resultDF <- rbind(resultDF,tempDF)
    }
    else if(type == 'acc')
    {
      res <- (TP+TN)/ (TP+TN+FP+FN)
      tempDF <- data.frame(type,res)
      names(tempDF) <- c("Metric","Score")
      resultDF <- rbind(resultDF,tempDF)
    }
    else{ # Capture if any invalid scoretypes were put in
      res <- NaN
      tempDF <- data.frame(type,res)
      names(tempDF) <- c("Metric","Score")
      resultDF <- rbind(resultDF,tempDF)
    }
    
    
    
  }
  return(resultDF)
  
}


