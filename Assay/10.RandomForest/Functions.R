# Functions might be used for random forest analyses

# Load library
library(plyr)
library(ggplot2)
library(reshape2)
library(sampling)
library(randomForest)
library(pROC)
#library(party)
library(parallel)
#library(cowplot)
library(ggpubr)
#library(DMwR)

# Self-defined functions

# Global setting
chaoTheme <- theme(axis.title=element_text(size=16,face="bold"),
                   text=element_text(size=14,face="bold"))

## Randomly select samples for training set
aveSam <- function(prf,phe,ID,grp,block=NULL,candy=NULL,size,pick=NULL,dw=NULL,ex.PID=NULL,sam=T,smote=F){
  ind <- which(colnames(phe)==ID)
  igrp <- which(colnames(phe)==grp)
  sphe <- phe[,c(ind,igrp)]
  colnames(sphe) <- c("ID","grp")
  id <- phe$ID
  tprf <- as.data.frame(t(prf))
  if(length(pick)>0){
    tprf <- tprf[,colnames(tprf)%in%pick]
  }
  tprf$ID <- rownames(tprf)
  if(length(candy)>0){
    cphe <- phe[,c(ind,which(colnames(phe)%in%candy))]
    colnames(cphe)[1] <- "ID"
    tprf <- merge(cphe,tprf,by="ID")
  }
  dat <- merge(sphe,tprf,by="ID")
  dat <- dat[which(!is.na(dat$grp)),]
  if(!is.null(block)){
    dats <- dat[which(dat$grp!=block),]
    dats$grp <- droplevels(dats$grp)
  }else{
    dats <- dat
  }
  ex.DNAID <- as.character(phe$DNAID[which(phe$PID%in%ex.PID)])
  if(!is.null(ex.PID)&&length(ex.DNAID)>0){
    dats <- dats[-which(dats$ID%in%ex.DNAID),]
  }
  #
  min.size <- min(table(dats$grp))
  if(min.size < size){
    size <- min.size -1
    warning(paste0("size large than possible. Choose the maximun available size: ",size,". "))
  }

  ###
  if(sam){
    avesampling <- strata(dats,stratanames=c("grp"),size=rep(size,length(table(dats$grp))), method="srswor")
    train<- getdata(dats,avesampling)[,1:ncol(dats)]
    train.id <- train$ID
    train <- train[,-1]
    rownames(train) <- train.id
    train.x <- train[,-ncol(train)]
    train.y <- train[,ncol(train)]
    test.id <- setdiff(dat$ID,train.id) #dats
    test <- dat[which(dat$ID%in%test.id),]
    
    rownames(test) <- test$ID
    test.y <- test$grp
    test.x <- test[,-c(1:2)]
  }else{
    train <- dats
    train.id <- train$ID
    rownames(train) <- train.id
    train.y <- train$grp
    train.x <- train[,-which(colnames(train)%in%c("grp","ID"))]
    test.x <-NULL
    test.y <- NULL
  }
  
  ###
  if(smote){
    sdat <- dats[which(dats$ID%in%train.id),]
    rownames(sdat) <- sdat$ID
    sdat <- sdat[,-which(colnames(sdat)=="ID")]
    grps <- sort(table(sdat$grp))
    # pOver <- floor(max(grps)/min(grps)) * 100
    # pUnder<- round((pOver + 100)/pOver,2) * 100
    gNames <- names(grps)
    dnew <- sdat
    for(i in 2:length(grps)-1){
      di <- sdat[which(sdat$grp%in%gNames[c(i,length(grps))]),]
      di$grp <- droplevels(di$grp)
      grpi <- sort(table(di$grp))
      pOver <- (round(max(grpi)/min(grpi),2) -1 ) * 100
      pEst <- floor(max(grpi)/min(grpi)) * 100
      refundNum <- min(grpi) * (pEst/100+1) - max(grpi) # because perc.over/100 cannot handle its fractional part
      keepIndex <- 1:(min(grpi) * pEst/100 - refundNum)
      newi <- SMOTE(grp ~., di, perc.over = pEst, k = 5, perc.under = 100 )
      dnew <- rbind(dnew,newi[which(rownames(newi)%in%keepIndex),])
    }
    train <- dnew
    train.id <- rownames(train)
    train.y <- train$grp
    train.x <- train[,-which(colnames(train)%in%c("grp","ID"))]
  }
  #exclude <- dat[which(dat$ID%in%ATB.0_2.DNAID),]
  #exclude.y <- exclude$grp
  #exclude.x <- exclude[,-c(1:2)]

  if(length(dw)>0){
    if(!dir.exists(dw)){dir.create(dw)}
    write.table(t(train.x),paste0(dw,grp,"_train_x.xls"),quote=F,col.names=NA,sep="\t")
    write.table(train.y,paste0(dw,grp,"_train_y.xls"),quote=F,col.names=NA,sep="\t")
    write.table(t(test.x),paste0(dw,grp,"_test_x.xls"),quote=F,col.names=NA,sep="\t")
    write.table(test.y,paste0(dw,grp,"_test_y.xls"),quote=F,col.names=NA,sep="\t")
  }
  return(list(train.x=train.x,train.y=train.y,
              test.x=test.x,test.y=test.y))
  #exclude.x=exclude.x,exclude.y=exclude.y))
}

# stratified sampling
stratiSam <- function(prf,phe,ID,grp,block=NULL,candy=NULL,
                      size=NULL,pct=NULL,strati=NULL,
                      pick=NULL,dw=NULL,ex.PID=NULL,sam=T,smote=F){
  ind <- which(colnames(phe)==ID)
  igrp <- which(colnames(phe)==grp)
  istra<- which(colnames(phe)==strati)
  sphe <- phe[,c(ind,igrp,istra)]
  colnames(sphe) <- c("ID","grp","strata")
  id <- phe$ID
  tprf <- as.data.frame(t(prf))
  if(length(pick)>0){
    tprf <- tprf[,colnames(tprf)%in%pick]
  }
  tprf$ID <- rownames(tprf)
  if(length(candy)>0){
    cphe <- phe[,c(ind,which(colnames(phe)%in%candy))]
    colnames(cphe)[1] <- "ID"
    tprf <- merge(cphe,tprf,by="ID")
  }
  dat <- merge(sphe,tprf,by="ID")
  dat <- dat[which(!is.na(dat$grp)),]
  if(!is.null(block)){
    dats <- dat[which(dat$grp!=block),]
    dats$grp <- droplevels(dats$grp)
  }else{
    dats <- dat
  }
  ex.DNAID <- as.character(phe$DNAID[which(phe$PID%in%ex.PID)])
  if(!is.null(ex.PID)&&length(ex.DNAID)>0){
    dats <- dats[-which(dats$ID%in%ex.DNAID),]
  }
  #
  
  if(sam){
    grpd <- dats[,c("ID","strata")]
    grpd <- grpd[order(grpd$strata),]
    ###
    if(!is.null(pct)){
      size=round(table(grpd$strata)*pct,0)
    }else{
      size=rep(size,length(table(grpd$strata)))
    }
    ###
    avesampling <- strata(grpd,stratanames=c("strata"),size=size, method="srswor")
    getID <- getdata(grpd,avesampling)[,1:2]
    train<- dats[which(dats$ID%in%getID$ID),]
    train.id <- train$ID
    train.y <- train$grp
    train.x <- train[,-c(1:3)]
    rownames(train.x) <- train.id
    
    test.id <- setdiff(dat$ID,train.id) #dats
    test <- dat[which(dat$ID%in%test.id),]
    
    rownames(test) <- test$ID
    test.y <- test$grp
    test.x <- test[,-c(1:3)]
  }else{
    train <- dats
    train.id <- train$ID
    rownames(train) <- train.id
    train.y <- train$grp
    train.x <- train[,-which(colnames(train)%in%c("grp","ID"))]
    test.x <-NULL
    test.y <- NULL
  }
  
  ###
  if(smote){
    sdat <- dats[which(dats$ID%in%train.id),]
    rownames(sdat) <- sdat$ID
    sdat <- sdat[,-which(colnames(sdat)=="ID")]
    grps <- sort(table(sdat$grp))
    # pOver <- floor(max(grps)/min(grps)) * 100
    # pUnder<- round((pOver + 100)/pOver,2) * 100
    gNames <- names(grps)
    dnew <- sdat
    for(i in 2:length(grps)-1){
      di <- sdat[which(sdat$grp%in%gNames[c(i,length(grps))]),]
      di$grp <- droplevels(di$grp)
      grpi <- sort(table(di$grp))
      pOver <- (round(max(grpi)/min(grpi),2) -1 ) * 100
      pEst <- floor(max(grpi)/min(grpi)) * 100
      refundNum <- min(grpi) * (pEst/100+1) - max(grpi) # because perc.over/100 cannot handle its fractional part
      keepIndex <- 1:(min(grpi) * pEst/100 - refundNum)
      newi <- SMOTE(grp ~., di, perc.over = pEst, k = 5, perc.under = 100 )
      dnew <- rbind(dnew,newi[which(rownames(newi)%in%keepIndex),])
    }
    train <- dnew
    train.id <- rownames(train)
    train.y <- train$grp
    train.x <- train[,-which(colnames(train)%in%c("grp","ID"))]
  }
  #exclude <- dat[which(dat$ID%in%ATB.0_2.DNAID),]
  #exclude.y <- exclude$grp
  #exclude.x <- exclude[,-c(1:2)]
  
  if(length(dw)>0){
    if(!dir.exists(dw)){dir.create(dw)}
    write.table(t(train.x),paste0(dw,grp,"_train_x.xls"),quote=F,col.names=NA,sep="\t")
    write.table(train.y,paste0(dw,grp,"_train_y.xls"),quote=F,col.names=NA,sep="\t")
    write.table(t(test.x),paste0(dw,grp,"_test_x.xls"),quote=F,col.names=NA,sep="\t")
    write.table(test.y,paste0(dw,grp,"_test_y.xls"),quote=F,col.names=NA,sep="\t")
  }
  return(list(train.x=train.x,train.y=train.y,
              test.x=test.x,test.y=test.y))
  #exclude.x=exclude.x,exclude.y=exclude.y))
}

# Cross-validation method to select candidate markers
rfcv1 <- function (trainx, trainy, cv.fold = 5, scale = "log", step = 0.5,ntree=500,
                   mtry = function(p) max(1, floor(sqrt(p))), recursive = FALSE,impGrp=NULL,
                   ...) {
  classRF <- is.factor(trainy)
  n <- nrow(trainx)
  p <- ncol(trainx)
  if (scale == "log") {
    k <- floor(log(p, base = 1/step))
    n.var <- round(p * step^(0:(k - 1)))
    n.var <- rev(sort(unique(c(n.var[1]-1,n.var[which(n.var>1)],1))))
  }else {
    n.var <- seq(from = p, to = 1, by = step)
  }
  k <- length(n.var)
  cv.pred <- vector(k, mode = "list")
  cv.vote <- vector(k, mode = "list")
  for (i in 1:k) {
    cv.pred[[i]] <- trainy
    cv.vote[[i]] <- matrix(NA,nrow=nrow(trainx),ncol=length(levels(trainy)),
                           dimnames=list(rownames(trainx),levels(trainy)))
  }
  if (classRF) {
    f <- trainy
  } else {
    f <- factor(rep(1:5, length = length(trainy))[order(order(trainy))])
  }
  nlvl <- table(f)
  idx <- numeric(n)
  for (i in 1:length(nlvl)) {
    idx[which(f == levels(f)[i])] <- sample(rep(1:cv.fold,
                                                length = nlvl[i]))
  }

  res = list()
  for (i in 1:cv.fold) {
    all.rf <- randomForest(trainx[idx != i, , drop = FALSE],
                           trainy[idx != i], 
                           trainx[idx == i, , drop = FALSE],
                           trainy[idx == i], mtry = mtry(p), importance = TRUE, norm.votes=T,ntree=ntree, ...)
    cv.pred[[1]][idx == i] <- all.rf$test$predicted
    cv.vote[[1]][idx == i,] <- all.rf$test$votes
    if(is.null(impGrp)){
      impvar <- (1:p)[order(importance(all.rf, type = 1),
                            decreasing = TRUE)]
    }else{
      impvar <- (1:p)[order(all.rf$importance[,impGrp]/all.rf$importanceSD[,impGrp],
                            decreasing = TRUE)]
    }
    
    res[[i]] <- impvar
    for (j in 2:k) {
      imp.idx <- impvar[1:n.var[j]]
      sub.rf <- randomForest(trainx[idx != i, imp.idx, drop = FALSE], 
                             trainy[idx != i], 
                             trainx[idx ==i, imp.idx, drop = FALSE], 
                             trainy[idx == i],
                             mtry = mtry(n.var[j]), importance = recursive, ...)
      cv.pred[[j]][idx == i] <- sub.rf$test$predicted
      cv.vote[[j]][idx == i,] <- sub.rf$test$votes
      if (recursive) {
        impvar <- (1:length(imp.idx))[order(importance(sub.rf,type = 1), decreasing = TRUE)]
      }
      NULL
    }
    NULL
  }
  if (classRF) {
    error.cv <- sapply(cv.pred, function(x) mean(trainy != x))
    #error.cv <- sapply(cv.pred, function(x) mean((as.numeric(trainy)-1 - x)^2))
  }
  else {
    error.cv <- sapply(cv.pred, function(x) mean((trainy - x)^2))
  }
  names(error.cv) <- names(cv.pred) <- n.var
  list(n.var = n.var, error.cv = error.cv, predicted = cv.pred, votes=cv.vote,res = res)
}

tt <- NULL
train.module <- function(prf,phe,ID,grp,size=NULL,max.marker=30,max.cv=0.4,
                         thread=5,rept=5,steps=2,cv.fold=5,rrf=10,block=NULL,
                         candy=NULL,pick=NULL,dw=NULL,ex.PID=NULL,mtry=3,
                         strati=NULL,pct=NULL){
  #step1
  min.cv <- 1
  marker.num=0
  # parallel method to speed up the process
  cl <- makeCluster(thread)
  max.try <- mtry
  tried = 0
  tmp.train <- NULL
  tmp.cv <- 1
  while((marker.num>max.marker|marker.num==0|min.cv > max.cv) & (tried <= max.try)){
    if(tried == max.try){
      tried <- "done"
      train.cv <- tmp.train
      
      error.cv <- sapply(train.cv, "[[", "error.cv")
      error.cv.rm <- rowMeans(error.cv)
      error.cv.sd <- apply(error.cv,1,sd)
      id <- max(which(error.cv.rm < min(error.cv.sd + error.cv.rm)))
      min.cv <- error.cv.rm[id]
      
      marker.num <- min(as.numeric(names(error.cv.rm)[id]))
      
      mode <- data.frame(num=train.cv[[1]]$n.var,error.cv)
      mmode <- melt(mode,id.vars = "num",variable.name = "times",value.name = "Err.cv")
      mmode$num  <- as.factor(mmode$num)
      
      smode <- data.frame(num=train.cv[[1]]$n.var,times="rm",Err.cv=error.cv.rm)
      smode$num <- as.factor(smode$num)
    }else{
      tried <- tried + 1
      
      if(is.null(strati)){
        tt <- aveSam(prf,phe,ID,grp,NULL,candy,
                     size,pick,ex.PID=ex.PID)
      }else{
        tt <- stratiSam(prf,phe,ID,grp,NULL,candy,
                        pct=pct,strati=strati,pick,ex.PID=ex.PID)
      }
      
      if(class(tt$train.y)=="factor"){ tt$train.y <- droplevels(tt$train.y) }
      
      #clusterEvalQ(cl,source("rfcv1.R"))
      clusterEvalQ(cl,library(randomForest))
      clusterExport(cl,c("tt","rfcv1"))
      train.cv <- parSapply(cl, 1:rept, function(i,...) {
        set.seed(i)
        rfcv1(tt$train.x, tt$train.y, cv.fold = cv.fold, step = steps, ...)
      },simplify = F)
      #
      error.cv <- sapply(train.cv, "[[", "error.cv")
      error.cv.rm <- rowMeans(error.cv)
      error.cv.sd <- apply(error.cv,1,sd)
      id <- max(which(error.cv.rm < min(error.cv.sd + error.cv.rm)))
      min.cv <- error.cv.rm[id]
      if(min.cv < tmp.cv){ 
        tmp.train <- train.cv
        tmp.cv <- min.cv
      }
      
      marker.num <- min(as.numeric(names(error.cv.rm)[id]))
      
      mode <- data.frame(num=train.cv[[1]]$n.var,error.cv)
      mmode <- melt(mode,id.vars = "num",variable.name = "times",value.name = "Err.cv")
      mmode$num  <- as.factor(mmode$num)
      
      smode <- data.frame(num=train.cv[[1]]$n.var,times="rm",Err.cv=error.cv.rm)
      smode$num <- as.factor(smode$num)
    }
    
    print(paste0("marker.num = ",marker.num," | min.cv = ",round(min.cv,2)," | tried: ",tried))
  }
  stopCluster(cl)

  #step 2
  marker.t <- table(unlist(lapply(train.cv, function(x) {
    lapply(x$res, "[", 1:marker.num)
  })))
  marker.t <- sort(marker.t, d = T)

  names(marker.t) <- colnames(tt$train.x)[as.numeric(names(marker.t))]
  marker.p <- names(marker.t)[1:marker.num]
  if(marker.num > 1){
    marker.anno <- as.data.frame(marker.t[1:marker.num])
  }else{
    marker.anno <- data.frame(id=names(marker.t[1:marker.num]),
                              marker.t[1:marker.num])
  }
  if(exists('mgs.V2.anno2')){
    colnames(marker.anno) <- c("mgs.V2","Freq")
    marker.anno <- merge(mgs.V2.anno2,marker.anno,by="mgs.V2",all.y=T)
  }else{
    colnames(marker.anno) <- c("ID","Freq")
  }


  rf.dat <- cbind(tt$train.x[, marker.p],grp=tt$train.y)
  train.rf <- randomForest(grp~.,data=rf.dat, importance = T,proximity=T)
  train.p <- NULL
  if(class(tt$train.y)=="factor"){
    for(i in 2:rrf){
      tmp.rf <- randomForest(grp~.,data=rf.dat, importance = T,proximity=T)
      if( mean(tmp.rf$err.rate[,1]) > mean(train.rf$err.rate[,1]) ){
        train.rf <- tmp.rf
      }
    }
    train.p <- predict(train.rf, type = "prob")
    train.dat <- data.frame(grp=tt$train.y,prob=train.p[,2])
  }else{
    train.p <- predict(train.rf)
    train.dat <- data.frame(grp=tt$train.y,prob=train.p)
  }
  
  
  #

  return(list(set=tt,cv=train.cv,mf=mmode,sf=smode,marker.num=marker.num,
              marker=marker.anno,rf=train.rf,train.pred=train.dat))
}

train.module2 <- function(prf,phe,ID,grp,size,
                         max.marker=30,max.cv=0.7,thread=5,rept=5,steps=2,cv.fold=5,
                         sam=F,rrf=10,smote=F,maxtry=10,impGrp=NULL,trees=500,
                         block=NULL,candy=NULL,pick=NULL,dw=NULL,ex.PID=NULL){
  #step1
  PD.min.cv <- 1
  min.cv <- 1
  marker.num=0
  tmp.min.cv <- 1
  ntree = trees
  # parallel method to speed up the process
  cl <- makeCluster(thread)
  clusterEvalQ(cl,library(randomForest))
  clusterExport(cl,c("tt","rfcv1"))
  tried <- 0
  while((marker.num>max.marker|marker.num==0|PD.min.cv > max.cv) && tried <= maxtry){
    #aveSam(prf,phe,ID,grp,block=NULL,candy=NULL,size,pick=NULL,dw=NULL,ex.ATB=F)
    tried = tried + 1
    tt <- aveSam(prf,phe,ID,grp,NULL,candy,
                 size,pick,ex.PID=ex.PID,sam=sam,smote=smote)
    
    tt$train.y <-droplevels(tt$train.y)
    
    if(tried > maxtry){
      train.cv <- min.train.cv
      tried <- "end"
    }else{
      train.cv <- parSapply(cl, 0:rept, function(i,...) {
        set.seed(i)
        rfcv1(tt$train.x, tt$train.y, cv.fold = cv.fold, step = steps, impGrp = impGrp, ntree = trees)
      },simplify = F)
    }
    #
    if(is.null(impGrp)){
      error.cv <- sapply(train.cv, "[[", "error.cv")
    }else{
      imp.err.cv <- t(matrix(unlist(lapply(train.cv,function(x){
        sapply(x$predicted,function(y) {
          mean((y!=tt$train.y)[which(tt$train.y==impGrp)])
        })
      })),ncol=length(train.cv[[1]]$n.var),byrow=T))
      rownames(imp.err.cv) <- names(train.cv[[1]]$predicted)
      error.cv <- imp.err.cv
    }
    
    error.cv.rm <- rowMeans(error.cv)
    min.cv <- min(error.cv.rm)
    id <- which(error.cv.rm == min(error.cv.rm))
    
    marker.num <- min(as.numeric(names(error.cv.rm)[id]))
    
    mode <- data.frame(num=train.cv[[1]]$n.var,error.cv)
    mmode <- melt(mode,id.vars = "num",variable.name = "times",value.name = "Err.cv")
    mmode$num  <- as.factor(mmode$num)
    
    smode <- data.frame(num=train.cv[[1]]$n.var,times="rm",Err.cv=error.cv.rm)
    smode$num <- as.factor(smode$num)
    if(length(which(grepl("PD",levels(tt$train.y))))>0){
      PD.cv <- reCalCv(train=NULL,train.cv=train.cv,grp="PD",train.mf=mmode,train.sf=smode,train.y=tt$train.y,rep=NULL)
      PD.min.cv <- min(PD.cv$Err.cv[which(PD.cv$grp=="PD")])

      print(sprintf("rfcv (cv=%.2f | PD.cv=%.2f | mk=%3d | tried %3s)",min.cv,PD.min.cv,marker.num,tried))
    }else{
      print(sprintf("rfcv (cv=%.2f | mk=%3d | tried %3s)",min.cv,marker.num,tried))
    }
    
    if(min.cv < tmp.min.cv){
      min.train.cv <- train.cv
      tmp.min.cv <- min.cv
    }
  }
  stopCluster(cl)
  #step 2
  marker.t <- table(unlist(lapply(train.cv, function(x) {
    lapply(x$res, "[", 1:marker.num)
  })))
  marker.t <- sort(marker.t, d = T)
  
  names(marker.t) <- colnames(tt$train.x)[as.numeric(names(marker.t))]
  marker.p <- names(marker.t)[1:marker.num]
  if(marker.num > 1){
    marker.anno <- as.data.frame(marker.t[1:marker.num])
  }else{
    marker.anno <- data.frame(id=names(marker.t[1:marker.num]),
                              marker.t[1:marker.num])
  }
  if(exists('mgs.V2.anno2')){
    colnames(marker.anno) <- c("mgs.V2","Freq")
    marker.anno <- merge(mgs.V2.anno2,marker.anno,by="mgs.V2",all.y=T)
  }else{
    colnames(marker.anno) <- c("ID","Freq")
  }
  
  
  rf.dat <- cbind(tt$train.x[, marker.p,F],grp=tt$train.y)
  train.rf <- randomForest(grp~.,data=rf.dat, importance = T,proximity=T,ntree=ntree)
  for(i in 2:rrf){
    tmp.rf <- randomForest(grp~.,data=rf.dat, importance = T,proximity=T,ntree=ntree)
    if( mean(tmp.rf$err.rate[,1]) > mean(train.rf$err.rate[,1]) ){
      train.rf <- tmp.rf
    }
  }
  train.p <- predict(train.rf, type = "prob")
  train.dat <- data.frame(grp=tt$train.y,prob=train.p[,2])
  #
  
  return(list(set=tt,cv=train.cv,mf=mmode,sf=smode,marker.num=marker.num,
              marker=marker.anno,rf=train.rf,train.pred=train.dat))
}

#for more than 2 groups
train.module3 <- function(prf,phe,ID,grp,size,
                          max.marker=30,max.cv=0.7,thread=5,rept=5,steps=2,cv.fold=5,
                          sam=F,rrf=10,smote=F,maxtry=10,impGrp=NULL,trees=500,
                          block=NULL,candy=NULL,pick=NULL,dw=NULL,ex.PID=NULL,...){
  #step1
  PD.min.cv <- 1
  min.cv <- 1
  marker.num=0
  tmp.min.cv <- 1
  ntree = trees
  # parallel method to speed up the process
  cl <- makeCluster(thread,...)
  clusterEvalQ(cl,library(randomForest))
  clusterExport(cl,c("tt","rfcv1"))
  tried <- 0
  while((marker.num>max.marker|marker.num==0|min.cv > max.cv) && tried <= maxtry){
    #aveSam(prf,phe,ID,grp,block=NULL,candy=NULL,size,pick=NULL,dw=NULL,ex.ATB=F)
    tried = tried + 1
    tt <- aveSam(prf,phe,ID,grp,NULL,candy,
                 size,pick,ex.PID=ex.PID,sam=sam,smote=smote)

    tt$train.y <-droplevels(tt$train.y)
    
    if(tried > maxtry){
      train.cv <- min.train.cv
      tried <- "end"
    }else{
      train.cv <- parSapply(cl, 0:rept, function(i,...) {
        set.seed(i)
        rfcv1(tt$train.x, tt$train.y, cv.fold = cv.fold, step = steps, impGrp = impGrp, ntree = trees,...)
      },simplify = F)
    }
    #
    if(is.null(impGrp)){
      err.cv <- t(matrix(unlist(lapply(train.cv,function(x){
        lv <- levels(tt$train.y)
        res <- sapply(x$predicted,function(y) {
          #mlv <- unlist(lapply(lv,FUN=function(x) mean((y!=tt$train.y)[which(tt$train.y==x)])))
          #res <- c(mlv,mean(mlv))
          #names(res) <- c(lv,"mean")
          res <- mean(unlist(lapply(lv,FUN=function(x) mean((y!=tt$train.y)[which(tt$train.y==x)]))))
          return(res)
        })
        return(res)
      })),ncol=length(train.cv[[1]]$n.var),byrow=T))
      rownames(err.cv) <- names(train.cv[[1]]$predicted)
      error.cv <- err.cv
    }else{
      imp.err.cv <- t(matrix(unlist(lapply(train.cv,function(x){
        sapply(x$predicted,function(y) {
          mean((y!=tt$train.y)[which(tt$train.y==impGrp)])
        })
      })),ncol=length(train.cv[[1]]$n.var),byrow=T))
      rownames(imp.err.cv) <- names(train.cv[[1]]$predicted)
      error.cv <- imp.err.cv
    }
    
    error.cv.rm <- rowMeans(error.cv)
    min.cv <- min(error.cv.rm)
    id <- which(error.cv.rm == min(error.cv.rm))
    
    marker.num <- min(as.numeric(names(error.cv.rm)[id]))
    
    mode <- data.frame(num=train.cv[[1]]$n.var,error.cv)
    mmode <- melt(mode,id.vars = "num",variable.name = "times",value.name = "Err.cv")
    mmode$num  <- as.factor(mmode$num)
    
    smode <- data.frame(num=train.cv[[1]]$n.var,times="rm",Err.cv=error.cv.rm)
    smode$num <- as.factor(smode$num)
    if(length(which(grepl("PD",levels(tt$train.y))))>0&&length(levels(tt$train.y))==3){
      PD.cv <- reCalCv(train=NULL,train.cv=train.cv,grp="PD",train.mf=mmode,train.sf=smode,train.y=tt$train.y,rep=NULL)
      PD.min.cv <- min(PD.cv$Err.cv[which(PD.cv$grp=="PD")])
      SD.cv <- reCalCv(train=NULL,train.cv=train.cv,grp="SD",train.mf=mmode,train.sf=smode,train.y=tt$train.y,rep=NULL)
      SD.min.cv <- min(SD.cv$Err.cv[which(SD.cv$grp=="SD")])
      PR.cv <- reCalCv(train=NULL,train.cv=train.cv,grp="PR",train.mf=mmode,train.sf=smode,train.y=tt$train.y,rep=NULL)
      PR.min.cv <- min(PR.cv$Err.cv[which(PR.cv$grp=="PR")])
      
      print(sprintf("rfcv (cv=%.2f | PD(%.2f)|SD(%.2f)|PR(%.2f) | mk=%3d | tried %3s)",
                    min.cv,PD.min.cv,SD.min.cv,PR.min.cv,marker.num,tried))
    }else{
      print(sprintf("rfcv (cv=%.2f | mk=%3d | tried %3s)",min.cv,marker.num,tried))
    }
    
    if(min.cv < tmp.min.cv){
      min.train.cv <- train.cv
      tmp.min.cv <- min.cv
    }
  }
  stopCluster(cl)
  #step 2
  marker.t <- table(unlist(lapply(train.cv, function(x) {
    lapply(x$res, "[", 1:marker.num)
  })))
  marker.t <- sort(marker.t, d = T)
  
  names(marker.t) <- colnames(tt$train.x)[as.numeric(names(marker.t))]
  marker.p <- names(marker.t)[1:marker.num]
  if(marker.num > 1){
    marker.anno <- as.data.frame(marker.t[1:marker.num])
  }else{
    marker.anno <- data.frame(id=names(marker.t[1:marker.num]),
                              marker.t[1:marker.num])
  }
  if(exists('mgs.V2.anno2')){
    colnames(marker.anno) <- c("mgs.V2","Freq")
    marker.anno <- merge(mgs.V2.anno2,marker.anno,by="mgs.V2",all.y=T)
  }else{
    colnames(marker.anno) <- c("ID","Freq")
  }
  
  
  rf.dat <- cbind(tt$train.x[, marker.p,F],grp=tt$train.y)
  train.rf <- randomForest(grp~.,data=rf.dat, importance = T,proximity=T,...)
  for(i in 2:rrf){
    tmp.rf <- randomForest(grp~.,data=rf.dat, importance = T,proximity=T,...)
    if( mean(tmp.rf$err.rate[,1]) < mean(train.rf$err.rate[,1]) ){
      train.rf <- tmp.rf
    }
  }
  train.p <- predict(train.rf, type = "prob")
  train.dat <- data.frame(grp=tt$train.y,train.p)
  #
  
  return(list(set=tt,cv=train.cv,mf=mmode,sf=smode,marker.num=marker.num,
              marker=marker.anno,rf=train.rf,train.pred=train.dat))
}

####
cv.fold=NULL
steps=NULL
train.direct <- function(set.a=NULL,set.k,set.n,set.s,phe,
                         max.marker=30,max.cv=0.7,thread=5,rept=5,steps=2,cv.fold=5,
                         block=NULL,candy=NULL,pick=NULL,dw=NULL,ex.PID=NULL){

  #step1
  min.cv <- 1
  marker.num=0
  # parallel method to speed up the process
  cl <- makeCluster(thread)
  max.try <- 0
  if(is.null(set.a)){
    tt <- list()
    tt$train.x <- rbind(set.k$train.x,set.n$train.x,set.s$train.x)
    tt$train.y <- as.factor(levels(set.k$train.y)[c(set.k$train.y,set.n$train.y,set.s$train.y)])
    tt$test.x <- rbind(set.k$test.x,set.n$test.x,set.s$test.x)
    tt$test.y <- as.factor(levels(set.k$test.y)[c(set.k$test.y,set.n$test.y,set.s$test.y)])
  }else{
    tt <- set.a
    set.a$train.x$DNAID<- rownames(set.a$train.x)
    set.a$test.x$DNAID<- rownames(set.a$test.x)

    train.rank <- rank(rownames(set.a$train.x))
    test.rank <- rank(rownames(set.a$test.x))
    tt$train.x <- merge(phe[which(!is.na(phe$DNAID)),c("DNAID",candy)],set.a$train.x,by="DNAID")
    tt$test.x <- merge(phe[,c("DNAID",candy)],set.a$test.x,by="DNAID")

    rownames(tt$train.x) <- tt$train.x$DNAID
    tt$train.x <- tt$train.x[train.rank,-1]
    rownames(tt$test.x) <- tt$test.x$DNAID
    tt$test.x <- tt$test.x[test.rank,-1]
  }
  #clusterEvalQ(cl,source("rfcv1.R"))
  clusterEvalQ(cl,library(randomForest))
  clusterExport(cl,c("tt","rfcv1","cv.fold","steps"))
  while(marker.num>max.marker|marker.num==0|min.cv > max.cv){
    print(paste0("CV start from marker.num=",marker.num,"|min.cv=",min.cv))
    train.cv <- parSapply(cl, 1:rept, function(i,...) {
      set.seed(i)
      rfcv1(tt$train.x, tt$train.y, cv.fold=cv.fold, step = steps)
    },simplify = F)
    #

    error.cv <- sapply(train.cv, "[[", "error.cv")
    error.cv.rm <- rowMeans(error.cv)
    min.cv <- min(error.cv.rm)
    id <- which(error.cv.rm == min(error.cv.rm))

    marker.num <- min(as.numeric(names(error.cv.rm)[id]))

    mode <- data.frame(num=train.cv[[1]]$n.var,error.cv)
    mmode <- melt(mode,id.vars = "num",variable.name = "times",value.name = "Err.cv")
    mmode$num  <- as.factor(mmode$num)

    smode <- data.frame(num=train.cv[[1]]$n.var,times="rm",Err.cv=error.cv.rm)
    smode$num <- as.factor(smode$num)
  }
  stopCluster(cl)
  print(paste0("CV end at marker.num=",marker.num,"|min.cv=",min.cv))
  #step 2
  marker.t <- table(unlist(lapply(train.cv, function(x) {
    lapply(x$res, "[", 1:marker.num)
  })))
  marker.t <- sort(marker.t, d = T)

  names(marker.t) <- colnames(tt$train.x)[as.numeric(names(marker.t))]
  marker.p <- names(marker.t)[1:marker.num]
  marker.anno <- as.data.frame(marker.t[1:marker.num])
  colnames(marker.anno) <- c("mgs.V2","Freq")
  marker.anno <- merge(marker.anno,mgs.V2.anno2,by="mgs.V2",all.x=T,sort=F)
  marker.anno <- marker.anno[rev(order(marker.anno$Freq)),]


  rf.dat <- cbind(tt$train.x[, marker.p],grp=tt$train.y)
  train.rf <- randomForest(grp~.,data=rf.dat, importance = T,proximity=T)
  train.p <- predict(train.rf, type = "prob")
  train.dat <- data.frame(grp=tt$train.y,prob=train.p[,2])
  #

  return(list(set=tt,cv=train.cv,mf=mmode,sf=smode,marker.num=marker.num,
              marker=marker.anno,rf=train.rf,train.pred=train.dat))
}


filter.repeat.marker <- function(model,rep.start,rep.end,num=NULL){
  train.cv <- model$cv[rep.start:rep.end]
  tt <- model$set
  
  if(is.null(num)){
    error.cv <- sapply(train.cv, "[[", "error.cv")
    error.cv.rm <- rowMeans(error.cv)
    min.cv <- min(error.cv.rm)
    id <- which(error.cv.rm == min(error.cv.rm))
    
    marker.num <- min(as.numeric(names(error.cv.rm)[id]))
  }else{
    marker.num <- num
  }
  marker.t <- table(unlist(lapply(train.cv, function(x) {
    lapply(x$res, "[", 1:marker.num)
  })))
  marker.t <- sort(marker.t, decreasing = T)

  names(marker.t) <- colnames(tt$train.x)[as.numeric(names(marker.t))]
  marker.p <- names(marker.t)[1:marker.num]
  marker.anno <- as.data.frame(marker.t[1:marker.num])
  colnames(marker.anno) <- c("mgs.V2","Freq")
  marker.anno <- merge(marker.anno,mgs.V2.anno2,by="mgs.V2",all.x=T,sort=F)
  marker.anno <- marker.anno[rev(order(marker.anno$Freq)),]
  return(marker.anno)
}

reCalCv <- function(train=NULL,train.cv=NULL,grp,train.mf=NULL,train.sf=NULL,train.y=NULL,rep=NULL,more=F){
  if(!is.null(train)){
    train.cv=train$cv
    train.mf=train$mf
    train.sf=train$sf
    train.y = train$set$train.y
  }
  rep <- ifelse(is.null(rep),length(train.cv),rep)
  mf.stat <- ddply(train.mf,"num",summarise,
                   Err.upper=quantile(Err.cv,.975),
                   Err.lower=quantile(Err.cv,.025),
                   Err.SD=sd(Err.cv))
  dat <- data.frame(grp="all",merge(train.sf,mf.stat,by="num"))
  for(i in grp){
    this.err.cv <- as.data.frame(matrix(unlist(lapply(train.cv,function(x){
      sapply(x$predicted,function(y) {
        mean((y!=train.y)[which(train.y==i)]) #[which(NSCLC.cln.train.cv$set$train.cv.y=="PID")]
      })
    })),ncol=length(train.cv[[1]]$n.var),byrow=T))
    this.err.cv <- this.err.cv[1:rep,]
    colnames(this.err.cv) <- train.cv[[1]]$n.var
    Emean <- colMeans(this.err.cv)
    Err.upper <- apply(this.err.cv,2,quantile,.975)
    Err.lower <- apply(this.err.cv,2,quantile,.025)
    Err.SD <- apply(this.err.cv,2,sd)
    this.err.cv.mean <- data.frame(grp=i,num=colnames(this.err.cv),times="mean",
                                 Err.cv=Emean,Err.upper=Err.upper,Err.lower=Err.lower,Err.SD=Err.SD)
    this.err.cv.mdat <- melt(data.frame(times=rownames(this.err.cv),this.err.cv),id.vars = "times",
                           variable.name = "num",value.name = "Err.cv")
    levels(this.err.cv.mdat$num) <- train.cv[[1]]$n.var
    dat <- rbind(dat,this.err.cv.mean)
    if(i=="PD"){the.this.err.cv <- this.err.cv}
  }
  if(more){
    return(list(cv=dat,detail=the.PD.err.cv))
  }else{
    return(dat)
  }
  #return(list(mdat=PD.err.cv.mdat,mean=PD.err.cv.mean))
}

optimal.num <- function(dat){
  num <- min(as.numeric(levels(dat$num)[dat$num[which(dat$Err.cv<=(min(dat$Err.cv)+dat$Err.SD[which(dat$Err.cv==min(dat$Err.cv))]))]]))
  pos <- list(x=num,y=dat$Err.cv[which(dat$num==num)])
  return(pos)
}


reCalFreq <- function(train,mk.num){
  dat <- sort(table(unlist(lapply(train$cv, function(x){lapply(x$res, "[", 1:mk.num)}))),d=T)[1:mk.num]
  dat <- data.frame(ID=colnames(train$set$train.x)[as.numeric(names(dat))],
                    Freq=as.numeric(dat))
  return(dat)
}

mkRep <- function(obs,pred){
  grp <- levels(obs)
  dat <- NULL
  lst <- list()
  auc <- NULL
  for( i in grp){
    ind <- which(colnames(pred)==i)
    iPred <- pred[,ind]
    iObs  <- ifelse(obs==i,1,0)
    iRoc  <- roc(iObs,iPred)
    iDat  <- data.frame(
      grp = i,
      specificity = iRoc$specificities,
      sensitivity = iRoc$sensitivities
    )
    lst <- c(lst,list(iRoc))
    auc <- c(auc,iRoc$auc)
    iDat <- iDat[rev(order(iDat$specificity)),]
    dat <- rbind(dat,iDat)
  }
  names(auc) <- grp
  names(lst) <- grp
  return(list(dat=dat,roc=lst,auc=auc))
}

plot.f1 <- function(model){
  model$mf$num <- factor(model$mf$num,levels=rev(levels(model$mf$num)))
  model$sf$num <- factor(model$sf$num,levels=rev(levels(model$sf$num)))
  ggplot(model$mf[which(model$mf$num!=0),],aes(x=num,y=Err.cv,group=times)) +
    #geom_line(aes(group=times),linetype=2,alpha=.1) +
    geom_boxplot(aes(group=num)) +
    geom_line(data=model$sf[which(model$sf$num!=0),],size=1,color="red") +
    geom_vline(xintercept = which(model$sf$num==model$marker.num),color="pink",size=1) +
    theme_bw()
}
#Fin.
