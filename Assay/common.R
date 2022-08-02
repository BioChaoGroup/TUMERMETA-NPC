#Common self defined functions

# check md5 of input data
### A function to validate file version by checking md5.
md5 <- function(f,m){
  if(tools::md5sum(f)==m){
    warning("MD5 check passed.\n")
    return(f)
  }else{
    stop("MD5 check failed. Please check the version and content!\n")
  }
}

# output control # Sometime we just need to render a Rmd with figures and mediate
                 # variables, there is no need to output the results again.
do.write.csv <- function(do=T,...){if(do){write.csv(...)}}
do.write.table <- function(do=T,...){if(do){write.table(...)}}

#
phaseCheckFun <- function(d,info,type=NULL){
  pid <- as.character(d$PID[1])
  tdat <- merge(
    data.frame(d[,c("Date","DNAID","PID")]),
    filter(info,PID==pid)[,c("Date","Times","tDay","Medicine")],
    by="Date",all=T
  )
  tdat <- tdat[order(tdat$Date),] #MAKE SURE it's ordered
  tdat$TTimes <- tdat$Times
  tdat$WTimes <- 0
  tdat$during <- "tail" # set default stage as tail
  tdat$duringDetail <- NA
  tdat$fix <- F
  t=0  #treatment phase aligned with treated times
  tt=0 #treatment phase aligned with weeks passed
  tD <- min(tdat$Date[which(tdat$Times==1)])
  date0 <- tD
  k=1
  for(i in 1:nrow(tdat)){
    theT <-  tdat$Times[i]
    if(!is.na(theT)){
      if(theT <= t + k){
        if(t>0&&tdat$Date[i]-tD<3){
          #treatment time too close, regard as once
          theT <- t
          k = k + 1
        }
        tD <- tdat$Date[i]
        tdat$duringDetail[which(tdat$TTimes==t & tD - tdat$Date<=7)] <- 
          tdat$Date[which(tdat$TTimes==t & tD - tdat$Date<=7)] - tD
        if(t==0){
          # extend previous 2weeks as baseline
          tdat$during[which(tdat$TTimes==t & tD - tdat$Date<=14)] <- "tail" 
        }else{
          tdat$during[which(tdat$TTimes==t & tD - tdat$Date<=7)] <- "tail"
        }
        t = theT
        if(tdat$Medicine[i]=="Keytruda"){
          if(tdat$tDay[i]>42&&tdat$tDay[i]<56){
            tt = 5
            #}else if(tdat$tDay[i]>84&&tdat$tDay[i]<98){
            #  tt = 8
          }else{
            tt = round(tdat$tDay[i]/14,0) + 1
          }
          
        }else{
          tt = round(tdat$tDay[i]/14,0) + 1
        }
      }else{
        print(tdat)
        stop(paste0("The ",pid," misordered with Times: ",tdat$Times[i]," at ",tdat$Date[i]))
      }
      tdat$TTimes[i] <- t
      tdat$WTimes[i] <- tt
      if(!is.na(tdat$DNAID[i])){
        tdat$TTimes[i] <- tdat$TTimes[i] - 1
        tdat$WTimes[i] <- tdat$WTimes[i] - 1
        
        tdat$during[i] <- "tail"
      }
    }else{
      tdat$TTimes[i] <- t
      tdat$WTimes[i] <- tt
      tt2 = round((tdat$Date[i] - date0)/14,2)
      if(tt2>tt+2){
        warning(paste0("The ",pid," got extended WTimes: ",tt2," from ",tt," at ",tdat$Date[i] - date0,
                       ". Ignored.\n"))
        tdat$WTimes[i] <- NA
      }else if(tt2>tt+0.7){
        warning(paste0("The ",pid," got extended WTimes: ",tt2," from ",tt," at ",tdat$Date[i] - date0,
                       ". fixed.\n"))
        tdat$WTimes[i] <- tt+1
      }
    }
    
    pDay <- tdat$Date[i]-tD
    if(pDay<=7){
      tdat$duringDetail[i] <- pDay
      if(pDay>0){tdat$during[i]<-"head"}
    }
    
  }
  if(is.null(type)){
    return(filter(tdat,!is.na(DNAID)))
  }else if(type=="treat"){
    return(filter(tdat,is.na(DNAID)))
  }
}

# choose the sampling stage according to the days after the first treatment
chooseCT <- function(x){
  x <- x[order(x$Day),]
  CTgrp <- rep(NA,nrow(x))
  
  # point selection function
  fun_Month <- function(d,l,CTgrp){
    t <- min(abs(x$Day-d))
    if(t<=14){
      CTgrp[which(abs(x$Day-d)==min(abs(x$Day-d)))[1]] <- l
    }
    return(CTgrp)
  }
  
  #choose 3 month(8 weeks) timepoint
  CTgrp <- fun_Month(84,"treat.3m",CTgrp)
  
  #choose 2 month(8 weeks) timepoint
  CTgrp <- fun_Month(56,"treat.2m",CTgrp)
  
  # choose 1 month(4 weeks) timepoint
  CTgrp <- fun_Month(28,"treat.1m",CTgrp)
  
  # choose baseline
  t <- which(x$Day<0)
  if(length(t)>0){
    CTgrp[which(x$Day==max(x$Day[t]))[1]] <- "baseline"
  }else if(min(x$Day) <= 3){
    CTgrp[which(x$Day==min(x$Day))[1]] <- "baseline"
  }
  
  ##merge
  x$CTgrp <- CTgrp
  return(x)
}


# define the best evaluation
best <- function(dat){
  dat <- dat[order(dat$Day),]
  if(length(grep("FD", dat$Evaluation))>0){
    pdd <- (dat$Day[which(dat$Evaluation=="FD")])[1]
  }else if(length(grep("PD", dat$Evaluation))>0){
    pdd <- (dat$Day[which(dat$Evaluation=="PD")])[1]
  }else{
    pdd <- NA
  }
  fb <- NULL
  bT <- NULL
  if(length(grep("PR", dat$Evaluation))>0){
    fb <- "PR"
    r3 <- "R"
    r6 <- "R"
    bT <- grep("PR",dat$Evaluation)
  }else if(length(grep("SD", dat$Evaluation))>0){
    if(length(grep("PD", dat$Evaluation))>0){
      PD1st <- (dat$Day[dat$Evaluation=="PD"])[1]
      if(PD1st>180){
        r3 <- 'R'
        r6 <- 'R'
      }else if (PD1st>90){
        r3 <- 'R'
        r6 <- 'NR'
      }else{
        r3 <- 'NR'
        r6 <- 'NR'
      }
    }else{
      ### Caution ### 4 patients are still missing PD records. 
      ### Following is a temporaty treat.
      r3 <- 'R'#"SD"
      r6 <- 'R'#"SD"
    }
    fb <- "SD"
    bT <- grep("SD",dat$Evaluation)
  }else if(length(grep("PD", dat$Evaluation))>0){
    fb <- "PD"
    r3 <- "NR"
    r6 <- "NR"
    bT <- grep("PD", dat$Evaluation)
  }else if(length(grep("FD", dat$Evaluation))>0){
    fb <- "FD"
    r3 <- "NR"
    r6 <- "NR"
    bT <- grep("FD", dat$Evaluation)
  } else{
    fb <- NA
    r3 <- NA
    r6 <- NA
    bT <- NA
  }
  dat$BestEvaluation <- fb
  dat$PFS3m <- r3
  dat$PFS6m <- r6
  dat$BestDay <- (dat$Day[bT])[1]
  dat$PD.Day <- pdd
  return(dat)
}

#FPS data generation (for 2 classed group)
makePFSdata <- function(dat,grp,maxDay){
  PFS <- dat
  PFS$grp <- PFS[,grp]

  ####
  pNum <- table(PFS$grp)
  PFSD <- data.frame(Day=seq(0,maxDay),G0=pNum[1],G1=pNum[2])
  PFSD <- melt(PFSD,id.vars = "Day", variable.name = "Group",value.name = "Num")
  PFSD$change <- F
  PFSD$PID <- NA
  for(i in 1:nrow(PFS)){
    if(PFS$grp[i]==names(pNum)[2]){
      PFSD$Num[which(PFSD$Group=="G1"&PFSD$Day>=PFS$PD.Day[i])] <- 
        PFSD$Num[which(PFSD$Group=="G1"&PFSD$Day>=PFS$PD.Day[i])] - 1
      PFSD$change[which(PFSD$Group=="G1"&PFSD$Day==PFS$PD.Day[i])] <- T
      PFSD$PID[which(PFSD$Group=="G1"&PFSD$Day==PFS$PD.Day[i])] <- as.character(PFS$PID[i])
    }else{
      PFSD$Num[which(PFSD$Group=="G0"&PFSD$Day>=PFS$PD.Day[i])] <- 
        PFSD$Num[which(PFSD$Group=="G0"&PFSD$Day>=PFS$PD.Day[i])] - 1
      PFSD$change[which(PFSD$Group=="G0"&PFSD$Day==PFS$PD.Day[i])] <- T
      PFSD$PID[which(PFSD$Group=="G0"&PFSD$Day==PFS$PD.Day[i])] <- as.character(PFS$PID[i])
    }
  }
  
  PFSD$PFS <- 0
  PFSD$PFS[which(PFSD$Group=="G0")] <- PFSD$Num[which(PFSD$Group=="G0")]/pNum[1]
  PFSD$PFS[which(PFSD$Group=="G1")] <- PFSD$Num[which(PFSD$Group=="G1")]/pNum[2]
  levels(PFSD$Group) <- names(pNum)
  colnames(PFSD)[which(colnames(PFSD)=="Group")] <- grp
  return(PFSD)
}

#function to calculate PERMANOVA 
# > * x: distance matrix
# > * y: phenotype
# > * z: Number of levels no more than level change as factor, such as: 0, 1

PermFun <- function(x, y, z){
  x.rn <- rownames(x)
  y.rn <- rownames(y)
  xy.rn <- intersect(x.rn, y.rn)
  x <- x[xy.rn, xy.rn]
  y <- y[xy.rn, ,F]
  y.cm <- colnames(y)
  y.nc <- ncol(y)
  out <- data.frame()
  for(i in 1:y.nc) {  
    phe <- y[, i]
    flag <- which(!is.na(phe))
    len <- length(flag)
    phe <- phe[flag]
    if (len == 0|length(unique(phe)) == 1) {
      tmp <- cbind(data.frame(phe=y.cm[i],len=len,Df=NA,
                              SumsOfSqs=NA, MeanSqs=NA, F.Model=NA, R2=NA, NA))
      colnames(tmp) <- colnames(out)
      out <- rbind(out,tmp)
      next
    }
    set.seed(0)
    dist <- as.dist(x[flag, flag])
    if(length(unique(phe)) <= 5) {
      phe <- as.factor(phe)
    }
    res <- adonis(dist ~ phe, permutations = 1000)
    tmp1 <- cbind(y.cm[i], len, res$aov.tab[1, ])
    out <- rbind(out, tmp1)
  }
  colnames(out) <- c("phenotype", "SampleNum", "Df", "SumsOfSqs", "MeanSqs", "F.Model", "R2", "Pr(>F)")
  rownames(out) <- NULL
  out$FDR <- p.adjust(out$`Pr(>F)`,method = "BH")
  return(out)
}


is_outlier <- function(x) {
  ind <- which(!is.na(x))
  x <- x[ind]
  out <- which(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
  return(ind[out])
}