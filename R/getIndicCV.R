getIndicCV = function(lp,lpnew,Surv.rsp,Surv.rsp.new,times.auc=seq(10,1000,10),times.prederr=1:500,train.fit,plot.it=FALSE,tmax.train=max(Surv.rsp[,"time"][ object$Surv.rsp[,"status"] == 1 ]),tmax.test=max(Surv.rsp.new[,"time"][ object$Surv.rsp.new[,"status"] == 1 ])){
  try(attachNamespace("survival"),silent=TRUE)
  #on.exit(try(unloadNamespace("survival"),silent=TRUE))
  try(attachNamespace("risksetROC"),silent=TRUE)
  on.exit(try(unloadNamespace("risksetROC"),silent=TRUE))
  try(attachNamespace("survcomp"),silent=TRUE)
  on.exit(try(unloadNamespace("survcomp"),silent=TRUE))
  #  library(survAUC)
  object <- NULL
  
  object$lp <- lp
  object$lpnew <- lpnew
  object$Surv.rsp <- Surv.rsp
  object$Surv.rsp.new <- Surv.rsp.new
  object$times.auc <- times.auc
  object$train.fit <- train.fit
  object$test.fit <- survival::coxph(object$Surv.rsp.new~object$lpnew, iter.max=0, init=1)
  
  object$nulltrain.fit <- survival::coxph(object$Surv.rsp~1)
  object$lp0 <- predict(object$nulltrain.fit)
  object$nulltest.fit <- survival::coxph(object$Surv.rsp.new~1)
  object$lp0new <- predict(object$nulltest.fit)
  
  object$tmax.train <- tmax.train
  object$tmax.test <- tmax.test
  
  # iAUC
  object$AUC_CD <- survAUC::AUC.cd(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.auc)  #CoxModel
  object$AUC_hc <- survAUC::AUC.hc(object$Surv.rsp, object$Surv.rsp.new, object$lpnew, object$times.auc)      #No model
  object$AUC_sh <- survAUC::AUC.sh(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.auc)  #CoxModel
  object$AUC_Uno <- list(auc=rep(0,length(object$times.auc)),times=times.auc,iauc=0)
  #class(object$AUC_Uno) <- "survAUC"
  #if(var(object$lpnew)>1e-8){try(
  object$AUC_Uno <- survAUC::AUC.uno(object$Surv.rsp, object$Surv.rsp.new, object$lpnew, object$times.auc)#  )}   #No model
  
  #AUC_CD$iauc
  #AUC_hc$iauc
  #AUC_sh$iauc
  #AUC_Uno$iauc
  
  if(plot.it){
    layout(matrix(1:4,nrow=2,byrow=TRUE))
    plot(object$AUC_CD)
    abline(h = 0.5)
    plot(object$AUC_hc)
    abline(h = 0.5)
    plot(object$AUC_sh)
    abline(h = 0.5)
    plot(object$AUC_Uno)
    abline(h = 0.5)
  }
  
  #iAUC HZ train.set
  #  library(risksetROC)
  ## first find the estimated survival probabilities at unique failure times
  object$surv.prob.train <- unique(survival::survfit(object$Surv.rsp~1)$surv)
  object$eta.train <- predict(object$train.fit)
  #model.score <- eta.train
  
  object$utimes.train <- unique( object$Surv.rsp[,"time"][ object$Surv.rsp[,"status"] == 1 ] )
  object$utimes.train <- object$utimes.train[ order(object$utimes.train) ]
  
  ## find AUC at unique failure times
  object$AUC_hz.train <- NULL
  object$AUC_hz.train$auc <- rep( NA, length(object$utimes.train) )
  for( j in 1:length(object$utimes.train) )
  {
    object$out.train <- risksetROC::CoxWeights(object$eta.train, object$Surv.rsp[,"time"], object$Surv.rsp[,"status"], object$utimes.train[j])
    object$AUC_hz.train$auc[j] <- object$out.train$AUC
  }
  ## integrated AUC to get concordance measure
  object$AUC_hz.train$times <- object$utimes.train
  object$AUC_hz.train$iauc <- risksetROC::IntegrateAUC( object$AUC_hz.train$auc, object$utimes.train, object$surv.prob.train, tmax=object$tmax.train)
  class(object$AUC_hz.train) <- "survAUC"
  if(plot.it){
    layout(matrix(1:4,nrow=2))
    plot(object$AUC_hz.train)
    abline(h = 0.5)
  }
  
  #iAUC HZ test.set
  #  library(risksetROC)
  ## first find the estimated survival probabilities at unique failure times
  object$surv.prob.test <- unique(survival::survfit(object$Surv.rsp.new~1)$surv)
  object$eta.test <- predict(object$test.fit)
  #model.score <- eta.test
  
  object$utimes.test <- unique( object$Surv.rsp.new[,"time"][ object$Surv.rsp.new[,"status"] == 1 ] )
  object$utimes.test <- object$utimes.test[ order(object$utimes.test) ]
  
  ## find AUC at unique failure times
  object$AUC_hz.test <- NULL
  object$AUC_hz.test$auc <- rep( NA, length(object$utimes.test) )
  for( j in 1:length(object$utimes.test) )
  {
    object$out.test <- risksetROC::CoxWeights( object$eta.test, object$Surv.rsp.new[,"time"], object$Surv.rsp.new[,"status"], object$utimes.test[j])
    object$AUC_hz.test$auc[j] <- object$out.test$AUC
  }
  ## integrated AUC to get concordance measure
  object$AUC_hz.test$times <- object$utimes.test
  object$AUC_hz.test$iauc <- NA
  if(length(object$utimes.test)>1){
    object$AUC_hz.test$iauc <- risksetROC::IntegrateAUC( object$AUC_hz.test$auc, object$utimes.test, object$surv.prob.test, tmax=object$tmax.test )
  } else {if(length(object$utimes.test)==1){object$AUC_hz.test$iauc<-object$AUC_hz.test$auc}}
  class(object$AUC_hz.test) <- "survAUC"
  if(plot.it){
    plot(object$AUC_hz.test)
    abline(h = 0.5)
  }
  
  ##time-dependent ROC curves
  #  library(survcomp)
  ##train
  mytdroc.train <- NULL
  mytdroc.train <- NULL
  object$AUC_survivalROC.train <- NULL
  object$AUC_survivalROC.train$auc <- rep(NA,length(object$utimes.train))
  object$AUC_survivalROC.train$iauc <- NA
  object$AUC_survivalROC.train$times <- object$utimes.train
  class(object$AUC_survivalROC.train) <- "survAUC"
  train.cc.ix <- complete.cases(object$lp, object$Surv.rsp[,"time"], object$Surv.rsp[,"status"], NULL)
  train.surv.event.cc.ix <- object$Surv.rsp.new[,"status"][train.cc.ix]
  if (all(sort(unique(train.surv.event.cc.ix)) == c(0, 1))) {
    for(i in 1:length(object$utimes.train)) {
      rr.train <- survcomp::tdrocc(x=object$lp, surv.time=object$Surv.rsp[,"time"], surv.event=object$Surv.rsp[,"status"], time=object$utimes.train[i], na.rm=TRUE, verbose=FALSE)
      mytdroc.train <- c(mytdroc.train, list(rr.train))
    }
    object$AUC_survivalROC.train$auc <- unlist(lapply(mytdroc.train, function(x) { return(x$AUC) }))
    cc.ix.train <- complete.cases(object$AUC_survivalROC.train$auc)
    auc.survivalROC.train.cc <- object$AUC_survivalROC.train$auc[cc.ix.train]
    time.train.cc <- object$utimes.train[cc.ix.train]
    if(length(time.train.cc)>0){
      diffs.train.cc <- c(time.train.cc[1], time.train.cc[2:length(time.train.cc)] - time.train.cc[1:(length(time.train.cc) - 1)])
      object$AUC_survivalROC.train$iauc <- sum(diffs.train.cc * auc.survivalROC.train.cc)/max(time.train.cc)
      if(plot.it){
        plot(object$AUC_survivalROC.train)
        abline(h = 0.5)
      }
    }
  }
  
  
  #  library(survcomp)
  ##test
  mytdroc.test <- NULL
  object$AUC_survivalROC.test <- NULL
  object$AUC_survivalROC.test$auc <- rep(NA,length(object$utimes.test))
  object$AUC_survivalROC.test$iauc <- NA
  object$AUC_survivalROC.test$times <- object$utimes.test
  class(object$AUC_survivalROC.test) <- "survAUC"
  test.cc.ix <- complete.cases(object$lpnew, object$Surv.rsp.new[,"time"], object$Surv.rsp.new[,"status"], NULL)
  test.surv.event.cc.ix <- object$Surv.rsp.new[,"status"][test.cc.ix]
  if (all(sort(unique(test.surv.event.cc.ix)) == c(0, 1))) {
    for(i in 1:length(object$utimes.test)) {
      rr.test <- survcomp::tdrocc(x=object$lpnew, surv.time=object$Surv.rsp.new[,"time"], surv.event=object$Surv.rsp.new[,"status"], time=object$utimes.test[i], na.rm=TRUE, verbose=FALSE)
      mytdroc.test <- c(mytdroc.test, list(rr.test))
    }
    object$AUC_survivalROC.test$auc <- unlist(lapply(mytdroc.test, function(x) { return(x$AUC) }))
    cc.ix.test <- complete.cases(object$AUC_survivalROC.test$auc)
    auc.survivalROC.test.cc <- object$AUC_survivalROC.test$auc[cc.ix.test]
    time.test.cc <- object$utimes.test[cc.ix.test]
    if(length(time.test.cc)>0){
      diffs.test.cc <- c(time.test.cc[1], time.test.cc[2:length(time.test.cc)] - time.test.cc[1:(length(time.test.cc) - 1)])
      object$AUC_survivalROC.test$iauc <- sum(diffs.test.cc * auc.survivalROC.test.cc)/max(time.test.cc)
      if(plot.it){
        plot(object$AUC_survivalROC.test)
        abline(h = 0.5)
      }
    }
  }
  #Surv.rsp A survival::Surv(.,.) object containing to the outcome of the test data.
  #lp The vector of predictors.
  #lp0 The vector of predictors obtained from the covariate-free null model.
  
  #prederr #ierror for iBrier
  object$prederr <- NULL
  object$times.prederr <- times.prederr[times.prederr>1]
  object$prederr$brier.unw <- survAUC::predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.prederr, type = "brier", int.type = "unweighted")
  object$prederr$robust.unw <- survAUC::predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.prederr, type = "robust", int.type = "unweighted")
  object$prederr$brier.w <- survAUC::predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.prederr, type = "brier", int.type = "weighted")
  object$prederr$robust.w <- survAUC::predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.prederr, type = "robust", int.type = "weighted")
  
  if(plot.it){
    layout(matrix(1:4,nrow=2))
    plot(object$prederr$brier.unw)
    abline(h = 0.25)
    plot(object$prederr$robust.unw)
    abline(h = 0.25)
    plot(object$prederr$brier.w)
    abline(h = 0.25)
    plot(object$prederr$robust.w)
    abline(h = 0.25)
  }
  
  return(object)
}