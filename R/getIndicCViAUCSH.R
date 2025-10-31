getIndicCViAUCSH = function(lp,lpnew,Surv.rsp,Surv.rsp.new,times.auc=seq(10,1000,10),times.prederr=1:500,train.fit,plot.it=FALSE,tmax.train=max(Surv.rsp[,"time"][ object$Surv.rsp[,"status"] == 1 ]),tmax.test=max(Surv.rsp.new[,"time"][ object$Surv.rsp.new[,"status"] == 1 ])){
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
  object$AUC_sh <- survAUC::AUC.sh(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.auc)  #CoxModel
  
  if(plot.it){
    plot(object$AUC_sh)
    abline(h = 0.5)
  }
  
  return(object)
}
