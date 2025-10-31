#' @title Internal bigPLScox functions
#' 
#' @name internal-bigPLScox
#' 
#' @description These are not to be called by the user.
#' 
#' @aliases ust spls.dv correctp.cox getIndic getIndicCV getIndicCViAUCSH 
#' getIndicCViAUCSurvROCTest NULL 
#' 
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@lecnam.net}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @references plsRcox, Cox-Models in a high dimensional setting in R, Frederic
#' Bertrand, Philippe Bastien, Nicolas Meyer and Myriam Maumy-Bertrand (2014).
#' Proceedings of User2014!, Los Angeles, page 152.\cr
#' 
#' Deviance residuals-based sparse PLS and sparse kernel PLS regression for
#' censored data, Philippe Bastien, Frederic Bertrand, Nicolas Meyer and Myriam
#' Maumy-Bertrand (2015), Bioinformatics, 31(3):397-404,
#' doi:10.1093/bioinformatics/btu660.
#' @keywords internal
NULL

logplik = function (x, time, status, b, method = c("breslow", "efron"), return.all = FALSE) {
  method <- match.arg(method)
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(b)) {
    b <- as.matrix(b)
  }
  res <- cox_partial_loglik_cpp(x, as.numeric(time), as.numeric(status), b, method, return.all)
  res
}

getmin2 = function (lambda, cvm, cvsd) {
  cvmin = min(cvm, na.rm = TRUE)
  idmin = cvm <= cvmin
  lambda.min = max(lambda[idmin], na.rm = TRUE)
  idminl = match(lambda.min, lambda)
  semin = (cvm + cvsd)[idminl]
  idmin2 = cvm >= semin
  #    if(lambda.min==-Inf){
  #    lambda.1se = -Inf
  #    } else {
  idmin2[idminl:length(idmin2)] = FALSE 
  lambda.1se = min(c(lambda[idmin2],min(lambda)), na.rm = TRUE)
  #    }
  list(lambda.min = lambda.min, lambda.1se = lambda.1se)
}


correctp.cox=function (x, y, eta, K, kappa, select, fit, verbose=FALSE)  {
  force(K)
  if (min(eta) < 0 | max(eta) >= 1) {
    if (max(eta) == 1) {
      stop("eta should be strictly less than 1!")
    }
    if (length(eta) == 1) {
      stop("eta should be between 0 and 1!")
    }
    else {
      stop("eta should be between 0 and 1! \n  Choose appropriate range of eta!")
    }
  }
  if (max(K) > ncol(x)) {
    stop("K cannot exceed the number of predictors! Pick up smaller K!")
  }
  if (max(K) >= nrow(x)) {
    stop("K cannot exceed the sample size! Pick up smaller K!")
  }
  if (min(K) <= 0 | !all(K%%1 == 0)) {
    if (length(K) == 1) {
      stop("K should be a positive integer!")
    }
    else {
      stop("K should be a positive integer! \n  Choose appropriate range of K!")
    }
  }
  if (kappa > 0.5 | kappa < 0) {
    if(verbose){cat("kappa should be between 0 and 0.5! kappa=0.5 is used. \n\n")}
    kappa <- 0.5
  }
  if (select != "pls2" & select != "simpls") {
    if(verbose){cat("Invalid PLS algorithm for variable selection.\n")}
    if(verbose){cat("pls2 algorithm is used. \n\n")}
    select <- "pls2"
  }
  fits <- c("regression", "canonical", "invariant", "classic")
  if (!any(fit == fits)) {
    if(verbose){cat("Invalid PLS algorithm for model fitting\n")}
    if(verbose){cat("regression algorithm is used. \n\n")}
    fit <- "regression"
  }
  list(K = K, eta = eta, kappa = kappa, select = select, fit = fit)
}

spls.dv <- function (Z, eta, kappa, eps, maxstep) 
{
  p <- nrow(Z)
  q <- ncol(Z)
  Znorm1 <- median(abs(Z))
  Z <- Z/Znorm1
  if (q == 1) {
    c <- ust(Z, eta)
  }
  if (q > 1) {
    M <- Z %*% t(Z)
    dis <- 10
    i <- 1
    if (kappa == 0.5) {
      c <- matrix(10, p, 1)
      c.old <- c
      while (dis > eps & i <= maxstep) {
        mcsvd <- svd(M %*% c)
        a <- mcsvd$u %*% t(mcsvd$v)
        c <- ust(M %*% a, eta)
        dis <- max(abs(c - c.old))
        c.old <- c
        i <- i + 1
      }
    }
    if (kappa > 0 & kappa < 0.5) {
      kappa2 <- (1 - kappa)/(1 - 2 * kappa)
      c <- matrix(10, p, 1)
      c.old <- c
      h <- function(lambda) {
        alpha <- solve(M + lambda * diag(p)) %*% M %*% 
          c
        obj <- t(alpha) %*% alpha - 1/kappa2^2
        return(obj)
      }
      if (h(eps) * h(1e+30) > 0) {
        while (h(eps) <= 1e+05) {
          M <- 2 * M
          c <- 2 * c
        }
      }
      while (dis > eps & i <= maxstep) {
        if (h(eps) * h(1e+30) > 0) {
          while (h(eps) <= 1e+05) {
            M <- 2 * M
            c <- 2 * c
          }
        }
        lambdas <- uniroot(h, c(eps, 1e+30))$root
        a <- kappa2 * solve(M + lambdas * diag(p)) %*% 
          M %*% c
        c <- ust(M %*% a, eta)
        dis <- max(abs(c - c.old))
        c.old <- c
        i <- i + 1
      }
    }
  }
  return(c)
}


# Helper to coerce prediction inputs while supporting big.matrix objects
.bigPLScox_coerce_newdata <- function(newdata, p) {
  if (inherits(newdata, "big.matrix")) {
    if (!requireNamespace("bigmemory", quietly = TRUE)) {
      stop("Package 'bigmemory' is required to handle big.matrix inputs")
    }
    newdata <- bigmemory::as.matrix(newdata[, , drop = FALSE])
  }
  if (is.data.frame(newdata)) {
    newdata <- data.matrix(newdata)
  }
  if (is.null(dim(newdata))) {
    if (length(newdata) != p)
      stop("'newdata' must be a numeric matrix with ncol = ",
           p, " or a vector of length = ", p, ".")
    dim(newdata) <- c(1, p)
  } else {
    newdata <- as.matrix(newdata)
    if (ncol(newdata) != p)
      stop("'newdata' must be a numeric matrix with ncol = ",
           p, " or a vector of length = ", p, ".")
  }
  storage.mode(newdata) <- "double"
  newdata
}

getIndic = function(lp,lpnew,Surv.rsp,Surv.rsp.new,times.auc=seq(10,1000,10),times.prederr=1:500,train.fit,train.fit.cph,tmax.train=365,tmax.test=365,TR,TE,plot.it=TRUE){
  try(attachNamespace("survival"),silent=TRUE)
  #on.exit(try(unloadNamespace("survival"),silent=TRUE))
  try(attachNamespace("risksetROC"),silent=TRUE)
  on.exit(try(unloadNamespace("risksetROC"),silent=TRUE))
  try(attachNamespace("survcomp"),silent=TRUE)
  on.exit(try(unloadNamespace("survcomp"),silent=TRUE))
  try(attachNamespace("rms"),silent=TRUE)
  on.exit(try(unloadNamespace("rms"),silent=TRUE))
  
  object <- NULL
  
  object$tmax.train <- tmax.train
  object$tmax.test <- tmax.test
  
  object$TR <- TR
  object$TE <- TE
  
  object$lp <- lp
  object$lpnew <- lpnew
  object$Surv.rsp <- Surv.rsp
  object$Surv.rsp.new <- Surv.rsp.new
  object$times.auc <- times.auc
  object$times.prederr <- times.prederr
  object$train.fit <- train.fit
  object$train.fit.cph <- train.fit.cph
  
  object$nulltrain.fit <- survival::coxph(object$Surv.rsp~1)
  object$lp0 <- predict(object$nulltrain.fit)
  object$nulltest.fit <- survival::coxph(object$Surv.rsp.new~1)
  object$lp0new <- predict(object$nulltest.fit)
  object$test.fit <- survival::coxph(object$Surv.rsp.new~object$lpnew, iter.max=0, init=1)
  
  object$Erestrain <- predict(train.fit, type='expected')
  object$Erestest <- predict(train.fit, newdata=TE, type='expected')
  
  #predicted Martingale resid from coxph
  object$mtrainresid <- Surv.rsp[,"status"] - object$Erestrain
  object$mtestresid <- Surv.rsp.new[,"status"] - object$Erestest
  #then var of them
  object$var_mtrainresid <- var(object$mtrainresid)
  object$var_mtestresid <- var(object$mtestresid)
  
  #LRT test stat and/or p-value
  
  #R^2 Neglekerke
  # cox$loglik[1] loglik with initial values of coef
  # cox$loglik[2] loglik with final values of coef
  # cox$n number of observations
  object$logtraintest <- -2 * (object$nulltrain.fit$loglik - object$train.fit$loglik[2])
  object$logtesttest <- -2 * (object$nulltest.fit$loglik - object$test.fit$loglik[2])
  
  object$rval <- NULL
  # Cox and Snell
  object$rval$train$rsq.cs <- c(rsq.cs = 1 - exp(-object$logtraintest/object$train.fit$n), maxrsq.cs = 1 -exp(2 * object$train.fit$loglik[1]/object$train.fit$n))
  object$rval$test$rsq.cs <- c(rsq.cs = 1 - exp(-object$logtesttest/object$test.fit$n), maxrsq.cs = 1 -exp(2 * object$nulltest.fit$loglik/object$test.fit$n))
  # Nagelkerke
  object$rval$train$rsq.nagel <- c(rsq.nagel = unname(object$rval$train$rsq.cs[1]/object$rval$train$rsq.cs[2]), maxrsq.nagel = 1)
  object$rval$test$rsq.nagel <- c(rsq.nagel = unname(object$rval$test$rsq.cs[1]/object$rval$test$rsq.cs[2]), maxrsq.nagel = 1)
  
  
  #  library(survAUC)
  # iAUC
  object$AUC_CD <- survAUC::AUC.cd(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.auc)  #CoxModel
  object$AUC_hc <- survAUC::AUC.hc(object$Surv.rsp, object$Surv.rsp.new, object$lpnew, object$times.auc)      #No model
  object$AUC_sh <- survAUC::AUC.sh(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.auc)  #CoxModel
  object$AUC_Uno <- survAUC::AUC.uno(object$Surv.rsp, object$Surv.rsp.new, object$lpnew, object$times.auc)    #No model
  
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
    layout(matrix(1:4,nrow=2,byrow=TRUE))
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
  
  
  object$HarrelC <- survAUC::BeggC(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew)  #CoxModel C-statistic by Begg et al.
  object$GonenHellerCI <- survAUC::GHCI(object$lpnew)  #CoxModel Gonen and Heller?s Concordance Index for Cox models
  
  object$rval$train$rsq.OXS <- survAUC::OXS(object$Surv.rsp, object$lp, object$lp0)
  object$rval$train$rsq.Nagelk <- survAUC::Nagelk(object$Surv.rsp, object$lp, object$lp0)
  object$rval$train$rsq.XO <- survAUC::XO(object$Surv.rsp, object$lp, object$lp0)
  
  object$rval$test$rsq.OXS <- survAUC::OXS(object$Surv.rsp.new, object$lpnew, object$lp0new)
  object$rval$test$rsq.Nagelk <- survAUC::Nagelk(object$Surv.rsp.new, object$lpnew, object$lp0new)
  object$rval$test$rsq.XO <- survAUC::XO(object$Surv.rsp.new, object$lpnew, object$lp0new)
  
  #Surv.rsp A survival::Surv(.,.) object containing to the outcome of the test data.
  #lp The vector of predictors.
  #lp0 The vector of predictors obtained from the covariate-free null model.
  
  #prederr #ierror for iBrier
  object$prederr <- NULL
  object$times.prederr <- times.prederr
  object$prederr$brier.unw <- survAUC::predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.prederr, type = "brier", int.type = "unweighted")
  object$prederr$robust.unw <- survAUC::predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.prederr, type = "robust", int.type = "unweighted")
  object$prederr$brier.w <- survAUC::predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.prederr, type = "brier", int.type = "weighted")
  object$prederr$robust.w <- survAUC::predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.prederr, type = "robust", int.type = "weighted")
  object$prederr$brier0.unw <- survAUC::predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp0, object$lp0new, object$times.prederr, type = "brier", int.type = "unweighted")
  object$prederr$robust0.unw <- survAUC::predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp0, object$lp0new, object$times.prederr, type = "robust", int.type = "unweighted")
  object$prederr$brier0.w <- survAUC::predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp0, object$lp0new, object$times.prederr, type = "brier", int.type = "weighted")
  object$prederr$robust0.w <- survAUC::predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp0, object$lp0new, object$times.prederr, type = "robust", int.type = "weighted")
  
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
  
  
  #integrated R2 from BS        max(T_i)=temps (censur? ou non) le plus grand=max(time). Fonction continue lin?aire par morceaux.
  # R2_{BS}(t)=1-BS(t)/BS_0(t)
  # iR2BS=1/max(T_i)int_0^{max(T_i)}R2_{BS}(t)dt.
  object$rval$test$R2.bs.unw <- 1-object$prederr$brier.unw$error/object$prederr$brier0.unw$error
  object$rval$test$R2.bs.unw[object$prederr$brier.unw$error==0] <- 0
  object$rval$test$R2.bs.unw[object$rval$test$R2.bs.unw<=0] <- 0
  object$rval$test$R2.rs.unw <- 1-object$prederr$robust.unw$error/object$prederr$robust0.unw$error
  object$rval$test$R2.rs.unw[object$prederr$robust.unw$error==0] <- 0
  object$rval$test$R2.rs.unw[object$rval$test$R2.rs.unw<=0] <- 0
  object$rval$test$R2.bs.w <- 1-object$prederr$brier.w$error/object$prederr$brier0.w$error
  object$rval$test$R2.bs.w[object$prederr$brier.w$error==0] <- 0
  object$rval$test$R2.bs.w[object$rval$test$R2.bs.w<=0] <- 0
  object$rval$test$R2.rs.w <- 1-object$prederr$robust.w$error/object$prederr$robust0.w$error
  object$rval$test$R2.rs.w[object$prederr$robust.w$error==0] <- 0
  object$rval$test$R2.rs.w[object$rval$test$R2.rs.w<=0] <- 0
  
  object$rval$test$iRbs.unw <- sum(object$rval$test$R2.bs.unw[-1]*diff(object$prederr$brier.unw$time))/max(object$prederr$brier.unw$time)
  object$rval$test$iRrs.unw <- sum(object$rval$test$R2.rs.unw[-1]*diff(object$prederr$robust.unw$time))/max(object$prederr$robust.unw$time)
  object$rval$test$iRbs.w <- sum(object$rval$test$R2.bs.w[-1]*diff(object$prederr$brier.w$time))/max(object$prederr$brier.w$time)
  object$rval$test$iRrs.w <- sum(object$rval$test$R2.rs.w[-1]*diff(object$prederr$robust.w$time))/max(object$prederr$robust.w$time)
  
  #UnoC C-statistic by Uno et al.
  object$Cstat <- survAUC::UnoC(object$Surv.rsp, object$Surv.rsp.new, object$lpnew) #no model
  
  #schemper Distance-based estimator of survival predictive accuracy proposed by Schemper and Henderson
  #Schemper and Henderson's estimator of the absolute deviation between survival functions
  #  library(rms)
  object$Schemper <- survAUC::schemper(object$train.fit.cph, object$TR, object$TE)
  return(invisible(object))
}

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


getIndicCViAUCSurvROCTest = function(lp,lpnew,Surv.rsp,Surv.rsp.new,times.auc=seq(10,1000,10),times.prederr=1:500,train.fit,plot.it=FALSE,tmax.train=max(Surv.rsp[,"time"][ object$Surv.rsp[,"status"] == 1 ]),tmax.test=max(Surv.rsp.new[,"time"][ object$Surv.rsp.new[,"status"] == 1 ])){
  try(attachNamespace("survival"),silent=TRUE)
  #on.exit(try(unloadNamespace("survival"),silent=TRUE))
  try(attachNamespace("survcomp"),silent=TRUE)
  on.exit(try(unloadNamespace("survcomp"),silent=TRUE))
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
  
  
  object$utimes.test <- unique( object$Surv.rsp.new[,"time"][ object$Surv.rsp.new[,"status"] == 1 ] )
  object$utimes.test <- object$utimes.test[ order(object$utimes.test) ]
  
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
  
  return(object)
}

