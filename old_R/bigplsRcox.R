#' Partial least squares Regression generalized linear models
#' 
#' This function implements an extension of Partial least squares Regression to
#' Cox Models.
#' 
#' 
#' A typical predictor has the form response ~ terms where response is the
#' (numeric) response vector and terms is a series of terms which specifies a
#' linear predictor for response. A terms specification of the form first +
#' second indicates all the terms in first together with all the terms in
#' second with any duplicates removed.
#' 
#' A specification of the form first:second indicates the the set of terms
#' obtained by taking the interactions of all terms in first with all terms in
#' second. The specification first*second indicates the cross of first and
#' second. This is the same as first + second + first:second.
#' 
#' The terms in the formula will be re-ordered so that main effects come first,
#' followed by the interactions, all second-order, all third-order and so on:
#' to avoid this pass a terms object as the formula.
#' 
#' Non-NULL weights can be used to indicate that different observations have
#' different dispersions (with the values in weights being inversely
#' proportional to the dispersions); or equivalently, when the elements of
#' weights are positive integers w_i, that each response y_i is the mean of w_i
#' unit-weight observations.
#' 
#' @aliases bigplsRcox bigplsRcoxmodel.default bigplsRcoxmodel.formula
#' @param Xplan a formula or a matrix with the eXplanatory variables (training)
#' dataset
#' @param time for right censored data, this is the follow up time. For
#' interval data, the first argument is the starting time for the interval.
#' @param time2 The status indicator, normally 0=alive, 1=dead. Other choices
#' are \code{TRUE/FALSE} (\code{TRUE} = death) or 1/2 (2=death). For interval
#' censored data, the status indicator is 0=right censored, 1=event at
#' \code{time}, 2=left censored, 3=interval censored. Although unusual, the
#' event indicator can be omitted, in which case all subjects are assumed to
#' have an event.
#' @param event ending time of the interval for interval censored or counting
#' process data only. Intervals are assumed to be open on the left and closed
#' on the right, \code{(start, end]}. For counting process data, event
#' indicates whether an event occurred at the end of the interval.
#' @param type character string specifying the type of censoring. Possible
#' values are \code{"right"}, \code{"left"}, \code{"counting"},
#' \code{"interval"}, or \code{"interval2"}. The default is \code{"right"} or
#' \code{"counting"} depending on whether the \code{time2} argument is absent
#' or present, respectively.
#' @param origin for counting process data, the hazard function origin. This
#' option was intended to be used in conjunction with a model containing time
#' dependent strata in order to align the subjects properly when they cross
#' over from one strata to another, but it has rarely proven useful.
#' @param typeres character string indicating the type of residual desired.
#' Possible values are \code{"martingale"}, \code{"deviance"}, \code{"score"},
#' \code{"schoenfeld"}, \code{"dfbeta"}, \code{"dfbetas"}, and
#' \code{"scaledsch"}. Only enough of the string to determine a unique match is
#' required.
#' @param collapse vector indicating which rows to collapse (sum) over. In
#' time-dependent models more than one row data can pertain to a single
#' individual. If there were 4 individuals represented by 3, 1, 2 and 4 rows of
#' data respectively, then \code{collapse=c(1,1,1,2,3,3,4,4,4,4)} could be used
#' to obtain per subject rather than per observation residuals.
#' @param weighted if \code{TRUE} and the model was fit with case weights, then
#' the weighted residuals are returned.
#' @param scaleX Should the \code{Xplan} columns be standardized ?
#' @param scaleY Should the \code{time} values be standardized ?
#' @param nt number of components to be extracted
#' @param limQ2set limit value for the Q2
#' @param dataPredictY predictor(s) (testing) dataset
#' @param pvals.expli should individual p-values be reported to tune model
#' selection ?
#' @param alpha.pvals.expli level of significance for predictors when
#' pvals.expli=TRUE
#' @param tol_Xi minimal value for Norm2(Xi) and \eqn{\mathrm{det}(pp' \times
#' pp)}{det(pp'*pp)} if there is any missing value in the \code{dataX}. It
#' defaults to \eqn{10^{-12}}{10^{-12}}
#' @param weights an optional vector of 'prior weights' to be used in the
#' fitting process. Should be \code{NULL} or a numeric vector.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param allres FALSE to return only the Cox model and TRUE for additionnal
#' results. See details. Defaults to FALSE.
#' @param dataXplan an optional data frame, list or environment (or object
#' coercible by \code{\link{as.data.frame}} to a data frame) containing the
#' variables in the model. If not found in \code{dataXplan}, the variables are
#' taken from \code{environment(Xplan)}, typically the environment from which
#' \code{coxDKplsDR} is called.
#' @param model_frame If \code{TRUE}, the model frame is returned.
#' @param method the method to be used in fitting the model. The default method
#' \code{"glm.fit"} uses iteratively reweighted least squares (IWLS).
#' User-supplied fitting functions can be supplied either as a function or a
#' character string naming a function, with a function which takes the same
#' arguments as \code{glm.fit}.
#' @param control a list of parameters for controlling the fitting process. For
#' \code{glm.fit} this is passed to \code{\link{glm.control}}.
#' @param sparse should the coefficients of non-significant predictors
#' (<\code{alpha.pvals.expli}) be set to 0
#' @param sparseStop should component extraction stop when no significant
#' predictors (<\code{alpha.pvals.expli}) are found
#' @param verbose Should some details be displayed ?
#' @param model_matrix If \code{TRUE}, the model matrix is returned.
#' @param contrasts.arg a list, whose entries are values (numeric matrices, 
#' functions or character strings naming functions) to be used as replacement 
#' values for the contrasts replacement function and whose names are the names 
#' of columns of data containing factors.
#' @param \dots arguments to pass to \code{plsRmodel.default} or to
#' \code{plsRmodel.formula}
#' @return Depends on the model that was used to fit the model.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@lecnam.net}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link[plsRglm]{plsR}} and \code{\link[plsRglm]{plsRglm}}
#' @references bigplsRcox, Cox-Models in a high dimensional setting in R, Frederic
#' Bertrand, Philippe Bastien, Nicolas Meyer and Myriam Maumy-Bertrand (2014).
#' Proceedings of User2014!, Los Angeles, page 152.\cr
#' 
#' Deviance residuals-based sparse PLS and sparse kernel PLS regression for
#' censored data, Philippe Bastien, Frederic Bertrand, Nicolas Meyer and Myriam
#' Maumy-Bertrand (2015), Bioinformatics, 31(3):397-404,
#' doi:10.1093/bioinformatics/btu660.
#' @keywords models regression
#' @examples
#' 
#' data(micro.censure)
#' data(Xmicro.censure_compl_imp)
#' 
#' X_train_micro <- apply((as.matrix(Xmicro.censure_compl_imp)),FUN="as.numeric",MARGIN=2)[1:80,]
#' X_train_micro_df <- data.frame(X_train_micro)
#' Y_train_micro <- micro.censure$survyear[1:80]
#' C_train_micro <- micro.censure$DC[1:80]
#' 
#' bigplsRcox(X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5)
#' bigplsRcox(~X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5)
#' 
#' bigplsRcox(Xplan=X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5,sparse=TRUE,
#' alpha.pvals.expli=.15)
#' bigplsRcox(Xplan=~X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5,sparse=TRUE,
#' alpha.pvals.expli=.15)
#' 
#' @export bigplsRcox
bigplsRcox <- function (formula, data, ...) UseMethod("bigplsRcoxmodel")

#' @rdname bigplsRcox
#' @aliases bigplsRcox
#' @export bigplsRcoxmodel
bigplsRcoxmodel <- bigplsRcox

#' @rdname bigplsRcox
#' @export

bigplsRcoxmodel.default <- function(formula = survival::Surv(time = time, status = status) ~ ., data, deepcopy_data=TRUE, scale.X=TRUE, scale.Y=TRUE, nt=2, type="double", allres=TRUE, verbose=TRUE, backingfile=NULL, backingpath=NULL, descriptorfile=NULL,...) {
  
  ##################################################
  #                                                #
  #    Initialization and formatting the inputs    #
  #                                                #
  ##################################################
  
  if(verbose){cat("____************************************************____\n")}

  modele <- "bigpls-cox"
  
  tempbackingpath <- backingpath
  if(is.null(descriptorfile)){
    if(verbose){cat("Using big matrix without backing and descriptor files\n")}
    big.data0.orig <- bigmemory::read.big.matrix(filename = data, sep = ",", 
                                                 skip = 0, header = TRUE, 
                                                 type=type)
    if(verbose){cat("Done\n")}
  } else {
    if(!file.exists(paste(backingpath,descriptorfile,sep=""))) {
      if(verbose){
        cat("Creating backing and descriptor files","\n")
        cat(paste(backingpath,descriptorfile,sep=""),"\n")
        cat(!file.exists(paste(backingpath,descriptorfile,sep="")),"\n")
        }
      big.data0.orig <- bigmemory::read.big.matrix(filename = data, sep = ",", 
                                                   skip = 0, header = TRUE, 
                                                   backingfile=backingfile,
                                                   backingpath=backingpath,
                                                   descriptorfile=descriptorfile,
                                                   type=type)
      if(verbose){cat("Done")}
    } else {
      if(verbose){
        cat("Attaching existing descriptor file","\n")
        cat(paste(backingpath,descriptorfile,sep=""),"\n")
        cat(!file.exists(paste(backingpath,descriptorfile,sep="")),"\n")
      }
      big.data0.orig <- bigmemory::attach.big.matrix(obj = descriptorfile,
                                                     backingpath=backingpath)
      if(verbose){cat("Done")}
    }
  }
  
  ind.col=1
  name.col.all <-(colnames(big.data0.orig)[-c(1,2)])
  name.col0 = name.col.all[ind.col]
  
  if(verbose){
    cat(name.col.all,"\n")
    cat(dim(big.data0.orig),"\n")
  }
  if(is.null(descriptorfile)){
    if(deepcopy_data){
      big.data0 <- bigmemory::deepcopy(big.data0.orig, type=type)
    } else {
      big.data0 <- big.data0.orig
    }
    }
  else {
    if(deepcopy_data){
      tempname <- sub(".", "", tempfile(pattern = "plsRcox", tmpdir = "", fileext = ""))
      if(verbose){
        cat(tempbackingpath,"\n")
        cat(paste(tempname,".bin",sep=""),"\n")
        cat(paste(tempname,".desc",sep=""),"\n")
      }
      big.data0 <- bigmemory::deepcopy(big.data0.orig, 
                                       type=type,
                                       backingfile = paste(tempname,"_temp.bin",sep=""),
                                       backingpath = tempbackingpath,
                                       descriptorfile = paste(tempname,"_temp.desc",sep=""))
    } else {
      big.data0 <- big.data0.orig
    }
  }
  
  if(verbose){
    cat(sum(abs(big.data0.orig[1:6,1:6]-big.data0[1:6,1:6])))
  }
  
  if(is.null(tempbackingpath)){
    residXX <- bigmemory::sub.big.matrix(
      big.data0,
      firstRow = 1,
      lastRow = NULL,
      firstCol = 3,
      lastCol = NULL
    )
  } else {
    residXX <- bigmemory::sub.big.matrix(
      big.data0,
      firstRow = 1,
      lastRow = NULL,
      firstCol = 3,
      lastCol = NULL,
      backingpath = tempbackingpath)
  }
  # Normalize data?
  
  cat("Compute scaler\n")  
  if(scale.X || scale.Y){
    resultsBigscale <- bigPLS::bigscale(formula=survival::Surv(time, status)~., data=big.data0, 
                                bigmemory.flag=TRUE, parallel.flag=TRUE, 
                                num.cores=2, norm.method = "standardize")
  } else {
    resultsBigscale <- bigPLS::bigscale(formula=survival::Surv(time, status)~., data=big.data0, 
                                bigmemory.flag=TRUE, parallel.flag=TRUE, 
                                num.cores=2, norm.method = "none")
  }
  

  lapply
  biganalytics::apply(big.data0, 2, mean)
  big.data0 <- 
  
    (macc[ind.col[j]][ind.row[i]] - mean[j]) / sd[j]
  
  x <- sweep(x, 2L, scale, "/", check.margin = FALSE)
  nt <- min(2,resultsBigscale$nc)

  res <- list(nr=NULL,nc=NULL,nt=nt,ww=NULL,wwnorm=NULL,wwetoile=NULL,tt=NULL,pp=NULL,CoeffC=NULL,uscores=NULL,YChapeau=NULL,residYChapeau=NULL) 
  res$name.col.all <- (resultsBigscale$col.names)[resultsBigscale$features.indices]
  res$nr <- resultsBigscale$nr
  res$nc <- resultsBigscale$nc-2
  res$temppred <- NULL
  
  
  ################################################
  ################################################
  ##                                            ##
  ##  Beginning of the loop for the components  ##
  ##                                            ##
  ################################################
  ################################################
  
  res$computed_nt <- 0
  
  for (kk in 1:nt) {
    cat(nt,"\n")  
    
    if(kk==1L){
      cat("Use scaler 1st loop\n")  
      resultsBigscale<-resultsBigscale  
    } else {
      cat("Do not use scaler additionnal loops\n")  
      resultsBigscale<-resultsBigscale  
    }
    
    tempww <- rep(0,res$nc)
    res$computed_nt <- kk
    
    ##############################################
    #                                            #
    #     Weight computation for each model      #
    #                                            #
    ##############################################
    
    ##############################################
    ######              PLS-COX             ######
    ##############################################
    if (modele %in% c("bigpls-cox")) {
      tempww <- simplify2array(lapply(name.col.all,bigPLS::partialbigSurvSGDv0, 
                                      datapath=big.data0, 
                                      resBigscale=resultsBigscale, 
                                      bigmemory.flag=TRUE,
                                      parallel.flag = TRUE))
    }
    stopifnot(length(tempww)==res$nc)

    ##############################################
    #                                            #
    # Computation of the components (model free) #
    #                                            #
    ##############################################

    tempwwnorm <- tempww/sqrt(drop(crossprod(tempww)))
#    return(tempwwnorm)
    
    BM2FBM <- function(bm) {
      bigstatsr::FBM(nrow = nrow(bm), ncol = ncol(bm), type = typeof(bm),
          backingfile = file.path(bigmemory::dir.name(bm), bigstatsr::sub_bk(bigmemory::file.name(bm))),
          create_bk = FALSE)
    }
    
    
    sub.big.data0 %*% tempwwnorm
    bigstatsr::big_prodVec(big.data0, tempwwnorm)
    
    test2 <- bigstatsr::big_prodMat(m2, a, ind.col = name.col.all)
    
    temptt <- XXwotNA%*%tempwwnorm/(XXNA%*%(tempwwnorm^2))
    
    temppp <- rep(0,res$nc)
    for (jj in 1:(res$nc)) {
      temppp[jj] <- crossprod(temptt,XXwotNA[,jj])/drop(crossprod(XXNA[,jj],temptt^2))
    }
    res$residXX <- XXwotNA-temptt%*%temppp
    
    if (na.miss.X & !na.miss.Y) {
      for (ii in 1:res$nr) {
        if(rcond(t(cbind(res$pp,temppp)[XXNA[ii,],,drop=FALSE])%*%cbind(res$pp,temppp)[XXNA[ii,],,drop=FALSE])<tol_Xi) {
          break_nt <- TRUE; res$computed_nt <- kk-1
          if(verbose){cat(paste("Warning : reciprocal condition number of t(cbind(res$pp,temppp)[XXNA[",ii,",],,drop=FALSE])%*%cbind(res$pp,temppp)[XXNA[",ii,",],,drop=FALSE] < 10^{-12}\n",sep=""))}
          if(verbose){cat(paste("Warning only ",res$computed_nt," components could thus be extracted\n",sep=""))}
          break
        }
      }
      rm(ii)
      if(break_nt==TRUE) {break}
    }
    
    if(!PredYisdataX){
      if (na.miss.PredictY & !na.miss.Y) {
        for (ii in 1:nrow(PredictYwotNA)) {
          if(rcond(t(cbind(res$pp,temppp)[PredictYNA[ii,],,drop=FALSE])%*%cbind(res$pp,temppp)[PredictYNA[ii,],,drop=FALSE])<tol_Xi) {
            break_nt <- TRUE; res$computed_nt <- kk-1
            if(verbose){cat(paste("Warning : reciprocal condition number of t(cbind(res$pp,temppp)[PredictYNA[",ii,",,drop=FALSE],])%*%cbind(res$pp,temppp)[PredictYNA[",ii,",,drop=FALSE],] < 10^{-12}\n",sep=""))}
            if(verbose){cat(paste("Warning only ",res$computed_nt," components could thus be extracted\n",sep=""))}
            break
          }
        }
        rm(ii)
        if(break_nt==TRUE) {break}
      }
    }
    
    
    res$ww <- cbind(res$ww,tempww)
    res$wwnorm <- cbind(res$wwnorm,tempwwnorm)
    res$tt <- cbind(res$tt,temptt)       
    res$pp <- cbind(res$pp,temppp)   
    
    
    ##############################################
    #                                            #
    #      Computation of the coefficients       #
    #      of the model with kk components       #
    #                                            #
    ##############################################
    
    ##############################################
    ######              PLS-COX            ######
    ##############################################
    if (modele %in% c("bigpls-cox")) {
      if (kk==1) {
        mf2b <- match.call(expand.dots = TRUE)
        m2b <- match(c(head(names(as.list(args(survival::coxph))),-2),head(names(as.list(args(survival::coxph.control))),-1)), names(mf2b), 0L)
        mf2b <- mf2b[c(1L, m2b)]
        mf2b$formula <- as.formula(YCsurv~1)
        mf2b[[1L]] <- as.name("coxph")
        tempconstcox <- eval(mf2b, parent.frame())
        res$AIC <- suppressWarnings(extractAIC(tempconstcox)[2])
        res$BIC <- suppressWarnings(extractAIC(tempconstcox, k = log(res$nr))[2])
        res$Coeffsmodel_vals <- rbind(rep(NA,5),matrix(rep(NA,5*nt),ncol=5))
        rm(tempconstcox)
        mf2b <- match.call(expand.dots = TRUE)
        m2b <- match(c(head(names(as.list(args(survival::coxph))),-2),head(names(as.list(args(survival::coxph.control))),-1)), names(mf2b), 0L)
        mf2b <- mf2b[c(1L, m2b)]
        mf2b$formula <- as.formula(YwotNA~.)
        #tt<-data.frame(res$tt); colnames(tt) <- paste("Comp_",1:length(tt),sep="")
        tttrain<-data.frame(YwotNA=YwotNA,tt=res$tt)
        #mf2b$data <- tt
        mf2b$data <- tttrain  
        mf2b$model <- TRUE
        mf2b[[1L]] <- as.name("coxph")
        tempregcox <- eval(mf2b, parent.frame())
        rm(tttrain)
        res$AIC <- cbind(res$AIC,extractAIC(tempregcox)[2])
        res$BIC <- cbind(res$BIC,extractAIC(tempregcox, k = log(res$nr))[2])
        res$Coeffsmodel_vals <- cbind(rbind(summary(tempregcox)$coefficients,matrix(rep(NA,5*(nt-kk)),ncol=5)))
        tempCoeffC <- as.vector(coef(tempregcox))
        res$CoeffCFull <- matrix(c(tempCoeffC,rep(NA,nt-kk)),ncol=1)
      } else {
        if (!(na.miss.X | na.miss.Y)) {
          mf2b <- match.call(expand.dots = TRUE)
          m2b <- match(c(head(names(as.list(args(survival::coxph))),-2),head(names(as.list(args(survival::coxph.control))),-1)), names(mf2b), 0L)
          mf2b <- mf2b[c(1L, m2b)]
          mf2b$formula <- as.formula(YwotNA~.)
          #tt<-data.frame(res$tt); colnames(tt) <- paste("Comp_",1:length(tt),sep="")
          tttrain<-data.frame(YwotNA=YwotNA,tt=res$tt)
          #mf2b$data <- tt
          mf2b$data <- tttrain  
          mf2b$model <- TRUE
          mf2b[[1L]] <- as.name("coxph")
          tempregcox <- eval(mf2b, parent.frame())
          rm(tttrain)
          res$AIC <- cbind(res$AIC,extractAIC(tempregcox)[2])
          res$BIC <- cbind(res$BIC,extractAIC(tempregcox, k = log(res$nr))[2])
          res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregcox)$coefficients,matrix(rep(NA,5*(nt-kk)),ncol=5)))
          tempCoeffC <- as.vector(coef(tempregcox))  
          res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffC,rep(NA,nt-kk)))
        }
        else
        {
          mf2b <- match.call(expand.dots = TRUE)
          m2b <- match(c(head(names(as.list(args(survival::coxph))),-2),head(names(as.list(args(survival::coxph.control))),-1)), names(mf2b), 0L)
          mf2b <- mf2b[c(1L, m2b)]
          mf2b$formula <- as.formula(YwotNA~.)
          #tt<-data.frame(res$tt); colnames(tt) <- paste("Comp_",1:length(tt),sep="")
          tttrain<-data.frame(YwotNA=YwotNA,tt=res$tt)
          #mf2b$data <- tt
          mf2b$data <- tttrain  
          mf2b$model <- TRUE
          mf2b[[1L]] <- as.name("coxph")
          tempregcox <- eval(mf2b, parent.frame())
          rm(tttrain)
          res$AIC <- cbind(res$AIC,extractAIC(tempregcox)[2])
          res$BIC <- cbind(res$BIC,extractAIC(tempregcox, k = log(res$nr))[2])
          res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregcox)$coefficients,matrix(rep(NA,5*(nt-kk)),ncol=5)))
          tempCoeffC <- as.vector(coef(tempregcox))  
          res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffC,rep(NA,nt-kk)))
        }
      }
      
      res$wwetoile <- (res$wwnorm)%*%solve(t(res$pp)%*%res$wwnorm)
      res$CoeffC <- tempCoeffC
      res$Std.Coeffs <- res$wwetoile%*%res$CoeffC
      rownames(res$Std.Coeffs) <- colnames(ExpliX)
    }
    
    
    ##############################################
    #                                            #
    #       Prediction of the components         #
    #     as if missing values (model free)      #
    #       For cross-validating the GLM         #
    #                                            #
    ##############################################
    
    
    
    
    
    if (!(na.miss.X | na.miss.Y)) {
      
      ##############################################
      #                                            #
      #             Cross validation               #
      #           without missing value            #
      #                                            #
      ##############################################
      
      ##############################################
      ######              PLS-COX            ######
      ##############################################
      if (modele %in% c("bigpls-cox")) {
        res$residYChapeau <- tempregcox$linear.predictors
        
        tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
        res$Coeffs <- tempCoeffs
        
        res$YChapeau <- as.matrix(predict(tempregcox, type='expected'))            
        res$Yresidus <- residuals(tempregcox,type="martingale")
      }
      
      ##############################################
    }
    
    else {
      if (na.miss.X & !na.miss.Y) {
        
        
        ##############################################
        #                                            #
        #             Cross validation               #
        #           with missing value(s)            #
        #                                            #
        ##############################################
        
        
        if (kk==1) {
          if(verbose){cat("____There are some NAs in X but not in Y____\n")}
        }
        
        ##############################################
        ######              PLS-COX            ######
        ##############################################
        if (modele %in% c("bigpls-cox")) {
          res$residYChapeau <- tempregcox$linear.predictors
          
          tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
          res$Coeffs <- tempCoeffs
          
          res$YChapeau <- as.matrix(predict(tempregcox, type='expected'))            
          res$Yresidus <- residuals(tempregcox,type="martingale")
        }
        
        ##############################################
      }
      
      else {
        if (kk==1) {
          if(verbose){cat("____There are some NAs both in X and Y____\n")}
        }
      }
    }
    
    
    
    
    ##############################################
    #                                            #
    #      Update and end of loop cleaning       #
    #        (Especially useful for PLS)         #
    #                                            #
    ##############################################
    
    
    
    ##############################################
    ######              PLS-COX            ######
    ##############################################
    if (modele %in% c("bigpls-cox")) {
      res$residY <- res$residY 
      res$residusY <- cbind(res$residusY,res$residY)
      
      rm(tempww)
      rm(tempwwnorm)
      rm(temptt)
      rm(temppp)
      rm(tempCoeffC)
    }
    
    if(verbose){cat("____Component____",kk,"____\n")}
  }
 
  
  ##############################################
  ##############################################
  ##                                          ##
  ##    End of the loop on the components     ##
  ##                                          ##
  ##############################################
  ##############################################
  
  if(res$computed_nt==0){
    if(verbose){cat("No component could be extracted please check the data for NA only lines or columns\n")}; stop()
  }
  
  
  if (pvals.expli&!(modele=="pls")) {
    res$Coeffsmodel_vals<-res$Coeffsmodel_vals[1:(dim(res$Coeffsmodel_vals)[1]-(nt-res$computed_nt)),]
  }
  
  
  ##############################################
  #                                            #
  #           Predicting components            #
  #                                            #
  ##############################################
  
  if (!(na.miss.PredictY | na.miss.Y)) {
    if(verbose){cat("____Predicting X without NA neither in X nor in Y____\n")}
    res$ttPredictY <- PredictYwotNA%*%res$wwetoile 
    colnames(res$ttPredictY) <- NULL
  }
  else {
    if (na.miss.PredictY & !na.miss.Y) {
      if(verbose){cat("____Predicting X with NA in X and not in Y____\n")}
      for (ii in 1:nrow(PredictYwotNA)) {  
        res$ttPredictY <- rbind(res$ttPredictY,t(solve(t(res$pp[PredictYNA[ii,],,drop=FALSE])%*%res$pp[PredictYNA[ii,],,drop=FALSE])%*%t(res$pp[PredictYNA[ii,],,drop=FALSE])%*%(PredictYwotNA[ii,])[PredictYNA[ii,]]))
      }
      colnames(res$ttPredictY) <- NULL
    }
    else {
      if(verbose){cat("____There are some NAs both in X and Y____\n")}
    }
  }

  
  ##############################################
  #                                            #
  #          Computing RSS, PRESS,             #
  #           Chi2, Q2 and Q2cum               #
  #                                            #
  ##############################################
  
  
  ##############################################
  ######              PLS-COX            ######
  ##############################################
  if (modele %in% c("bigpls-cox")) {
    res$InfCrit <- t(rbind(res$AIC, res$BIC))
    dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt,sep=""), c("AIC", "BIC"))
  }
  
  ##########################################
  #                                        #
  #          Predicting responses          #
  #                                        #
  ##########################################
  
  
  
  ##############################################
  ######              PLS-COX            ######
  ##############################################
  if (modele %in% c("bigpls-cox")) {
    res$YChapeau <- as.matrix(predict(tempregcox, type='expected'))            
    rownames(res$YChapeau) <- rownames(ExpliX)
    
    ttpred <- data.frame(tt=res$ttPredictY)
    res$Std.ValsPredictY <- predict(tempregcox,newdata=ttpred, type = "lp")
    res$ValsPredictY <- predict(tempregcox,newdata=ttpred,type = "risk")
    
    res$Std.XChapeau <- res$tt%*%t(res$pp)
    rownames(res$Std.XChapeau) <- rownames(ExpliX)
    names(res$CoeffC) <- paste("Coeff_Comp_Reg",1:res$computed_nt)
    rownames(res$Coeffs) <- colnames(ExpliX)
    res$FinalModel <- tempregcox
  }
  
  rownames(res$pp) <- colnames(ExpliX)
  colnames(res$pp) <- paste("Comp_",1:res$computed_nt,sep="")
  rownames(res$ww) <- colnames(ExpliX)
  colnames(res$ww) <- paste("Comp_",1:res$computed_nt,sep="")
  rownames(res$wwnorm) <- colnames(ExpliX)
  colnames(res$wwnorm) <- paste("Comp_",1:res$computed_nt,sep="")
  rownames(res$wwetoile) <- colnames(ExpliX)
  colnames(res$wwetoile) <- paste("Coord_Comp_",1:res$computed_nt,sep="")
  rownames(res$tt) <- rownames(ExpliX)
  colnames(res$tt) <- paste("Comp_",1:res$computed_nt,sep="")
  res$XXwotNA <- XXwotNA
  if(verbose){cat("****________________________________________________****\n")}
  if(verbose){cat("\n")}
  class(res) <- "bigplsRcoxmodel"
  return(res)
}

