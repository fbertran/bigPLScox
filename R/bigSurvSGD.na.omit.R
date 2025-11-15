#' Fit Survival Models with Stochastic Gradient Descent
#'
#' Performs stochastic gradient descent optimisation for large-scale survival
#' models after removing observations with missing values.
#'
#' @param formula Model formula describing the survival outcome and the set of
#' predictors to include in the optimisation.
#' @param data Input data set or connection to a big-memory backed design
#' matrix that contains the variables referenced in \code{formula}.
#' @param norm.method Normalization strategy applied to the feature matrix
#' before optimisation, for example centring or standardising columns.
#' @param features.mean Optional pre-computed column means used when
#' normalising the features so that repeated fits can reuse shared
#' statistics.
#' @param features.sd Optional pre-computed column standard deviations used in
#' concert with \code{features.mean} for scaling the predictors.
#' @param opt.method Gradient based optimisation routine to employ, such as
#' vanilla SGD or adaptive methods like Adam.
#' @param beta.init Vector of starting values for the regression coefficients
#' supplied when warm-starting the optimisation.
#' @param beta.type Indicator controlling how \code{beta.init} is interpreted,
#' for example whether the coefficients correspond to the original or
#' normalised scale.
#' @param lr.const Base learning-rate constant used by the stochastic
#' gradient descent routine.
#' @param lr.tau Learning-rate decay horizon or damping factor that moderates
#' the step size schedule.
#' @param strata.size Number of observations drawn per stratum when building
#' mini-batches for the optimisation loop.
#' @param batch.size Total number of observations assembled into each
#' stochastic gradient batch.
#' @param num.epoch Number of passes over the training data used during the
#' optimisation.
#' @param b1 First exponential moving-average rate used by adaptive methods
#' such as Adam to smooth gradients.
#' @param b2 Second exponential moving-average rate used by adaptive methods
#' to smooth squared gradients.
#' @param eps Numerical stabilisation constant added to denominators when
#' updating the adaptive moments.
#' @param inference.method Inference approach requested after fitting, for
#' example naive asymptotics or bootstrap resampling.
#' @param num.boot Number of bootstrap replicates to draw when
#' \code{inference.method} relies on resampling.
#' @param num.epoch.boot Number of optimisation epochs to run within each
#' bootstrap replicate.
#' @param boot.method Type of bootstrap scheme to apply, such as ordinary or
#' stratified resampling.
#' @param lr.const.boot Learning-rate constant used during bootstrap refits.
#' @param lr.tau.boot Learning-rate decay factor applied during bootstrap
#' refits.
#' @param num.sample.strata Number of strata sampled without replacement during
#' each bootstrap iteration when stratified resampling is selected.
#' @param sig.level Significance level used when constructing confidence
#' intervals or hypothesis tests.
#' @param beta0 Optional vector of coefficients under the null hypothesis when
#' performing hypothesis tests.
#' @param alpha Elastic-net mixing parameter controlling the relative weight of
#' \eqn{\ell_1} and \eqn{\ell_2} regularisation penalties.
#' @param lambda Sequence of regularisation strengths supplied explicitly for
#' penalised estimation.
#' @param nlambda Number of automatically generated \code{lambda} values when a
#' grid is produced internally.
#' @param num.strata.lambda Number of strata used when tuning \code{lambda} via
#' cross-validation or other search procedures.
#' @param lambda.scale Scale on which the \code{lambda} grid is generated, for
#' example logarithmic or linear spacing.
#' @param parallel.flag Logical flag enabling parallel computation of
#' gradients or bootstrap replicates.
#' @param num.cores Number of processing cores to use when parallel execution
#' is enabled.
#' @param bigmemory.flag Logical flag indicating whether intermediate matrices
#' should be stored using \pkg{bigmemory} backed objects.
#' @param num.rows.chunk Row chunk size to use when streaming data from an
#' on-disk matrix representation.
#' @param col.names Optional character vector of column names associated with
#' the feature matrix.
#' @param type Type of survival model to fit, for example Cox proportional
#' hazards or accelerated failure time variants.
#'
#' @return A fitted model object storing the learned coefficients, optimisation
#' metadata, and any requested inference summaries.
#' coef: Log of hazards ratio. If no inference is used, it returns a vector for 
#' estimated coefficients: If inference is used, it returns a matrix including 
#' estimates and confidence intervals of coefficients. In case of penalization, 
#' it resturns a matrix with columns corresponding to lambdas.
#' coef.exp: Exponentiated version of coef (hazards ratio).
#' lambda: Returns lambda(s) used for penalizarion.
#' alpha: Returns alpha used for penalizarion.
#' features.mean: Returns means of features, if given or calculated
#' features.sd: Returns standard deviations of features, if given or calculated.
#' @export
#' @seealso See Also \code{\link[bigSurvSGD]{bigSurvSGD}},
#' \code{\link{bigscale}} for constructing normalised design matrices and 
#' \code{\link{partialbigSurvSGDv0}} for partial fitting pipelines.
#'
#' @examples
#' \donttest{
#' data(micro.censure, package = "bigPLScox")
#' surv_data <- stats::na.omit(micro.censure[, c("survyear", "DC", "sexe", "Agediag")])
#' # Increase num.epoch and num.boot for real use
#' fit <- bigSurvSGD.na.omit(
#'    survival::Surv(survyear, DC) ~ .,
#'    data = surv_data,
#'    norm.method = "standardize",
#'    opt.method = "adam",
#'    batch.size = 16,
#'    num.epoch = 2,
#'  )
#' }
#' 
bigSurvSGD.na.omit <- function (formula = survival::Surv(time = time, status = status) ~ ., data, 
          norm.method = "standardize", features.mean = NULL, features.sd = NULL, 
          opt.method = "AMSGrad", beta.init = NULL, beta.type = "averaged", 
          lr.const = 0.12, lr.tau = 0.5, strata.size = 20, batch.size = 1, 
          num.epoch = 100, b1 = 0.9, b2 = 0.99, eps = 1e-08, inference.method = "plugin", 
          num.boot = 1000, num.epoch.boot = 100, boot.method = "SGD", 
          lr.const.boot = 0.12, lr.tau.boot = 0.5, num.sample.strata = 1000, 
          sig.level = 0.05, beta0 = 0, alpha = NULL, lambda = NULL, 
          nlambda = 100, num.strata.lambda = 10, lambda.scale = 1, 
          parallel.flag = FALSE, num.cores = NULL, bigmemory.flag = FALSE, 
          num.rows.chunk = 1e+06, col.names = NULL, type="float") 
{
  if (!bigmemory.flag) {
    if (!is.data.frame(data)) {
      big.data <- read.csv(file = data)
    }
    else {
      big.data <- data
    }
  }
  else {
    if(bigmemory::is.big.matrix(data)){big.data <- data} else{
      if (!is.null(col.names)) {
      big.data <- bigmemory::read.big.matrix(filename = data, sep = ",", 
                                  skip = 0, header = TRUE, 
                                  col.names = col.names, type=type)
    }
    else {
      big.data <- bigmemory::read.big.matrix(filename = data, sep = ",", 
                                  skip = 0, header = TRUE, type=type)
    }
    }
  }
  num.rows.big <- nrow(big.data)
  if (ncol(big.data) < 3) {
    stop("data must have 3 or more columns: time, status, and at least one feature")
  }
  if (nrow(big.data) < 2) {
    stop("Sample size is too small (with size less than 2)")
  }
  all.variables <- all.vars(formula)
  surv.indices <- match(all.variables[1:2], colnames(big.data))
  if (length(all.variables) == 3 & all.variables[3] == ".") {
    features.indices <- setdiff(1:NCOL(big.data), surv.indices)
    sub.col.names <- colnames(big.data)[features.indices]
  }
  else {
    features.indices <- match(all.variables[3:length(all.variables)], 
                              colnames(big.data))
    sub.col.names <- all.variables[3:length(all.variables)]
  }
  chengeStrataBatch <- (strata.size > floor(num.rows.big/batch.size)) & 
    (strata.size > 2)
  while ((strata.size > floor(num.rows.big/batch.size)) & (strata.size > 
                                                           2)) {
    if (batch.size > 1) {
      batch.size <- max(floor(batch.size/2), 1)
    }
    else {
      strata.size <- max(floor(num.rows.big/batch.size), 
                         2)
    }
  }
  if (chengeStrataBatch) {
    warning(paste0("Strata size times batch size is greater than number of observations.\n This package resizes them to strata size = ", 
                   strata.size, " and batch size = ", batch.size))
  }
  num.sub.sample <- floor(num.rows.big/num.rows.chunk)
  chunks.length <- c(0, rep(num.rows.chunk, floor(num.rows.big/num.rows.chunk)), 
                     if (num.rows.big%%num.rows.chunk != 0) {
                       num.rows.big%%num.rows.chunk
                     })
  if (is.null(features.mean) & is.null(features.sd)) {
    if (norm.method == "center") {
      n2 <- 0
      features.mean <- 0
      for (i in 1:(length(chunks.length) - 1)) {
        indices.chunk <- (sum(chunks.length[1:i]) + 1):(sum(chunks.length[1:(i + 
                                                                               1)]))
        sub.data <- big.data[indices.chunk, features.indices]
        n1 <- NROW(sub.data)
        if (NCOL(sub.data) > 1) {
          M1 <- colMeans(sub.data, na.rm = TRUE)
        }
        else {
          M1 <- mean(sub.data, na.rm = TRUE)
        }
        features.mean <- (n1 * M1 + n2 * features.mean)/(n1 + 
                                                           n2)
        n2 <- n1 + n2
      }
      features.sd <- rep(1, NCOL(sub.data))
    }
    else if (norm.method == "scale" || norm.method == "standardize") {
      n2 <- 0
      features.mean <- 0
      features.sd <- 0
      for (i in 1:(length(chunks.length) - 1)) {
        indices.chunk <- (sum(chunks.length[1:i]) + 1):(sum(chunks.length[1:(i + 
                                                                               1)]))
        sub.data <- big.data[indices.chunk, features.indices]
        n1 <- NROW(sub.data)
        if (NCOL(sub.data) > 1) {
          M1 <- colMeans(sub.data, na.rm = TRUE)
        }
        else {
          M1 <- mean(sub.data, na.rm = TRUE)
        }
        if (NCOL(sub.data) > 1) {
          S1 <- colMeans(sub.data^2, na.rm = TRUE) - 
            M1^2
        }
        else {
          S1 <- mean(sub.data^2, na.rm = TRUE) - M1^2
        }
        M2 <- features.mean
        S2 <- features.sd
        features.mean <- (n1 * M1 + n2 * M2)/(n1 + n2)
        features.sd <- 1/(n1 + n2) * (n1 * S1 + n2 * 
                                        S2 + (n1 * n2)/(n1 + n2) * (M1 - M2)^2)
        n2 <- n1 + n2
      }
      features.sd <- sqrt(features.sd)
    }
    else {
      features.mean <- rep(0, NCOL(sub.data))
      features.sd <- rep(1, NCOL(sub.data))
    }
  }
  if (sum(features.sd == 0) > 0) {
    stop(paste0("feature(s) ", colnames(big.data)[features.indices][which(features.sd == 
                                                                            0)], " is/are constant without any variability"))
  }
  lambda.max <- function(big.data, strata.size, num.rows.big, 
                         num.rows.chunk, num.strata.lambda, features.indices, 
                         surv.indices) {
    num.rows.chunk <- strata.size * floor(num.rows.big/strata.size)
    num.round <- num.strata.lambda/floor(num.rows.big/strata.size)
    num.round.vec <- c(rep(floor(num.rows.big/strata.size), 
                           floor(num.round)), if ((num.round%%1) > 0) {
                             ceiling((num.round%%1) * floor(num.rows.big/strata.size))
                           })
    indices.stratas <- NULL
    for (i in 1:length(num.round.vec)) {
      indices.stratas <- c(indices.stratas, sample(1:num.rows.big, 
                                                   num.round.vec[i] * strata.size, replace = FALSE))
    }
    gt.sum <- 0
    num.round <- length(indices.stratas)/num.rows.chunk
    num.round.vec <- cumsum(c(0, rep(num.rows.chunk, floor(num.round)), 
                              if ((num.round%%1) > 0) {
                                ceiling((num.round%%1) * num.rows.chunk)
                              }))
    for (i.chunk in 1:(length(num.round.vec) - 1)) {
      indices.chosen <- indices.stratas[(num.round.vec[i.chunk] + 
                                           1):(num.round.vec[i.chunk + 1])]
      sub.data <- na.omit(as.matrix(big.data[indices.chosen, c(surv.indices, 
                                                       features.indices)]))
      
      gt.sum <- gt.sum + bigSurvSGD::lambdaMaxC(sub.data, strata.size, 
                                    norm.method, features.mean, features.sd)
    }
    return(max(abs(gt.sum))/num.strata.lambda)
  }
  do.one <- function(big.data = big.data, num.rows.chunk, num.rows.big, 
                     norm.method, opt.method, beta.init, beta.type, lr.const, 
                     lr.tau, strata.size, batch.size, num.epoch, b1, b2, eps, 
                     lambda, alpha, bootstrap, features.indices, surv.indices, 
                     features.mean, features.sd) {
    beta.norm <- beta.init
    t <- 0
    m <- rep(0, length(beta.init))
    v <- rep(0, length(beta.init))
    vHat <- rep(0, length(beta.init))
    if (bootstrap) {
      sample.indices.all <- sample(1:num.rows.big, num.rows.big, 
                                   replace = TRUE)
    }
    else {
      sample.indices.all <- sample(1:num.rows.big, num.rows.big, 
                                   replace = FALSE)
    }
    num.rows.chunk <- strata.size * floor(num.rows.chunk/strata.size)
    num.round <- num.rows.big/num.rows.chunk
    num.round.vec <- cumsum(c(0, rep(num.rows.chunk, floor(num.round)), 
                              if ((num.round%%1) > 0) {
                                ceiling((num.round%%1) * num.rows.chunk)
                              }))
    for (n_e in 1:num.epoch) {
      sample.indices <- sample(sample.indices.all, num.rows.big, 
                               replace = FALSE)
      for (i.chunk in 1:(length(num.round.vec) - 1)) {
        indices.chosen <- sample.indices[(num.round.vec[i.chunk] + 
                                            1):(num.round.vec[i.chunk + 1])]
        sub.data <- na.omit(as.matrix(big.data[indices.chosen, 
                                       c(surv.indices, features.indices)]))
        time.unique.sorted <- sort(unique(sub.data[, 
                                                   1]))
        sd.jitter <- 0.1 * min(abs(time.unique.sorted[2:length(time.unique.sorted)] - 
                                     time.unique.sorted[1:(length(time.unique.sorted) - 
                                                             1)]))
        sub.data[, 1] <- sub.data[, 1] + rnorm(NROW(sub.data), 
                                               mean = 0, sd = sd.jitter)
        oneChunkResult <- bigSurvSGD::oneChunkC(sub.data, beta.init, 
                                    beta.type, strata.size, batch.size, t, m, v, 
                                    vHat, lr.const, lr.tau, opt.method, norm.method, 
                                    b1, b2, eps, lambda, alpha, features.mean, 
                                    features.sd)
        beta.init <- oneChunkResult$beta
        beta.ave <- oneChunkResult$betaAve
        t <- oneChunkResult$t
        tAve <- oneChunkResult$tAve
        m <- oneChunkResult$m
        v <- oneChunkResult$v
        vHat <- oneChunkResult$vHat
        if (beta.type == "averaged") {
          beta.norm <- ((t - tAve) * beta.norm + tAve * 
                          oneChunkResult$betaAve)/t
        }
        else {
          beta.norm <- beta.init
        }
      }
    }
    return(beta.norm)
  }
  do.plugin <- function(k, big.data = big.data, num.rows.chunk, 
                        num.rows.big, norm.method, beta.hat, strata.size, num.sample.strata, 
                        surv.indices, features.indices, features.mean, features.sd) {
    num.round <- num.sample.strata/floor((num.rows.big - 
                                            1)/(strata.size - 1))
    num.round.vec <- c(rep(floor((num.rows.big - 1)/(strata.size - 
                                                       1)), floor(num.round)), if ((num.round%%1) > 0) {
                                                         ceiling((num.round%%1) * floor((num.rows.big - 1)/(strata.size - 
                                                                                                              1)))
                                                       })
    indices.Not.k <- (1:num.rows.big)[-k]
    indices.stratas <- NULL
    for (i in 1:length(num.round.vec)) {
      indices.stratas <- c(indices.stratas, sample(indices.Not.k, 
                                                   num.round.vec[i] * (strata.size - 1), replace = FALSE))
    }
    r.cum <- 0
    h.cum <- 0
    num.rows.chunk <- (strata.size - 1) * floor(num.rows.chunk/(strata.size - 
                                                                  1))
    num.round <- length(indices.stratas)/num.rows.chunk
    num.round.vec <- cumsum(c(0, rep(num.rows.chunk, floor(num.round)), 
                              if ((num.round%%1) > 0) {
                                ceiling((num.round%%1) * num.rows.chunk)
                              }))
    for (i.chunk in 1:(length(num.round.vec) - 1)) {
      indices.chosen <- indices.stratas[(num.round.vec[i.chunk] + 
                                           1):(num.round.vec[i.chunk + 1])]
      sub.data <- na.omit(as.matrix(big.data[c(k, indices.chosen), 
                                     c(surv.indices, features.indices)]))
      oneObsResults <- bigSurvSGD::oneObsPlugingC(sub.data, beta.hat, 
                                      strata.size, norm.method, features.mean, features.sd)
      r.cum <- r.cum + oneObsResults$Grad
      h.cum <- h.cum + oneObsResults$Hessian
    }
    r.cum <- r.cum/num.sample.strata
    h.cum <- h.cum/num.sample.strata
    r.cum <- matrix(r.cum, ncol = 1) %*% matrix(r.cum, nrow = 1)
    list(r.cum = r.cum, h.cum = h.cum)
  }
  if (is.null(beta.init)) {
    beta.init <- matrix(0, length(features.indices), 1)
  }
  if (!is.null(alpha) & is.null(lambda)) {
    if (alpha < 0) {
      warning("alpha < 0; set to 0")
      alpha <- 0
    }
    if (alpha > 1) {
      warning("alpha > 1; set to 1")
      alpha <- 1
    }
    lam.max <- lambda.max(big.data, strata.size, num.rows.big, 
                          norm.method, num.strata.lambda, features.indices, 
                          surv.indices) * lambda.scale
    if (length(features.indices) < num.rows.big) {
      lam.min <- 0.01 * lam.max
    }
    else {
      lam.min <- 1e-04 * lam.max
    }
    lambdaAll <- round(lam.max * (lam.min/lam.max)^seq(0, 
                                                       1, l = nlambda), 10)
    beta.hat <- matrix(0, length(beta.init), length(lambdaAll))
    colnames(beta.hat) <- paste0("lambda=", lambdaAll)
    rownames(beta.hat) <- sub.col.names
    for (l in 1:length(lambdaAll)) {
      results.hat <- do.one(big.data = big.data, num.rows.chunk,
                            num.rows.big, norm.method, opt.method, beta.init,
                            beta.type = "single", lr.const, lr.tau, strata.size,
                            batch.size, num.epoch, b1, b2, eps, lambda = lambdaAll[l],
                            alpha, bootstrap = FALSE, features.indices, surv.indices,
                            features.mean, features.sd)
      beta.init <- results.hat
      beta.hat[, l] <- results.hat
    }
    beta.hat[features.sd != 0, ] <- beta.hat[features.sd != 
                                               0, ]/matrix(rep(features.sd[features.sd != 0], length(lambdaAll)), 
                                                           nrow(beta.hat), length(lambdaAll), byrow = FALSE)
    beta.hat.exp <- exp(beta.hat)
    lambda <- lambdaAll
  }
  else if (!is.null(alpha) & !is.null(lambda)) {
    if (alpha < 0) {
      warning("alpha < 0; set to 0")
      alpha <- 0
    }
    if (alpha > 1) {
      warning("alpha > 1; set to 1")
      alpha <- 1
    }
    if (lambda < 0) {
      warning("lambda < 0; set to 0")
      lambda <- 0
    }
    beta.hat <- do.one(big.data = big.data, num.rows.chunk, 
                       num.rows.big, norm.method, opt.method, beta.init, 
                       beta.type = "single", lr.const, lr.tau, strata.size, 
                       batch.size, num.epoch, b1, b2, eps, lambda, alpha, 
                       bootstrap = FALSE, features.indices, surv.indices, 
                       features.mean, features.sd)
    beta.hat[features.sd != 0] <- beta.hat[features.sd != 
                                             0]/features.sd[features.sd != 0]
    names(beta.hat) <- sub.col.names
    beta.hat.exp <- exp(beta.hat)
  }
  else {
    alpha <- 0
    lambda <- 0
    beta.hat <- do.one(big.data = big.data, num.rows.chunk, 
                       num.rows.big, norm.method, opt.method, beta.init, 
                       beta.type, lr.const, lr.tau, strata.size, batch.size, 
                       num.epoch, b1, b2, eps, lambda = 0, alpha = 0, bootstrap = FALSE, 
                       features.indices, surv.indices, features.mean, features.sd)
    if (inference.method == "bootstrap") {
      if (!parallel.flag) {
        beta.boot <- matrix(0, length(beta.hat), num.boot)
        for (i in 1:num.boot) {
          beta.boot[, i] <- do.one(big.data = big.data, 
                                   num.rows.chunk, num.rows.big, norm.method, 
                                   opt.method = boot.method, beta.hat, beta.type, 
                                   lr.const = lr.const.boot, lr.tau = lr.tau.boot, 
                                   strata.size, batch.size, num.epoch.boot, 
                                   b1, b2, eps, lambda = 0, alpha = 0, bootstrap = TRUE, 
                                   features.indices, surv.indices, features.mean, 
                                   features.sd)
        }
      }
      else {
        if (is.null(num.cores)) {
          num.cores <- parallel::detectCores()
        }
        doParallel::registerDoParallel(cores = num.cores)
        beta.boot <- foreach::foreach(i = 1:num.boot, .combine = cbind) %dopar%
          {
            do.one(big.data = big.data, num.rows.chunk, 
                   num.rows.big, norm.method, opt.method = boot.method, 
                   beta.hat, beta.type, lr.const = lr.const.boot, 
                   lr.tau = lr.tau.boot, strata.size, batch.size, 
                   num.epoch.boot, b1, b2, eps, lambda = 0, 
                   alpha = 0, bootstrap = TRUE, features.indices, 
                   surv.indices, features.mean, features.sd)
          }
      }
      quantiles.boot <- apply((beta.boot - matrix(rep(beta.hat, 
                                                      num.boot), length(beta.hat), num.boot, byrow = FALSE)), 
                              1, function(x) quantile(x, probs = c((sig.level/2), 
                                                                   (1 - sig.level/2)), na.rm = TRUE))
      beta.hat <- cbind(beta.hat, beta.hat - quantiles.boot[2, 
      ], beta.hat - quantiles.boot[1, ])
      rownames(beta.hat) <- sub.col.names
      colnames(beta.hat) <- c("estimate", paste0("lower ", 
                                                 (1 - sig.level), "%CI"), paste0("upper ", (1 - 
                                                                                              sig.level), "%CI"))
      beta.hat[features.sd != 0, 1:3] <- beta.hat[features.sd != 
                                                    0, 1:3]/matrix(rep(features.sd[features.sd != 
                                                                                     0], 3), nrow(beta.hat), 3, byrow = FALSE)
      beta.hat.exp <- exp(beta.hat)
    }
    else if (inference.method == "plugin") {
      v.hat <- h.hat <- 0
      if (!parallel.flag) {
        for (k in 1:num.rows.big) {
          results.vh <- do.plugin(k, big.data = big.data, 
                                  num.rows.chunk, num.rows.big, norm.method, 
                                  beta.hat, strata.size, num.sample.strata, 
                                  surv.indices, features.indices, features.mean, 
                                  features.sd)
          v.hat <- v.hat + results.vh$r.cum
          h.hat <- h.hat + results.vh$h.cum
        }
      }
      else {
        if (is.null(num.cores)) {
          num.cores <- parallel::detectCores()
        }
        doParallel::registerDoParallel(cores = num.cores)
        grad.hes <- foreach::foreach(k = 1:num.rows.big) %dopar% 
          {
            do.plugin(k, big.data = big.data, num.rows.chunk, 
                      num.rows.big, norm.method, beta.hat, strata.size, 
                      num.sample.strata, surv.indices, features.indices, 
                      features.mean, features.sd)
          }
        for (i in 1:length(grad.hes)) {
          v.hat <- v.hat + grad.hes[[i]]$r.cum
          h.hat <- h.hat + grad.hes[[i]]$h.cum
        }
      }
      v.hat <- strata.size^2 * v.hat/num.rows.big
      h.hat <- h.hat/num.rows.big
      Sigma.hat <- solve(h.hat) %*% v.hat %*% solve(h.hat)
      se.hat <- sqrt(diag(Sigma.hat)/num.rows.big)
      t.stat <- (beta.hat - beta0)/se.hat
      p.value <- 2 * pt(-abs(t.stat), df = num.rows.big - 
                          1)
      beta.hat <- cbind(beta.hat, beta.hat - qt(1 - (sig.level/2), 
                                                df = num.rows.big) * se.hat, beta.hat + qt(1 - 
                                                                                             (sig.level/2), df = num.rows.big) * se.hat, t.stat, 
                        p.value)
      colnames(beta.hat) <- c("estimate", paste0("lower ", 
                                                 (1 - sig.level), "%CI"), paste0("upper ", (1 - 
                                                                                              sig.level), "%CI"), "z", "p-value")
      rownames(beta.hat) <- sub.col.names
      beta.hat[features.sd != 0, 1:3] <- beta.hat[features.sd != 
                                                    0, 1:3]/matrix(rep(features.sd[features.sd != 
                                                                                     0], 3), nrow(beta.hat), 3, byrow = FALSE)
      beta.hat.exp <- beta.hat
      beta.hat.exp[, c(1, 2, 3)] <- exp(beta.hat[, c(1, 
                                                     2, 3)])
    }
    else {
      names(beta.hat) <- sub.col.names
      beta.hat[features.sd != 0] <- beta.hat[features.sd != 
                                               0]/features.sd[features.sd != 0]
      beta.hat.exp <- exp(beta.hat)
    }
  }
  out <- NULL
  out$coef <- beta.hat
  out$coef.exp <- beta.hat.exp
  out$lambda <- lambda
  out$alpha <- alpha
  out$features.mean <- features.mean
  out$features.sd <- features.sd
  out$call <- match.call()
  class(out) <- "bigSurvSGD"
  out
}


#' Incremental Survival Model Fitting with Pre-Scaled Data 
#' 
#' Loads a previously scaled design matrix and continues the stochastic
#' gradient optimisation for a subset of variables.
#' 
#' @param name.col Character vector containing the column names that should be 
#' included in the partial fit. 
#' @param datapath File system path or connection where the big-memory backing 
#' file for the scaled design matrix is stored.
#' @param ncores Number of processor cores allocated to the partial fitting
#' procedure. Defaults to \code{1}.
#' @param resBigscale Result object returned by \code{\link{bigscale}}
#' containing scaling statistics to be reused. By default the helper reuses the
#' globally cached \code{resultsBigscale} object created by
#' \code{\link{bigscale}}.
#' @param bigmemory.flag Logical flag determining whether big-memory backed
#' matrices are used when loading and updating the design matrix. Defaults to
#' \code{FALSE}.
#' @param parallel.flag Logical flag toggling the use of parallelised
#' stochastic gradient updates. Defaults to \code{FALSE}.
#' @param inf.mth Inference method requested for the partial fit, such as
#' \code{"none"}, \code{"asymptotic"}, or bootstrap summaries. Defaults to
#' \code{"none"}.
#' 
#' @return Either a numeric vector of log hazard-ratio coefficients or, when
#' inference is requested, a matrix whose columns correspond to the inferred
#' coefficient summaries for each penalisation setting.
#' 
#' @export
#' @seealso [bigscale()], [bigSurvSGD.na.omit()] and \link[bigSurvSGD]{bigSurvSGD}.
#' 
#' @examples
#' \donttest{
#' data(micro.censure, package = "bigPLScox")
#' surv_data <- stats::na.omit(
#'   micro.censure[, c("survyear", "DC", "sexe", "Agediag")]
#' )
#' scaled <- bigscale(
#'   survival::Surv(survyear, DC) ~ .,
#'   data = surv_data,
#'   norm.method = "standardize",
#'   batch.size = 16
#' )
#' datapath <- tempfile(fileext = ".csv")
#' utils::write.csv(surv_data, datapath, row.names = FALSE)
#' 
#' continued <- partialbigSurvSGDv0(
#'   name.col = c("Agediag", "sexe"),
#'   datapath = datapath,
#'   ncores = 1,
#'   resBigscale = scaled,
#'   bigmemory.flag = FALSE,
#'   parallel.flag = FALSE,
#'   inf.mth = "none"
#' )
#' }
# 
#' 
partialbigSurvSGDv0 <-
  function(name.col,
           datapath,
           ncores = 1,
           resBigscale,
           bigmemory.flag = FALSE,
           parallel.flag = FALSE,
           inf.mth = "none") {
    time.col <- NULL
    status.col <- NULL
    if (!is.null(resBigscale$col.names)) {
      if (!is.null(resBigscale$time.indices) &&
          length(resBigscale$time.indices) >= 1) {
        time.col <- resBigscale$col.names[resBigscale$time.indices[1]]
      }
      if (!is.null(resBigscale$cens.indices) &&
          length(resBigscale$cens.indices) >= 1) {
        status.col <- resBigscale$col.names[resBigscale$cens.indices[1]]
      }
    }
    if (is.null(time.col) || is.na(time.col)) {
      time.col <- "time"
    }
    if (is.null(status.col) || is.na(status.col)) {
      status.col <- "status"
    }
    feature.terms <- paste(sprintf("`%s`", name.col), collapse = " + ")
    surv.term <- sprintf("survival::Surv(`%s`, `%s`)", time.col, status.col)
    
    form.name.col <-
      stats::as.formula(sprintf("%s ~ %s", surv.term, feature.terms))
    coefres <-
      bigSurvSGD.na.omit(
        formula = form.name.col,
        data = datapath,
        inference.method = inf.mth,
        parallel.flag = parallel.flag,
        num.cores = ncores,
        features.mean = resBigscale$features.mean[name.col],
        features.sd = resBigscale$features.sd[name.col],
        bigmemory.flag = bigmemory.flag
      )$coef
    #  coefres <- coefres[!is.na(coefres)]
    return(coefres)
  }

