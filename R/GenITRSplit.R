GenITRSplit <- function(data = list(predictor, treatment, outcome),
           propensityEst = NULL,
           outcomeEst = NULL,
           sampleSplitIndex = NULL,
           outcomeModel = c('lm', 'glmnet', 'kernel', 'others'),
           outcomeFormula = outcomeFormula,
           propensityModel = c('lm', 'glmnet', 'kernel'),
           propensityFormula = NULL,
           screeningMethod = "SIRS",
           outcomeScreeningFamily = "Gaussian") {
    ## fit propnesity and outcome model
    size <- dim(data$predictor)[1]
    if (is.null(sampleSplitIndex)) {
      sampleSplitIndex <- (rnorm(size) > 0)
    }
    if (is.null(outcomeEst)) {
      predictedOutcomeAll <-
        getOutcomeModel(
          data,
          method = outcomeModel,
          sampleSplitIndex = sampleSplitIndex,
          Formula = outcomeFormula,
          predictAll = TRUE,
          screeningMethod = screeningMethod,
          outcomeScreeningFamily = outcomeScreeningFamily
        )
    } else {
      predictedOutcomeAll <- NULL
      predictedOutcomeAll$control <- outcomeEst$control
      predictedOutcomeAll$treatment <- outcomeEst$treatment
    }
    predictedOutcome <- NULL
    predictedOutcome$control <-
      predictedOutcomeAll$control[sampleSplitIndex]
    predictedOutcome$treatment <-
      predictedOutcomeAll$treatment[sampleSplitIndex]
    if (is.null(propensityEst)) {
      predictedPropensityAll <-
        getPropensityModel(
          data,
          method = propensityModel,
          sampleSplitIndex = sampleSplitIndex,
          Formula = propensityFormula,
          predictAll = TRUE,
          screeningMethod = screeningMethod
        )
    } else {
      predictedPropensityAll <- propensityEst
    }
    predictedPropensity <- predictedPropensityAll[sampleSplitIndex]

    ### estimate coefficient by concordance-assisted learning
    ndim <- NCOL(data$predictor)
    weight_positive <-
      (data$outcome[sampleSplitIndex] - predictedOutcome$treatment) * (data$treatment[sampleSplitIndex] == 1) /
      predictedPropensity + predictedOutcome$treatment
    weight_negative <-
      (data$outcome[sampleSplitIndex] - predictedOutcome$control) * (data$treatment[sampleSplitIndex] == -1) /
      (1 - predictedPropensity) + predictedOutcome$control
    weight_diff <- weight_positive - weight_negative
    pair_diff <- outer(weight_diff, weight_diff, '-')
    targetFun <- function(coef) {
      link <- as.vector(data$predictor[sampleSplitIndex, ] %*% c(1,coef))
      sgn_diff <- (outer(link, link, '-') >= 0)
      - mean(as.vector(pair_diff * sgn_diff))
    }

    msg <- try(fit_coef <- optim(par = rep(0, ndim-1), targetFun, method="BFGS"), silent = TRUE)
    if("try-error" %in% class(msg)){
      fit_coef <- optim(par = rep(0, ndim-1), targetFun)
    }
    diff <- max(abs(fit_coef$par))
    iter <- 1
    while (diff > 10 ^ -3 & iter < 10) {
      betaPast <- fit_coef$par
      fit_coef <- optim(par = betaPast, targetFun)
      diff <-
        max(abs(fit_coef$par - betaPast))
      iter <- iter + 1
    }
    coef <- fit_coef$par

    ### estimate threshold
    link <- as.vector(data$predictor[sampleSplitIndex, ] %*% c(1,coef))
    targetFun <- function(thresh) {
      sgn_est <- (link >= thresh)
      - mean(weight_positive * (sgn_est == TRUE) + weight_negative * (sgn_est ==
                                                                        FALSE))
    }
    fit_thresh <-
      optimize(targetFun, interval = c(min(link), max(link)))
    thresh <- fit_thresh$minimum

    ### estimate sd of coef
    # normalize coef
    sampleSize <- sum(sampleSplitIndex)
    link_diff <- outer(link, link, '-')
    V <- array(0, c(ndim-1, ndim-1))
    Delta <- array(0, c(ndim-1, ndim-1))
    for (i in 1:sampleSize) {
      Deltatmp <- array(0, c(ndim-1, ndim-1))
      for (j in 1:sampleSize) {
        if (i != j) {
          hopt <- 2*sqrt(sum((data$predictor[sampleSplitIndex, ][i, ] - data$predictor[sampleSplitIndex, ][j, ])^2))/sqrt(sampleSize)
          weight_v <- -pair_diff[i, j] * exp(-0.5 * (link_diff[j, i]^2) / (hopt ^ 2)) * link_diff[j, i] /
            sqrt(2 * pi)
          weight_delta <- pair_diff[i, j] * exp(-0.5 * (link_diff[j, i])^2 / (hopt ^ 2)) / sqrt(2 *
                                                                                    pi)
          V <- V + weight_v * (data$predictor[sampleSplitIndex, -1][i, ] - data$predictor[sampleSplitIndex, -1][j, ]) %*% t(data$predictor[sampleSplitIndex, -1][i, ] -
                                                                                                                       data$predictor[sampleSplitIndex, -1][j, ]) / hopt ^ 2
          Deltatmp <- Deltatmp + weight_delta * (data$predictor[sampleSplitIndex, -1][i, ] - data$predictor[sampleSplitIndex, -1][j, ]) / hopt
        }
      }
      Delta <- Delta + Deltatmp %*% t(Deltatmp)
    }
    V <- V / (sampleSize * (sampleSize - 1))
    Delta <- 4 * Delta / (sampleSize^3)
    list(
      coef = c(1,coef),
      thresh = thresh,
      cov = solve(V+0.001*diag(1, ndim-1, ndim-1)) %*% Delta %*% solve(V+0.001*diag(1, ndim-1, ndim-1)),
      D = weight_diff
    )
  }
