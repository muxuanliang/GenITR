GenValueInfer <- function(data = list(predictor, treatment, outcome),
           dataRef = NULL,
           compareFun = function(y, yr) {
             as.numeric(y >= yr)
           },
           method = c("concordance", "kernel", "song"),
           propensityEst = NULL,
           outcomeEst = NULL,
           propensityModel = 'kernel',
           outcomeModel = 'kernel',
           outcomeFormula = NULL,
           propensityFormula = NULL,
           screeningMethod = "SIRS",
           outcomeScreeningFamily = 'Gaussian',
           standardize = TRUE, alpha=0.8) {
    totalSampleSize <- NROW(data$predictor)
    numberPredictor <- NCOL(data$predictor)

    sampleSplitIndex <- (rbinom(totalSampleSize,1,alpha) > 0)
    # standardize or not
    if (standardize) {
      data$predictor <- scale(data$predictor)
    }

    # set reference
    if (is.null(dataRef)) {
      dataRef <-
        list(predictor = data$predictor[data$treatment == 0, ],
             outcome = data$outcome[data$treatment == 0])
    }

    # fit kernel regression
    # fit kernel regression
    if (is.null(compareFun)) {
      QEst <- data$outcome
    } else {
      fitMAVE <- MAVE::mave(outcome ~ predictor, data = dataRef, method="csMAVE")
      selectDim <- MAVE::mave.dim(fitMAVE)
      reducedDim <- fitMAVE$dir[[selectDim$dim.min]]
      QEst <-
        apply(cbind(data$outcome, data$predictor %*% reducedDim), 1, function(t) {
          ks(dataRef$predictor %*% reducedDim,
             compareFun(t[1], dataRef$outcome),
             t[-1])
        })
    }

    # infer on split data
    fit <- NULL
    data$outcome <- QEst
    fit <-
      GenValueInferSplit(
        data = data,
        dataRef = NULL,
        sampleSplitIndex = sampleSplitIndex,
        compareFun = NULL,
        method = method,
        propensityEst = propensityEst,
        outcomeEst = outcomeEst,
        propensityModel = propensityModel,
        outcomeModel = outcomeModel,
        outcomeFormula = outcomeFormula,
        propensityFormula = propensityFormula,
        screeningMethod = screeningMethod,
        outcomeScreeningFamily = outcomeScreeningFamily,
        standardize = FALSE
      )
    value <- fit$value
    se <- fit$sd / sqrt(sum(!sampleSplitIndex))
    return(list(
      value = value,
      se = se,
      upper = value + 1.96 * se,
      lower = value - 1.96 * se))
  }
