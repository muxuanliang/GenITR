GenValueInferSplit <-
  function(data = list(predictor, treatment, outcome),
           dataRef = NULL,
           sampleSplitIndex = NULL,
           compareFun = function(y, yr) {
             as.numeric(y >= yr)
           },
           method = c("concordance", "kernel", "kernel_aug"),
           propensityEst = NULL,
           outcomeEst = NULL,
           propensityModel = 'kernel',
           outcomeModel = 'kernel',
           outcomeFormula = NULL,
           propensityFormula = NULL,
           screeningMethod = "SIRS",
           outcomeScreeningFamily = 'Gaussian',
           standardize = TRUE) {
    totalSampleSize <- NROW(data$predictor)
    numberPredictor <- NCOL(data$predictor)

    if (is.null(sampleSplitIndex)) {
      sampleSplitIndex <- (rnorm(totalSampleSize) > 0)
    }

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
    if (is.null(compareFun)) {
      QEst <- data$outcome
    } else {
      fitMAVE <- MAVE::mave(outcome ~ predictor, data = dataRef)
      selectDim <- MAVE::mave.dim(fitMAVE)
      reducedDim <- fitMAVE$dir[[selectDim$dim.min]]
      QEst <-
        apply(cbind(data$outcome, data$predictor %*% reducedDim), 1, function(t) {
          ks(dataRef$predictor %*% reducedDim,
             compareFun(t[1], dataRef$outcome),
             t[-1])
        })
    }

    if (method == 'concordance') {
      fit <- GenITR(
          data = list(
            predictor = data$predictor[sampleSplitIndex, ],
            treatment = data$treatment[sampleSplitIndex],
            outcome = QEst[sampleSplitIndex]
          ),
          dataRef = dataRef,
          compareFun = compareFun,
          propensityEst = propensityEst[sampleSplitIndex],
          outcomeEst = switch(
            is.null(outcomeEst),
            'TRUE'=NULL,
            'FALSE'=list(
              control = outcomeEst$control[sampleSplitIndex],
              treatment = outcomeEst$treatment[sampleSplitIndex]
            )
          ),
          propensityModel = propensityModel,
          outcomeModel = outcomeModel,
          outcomeFormula = outcomeFormula,
          propensityFormula = propensityFormula,
          screeningMethod = screeningMethod,
          outcomeScreeningFamily = outcomeScreeningFamily,
          standardize = FALSE
        )
      predictedTreatment <-
        as.numeric(data$predictor[!sampleSplitIndex, ] %*% fit$coef >= fit$thresh)
    } else if (method == 'kernel') {
      fit <-
        ks(data$predictor[sampleSplitIndex &
                            (data$treatment == 1), ], QEst[sampleSplitIndex &
                                                             (data$treatment == 1)], data$predictor[!sampleSplitIndex, ]) - ks(data$predictor[sampleSplitIndex &
                                                                                                                                                (data$treatment == 0), ], QEst[sampleSplitIndex &
                                                                                                                                                                                 (data$treatment == 0)], data$predictor[!sampleSplitIndex, ])
      predictedTreatment <- as.numeric(fit >= 0)
    } else if (method == 'kernel_aug') {
      fit_pre <-GenITR(
        data = list(
          predictor = data$predictor[sampleSplitIndex, ],
          treatment = data$treatment[sampleSplitIndex],
          outcome = QEst[sampleSplitIndex]
        ),
        dataRef = dataRef,
        compareFun = compareFun,
        propensityEst = propensityEst[sampleSplitIndex],
        outcomeEst = switch(
          is.null(outcomeEst),
          'TRUE'=NULL,
          'FALSE'=list(
            control = outcomeEst$control[sampleSplitIndex],
            treatment = outcomeEst$treatment[sampleSplitIndex]
          )
        ),
        propensityModel = propensityModel,
        outcomeModel = outcomeModel,
        outcomeFormula = outcomeFormula,
        propensityFormula = propensityFormula,
        screeningMethod = screeningMethod,
        outcomeScreeningFamily = outcomeScreeningFamily,
        standardize = FALSE
      )
      fit <- ks(fit_pre$predictor[[1]], fit_pre$D[[1]], data$predictor[!sampleSplitIndex, ]) + ks(fit_pre$predictor[[2]], fit_pre$D[[2]], data$predictor[!sampleSplitIndex, ])
      predictedTreatment <- as.numeric(fit >= 0)
    }

    if (is.null(outcomeEst)) {
      predictedOutcomeAll <- getOutcomeModel(
          data = list(
            predictor = data$predictor,
            treatment = data$treatment,
            outcome = QEst
          ),
          method = 'kernel',
          sampleSplitIndex = (!sampleSplitIndex),
          predictAll = TRUE,
          screeningMethod = "SIRS"
        )
      predictedOutcome <- NULL
      predictedOutcome$control <-
        predictedOutcomeAll$control[!sampleSplitIndex]
      predictedOutcome$treatment <-
        predictedOutcomeAll$treatment[!sampleSplitIndex]
    } else {
      predictedOutcome <- NULL
      predictedOutcome$control <-
        outcomeEst$control[!sampleSplitIndex]
      predictedOutcome$treatment <-
        outcomeEst$treatment[!sampleSplitIndex]
    }
    if (is.null(propensityEst)) {
      predictedPropensityAll <-
        getPropensityModel(
          data = list(
            predictor = data$predictor,
            treatment = data$treatment,
            outcome = QEst
          ),
          method = 'kernel',
          sampleSplitIndex = (!sampleSplitIndex),
          predictAll = TRUE,
          screeningMethod = "SIRS"
        )
      predictedPropensity <- predictedPropensityAll[!sampleSplitIndex]
    } else {
      predictedPropensity <- propensityEst[!sampleSplitIndex]
    }

    aipw <-
      (predictedTreatment == data$treatment[!sampleSplitIndex]) / (
        predictedPropensity * data$treatment[!sampleSplitIndex] + (1 - predictedPropensity) *
          (1 - data$treatment[!sampleSplitIndex])
      ) * (
        QEst[!sampleSplitIndex] - data$treatment[!sampleSplitIndex] * predictedOutcome$treatment -
          (1 - data$treatment[!sampleSplitIndex]) * predictedOutcome$control
      ) +
      (predictedTreatment > 0) * predictedOutcome$treatment + (predictedTreatment <=
                                                                 0) * predictedOutcome$control

    valueMean <- mean(aipw, na.rm = TRUE)
    valueSe <- sd(aipw, na.rm = TRUE) / sqrt(sum(!sampleSplitIndex))

    return(
      list(
        value = valueMean,
        se = valueSe,
        upper = valueMean + 1.96 * valueSe,
        lower = valueMean - 1.96 * valueSe,
        sd = sd(aipw, na.rm = TRUE)
      )
    )
  }
