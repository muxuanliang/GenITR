#' Internal value estimation on a single sample split
#'
#' Internal helper function used by `GenValueInfer()` to estimate the value of
#' an individualized treatment rule on one sample split. The function optionally
#' standardizes predictors, constructs a reference dataset, estimates a
#' transformed outcome, fits a treatment rule using one of several methods, and
#' then evaluates the rule on the held-out subsample using an augmented inverse
#' probability weighted estimator.
#'
#' @param data A list containing:
#' \describe{
#'   \item{predictor}{A numeric matrix or data frame of predictors/covariates.}
#'   \item{treatment}{A treatment assignment vector.}
#'   \item{outcome}{An outcome vector.}
#' }
#' @param dataRef Optional reference dataset used to estimate the comparison-based
#' transformed outcome. If `NULL`, the untreated group (`treatment == 0`) from
#' `data` is used as the reference sample.
#' @param sampleSplitIndex Optional logical vector indicating the training split.
#' If `NULL`, a random split is generated.
#' @param compareFun A function comparing an observed outcome `y` with a
#' reference outcome `yr`. The default is
#' `function(y, yr) as.numeric(y >= yr)`. If `NULL`, the observed outcome is
#' used directly.
#' @param method Character string specifying the method used to estimate the
#' treatment rule. Supported options are `"concordance"`, `"kernel"`,
#' `"kernel_aug"`, and `"song"`.
#' @param propensityEst Optional pre-estimated propensity scores.
#' @param outcomeEst Optional pre-estimated outcome regression values. If
#' provided, it should be a list with components `control` and `treatment`.
#' @param propensityModel Character string specifying the propensity score model.
#' Default is `"kernel"`.
#' @param outcomeModel Character string specifying the outcome regression model.
#' Default is `"kernel"`.
#' @param outcomeFormula Optional formula used for the outcome model.
#' @param propensityFormula Optional formula used for the propensity model.
#' @param screeningMethod Variable screening method used in nuisance model
#' fitting. Default is `"SIRS"`.
#' @param outcomeScreeningFamily Family used for outcome screening. Default is
#' `"Gaussian"`.
#' @param standardize Logical; whether to standardize predictors before fitting.
#' Default is `TRUE`.
#'
#' @return A list with components:
#' \describe{
#'   \item{value}{Estimated value of the treatment rule.}
#'   \item{se}{Estimated standard error of the value estimator.}
#'   \item{upper}{Upper bound of the approximate 95\% confidence interval.}
#'   \item{lower}{Lower bound of the approximate 95\% confidence interval.}
#'   \item{sd}{Estimated standard deviation of the influence-function-based
#'   pseudo-outcome.}
#' }
#'
#' @details
#' This function is not intended for direct use by package users. It is the
#' computational engine behind `GenValueInfer()`. Depending on `method`, it
#' estimates the treatment rule using:
#' \enumerate{
#'   \item `"song"`: a rule estimated by `GenITR()` using the observed or
#'   transformed outcome;
#'   \item `"concordance"`: a concordance-based rule estimated by `GenITR()`;
#'   \item `"kernel"`: a kernel-smoothed contrast between treatment groups;
#'   \item `"kernel_aug"`: an augmented kernel rule using fitted contrast
#'   quantities from `GenITR()`.
#' }
#'
#' The final value is computed on the held-out subsample using an augmented
#' inverse probability weighted estimator.
#'
#' @keywords internal
#' @noRd
GenValueInferSplit <-
  function(data = list(predictor, treatment, outcome),
           dataRef = NULL,
           sampleSplitIndex = NULL,
           compareFun = function(y, yr) {
             as.numeric(y >= yr)
           },
           method = c("concordance", "kernel", "kernel_aug", "song"),
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

    # if song's method is implemented
    if (method == 'song'){
      fit <- GenITR(
        data = list(
          predictor = data$predictor[sampleSplitIndex, ],
          treatment = data$treatment[sampleSplitIndex],
          outcome = data$outcome[sampleSplitIndex]
        ),
        dataRef = NULL,
        compareFun = NULL,
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
    } else if (method == 'concordance') {
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
