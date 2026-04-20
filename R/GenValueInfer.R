#' Estimate the value of an individualized treatment rule
#'
#' This function estimates the value of an individualized treatment rule using
#' a sample-splitting procedure. It optionally standardizes predictors, constructs
#' a reference sample, estimates a transformed outcome `QEst`, and then calls
#' `GenValueInferSplit()` on the held-out split to estimate the rule value and
#' its standard error.
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
#' @param compareFun A function comparing an observed outcome `y` with a reference
#' outcome `yr`. The default is `function(y, yr) as.numeric(y >= yr)`. If `NULL`,
#' the observed outcome is used directly.
#' @param method Character string specifying the inference method passed to
#' `GenValueInferSplit()`. Options include `"concordance"`, `"kernel"`, and `"song"`.
#' @param propensityEst Optional pre-estimated propensity scores.
#' @param outcomeEst Optional pre-estimated outcome regression values.
#' @param propensityModel Character string specifying the propensity score model.
#' Default is `"kernel"`.
#' @param outcomeModel Character string specifying the outcome regression model.
#' Default is `"kernel"`.
#' @param outcomeFormula Optional formula used for the outcome model.
#' @param propensityFormula Optional formula used for the propensity model.
#' @param screeningMethod Variable screening method used in nuisance model fitting.
#' Default is `"SIRS"`.
#' @param outcomeScreeningFamily Family used for outcome screening. Default is
#' `"Gaussian"`.
#' @param standardize Logical; whether to standardize predictors before fitting.
#' Default is `TRUE`.
#' @param alpha Proportion of the sample used for the training split. The remaining
#' observations are used for inference. Default is `0.8`.
#'
#' @return A list with components:
#' \describe{
#'   \item{value}{Estimated value of the individualized treatment rule.}
#'   \item{se}{Estimated standard error of the value estimator.}
#'   \item{upper}{Upper bound of the approximate 95\% confidence interval.}
#'   \item{lower}{Lower bound of the approximate 95\% confidence interval.}
#' }
#'
#' @details
#' The procedure first randomly splits the data into a training subset and an
#' inference subset according to `alpha`. If `compareFun` is not `NULL`, the
#' function estimates a transformed outcome `QEst` by comparing each observed
#' outcome to the reference outcomes through sufficient dimension reduction
#' (`MAVE`) and kernel smoothing. The transformed outcome is then passed to
#' `GenValueInferSplit()` to estimate the rule value on the inference split.
#'
#' The reported confidence interval is based on a normal approximation:
#' `value +/- 1.96 * se`.
#'
#' @export

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
