#' Internal fitting routine for one sample split in GenITR
#'
#' Internal helper function used by `GenITR()` to fit the individualized
#' treatment rule on a single sample split. This function estimates nuisance
#' quantities (outcome regression and propensity score), computes concordance-based
#' weights, estimates the coefficient vector and decision threshold, and returns
#' an estimated covariance matrix for the coefficient estimator.
#'
#' @param data A list containing:
#' \describe{
#'   \item{predictor}{A numeric matrix or data frame of predictors/covariates.}
#'   \item{treatment}{A treatment assignment vector, typically coded as 1 and -1.}
#'   \item{outcome}{An outcome vector.}
#' }
#' @param propensityEst Optional pre-estimated propensity scores for all subjects.
#' If provided, `getPropensityModel()` is not called.
#' @param outcomeEst Optional pre-estimated outcome model values. Should be a list
#' with components `control` and `treatment`. If provided, `getOutcomeModel()`
#' is not called.
#' @param sampleSplitIndex Optional logical vector indicating the subsample used
#' for estimation. If `NULL`, a random split is generated.
#' @param outcomeModel Character string specifying the outcome regression model.
#' Allowed options include `"lm"`, `"glmnet"`, `"kernel"`, and `"others"`.
#' @param outcomeFormula Optional formula used for the outcome model.
#' @param propensityModel Character string specifying the propensity score model.
#' Allowed options include `"lm"`, `"glmnet"`, and `"kernel"`.
#' @param propensityFormula Optional formula used for the propensity model.
#' @param screeningMethod Variable screening method passed to nuisance model fitting.
#' Default is `"SIRS"`.
#' @param outcomeScreeningFamily Family used in outcome screening. Default is
#' `"Gaussian"`.
#' @param thresholdMethod Method for estimating the treatment threshold. Intended
#' options include optimization-based and grid-based approaches.
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Estimated coefficient vector, including the intercept fixed at 1.}
#'   \item{thresh}{Estimated decision threshold.}
#'   \item{cov}{Estimated covariance matrix for the non-intercept coefficients.}
#'   \item{D}{Estimated subject-level contrast quantity `weight_diff`.}
#' }
#'
#' @details
#' This function is an internal computational engine for `GenITR()`. It is not
#' intended to be called directly by package users.
#'
#' @keywords internal
#' @noRd
GenITRSplit <- function(data = list(predictor, treatment, outcome),
           propensityEst = NULL,
           outcomeEst = NULL,
           sampleSplitIndex = NULL,
           outcomeModel = c('lm', 'glmnet', 'kernel', 'others'),
           outcomeFormula = outcomeFormula,
           propensityModel = c('lm', 'glmnet', 'kernel'),
           propensityFormula = NULL,
           screeningMethod = "SIRS",
           outcomeScreeningFamily = "Gaussian", thresholdMethod="optim") {
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
# <<<<<<< HEAD
#     if (thresholdMethod=="optim"){
#       fit_thresh <-
#         optimize(targetFun, interval = c(min(link), max(link)))
#       thresh <- fit_thresh$minimum
#     } else {
#       grid <- link
#       value <- sapply(grid, targetFun)
#       thresh <- grid[which.min(value)]
#     }
#
# =======
    fit_thresh <- sapply(link, targetFun)
    thresh <- link[which.min(fit_thresh)]
#>>>>>>> 914dacf (minor changes)

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
      cov = solve(V+0.001*diag(1, ndim-1, ndim-1)) %*% Delta %*% solve(V+0.001*diag(1, ndim-1, ndim-1)), #0.01
      D = weight_diff
    )
  }
