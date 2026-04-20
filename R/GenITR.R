#' Fit the GenITR estimator with sample splitting and aggregation
#'
#' This function fits the GenITR estimator for individualized treatment rules.
#' It optionally standardizes predictors, constructs a reference group, estimates
#' a transformed outcome quantity `QEst`, applies the `GenITRSplit()` procedure
#' on two complementary sample splits, and then aggregates the results to produce
#' a final normalized coefficient vector, threshold, and standard error estimates.
#'
#' @param data A list containing:
#' \describe{
#'   \item{predictor}{A numeric matrix or data frame of predictors/covariates.}
#'   \item{treatment}{A binary treatment assignment vector.}
#'   \item{outcome}{An outcome vector.}
#' }
#' @param dataRef Optional reference dataset used to estimate the comparison-based
#' transformed outcome. If `NULL`, the untreated group (`treatment == 0`) from `data`
#' is used as the reference.
#' @param compareFun A function comparing an observed outcome `y` with a reference
#' outcome `yr`. The default is `function(y, yr) as.numeric(y >= yr)`. If `NULL`,
#' the raw outcome is used directly as `QEst`.
#' @param sampleSplitIndex Optional logical vector indicating the first sample split.
#' If `NULL`, a random split is generated.
#' @param propensityEst Optional pre-estimated propensity scores.
#' @param outcomeEst Optional pre-estimated outcome regression values.
#' @param propensityModel Character string specifying the propensity score model.
#' Default is `"kernel"`.
#' @param outcomeModel Character string specifying the outcome model.
#' Default is `"kernel"`.
#' @param outcomeFormula Optional formula for fitting the outcome model.
#' @param propensityFormula Optional formula for fitting the propensity model.
#' @param screeningMethod Variable screening method used in `GenITRSplit()`.
#' Default is `"SIRS"`.
#' @param outcomeScreeningFamily Family used for screening the outcome model.
#' Default is `"Gaussian"`.
#' @param standardize Logical; whether to standardize predictors before fitting.
#' Default is `TRUE`.
#' @param thresholdMethod Method used to estimate the treatment threshold in
#' `GenITRSplit()`. Default is `"optim"`.
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Normalized estimated coefficient vector.}
#'   \item{thresh}{Estimated treatment threshold, averaged over the two splits.}
#'   \item{QEst}{Estimated transformed outcome used in the fitting procedure.}
#'   \item{se}{Estimated standard errors of the normalized coefficient vector.}
#'   \item{D}{A list containing the estimated reduced-dimension directions from
#'   the two sample splits.}
#'   \item{predictor}{A list containing the predictor matrices in the two splits.}
#' }
#'
#' @details
#' The procedure works as follows:
#' \enumerate{
#'   \item Standardize predictors if requested.
#'   \item Construct a reference sample if not provided.
#'   \item If `compareFun` is not `NULL`, estimate a comparison-based transformed
#'   outcome `QEst` using sufficient dimension reduction (`MAVE`) and kernel smoothing.
#'   Otherwise, use the observed outcome directly.
#'   \item Apply `GenITRSplit()` on two complementary sample splits.
#'   \item Average the coefficient, threshold, and covariance estimates from the
#'   two splits.
#'   \item Normalize the coefficient vector and compute standard errors using the
#'   delta method transformation.
#' }
#'
#' @export
GenITR <- function(data=list(predictor, treatment, outcome), dataRef=NULL, compareFun = function(y, yr){as.numeric(y >= yr)}, sampleSplitIndex = NULL, propensityEst = NULL, outcomeEst = NULL, propensityModel = 'kernel', outcomeModel = 'kernel', outcomeFormula = NULL, propensityFormula = NULL,
                   screeningMethod="SIRS", outcomeScreeningFamily = 'Gaussian', standardize = TRUE, thresholdMethod = "optim"){
  size <- dim(data$predictor)[1]
  p <- dim(data$predictor)[2]
  data$predictor <- as.matrix(data$predictor)
  if (is.null(sampleSplitIndex)){
    sampleSplitIndex <- (rnorm(size) > 0)
  }

  # standardize or not
  if(standardize){
    data$predictor <- scale(data$predictor)
  }

  # set reference
  if (is.null(dataRef)){
    dataRef <- list(predictor = data$predictor[data$treatment == 0,], outcome = data$outcome[data$treatment == 0])
  }

  # fit kernel regression
  if (is.null(compareFun)){
    QEst <- data$outcome
  } else {
    dataRef$predictor <- as.matrix(dataRef$predictor)
    fitMAVE <- MAVE::mave(outcome~predictor, data=dataRef, method="csMAVE")
    selectDim <- MAVE::mave.dim(fitMAVE)
    reducedDim <- fitMAVE$dir[[selectDim$dim.min]]
    QEst <- apply(cbind(data$outcome, data$predictor %*% reducedDim), 1, function(t){ks(dataRef$predictor %*% reducedDim, compareFun(t[1], dataRef$outcome), array(t[-1], c(1,selectDim$dim.min)))})
  }


  # fit proposed
  fit <- NULL
  fit[[1]] <- GenITRSplit(data = list(predictor = data$predictor, treatment = data$treatment, outcome = QEst), propensityEst = propensityEst, outcomeEst = outcomeEst, sampleSplitIndex = sampleSplitIndex,
                     outcomeModel = outcomeModel, outcomeFormula = outcomeFormula, propensityModel = propensityModel,
                     propensityFormula = propensityFormula, screeningMethod = screeningMethod,
                     outcomeScreeningFamily = outcomeScreeningFamily, thresholdMethod = thresholdMethod)
  fit[[2]] <- GenITRSplit(data = list(predictor = data$predictor, treatment = data$treatment, outcome = QEst), propensityEst = propensityEst, outcomeEst = outcomeEst, sampleSplitIndex = (!sampleSplitIndex),
                          outcomeModel = outcomeModel, outcomeFormula = outcomeFormula, propensityModel = propensityModel,
                          propensityFormula = propensityFormula, screeningMethod = screeningMethod,
                          outcomeScreeningFamily = outcomeScreeningFamily, thresholdMethod = thresholdMethod)

  # aggregate
  coef <- (fit[[1]]$coef+fit[[2]]$coef)/2
  thresh <- (fit[[1]]$thresh+fit[[2]]$thresh)/2
  cov <- (fit[[1]]$cov+fit[[2]]$cov)/2
  cov_full <- matrix(0, p, p)
  cov_full[2:p, 2:p] <- cov

  # normalized
  coef_norm <- coef/sqrt(sum(coef^2))
  trans_matrix <- 1/sqrt(sum(coef^2)) * (pracma::eye(p)-coef %*% t(coef)/sum(coef^2))
  cov_trans <- trans_matrix %*% cov_full %*% t(trans_matrix)

  list(coef = coef_norm, thresh = thresh, QEst = QEst, se = sqrt(diag(cov_trans)/size), D = list(fit[[1]]$D, fit[[2]]$D), predictor = list(data$predictor[sampleSplitIndex,], data$predictor[!sampleSplitIndex,]))
}
