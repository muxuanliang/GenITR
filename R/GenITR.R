GenITR <- function(data=list(predictor, treatment, outcome), dataRef=NULL, compareFun = function(y, yr){as.numeric(y >= yr)}, propensityEst = NULL, outcomeEst = NULL, propensityModel = 'kernel', outcomeModel = 'kernel', outcomeFormula = NULL, propensityFormula = NULL,
                   screeningMethod="SIRS", outcomeScreeningFamily = 'Gaussian', standardize = TRUE){
  size <- dim(data$predictor)[1]
  sampleSplitIndex <- (rnorm(size) > 0)

  # standardize or not
  if(standardize){
    data$predictor <- scale(data$predictor)
  }

  # set reference
  if (is.null(dataRef)){
    dataRef <- list(predictor = data$predictor[data$treatment == 0,], outcome = data$outcome[data$treatment == 0])
  }

  # fit kernel regression
  QEst <- apply(cbind(data$outcome, data$predictor), 1, function(t){ks(dataRef$predictor, compareFun(t[1], dataRef$outcome), t[-1])})

  # fit proposed
  fit <- NULL
  fit[[1]] <- GenITRSplit(data = list(predictor = data$predictor, treatment = data$treatment, outcome = QEst), propensityEst = propensityEst, outcomeEst = outcomeEst, sampleSplitIndex = sampleSplitIndex,
                     outcomeModel = outcomeModel, outcomeFormula = outcomeFormula, propensityModel = propensityModel,
                     propensityFormula = propensityFormula, screeningMethod = screeningMethod,
                     outcomeScreeningFamily = outcomeScreeningFamily)
  fit[[2]] <- GenITRSplit(data = list(predictor = data$predictor, treatment = data$treatment, outcome = QEst), propensityEst = propensityEst, outcomeEst = outcomeEst, sampleSplitIndex = (!sampleSplitIndex),
                          outcomeModel = outcomeModel, outcomeFormula = outcomeFormula, propensityModel = propensityModel,
                          propensityFormula = propensityFormula, screeningMethod = screeningMethod,
                          outcomeScreeningFamily = outcomeScreeningFamily)

  # aggregate
  coef <- (fit[[1]]$coef+fit[[2]]$coef)/2
  thresh <- (fit[[1]]$thresh+fit[[2]]$thresh)/2
  cov <- (fit[[1]]$cov+fit[[2]]$cov)/2

  list(coef = coef, thresh = thresh, QEst = QEst, se = sqrt(diag(cov)/size))
}
