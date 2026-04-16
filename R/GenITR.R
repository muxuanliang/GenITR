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
