GenValueInfer <- function(data=list(predictor, treatment, outcome), dataRef=NULL, compareFun = function(y, yr){as.numeric(y >= yr)}, method = c("concordance", "kernel"), propensityEst = NULL, outcomeEst = NULL, propensityModel = 'kernel', outcomeModel = 'kernel', outcomeFormula = NULL, propensityFormula = NULL,
                                  screeningMethod="SIRS", outcomeScreeningFamily = 'Gaussian', standardize = TRUE){
  totalSampleSize <- NROW(data$predictor)
  numberPredictor <- NCOL(data$predictor)
  index_seq <- c(1:totalSampleSize)

  sampleSplitIndex <- (rnorm(totalSampleSize) > 0)
  # standardize or not
  if(standardize){
    data$predictor <- scale(data$predictor)
  }

  # set reference
  if (is.null(dataRef)){
    dataRef <- list(predictor = data$predictor[data$treatment == 0,], outcome = data$outcome[data$treatment == 0])
  }

  # fit kernel regression
  # fit kernel regression
  if (is.null(compareFun)){
    QEst <- data$outcome
  } else {
    QEst <- apply(cbind(data$outcome, data$predictor), 1, function(t){ks(dataRef$predictor, compareFun(t[1], dataRef$outcome), t[-1])})
  }

  # infer on split data
  fit <- NULL
  data$outcome <- QEst
  fit[[1]] <- GenValueInferSplit(data=data, dataRef=NULL, sampleSplitIndex = sampleSplitIndex, compareFun = NULL, method = method, propensityEst = propensityEst, outcomeEst = outcomeEst, propensityModel = propensityModel, outcomeModel = outcomeModel, outcomeFormula = outcomeFormula, propensityFormula = propensityFormula,
                                 screeningMethod=screeningMethod, outcomeScreeningFamily = outcomeScreeningFamily, standardize = FALSE)
  fit[[2]] <- GenValueInferSplit(data=data, dataRef=NULL, sampleSplitIndex = (!sampleSplitIndex), compareFun = NULL, method = method, propensityEst = propensityEst, outcomeEst = outcomeEst, propensityModel = propensityModel, outcomeModel = outcomeModel, outcomeFormula = outcomeFormula, propensityFormula = propensityFormula,
                                 screeningMethod=screeningMethod, outcomeScreeningFamily = outcomeScreeningFamily, standardize = FALSE)
  value <- (fit[[1]]$value+fit[[2]]$value)/2
  se <- sqrt((fit[[1]]$sd^2+fit[[2]]$sd^2)/2)/sqrt(totalSampleSize)
  return(list(value = value, se = se, upper = value+1.96*se, lower = value-1.96*se))
}
