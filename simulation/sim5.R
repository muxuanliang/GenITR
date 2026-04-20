library(GenITR)

args <- (commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}

#n <- 500
#N <- 800
p <- 4

beta <- c(1,1,-1,1)*0.5
gamma1 <- c(1,-1,1,1)*0.5
gamma2 <- c(1,0,-1,0)*0.5
delta <- delta/10

# commandline par
#case <- 1
#index <- 1

# set seed
set.seed(index)
#set rho
rho <- 0.5
variance <- array(0, c(p,p))
for(i in 1:p){
  for(j in 1:p){
    variance[i,j]=rho^{abs(i-j)}
  }
}
discont_effect <- function(eff){
  (eff>0.2)*(eff-0.2)+(eff< -0.2)*(eff+0.2)
}

# sample predictor
generateData <- function(sampleSize, case, reference = FALSE, delta=0){
  x <- MASS::mvrnorm(sampleSize, mu = rep(0, times=p), Sigma = variance) #(array(runif(sampleSize*p), c(sampleSize, p))-0.5) %*% variance
  x <- (x<=-3)*(-3)+x+(x>=3)*3
  x[,1] <- as.numeric(x[,1]>0)
  prob <- apply(x, 1, function(t){exp(0.25*(t[1]^2+t[2]^2+t[1]*t[2]))/(1+exp(0.25*(t[1]^2+t[2]^2+t[1]*t[2])))})
  probRef <- apply(x, 1, function(t){exp(t%*%beta)/(1+exp(t%*%beta))})
  trt <- sapply(prob, function(t){ rbinom(1, 1, prob = t)})
  if (reference){
    x <- 6*(array(runif(sampleSize*p), c(sampleSize, p))-0.5) %*% variance #(array(runif(sampleSize*p), c(sampleSize, p))-0.5) %*% variance
    x[,1] <- as.numeric(x[,1]>0)
    mix <- rbinom(sampleSize, 1, prob = delta)
    probRef <- apply(x, 1, function(t){exp(t%*%beta)/(1+exp(t%*%beta))})
    trt <- mix*sapply(probRef, function(t){ rbinom(1, 1, prob = as.numeric(probRef>0.5))})
  }
  # sample outcome
  mu <- switch(case,
               '3' = 1+sin(x %*% gamma1)+0.5*(x %*% gamma2)^2,
               '4' = 1+x[,1]*x[,2]+0.5*(x[,3])^2,
               '5' = 1+sin(x %*% gamma1)+0.5*(x %*% gamma2)^2)
  d <- switch(case,
              '3' = (x %*% beta)^3,
              '4' = discont_effect((x %*% beta)^3),
              '5' = 0.2*(x[,2])^2+2 * x %*% beta)
  y <- mu + trt * d + rnorm(sampleSize, 0, 0.5-trt*0.3)
  QEst <- pnorm(y-mu, 0, 0.5) * (1-delta * (probRef>0.5))+pnorm(y-mu-d, 0, 0.2) * delta * (probRef>0.5)
  list(predictor=x, treatment=trt, outcome=y, propensityTRUE = prob, outcomeTRUE = list(control=mu+rnorm(sampleSize, 0, 0.5), treatment=mu+d+rnorm(sampleSize, 0, 0.2)), QEst = QEst, d=d)
}
data <- generateData(n, case)
dataRef <- generateData(N, case, reference = TRUE, delta = delta)

# fit proposed
sampleSplitIndex <- (rnorm(n) > 0)
dataQ <- data
dataQ$outcome <- data$QEst
dataRefQ <- dataRef
dataRefQ$outcome <- dataRef$QEst
fit_mcid1_ccd <- GenITR(data=dataQ, dataRef=dataRefQ, compareFun = NULL, sampleSplitIndex = sampleSplitIndex)
fit_mcid2_ccd <- GenITR(data=data, dataRef=dataRef, compareFun = function(y, yr){as.numeric(y >= yr)}, sampleSplitIndex = sampleSplitIndex)#GenITR(data=data, dataRef=dataRef, compareFun = function(y, yr){as.numeric(y >= (yr+1))}, sampleSplitIndex = sampleSplitIndex, propensityEst = data$propensityTRUE)
fit_mcid3_ccd <- GenITR(data=data, dataRef=NULL, compareFun = NULL, sampleSplitIndex = sampleSplitIndex)

# infer value
infer_mcid1_ccd <- GenValueInfer(data=dataQ, dataRef=dataRefQ, compareFun = NULL, method = "concordance")
infer_mcid2_ccd <- GenValueInfer(data=data, dataRef=dataRef, compareFun = function(y, yr){as.numeric(y >= yr)}, method = "concordance")
infer_mcid1_kernel <- GenValueInfer(data=dataQ, dataRef=dataRefQ, compareFun = NULL, method = "kernel", propensityEst = dataQ$propensityTRUE)
infer_mcid2_kernel <- GenValueInfer(data=data, dataRef=dataRef, compareFun = function(y, yr){as.numeric(y >= yr)}, method = "kernel")
infer_mcid1_kernel_aug <- GenValueInfer(data=dataQ, dataRef=dataRefQ, compareFun = NULL, method = "kernel_aug", propensityEst = dataQ$propensityTRUE)
infer_mcid2_kernel_aug <- GenValueInfer(data=data, dataRef=dataRef, compareFun = function(y, yr){as.numeric(y >= yr)}, method = "kernel_aug")
infer_mcid3_ccd <- GenValueInfer(data=data, dataRef=NULL, compareFun = function(y, yr){as.numeric(y >= yr)}, method = "song")

# testing
data_test <- generateData(10^4, case, reference = TRUE, delta = delta)
test_func <- function(itr){
  mean((itr * data_test$outcomeTRUE$treatment +
          (1-itr) * data_test$outcomeTRUE$control) >= data_test$outcome)
}
pcid_func <- function(itr){
  1-mean(abs(as.numeric(itr)-(data_test$d>=0)))
}
value_mcid1_ccd <- test_func(data_test$predictor %*% fit_mcid1_ccd$coef>=fit_mcid1_ccd$thresh)
value_mcid1_kernel <- test_func(ks(data$predictor[(data$treatment==1),], fit_mcid1_ccd$QEst[(data$treatment==1)], data_test$predictor) - ks(data$predictor[(data$treatment==0),], fit_mcid1_ccd$QEst[(data$treatment==0)], data_test$predictor)>0)
value_mcid1_kernel_aug <- test_func(ks(data$predictor[sampleSplitIndex,], fit_mcid1_ccd$D[[1]], data_test$predictor) + ks(data$predictor[!sampleSplitIndex,], fit_mcid1_ccd$D[[2]], data_test$predictor)>0)
value_mcid2_ccd <- test_func(data_test$predictor %*% fit_mcid2_ccd$coef>=fit_mcid2_ccd$thresh)
value_mcid2_kernel <- test_func(ks(data$predictor[(data$treatment==1),], fit_mcid2_ccd$QEst[(data$treatment==1)], data_test$predictor) - ks(data$predictor[(data$treatment==0),], fit_mcid2_ccd$QEst[(data$treatment==0)], data_test$predictor)>0)
value_mcid2_kernel_aug <- test_func(ks(data$predictor[sampleSplitIndex,], fit_mcid2_ccd$D[[1]], data_test$predictor) + ks(data$predictor[!sampleSplitIndex,], fit_mcid2_ccd$D[[2]], data_test$predictor)>0)
value_mcid3_ccd <- test_func(data_test$predictor %*% fit_mcid3_ccd$coef>=fit_mcid3_ccd$thresh)

pcid_mcid1_ccd <- pcid_func(data_test$predictor %*% fit_mcid1_ccd$coef>=fit_mcid1_ccd$thresh)
pcid_mcid1_kernel <- pcid_func(ks(data$predictor[(data$treatment==1),], fit_mcid1_ccd$QEst[(data$treatment==1)], data_test$predictor) - ks(data$predictor[(data$treatment==0),], fit_mcid1_ccd$QEst[(data$treatment==0)], data_test$predictor)>0)
pcid_mcid1_kernel_aug <- pcid_func(ks(data$predictor[sampleSplitIndex,], fit_mcid1_ccd$D[[1]], data_test$predictor) + ks(data$predictor[!sampleSplitIndex,], fit_mcid1_ccd$D[[2]], data_test$predictor)>0)
pcid_mcid2_ccd <- pcid_func(data_test$predictor %*% fit_mcid2_ccd$coef>=fit_mcid2_ccd$thresh)
pcid_mcid2_kernel <- pcid_func(ks(data$predictor[(data$treatment==1),], fit_mcid2_ccd$QEst[(data$treatment==1)], data_test$predictor) - ks(data$predictor[(data$treatment==0),], fit_mcid2_ccd$QEst[(data$treatment==0)], data_test$predictor)>0)
pcid_mcid2_kernel_aug <- pcid_func(ks(data$predictor[sampleSplitIndex,], fit_mcid2_ccd$D[[1]], data_test$predictor) + ks(data$predictor[!sampleSplitIndex,], fit_mcid2_ccd$D[[2]], data_test$predictor)>0)
pcid_mcid3_ccd <- pcid_func(data_test$predictor %*% fit_mcid3_ccd$coef>=fit_mcid3_ccd$thresh)

value <- data.frame(value=c(value_mcid1_ccd, value_mcid1_kernel, value_mcid1_kernel_aug, value_mcid2_ccd, value_mcid2_kernel, value_mcid2_kernel_aug), compareFunc = rep(1:2, each=3), method = rep(c('proposed', 'kernel', 'kernel_aug'), by =2))
value <- rbind(value, data.frame(value=value_mcid3_ccd, compareFunc=3, method="song"))
pcid <- data.frame(pcid=c(pcid_mcid1_ccd, pcid_mcid1_kernel, pcid_mcid1_kernel_aug, pcid_mcid2_ccd, pcid_mcid2_kernel, pcid_mcid2_kernel_aug), compareFunc = rep(1:2, each=3), method = rep(c('proposed', 'kernel', 'kernel_aug'), by =2))
pcid <- rbind(pcid, data.frame(pcid=pcid_mcid3_ccd, compareFunc=3, method="song"))
est <- data.frame(est=c(fit_mcid1_ccd$coef, fit_mcid2_ccd$coef, fit_mcid3_ccd$coef),
                  se=c(fit_mcid1_ccd$se, fit_mcid2_ccd$se, fit_mcid3_ccd$se),
                  compareFunc = rep(c(1,2,3), each=4), index = sapply(rep(1:4, by = 3), function(t){paste0("beta", t)}))
infer <- data.frame(upper=c(infer_mcid1_ccd$upper, infer_mcid1_kernel$upper, infer_mcid1_kernel_aug$upper, infer_mcid2_ccd$upper, infer_mcid2_kernel$upper, infer_mcid2_kernel_aug$upper, infer_mcid3_ccd$upper),
                    lower=c(infer_mcid1_ccd$lower, infer_mcid1_kernel$lower, infer_mcid1_kernel_aug$lower, infer_mcid2_ccd$lower, infer_mcid2_kernel$lower, infer_mcid2_kernel_aug$lower, infer_mcid3_ccd$lower),
                    compareFunc = c(rep(1:2, each=3),3), method = c(rep(c('proposed', 'kernel', 'kernel_aug'), times=2), "song"))


# write.csv(value, file = paste0("result_change_setting_cor/val_n", n,"_N", N,"_index", index,"_case",case,".csv"))
# write.csv(pcid, file = paste0("result_change_setting_cor/pcid_n", n,"_N", N,"_index", index,"_case",case,".csv"))
# write.csv(est, file = paste0("result_change_setting_cor/est_n", n,"_N", N,"_index", index,"_case",case,".csv"))
# write.csv(infer, file = paste0("result_change_setting_cor/infer_n", n,"_N", N,"_index", index,"_case",case,".csv"))

 write.csv(value, file = paste0("result_change_setting_ehr_cor/val_delta", delta,"_n", n,"_N", N,"_index", index,"_case",case,".csv"))
 write.csv(pcid, file = paste0("result_change_setting_ehr_cor/pcid_delta", delta,"_n", n,"_N", N,"_index", index,"_case",case,".csv"))
 write.csv(est, file = paste0("result_change_setting_ehr_cor/est_delta", delta,"_n", n,"_N", N,"_index", index,"_case",case,".csv"))
 write.csv(infer, file = paste0("result_change_setting_ehr_cor/infer_delta", delta,"_n", n,"_N", N,"_index", index,"_case",case,".csv"))
