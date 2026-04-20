library(GenITR)

args <- (commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}

#n <- 500
#N <- 800
p <- 4

beta <- c(1,1,-1,1) * 0.5
gamma1 <- c(1,-1,1,1) * 0.5
gamma2 <- c(1,0,-1,0) * 0.5

# commandline par
#case <- 1
#index <- 1

# set seed
set.seed(index)
#set rho
rho <- 0
variance <- array(0, c(p,p))
for(i in 1:p){
  for(j in 1:p){
    variance[i,j]=rho^{abs(i-j)}
  }
}
discont_effect <- function(eff){
  (eff>0)*(eff-0)+(eff< -0)*(eff+0)
}

# sample predictor
generateData <- function(sampleSize, case, onlyControl = FALSE){
  x <- MASS::mvrnorm(sampleSize, mu = rep(0, times=p), Sigma = variance) #(array(runif(sampleSize*p), c(sampleSize, p))-0.5) %*% variance
  x <- (x<=-3)*(-3)+x+(x>=3)*3
  prob <- apply(x, 1, function(t){exp(0.25*(t[1]^2+t[2]^2+t[1]*t[2]))/(1+exp(0.25*(t[1]^2+t[2]^2+t[1]*t[2])))})
  trt <- sapply(prob, function(t){ rbinom(1, 1, prob = t)})
  # sample outcome
  mu <- switch(case,
               '1' = 1+x %*% gamma1,
               '2' = 1+x %*% gamma1,
               '3' = 1+sin(x %*% gamma1)+0.5*(x %*% gamma2)^2,
               '4' = 1+x[,1]*x[,2]+0.5*(x[,3])^2)
  d <- switch(case,
              '1' = 2 * x %*% beta,
              '2' = exp(x %*% beta)-1,
              '3' = (x %*% beta)^3,
              '4' = (x %*% beta)^3)
  y <- mu + trt * d + rnorm(sampleSize, 0, 0.5)
  if (onlyControl) {
    y <- mu + rnorm(sampleSize, 0, 0.5)
    trt <- trt * 0
  }
  list(predictor=x, treatment=trt, outcome=y, propensityTRUE = prob, outcomeTRUE = list(control=mu, treatment=mu+d), QEst = pnorm(y-mu-1, 0, 0.5), d=d)
}
data <- generateData(n, case)
dataRef <- generateData(N, case, onlyControl = TRUE)

# fit proposed
sampleSplitIndex <- (rnorm(n) > 0)
dataQ <- data
dataQ$outcome <- data$QEst
dataRefQ <- dataRef
dataRefQ$outcome <- dataRef$QEst
fit_mcid1_ccd <- GenITR(data=dataQ, dataRef=dataRefQ, compareFun = NULL, sampleSplitIndex = sampleSplitIndex)
fit_mcid2_ccd <- GenITR(data=data, dataRef=dataRef, compareFun = function(y, yr){as.numeric(y >= (yr+1))}, sampleSplitIndex = sampleSplitIndex)#GenITR(data=data, dataRef=dataRef, compareFun = function(y, yr){as.numeric(y >= (yr+1))}, sampleSplitIndex = sampleSplitIndex, propensityEst = data$propensityTRUE)

# infer value
infer_mcid1_ccd <- GenValueInfer(data=dataQ, dataRef=dataRefQ, compareFun = NULL, method = "concordance")
infer_mcid2_ccd <- GenValueInfer(data=data, dataRef=dataRef, compareFun = function(y, yr){as.numeric(y >= (yr+1))}, method = "concordance")
infer_mcid1_kernel <- GenValueInfer(data=dataQ, dataRef=dataRefQ, compareFun = NULL, method = "kernel", propensityEst = dataQ$propensityTRUE)
infer_mcid2_kernel <- GenValueInfer(data=data, dataRef=dataRef, compareFun = function(y, yr){as.numeric(y >= (yr+1))}, method = "kernel")
infer_mcid1_kernel_aug <- GenValueInfer(data=dataQ, dataRef=dataRefQ, compareFun = NULL, method = "kernel_aug", propensityEst = dataQ$propensityTRUE)
infer_mcid2_kernel_aug <- GenValueInfer(data=data, dataRef=dataRef, compareFun = function(y, yr){as.numeric(y >= (yr+1))}, method = "kernel_aug")

# testing
x_test <- MASS::mvrnorm(10^4, mu = rep(0, times=p), Sigma = variance)
x_test <- (x_test<=-3)*(-3)+x_test+(x_test>=3)*3
mu_test <- switch(case,
                  '1' = 1+x_test %*% gamma1,
                  '2' = 1+x_test %*% gamma1,
                  '3' = 1+sin(x_test %*% gamma1)+0.5*(x_test %*% gamma2)^2,
                  '4' = 1+x_test[,1]*x_test[,2]+0.5*(x_test[,3])^2)
d_test <- switch(case,
                 '1' = 2 * x_test %*% beta,
                 '2' = exp(x_test %*% beta)-1,
                 '3' = (x_test %*% beta)^3,
                 '4' = (x_test %*% beta)^3)
value_mcid1_ccd <- mean(((x_test %*% fit_mcid1_ccd$coef>=fit_mcid1_ccd$thresh) * d_test - 1 + rnorm(10^4, 0, 0.5)) > rnorm(10^4, 0,0.5))
value_mcid1_kernel <- mean(((ks(data$predictor[(data$treatment==1),], fit_mcid1_ccd$QEst[(data$treatment==1)], x_test) - ks(data$predictor[(data$treatment==0),], fit_mcid1_ccd$QEst[(data$treatment==0)], x_test)>0) * d_test - 1+ rnorm(10^4, 0,0.5)) > rnorm(10^4, 0,0.5))
value_mcid1_kernel_aug <- mean(((ks(data$predictor[sampleSplitIndex,], fit_mcid1_ccd$D[[1]], x_test) + ks(data$predictor[!sampleSplitIndex,], fit_mcid1_ccd$D[[2]], x_test)>0) * d_test - 1 + rnorm(10^4, 0,0.5)) > rnorm(10^4, 0,0.5))
value_mcid2_ccd <- mean(((x_test %*% fit_mcid2_ccd$coef>=fit_mcid2_ccd$thresh) * d_test - 1 + rnorm(10^4, 0,0.5)) > (rnorm(10^4, 0,0.5)))
value_mcid2_kernel <- mean(((ks(data$predictor[(data$treatment==1),], fit_mcid2_ccd$QEst[(data$treatment==1)], x_test) - ks(data$predictor[(data$treatment==0),], fit_mcid2_ccd$QEst[(data$treatment==0)], x_test)>0) * d_test - 1+ rnorm(10^4, 0,0.5)) > (rnorm(10^4, 0,0.5)))
value_mcid2_kernel_aug <- mean(((ks(data$predictor[sampleSplitIndex,], fit_mcid2_ccd$D[[1]], x_test) + ks(data$predictor[!sampleSplitIndex,], fit_mcid2_ccd$D[[2]], x_test)>0) * d_test - 1 + rnorm(10^4, 0,0.5)) > (rnorm(10^4, 0,0.5)))

value <- data.frame(value=c(value_mcid1_ccd, value_mcid1_kernel, value_mcid1_kernel_aug, value_mcid2_ccd, value_mcid2_kernel, value_mcid2_kernel_aug), compareFunc = rep(1:2, each=3), method = rep(c('proposed', 'kernel', 'kernel_aug'), by =2))
est <- data.frame(est=c(fit_mcid1_ccd$coef, fit_mcid2_ccd$coef),
                  se=c(fit_mcid1_ccd$se, fit_mcid2_ccd$se),
                  compareFunc = rep(1:2, each=4), index = sapply(rep(1:4, by = 2), function(t){paste0("beta", t)}))
infer <- data.frame(upper=c(infer_mcid1_ccd$upper, infer_mcid1_kernel$upper, infer_mcid1_kernel_aug$upper, infer_mcid2_ccd$upper, infer_mcid2_kernel$upper, infer_mcid2_kernel_aug$upper),
                    lower=c(infer_mcid1_ccd$lower, infer_mcid1_kernel$lower, infer_mcid1_kernel_aug$lower, infer_mcid2_ccd$lower, infer_mcid2_kernel$lower, infer_mcid2_kernel_aug$lower),
                    compareFunc = rep(1:2, each=3), method = rep(c('proposed', 'kernel', 'kernel_aug'), by =2))
inferSplit1 <- data.frame(upper=c(infer_mcid1_ccd$split$upper[1], infer_mcid1_kernel$split$upper[1], infer_mcid1_kernel_aug$split$upper[1], infer_mcid2_ccd$split$upper[1], infer_mcid2_kernel$split$upper[1], infer_mcid2_kernel_aug$split$upper[1]),
                          lower=c(infer_mcid1_ccd$split$lower[1], infer_mcid1_kernel$split$lower[1], infer_mcid1_kernel_aug$split$lower[1], infer_mcid2_ccd$split$lower[1], infer_mcid2_kernel$split$lower[1], infer_mcid2_kernel_aug$split$lower[1]),
                          compareFunc = rep(1:2, each=3), method = rep(c('proposed', 'kernel', 'kernel_aug'), by =2))
inferSplit2 <- data.frame(upper=c(infer_mcid1_ccd$split$upper[2], infer_mcid1_kernel$split$upper[2], infer_mcid1_kernel_aug$split$upper[2], infer_mcid2_ccd$split$upper[2], infer_mcid2_kernel$split$upper[2], infer_mcid2_kernel_aug$split$upper[2]),
                          lower=c(infer_mcid1_ccd$split$lower[2], infer_mcid1_kernel$split$lower[2], infer_mcid1_kernel_aug$split$lower[2], infer_mcid2_ccd$split$lower[2], infer_mcid2_kernel$split$lower[2], infer_mcid2_kernel_aug$split$lower[2]),
                          compareFunc = rep(1:2, each=3), method = rep(c('proposed', 'kernel', 'kernel_aug'), by =2))

write.csv(value, file = paste0("result2_large_N/val_n", n,"_N", N,"_index", index,"_case",case,".csv"))
write.csv(est, file = paste0("result2_large_N/est_n", n,"_N", N,"_index", index,"_case",case,".csv"))
write.csv(infer, file = paste0("result2_large_N/infer_n", n,"_N", N,"_index", index,"_case",case,".csv"))
write.csv(inferSplit1, file = paste0("result2_large_N/infer1_n", n,"_N", N,"_index", index,"_case",case,".csv"))
write.csv(inferSplit2, file = paste0("result2_large_N/infer2_n", n,"_N", N,"_index", index,"_case",case,".csv"))

