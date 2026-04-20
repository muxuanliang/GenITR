library(GenITR)

args <- (commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}

#n <- 500
N <- 3200
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
fit_mcid2_ccd <- GenITR(data=data, dataRef=NULL, compareFun = function(y, yr){as.numeric(y >= yr)}, sampleSplitIndex = sampleSplitIndex)#GenITR(data=data, dataRef=dataRef, compareFun = function(y, yr){as.numeric(y >= (yr+1))}, sampleSplitIndex = sampleSplitIndex, propensityEst = data$propensityTRUE)
fit_mcid3_ccd <- GenITR(data=data, dataRef=NULL, compareFun = NULL, sampleSplitIndex = sampleSplitIndex)

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

write.csv(value, file = paste0("result_small_N_cor/val_n", n,"_index", index,"_case",case,".csv"))
write.csv(pcid, file = paste0("result_small_N_cor/pcid_n", n,"_index", index,"_case",case,".csv"))

