library(doParallel)
n.cores <- detectCores()
cl <- makeCluster(n.cores)
registerDoParallel(cl)

load("echo_data_imputed.rda")
result <- foreach(index=1:500, .combine = rbind, .packages = c("GenITR", "MAVE"), .errorhandling = "remove") %dopar% {
# set seed
set.seed(index)
# import data

data <- list(outcome=echo_data_imputed$icu_los_day*(1-echo_data_imputed$mort_28_day_int)+200*echo_data_imputed$mort_28_day_int,
             treatment=echo_data_imputed$echo_int,
             predictor=cbind(age=echo_data_imputed$age, gender=echo_data_imputed$gender, weight=echo_data_imputed$weight, saps=echo_data_imputed$saps,
                             sofa=echo_data_imputed$sofa, elix_score=echo_data_imputed$elix_score,
                             echo_data_imputed[,26:110])) # 110
data$predictor <- data$predictor[,!(endsWith(colnames(data$predictor), "max")|endsWith(colnames(data$predictor), "min"))]
data$predictor<- data$predictor[,(!endsWith(colnames(data$predictor), "first")|startsWith(colnames(data$predictor), "vs"))]
conti_var_index<- colnames(data$predictor) %in% colnames(data$predictor)[c(4,5,which(startsWith(colnames(data$predictor), "vs")))]
data$predictor[,conti_var_index] <- apply(data$predictor[,conti_var_index], 2, scale)
data$predictor[,c(1,3)] <- apply(data$predictor[,c(1,3)],2,function(t){
  quan <- quantile(t, probs = c(0.5))
  sapply(t, function(s){
    sum((s-quan)>0)
  })
})
data$predictor[,6] <- sapply(data$predictor[,6], function(s){
  quan <- c(0,4)
  sum((s-quan)>0)
})
data$predictor[,c(1,2,3,6,7)] <- apply(data$predictor[,c(1,2,3,6,7)], 2, as.factor)
expr <- "data$predictor <- model.matrix(~("
for (i in colnames(data$predictor)){
  expr <- paste0(expr, "data$predictor$", i,"+")
}
expr <- paste0(expr, "0))")
eval(expr=parse(text=expr))
data$predictor <- data$predictor[,-1]
conti_var_index <- which(apply(data$predictor,2,function(t){sum((t!=0)&(t!=1))})>0)
predictor_conti <- as.data.frame(data$predictor[, conti_var_index])
predictor_discrete <- as.data.frame(data$predictor[, !(c(1:ncol(data$predictor)) %in% conti_var_index)])
expr <- "predictor_conti <- model.matrix(~("
for (i in colnames(predictor_conti)){
  expr <- paste0(expr, "predictor_conti$'", i,"'+")
}
expr <- paste0(expr, "0)*(")
for (i in colnames(predictor_conti)){
  expr <- paste0(expr, "predictor_conti$'", i,"'+")
}
expr <- paste0(expr, "0))")
eval(expr=parse(text=expr))

sd_discrete <- apply(predictor_discrete,2,sd)
predictor_discrete <- predictor_discrete[, order(sd_discrete, decreasing = TRUE)[1:10]]

expr <- "predictor_discrete <- model.matrix(~("
for (i in colnames(predictor_discrete)){
  expr <- paste0(expr, "predictor_discrete$'", i,"'+")
}
expr <- paste0(expr, "-1)*(")
for (i in colnames(predictor_discrete)){
  expr <- paste0(expr, "predictor_discrete$'", i,"'+")
}
expr <- paste0(expr, "-1))")
eval(expr=parse(text=expr))
#sd_discrete <- apply(predictor_discrete,2,sd)
#predictor_discrete <- predictor_discrete[, order(sd_discrete, decreasing = TRUE)[1:5]]
data$predictor <- cbind(predictor_conti, predictor_discrete)
vars <- apply(data$predictor,2,sd)
data$predictor <- data$predictor[,vars>0.1]
print(colnames(data$predictor))

total_sample <- NROW(data$predictor)
external_index <- (rnorm(total_sample)>0)
training_index <- (rnorm(total_sample)>0)&(rnorm(total_sample)>0)
data_train <- list(outcome=data$outcome[training_index & !external_index], treatment=data$treatment[training_index & !external_index], predictor=data$predictor[training_index & !external_index,])
data_test <- list(outcome=data$outcome[!training_index & !external_index], treatment=data$treatment[!training_index & !external_index], predictor=data$predictor[!training_index & !external_index,])
data_ext <- list(outcome=data$outcome[training_index|external_index], treatment=data$treatment[training_index|external_index], predictor=data$predictor[training_index|external_index,])
data_control <- list(outcome=data$outcome[data$treatment==0], treatment=data$treatment[data$treatment==0], predictor=data$predictor[data$treatment==0,])


# fit proposed
sampleSplitIndex <- (rnorm(length(data_train$treatment)) > 0)
weight <- rep(0, length(data_ext$treatment))
weight[data_ext$treatment==1] <- runif(sum(data_ext$treatment==1))
weight[data_ext$treatment==0] <- runif(sum(data_ext$treatment==0), min = 1, max = 2)

s_seq <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
for (s in s_seq){
  expr <- paste0("dataRef",s,"<-list(outcome=data_ext$outcome[weight<=",s," | weight>=1+",s,"],
                 treatment=data_ext$treatment[weight<=",s," | weight>=1+",s,"], predictor=data_ext$predictor[weight<=",s," | weight>=1+",s,",])")
  eval(expr=parse(text=expr))
  expr <- paste0("fit_mcid1_ccd_",s," <- GenITR(data=data_train, dataRef=dataRef",s,", compareFun = function(y, yr){as.numeric(y >= yr)}, sampleSplitIndex = sampleSplitIndex, standardize = FALSE, outcomeModel = 'dr', propensityModel = 'dr')")
  eval(expr=parse(text=expr))
}

fit_mcid1_ccd <- fit_mcid1_ccd_0
fit_mcid2_ccd <- GenITR(data=data_train, dataRef=NULL, compareFun = function(y, yr){as.numeric(y >= yr)}, sampleSplitIndex = sampleSplitIndex, standardize = FALSE)#GenITR(data=data, dataRef=dataRef, compareFun = function(y, yr){as.numeric(y >= (yr+1))}, sampleSplitIndex = sampleSplitIndex, propensityEst = data$propensityTRUE)
fit_mcid3_ccd <- GenITR(data=data_train, dataRef=NULL, compareFun = NULL, sampleSplitIndex = sampleSplitIndex, standardize = FALSE)

# testing
compareFun = function(y, yr){as.numeric(y >= yr)}
#estimate Q
fitMAVE <- MAVE::mave(outcome~predictor, data=data_control, method="csMAVE")
selectDim <- MAVE::mave.dim(fitMAVE)
reducedDim <- fitMAVE$dir[[selectDim$dim.min]]
QEst <- apply(cbind(data$outcome, data$predictor %*% reducedDim), 1, function(t){ks(data_control$predictor %*% reducedDim, compareFun(t[1], data_control$outcome), array(t[-1], c(1,selectDim$dim.min)))})
# propensity
predictedPropensityAll <-
  getPropensityModel(
    data,
    method = "dr",
    sampleSplitIndex = !training_index & !external_index
  )
predictedPropensity<- predictedPropensityAll[!training_index & !external_index]
# outcome model
predictedOutcomeAll <- getOutcomeModel(
  data = list(
    predictor = data$predictor,
    treatment = data$treatment,
    outcome = QEst
  ),
  method = 'dr',
  sampleSplitIndex = !training_index & !external_index
)
predictedOutcome <- list(control=predictedOutcomeAll$control[!training_index & !external_index], treatment=predictedOutcomeAll$treatment[!training_index & !external_index])

calValue <- function(predictedTreatment){
  aipw <-
    (predictedTreatment == data$treatment[!training_index & !external_index]) / (
      predictedPropensity * data$treatment[!training_index & !external_index] + (1 - predictedPropensity) *
        (1 - data$treatment[!training_index & !external_index])
    ) * (
      QEst[!training_index & !external_index] - data$treatment[!training_index & !external_index] * predictedOutcome$treatment -
        (1 - data$treatment[!training_index & !external_index]) * predictedOutcome$control
    ) +
    (predictedTreatment > 0) * predictedOutcome$treatment + (predictedTreatment <=
                                                               0) * predictedOutcome$control
  mean(aipw, na.rm = TRUE)
}

value_mcid1_ccd <- calValue(data_test$predictor %*% fit_mcid1_ccd$coef>=fit_mcid1_ccd$thresh)
value_mcid1_kernel <- calValue(ks(data_train$predictor[(data_train$treatment==1),], fit_mcid1_ccd$QEst[(data_train$treatment==1)], data_test$predictor) - ks(data_train$predictor[(data_train$treatment==0),], fit_mcid1_ccd$QEst[(data_train$treatment==0)], data_test$predictor)>0)
value_mcid1_kernel_aug <- calValue(ks(data_train$predictor[sampleSplitIndex,], fit_mcid1_ccd$D[[1]], data_test$predictor) + ks(data_train$predictor[!sampleSplitIndex,], fit_mcid1_ccd$D[[2]], data_test$predictor)>0)
value_mcid2_ccd <- calValue(data_test$predictor %*% fit_mcid2_ccd$coef>=fit_mcid2_ccd$thresh)
value_mcid2_kernel <- calValue(ks(data_train$predictor[(data_train$treatment==1),], fit_mcid2_ccd$QEst[(data_train$treatment==1)], data_test$predictor) - ks(data_train$predictor[(data_train$treatment==0),], fit_mcid2_ccd$QEst[(data_train$treatment==0)], data_test$predictor)>0)
#value_mcid2_kernel_aug <- calValue(ks(data_train$predictor[sampleSplitIndex,], fit_mcid2_ccd$D[[1]][sampleSplitIndex], data_test$predictor) + ks(data_train$predictor[!sampleSplitIndex,], fit_mcid2_ccd$D[[2]], data_test$predictor)>0)
value_mcid3_ccd <- calValue(data_test$predictor %*% fit_mcid3_ccd$coef>=fit_mcid3_ccd$thresh)

#value <- data.frame(value=c(value_mcid1_ccd, value_mcid1_kernel, value_mcid1_kernel_aug, value_mcid2_ccd, value_mcid2_kernel, value_mcid2_kernel_aug), compareFunc = rep(1:2, each=3), method = rep(c('proposed', 'kernel', 'kernel_aug'), by =2))
value <- data.frame(value=c(value_mcid1_ccd, value_mcid1_kernel, value_mcid2_ccd, value_mcid2_kernel), compareFunc = rep(1:2, each=2), method = rep(c('proposed', 'kernel'), by =2))
value <- rbind(value, data.frame(value=value_mcid3_ccd, compareFunc=3, method="song"))


value_mcid1_ccd_1 <- calValue(data_test$predictor %*% fit_mcid1_ccd_1$coef>=fit_mcid1_ccd_1$thresh)
value_mcid1_ccd_2 <- calValue(data_test$predictor %*% fit_mcid1_ccd_0.2$coef>=fit_mcid1_ccd_0.2$thresh)
value_mcid1_ccd_3 <- calValue(data_test$predictor %*% fit_mcid1_ccd_0.4$coef>=fit_mcid1_ccd_0.4$thresh)
value_mcid1_ccd_4 <- calValue(data_test$predictor %*% fit_mcid1_ccd_0.6$coef>=fit_mcid1_ccd_0.6$thresh)
value_mcid1_ccd_5 <- calValue(data_test$predictor %*% fit_mcid1_ccd_0.8$coef>=fit_mcid1_ccd_0.8$thresh)
value_mcid1_ccd_0<- calValue(data_test$predictor %*% fit_mcid1_ccd_0$coef>=fit_mcid1_ccd_0$thresh)

value <- rbind(value, data.frame(value=c(value_mcid1_ccd_0, value_mcid1_ccd_1, value_mcid1_ccd_2,
                                         value_mcid1_ccd_3, value_mcid1_ccd_4, value_mcid1_ccd_5), compareFunc = rep(1, each=6), method = rep(c('proposed'), by =6)))


write.csv(value, file = paste0("result_real_data/val_index", index,"_more_ref.csv"))
}
stopCluster(cl)
