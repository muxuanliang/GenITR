library(GenITR)

load("echo_data_imputed.rda")

# commandline par
#index <- 1

# set seed
set.seed(0)
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
training_index <- rbinom(total_sample, 1, prob=0.25)
data_train <- list(outcome=data$outcome[training_index], treatment=data$treatment[training_index], predictor=data$predictor[training_index,])
data_ext <- list(outcome=data$outcome[!training_index], treatment=data$treatment[!training_index], predictor=data$predictor[!training_index,])

# fit proposed
sampleSplitIndex <- (rnorm(length(data_train$treatment)) > 0)
dataRef <- list(outcome=data_ext$outcome, treatment=data_ext$treatment, predictor=data_ext$predictor)
fit_mcid1_ccd <- GenITR(data=data_train, dataRef=dataRef, compareFun = function(y, yr){as.numeric(y >= yr)}, sampleSplitIndex = sampleSplitIndex, standardize = FALSE, outcomeModel = 'dr', propensityModel = 'dr')
fit_mcid2_ccd <- GenITR(data=data_train, dataRef=NULL, compareFun = function(y, yr){as.numeric(y >= yr)}, sampleSplitIndex = sampleSplitIndex, standardize = FALSE)#GenITR(data=data, dataRef=dataRef, compareFun = function(y, yr){as.numeric(y >= (yr+1))}, sampleSplitIndex = sampleSplitIndex, propensityEst = data$propensityTRUE)
fit_mcid3_ccd <- GenITR(data=data_train, dataRef=NULL, compareFun = NULL, sampleSplitIndex = sampleSplitIndex, standardize = FALSE)

colnames(data_train$predictor)[abs(fit_mcid1_ccd$coef)>=1.96*fit_mcid1_ccd$se]
fit_mcid1_ccd$coef[abs(fit_mcid1_ccd$coef)>=1.96*fit_mcid1_ccd$se]
fit_mcid1_ccd$coef[abs(fit_mcid1_ccd$coef)>=1.96*fit_mcid1_ccd$se]-1.96*fit_mcid1_ccd$se[abs(fit_mcid1_ccd$coef)>=1.96*fit_mcid1_ccd$se]
fit_mcid1_ccd$coef[abs(fit_mcid1_ccd$coef)>=1.96*fit_mcid1_ccd$se]+1.96*fit_mcid1_ccd$se[abs(fit_mcid1_ccd$coef)>=1.96*fit_mcid1_ccd$se]

