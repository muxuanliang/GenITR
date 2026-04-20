library(plyr)
library(reshape2)
library(ggplot2)
library(grid)
source('multiplot.R')
delta.seq <- 0.1*c(0:5)
case.seq <- c(1:5)
result_cov <- result_val <- NULL
for (delta in delta.seq){
    for (case in case.seq){
      res.value <- res.pcid <- res.est <- NULL
      rep <- 0
      for (index in c(1:500)){
        name <- paste0("result_change_setting_ehr/val_delta", delta, "_n", 350,"_N",3200, "_index", index, "_case", case,".csv")
        if(file.exists(name)){
          tmp <-read.csv(name)
          res.value <- rbind(res.value, tmp)
        }
        name <- paste0("result_change_setting_ehr/pcid_delta", delta, "_n", 350, "_N", 3200, "_index", index, "_case", case,".csv")
        if(file.exists(name)){
          tmp <-read.csv(name)
          res.pcid <- rbind(res.pcid, tmp)
        }
        name <- paste0("result_change_setting_ehr/est_delta", delta, "_n", 350,"_N",3200, "_index", index, "_case", case,".csv")
        if(file.exists(name)){
          tmp <-read.csv(name)
          res.est <- rbind(res.est, tmp)
          rep <- rep+1
        }
      }
      if(!is.null(res.est)){
        info <- data.frame(delta = rep(delta, 7), case = rep(case, 7), metric = rep("value", 7))
        result_val <- rbind(result_val, cbind(aggregate(res.value$value, by=list(method=res.value$method, compareFun=res.value$compareFunc), FUN=function(t){mean(t, rm.na=TRUE)}), info))
        info <- data.frame(delta = rep(delta, 7), case = rep(case, 7), metric = rep("pcid", 7))
        result_val <- rbind(result_val, cbind(aggregate(res.pcid$pcid, by=list(method=res.pcid$method, compareFun=res.pcid$compareFunc), FUN=function(t){mean(t, rm.na=TRUE)}), info))

        beta_opt <- aggregate(res.est$est, by=list(index=res.est$index, compareFun=res.est$compareFunc), FUN=function(t){median(t, rm.na=TRUE)})
        sd_beta_opt <- aggregate(res.est$est, by=list(index=res.est$index,compareFun=res.est$compareFunc), FUN=function(t){sd(t)})
        sd_beta_est <- aggregate(res.est$se, by=list(method=res.est$index,compareFun=res.est$compareFunc), FUN=function(t){median(t, rm.na=TRUE)})
        coverage <- (abs(res.est$est - rep(beta_opt$x, times=rep)) < 1.96*res.est$se)
        result_tmp <- data.frame(coverage=coverage, index = rep(beta_opt$index, times=rep), compareFun=rep(beta_opt$compareFun, times=rep))
        result_cov <- rbind(result_cov, cbind(aggregate(result_tmp$coverage, by=list(index=result_tmp$index, compareFun=result_tmp$compareFun), FUN=function(t){mean(t, rm.na=TRUE)}), delta = rep(delta, 12), case = rep(case, 12)))
      }
    }
}
result_infer <- NULL
for (delta in delta.seq){
    for (case in case.seq){
      res <- NULL
      rep <- 0
      for (index in c(1:500)){
        name <- paste0("result_change_setting_ehr/infer_delta", delta, "_n", 350,"_N",3200, "_index", index, "_case", case,".csv")
        if(file.exists(name)){
          tmp <-read.csv(name)
          res <- rbind(res, tmp)
          rep <- rep+1
        }
      }
      if(!is.null(res)){
        info <- data.frame(delta = rep(delta, 7), case = rep(case, 7))

        value_opt <- aggregate((res$upper+res$lower)/2, by=list(method=res$method,compareFun=res$compareFunc), FUN=function(t){median(t, rm.na=TRUE)})
        res$mean <- (res$upper+res$lower)/2
        res$se <- (res$upper-res$lower)/2
        res$coverage <- abs(res$mean- rep(value_opt$x, times=rep))<=res$se
        result_tmp <- data.frame(coverage=res$coverage, method = res$method, compareFun=res$compareFunc)
        result_infer <- rbind(result_infer, cbind(aggregate(result_tmp$coverage, by=list(method=result_tmp$method, compareFun=result_tmp$compareFun), FUN=function(t){mean(t, rm.na=TRUE)}), info))
      }
    }
}
result_val$compareFun <- apply(result_val,1,
                               function(t){
                                 if(t[2]==1){
                                   "Oracle"} else if (t[2]==2){
                                     "EstimatedQ"
                                   } else if (t[2]==3){
                                     "EstimatedQ"
                                   }})
result_val$compareFun <- as.factor(result_val$compareFun)
colnames(result_val)[2] <-"QFun"
result_val$case <- apply(result_val,1, function(t){paste0("Case ",t[5])})
result_val$method <- factor(result_val$method, levels = c('kernel', 'kernel_aug', 'song', 'proposed'))
result_val$method <- revalue(result_val$method,c('kernel'='Kernel', 'kernel_aug'="Kernel Aug","song"="CAL-DR", "proposed"='Proposed'))
result_val$QFun <- factor(result_val$QFun, levels = c('Oracle', 'EstimatedQ', 'Contrast'))
g2.value<- ggplot(result_val[result_val$metric=='value',], aes(delta, x, color=method, linetype=QFun)) + facet_grid(case~., scales = 'free_y')+
  geom_line() + geom_point() + ylab('Generalized value (external reference)') + xlab("Percentage of patients assigned optimal treatment")+
  scale_color_manual(values=c('#6495ed','#ffc66b', '#ff966c', '#ff6666')) + scale_linetype_manual(values=c('dashed', 'solid'))
g2.pcid<- ggplot(result_val[result_val$metric=='pcid',], aes(delta, x, color=method, linetype=QFun)) + facet_grid(case~., scales = 'free_y')+
  geom_line() + geom_point() + ylab('Percentage of making correct decisions (external reference)') + xlab("Percentage of patients assigned optimal treatment")+
  scale_color_manual(values=c('#6495ed','#ffc66b', '#ff966c', '#ff6666')) + scale_linetype_manual(values=c('dashed', 'solid'))
grid_arrange_shared_legend(g2.value,g2.pcid,nrow = 1, ncol = 2, position = 'bottom')
