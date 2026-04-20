library(plyr)
library(reshape2)
library(ggplot2)
library(grid)
source('multiplot.R')
n.seq <- c(350, 500, 800)
case.seq <- c(4:5)
result_infer <- result_cov <- result_val <- NULL
for (n in n.seq){
    for (case in case.seq){
      res <- res.value <- res.pcid <- NULL
      rep <- 0
      for (index in c(1:500)){
        name <- paste0("result_small_N_cor/val_n", n,"_index", index, "_case", case,".csv")
        if(file.exists(name)){
          tmp <-read.csv(name)
          tmp$value <- 2*tmp$value-1
          res.value <- rbind(res.value, tmp)
        }
        name <- paste0("result_small_N_cor/pcid_n", n,"_index", index, "_case", case,".csv")
        if(file.exists(name)){
          tmp <-read.csv(name)
          res.pcid <- rbind(res.pcid, tmp)
        }
      }
      if(!is.null(res.value)){
        info <- data.frame(sample_size = rep(n, 7), case = rep(case, 7), metric = rep("value", 7))
        result_val <- rbind(result_val, cbind(aggregate(res.value$value, by=list(method=res.value$method, compareFun=res.value$compareFunc), FUN=function(t){mean(t, rm.na=TRUE)}), info))
        info <- data.frame(sample_size = rep(n, 7), case = rep(case, 7), metric = rep("pcid", 7))
        result_val <- rbind(result_val, cbind(aggregate(res.pcid$pcid, by=list(method=res.pcid$method, compareFun=res.pcid$compareFunc), FUN=function(t){mean(t, rm.na=TRUE)}), info))
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
#result_val$case <- apply(result_val,1, function(t){paste0("Case ",t[5])})
result_val$case <- apply(result_val,1, function(t){paste0("Case ",as.numeric(t[5])+2)})
result_val$method <- factor(result_val$method, levels = c('kernel', 'kernel_aug', 'song', 'proposed'))
result_val$method <- revalue(result_val$method,c('kernel'='Kernel', 'kernel_aug'="Kernel Aug","song"="CAL-DR", "proposed"='Proposed'))
result_val$QFun <- factor(result_val$QFun, levels = c('Oracle', 'EstimatedQ'))
g1.value<- ggplot(result_val[result_val$metric=='value',], aes(sample_size, x, color=method, linetype=QFun)) + facet_grid(case~., scales = 'free_y')+
  geom_line() + geom_point() + ylab('Generalized value (internal reference)') + xlab("Training sample size")+
  scale_color_manual(values=c('#6495ed','#ffc66b', '#ff966c', '#ff6666')) + scale_linetype_manual(values=c('dashed', 'solid'))
g1.pcid<- ggplot(result_val[result_val$metric=='pcid',], aes(sample_size, x, color=method, linetype=QFun)) + facet_grid(case~., scales = 'free_y')+
  geom_line() + geom_point() + ylab('Percentage of making correct decisions (internal reference)') + xlab("Training sample size")+
  scale_color_manual(values=c('#6495ed','#ffc66b', '#ff966c', '#ff6666')) + scale_linetype_manual(values=c('dashed', 'solid'))
grid_arrange_shared_legend(g1.value,g1.pcid,g2.value,g2.pcid, nrow = 2, ncol = 2, position = 'bottom')

