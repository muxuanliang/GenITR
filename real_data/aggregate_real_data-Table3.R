library(plyr)
library(reshape2)
library(ggplot2)
library(grid)
source('hutch/GenITR/multiplot.R')

result_val <- NULL
res <- res.value <- NULL
rep <- 0
ref <- c(0,0.2,0.4,0.6, 0.8,1.0)
for (index in c(1:500)){
  name <- paste0("real_data/result_real_data/val_index", index,"_more_ref.csv")
  if(file.exists(name)){
    tmp <-read.csv(name)
    res.value <- rbind(res.value, cbind(tmp[6:11,], ref=ref))
  }
}
res.value <- res.value[,-1]
if(!is.null(res.value)){
  result_val <- aggregate(res.value$value, by=list(ref=res.value$ref,method=res.value$method, compareFun=res.value$compareFunc), FUN=function(t){mean(t, rm.na=TRUE)})
  result_sd <- aggregate(res.value$value, by=list(ref=res.value$ref,method=res.value$method, compareFun=res.value$compareFunc), FUN=function(t){sd(t)})
}

res <- res.value <- NULL
rep <- 0
for (index in c(1:500)){
  name <- paste0("hutch/GenITR/result_real_data/val_index", index,".csv")
  if(file.exists(name)){
    tmp <-read.csv(name)
    res.value <- rbind(res.value, tmp[-1,])
  }
}
result_val <- rbind(result_val,
                    cbind(ref=0, aggregate(res.value$value, by=list(method=res.value$method, compareFun=res.value$compareFunc),
                                           FUN=function(t){mean(t, rm.na=TRUE)})))
result_sd <- rbind(result_sd,
                   cbind(ref=0, aggregate(res.value$value, by=list(method=res.value$method, compareFun=res.value$compareFunc),
                                          FUN=function(t){sd(t)})))

result_val$compareFun <- apply(result_val,1,
                               function(t){
                                 if(t[3]==1){
                                   "External"} else if (t[3]==2){
                                     "Internal"
                                   } else if (t[3]==3){
                                     "Contrast"
                                   }})
result_sd$compareFun <- apply(result_sd,1,
                               function(t){
                                 if(t[3]==1){
                                   "External"} else if (t[3]==2){
                                     "Internal"
                                   } else if (t[3]==3){
                                     "Contrast"
                                   }})
# Rename the aggregated columns from 'x' to something descriptive
names(result_val)[names(result_val) == "x"] <- "mean_val"
names(result_sd)[names(result_sd) == "x"] <- "sd_val"

# Merge the mean and standard deviation data frames
plot_data <- merge(result_val, result_sd, by = c("ref", "method", "compareFun"))

# Filter the data for 'proposed' method and 'External' compare function
filtered_data <- subset(plot_data, method == "proposed" & compareFun == "External")

# 1. Extract the specific y-intercept values from your main dataset
val_song <- plot_data$mean_val[plot_data$method == "song"]
val_proposed_internal <- plot_data$mean_val[plot_data$method == "proposed" & plot_data$compareFun == "Internal" & plot_data$ref == 0.0]

# 2. Build the plot with a unified legend
ggplot(filtered_data, aes(x = ref, y = mean_val)) +

  # Main line and points: Map to "Proposed (External)"
  geom_line(aes(color = "Proposed (External)", linetype = "Proposed (External)"), linewidth = 0.5) +
  geom_point(aes(color = "Proposed (External)"), size = 1.5) +

  # Horizontal line for Song / CAL-DR: Map to "CAL-DR"
  geom_hline(aes(yintercept = val_song, color = "CAL-DR", linetype = "CAL-DR"), linewidth = 0.5) +

  # Horizontal line for Proposed Internal: Map to "Proposed (Internal)"
  geom_hline(aes(yintercept = val_proposed_internal, color = "Proposed (Internal)", linetype = "Proposed (Internal)"), linewidth = 0.5) +

  # Define the exact colors for the legend
  scale_color_manual(name = "Method",
                     breaks = c("Proposed (External)", "Proposed (Internal)", "CAL-DR"),
                     values = c("Proposed (External)" = "blue",
                                "Proposed (Internal)" = "red",
                                "CAL-DR" = "red")) +

  # Define the exact linetypes for the legend
  scale_linetype_manual(name = "Method",
                        breaks = c("Proposed (External)", "Proposed (Internal)", "CAL-DR"),
                        values = c("Proposed (External)" = "solid",
                                   "Proposed (Internal)" = "solid",
                                   "CAL-DR" = "dashed")) +

  # Clean up the theme and labels
  theme_bw() +
  labs(x = "Mixture of Treatments vs. Controls",
       y = "Mean Generalized Value") +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",         # Moves the legend to the bottom
    legend.title = element_text(face = "bold")
  )
