# @Author: Yu Zhu
# @Email: 1830416012@stu.suda.edu.cn
# @Last Modified by: Yu Zhu
# @Address: Department of Bioinformatics, Medical College, Soochow University

## NOTE: R 4.0.2 code

library(pROC)
library(ggplot2)

###########################################################################
#=======================---------------------------=======================#
#                              ROC and Plot                               #
#=======================---------------------------=======================#
###########################################################################

data <- read.csv("/Users/chnzhuyu/Desktop/EvaluateResult.csv")

roc1 <- roc(data$label, data$probability)
# plot(roc1, print.auc = TRUE, plot = TRUE, print.thres = TRUE)
auc(roc1)
ci.auc(roc1)

# Plot ROC Curve
g <- ggroc(roc1, 
           legacy.axes = TRUE, 
           colour = "darkblue", 
           size = 0.8)
g
g + theme_minimal() + 
    xlab("1 - Specificity") + 
    ylab("Sensitivity") +
    ggtitle("ROC Curve") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    #scale_y_continuous(expand = c(0, 0)) +     
    #scale_x_continuous(expand = c(0, 0)) +
    #geom_ribbon(aes(ymin = 0, ymax = sensitivity), 
    #            fill = "red", 
    #            alpha = 0.3) +
    #geom_area(x = roc$specificities, 
    #          y = roc$sensitivities, 
    #          fill = "blue", 
    #          alpha = 0.2) +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color="darkgrey", 
                 linetype="dashed") +
    annotate("text", 
             x = 0.4, 
             y = 0.7, 
             label = paste("AUC = ", round(auc(roc1), digits = 3)), 
             size = 5, 
             fontface = "bold")

ggsave("ROC_Plot.pdf")
