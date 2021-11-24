# @Author: Yu Zhu
# @Email: 1830416012@stu.suda.edu.cn
# @Last Modified by: Yu Zhu

########################################################################################################################
### A volcano plot to show Differentially Expressed Gene(DEG).
### @ Packages in need:
###   limma
###   ggplot2
########################################################################################################################

########################################################################################################################
#
# step0: Load the needed packages and set the working directory.
#
########################################################################################################################

library("limma")
library("ggplot2")

setwd("/Users/chnzhuyu/R_wd")

########################################################################################################################
#
# step1: Read data and calculate p values and logFC with limma.
#
########################################################################################################################

# Read gene expression data from file and set the rownames with first column.
eps_data <- read.table("expression.txt", header = T, sep = "\t", row.names = 1)
dim(eps_data)
#[1] 27980     6

# Set group information.
group <- rep(c("C", "T"), each = 3)
design <- model.matrix(~group)

# Filter genes which are less expressed.
keep <- rowSums(eps_data > 0) >= 3
eps_data.filtered <- eps_data[keep,]
dim(eps_data)
#[1] 27980     6

# Utilize Empirical Bayesian model for differential gene expression analysis.
fit <- lmFit(eps_data.filtered, design)
fit <- eBayes(fit, trend = T)
deg <- topTable(fit, coef = ncol(design), n = Inf, sort = "p", lfc = 0)
head(deg)
#> head(deg)
#                logFC  AveExpr          t      P.Value    adj.P.Val        B
#MIR4777   -0.20597272 2.947778 -187.80656 3.532394e-10 9.883638e-06 9.170750
#CACNG4    -0.01105989 2.850322  -77.93669 2.081970e-08 1.719362e-04 8.561121
#EHF        3.22635088 5.120782   76.47008 2.273531e-08 1.719362e-04 8.535608
#SHISAL2A  -0.01033358 2.849958  -73.90145 2.663484e-08 1.719362e-04 8.487839
#ARMC2-AS1  0.05145767 2.870520   71.65755 3.072484e-08 1.719362e-04 8.442595
#PROK1      0.02712044 2.858352   68.17150 3.871043e-08 1.756882e-04 8.364969

# Create a new column with all gene names.
gene <- rownames(deg)

# Get a data frame composed of gene names, logFC and p values.
dataset <- cbind(gene, deg["logFC"], deg["P.Value"])
head(dataset)
#> head(dataset)
#               gene       logFC      P.Value
#MIR4777     MIR4777 -0.20597272 3.532394e-10
#CACNG4       CACNG4 -0.01105989 2.081970e-08
#EHF             EHF  3.22635088 2.273531e-08
#SHISAL2A   SHISAL2A -0.01033358 2.663484e-08
#ARMC2-AS1 ARMC2-AS1  0.05145767 3.072484e-08
#PROK1         PROK1  0.02712044 3.871043e-08

########################################################################################################################
#
# step2: Plot the graph into a pdf file with ggplot2.
#
########################################################################################################################

# Create a pdf file for plotting.
pdf(file = "VolcanoPlot.pdf")
# Set cut offs for plotting.
cut_off_pvalue = 0.05
cut_off_logFC = 0
# Set colors to represent up regulated genes, down regulated genes and stable genes.
dataset$col = ifelse(dataset$P.Value < cut_off_pvalue & abs(dataset$logFC) >= cut_off_logFC, ifelse(dataset$logFC> cut_off_logFC, 'Up', 'Down'), 'Stable')
# Make a volcano plot with ggplot().
ggplot(
  # Set data and point colors for plotting.
  dataset, aes(x = logFC, y = -log10(P.Value), colour = col)) +
  # Set point size.
  geom_point(alpha = 0.4, size = 3.5) +
  # Set different colors for different expression.
  scale_color_manual(values = c("#546de5", "#d2dae2","#ff4757")) +
  # Set plotting range.
  xlim(c(-4, 4)) +
  # Add a title to the plot.
  ggtitle("Volcano Plot") +
  
  # Add grid lines.
  geom_vline(xintercept = c(-1,1),lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue), lty = 4, col = "black", lwd = 0.8) +
  
  # Set axes labels.
  labs(x = "log2 fold change", y = "-log10 p-value")+
  # Set a white background.
  theme_bw()+
  
  # Add legends to the plot.
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right", legend.title = element_blank())

# Shut down the current device and save the pdf file.
dev.off()
