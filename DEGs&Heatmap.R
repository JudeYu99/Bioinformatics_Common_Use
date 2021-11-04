# @Author: Yu Zhu
# @Email: 1830416012@stu.suda.edu.cn
# @Timestamp for Creation: 2020-10-29 9:40:12
# @Last Modified by: Yu Zhu
# @Timestamp for Last Modification: 2021-10-09 13:52:17

########################################################################################################################
#
# step0: Set working directory and load the needed packages.
#
########################################################################################################################

setwd("~/test")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("org.Hs.eg.db")


library("limma")
library("org.Hs.eg.db")
library("tidyr")
library("pheatmap")

########################################################################################################################
#
# step1: mRNA process
#
########################################################################################################################

mRNA_input <- read.table("GSE10072_series_matrix.txt", header = T, sep = "\t", comment.char = "!")

# Convert ID_REF to gene symbol
mRNA_input <- separate(mRNA_input, ID_REF, c("ENTREZID", "at"), sep = "_")
mRNA_input <- mRNA_input[,-2]
mRNA_gene <- mRNA_input[, 1]
mRNA_gene_modified <- bitr(mRNA_gene, fromType = "ENTREZID", toType = c("SYMBOL"), OrgDb = org.Hs.eg.db)
mRNA_output <- merge(mRNA_gene_modified, mRNA_input, by = "ENTREZID")

# Get mRNA expression matrix and all mRNA gene symbol profile
write.table(mRNA_output[, -1], file = "mRNA_processed.txt", quote = F, row.names = F, sep = "\t")

# DEG analysis with limma
mRNA_data <- read.table("mRNA_processed.txt", header = TRUE, sep = "\t", row.names = 1)
m_sample <- data.frame(sample = rep(c("C", "T"), c(3, 3)))
rownames(m_sample) <- colnames(mRNA_data)
groups <- m_sample$sample
f <- factor(groups,levels = c("C", "T"))
design <- model.matrix(~ 0 + f)
colnames(design) <- c("C","T")
rownames(design) <- colnames(mRNA_data)
data.fit <- lmFit(mRNA_data, design)
contrast.matrix <- makeContrasts(C-T, levels = design)
data.fit.con <- contrasts.fit(data.fit, contrast.matrix)
data.fit.eb <- eBayes(data.fit.con)
options(digits = 2)
top <- topTable(data.fit.eb, coef = 1, number = nrow(mRNA_data), adjust.method = "BH", lfc = 0)
top <- cbind(gene = rownames(top), top)
mRNA_data <- cbind(gene = rownames(mRNA_data), mRNA_data)

# Filter up-regulated and down-regulated mRNA, p value threshold: adjusted 0.01, logFC threshold: 1
up <- subset(top, adj.P.Val < 0.01 & logFC > 1)
up_mRNA <- merge(mRNA_data, up, by = "gene")
down <- subset(top, adj.P.Val < 0.01 & logFC< -1)
down_mRNA <- merge(mRNA_data, down, by = "gene")

# Save mRNA's gene symbol to files
write.table(up_mRNA[, 1], "up_mRNA.txt", col.names = F, quote = F, row.names = F, sep = "\t")
write.table(down_mRNA[, 1], "down_mRNA.txt", col.names = F, quote = F, row.names = F, sep = "\t")
#write.table(up_mRNA[, c("gene", "logFC")], "up_mRNA_GSEA.txt", col.names = F, quote = F, row.names = F, sep = "\t")
#write.table(down_mRNA[, c("gene", "logFC")], "down_mRNA_GSEA.txt", col.names = F, quote = F, row.names = F, sep = "\t")


########################################################################################################################
#
# step2: Plot heatmap
#
########################################################################################################################

mRNA <- rbind(up_mRNA, down_mRNA)
mRNA <- mRNA[, 1:7]
rownames(mRNA) <- mRNA[, 1]
mRNA <- mRNA[, -1]
annotation <- data.frame(Group = rep(c("C", "T"), c(3, 3)))
rownames(annotation) <- colnames(mRNA)

pdf("mRNA_heatmap.pdf")
pheatmap(mRNA, 
         annotation_col = annotation,
         show_rownames = F,
         annotation_legend = T,
         color = colorRampPalette(c("green", "black", "red"))(300),
         scale="row",
         main = "Heatmap for Differentially Expressed Genes",
         cluster_rows=T, 
         cluster_cols=T)
dev.off()
