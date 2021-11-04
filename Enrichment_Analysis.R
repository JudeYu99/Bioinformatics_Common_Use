# @Author: Yu Zhu
# @Email: 1830416012@stu.suda.edu.cn
# @Timestamp for Creation: 2020-10-29 9:40:12
# @Last Modified by: Yu Zhu
# @Timestamp for Last Modification: 2021-11-03 10:39:04


########################################################################################################################
#
# step0: Load the needed package.
#
########################################################################################################################

library("limma")
library("org.Hs.eg.db")
library("clusterProfiler")
library("tidyr")
library("ggplot2")

########################################################################################################################
#
#  clusterProfiler Enrichment Analysis
#
########################################################################################################################

input <- read.table("GSE39791.txt", header = T, sep = "\t")

rownames(input) <- input[, 1]
data <- input[, -1]
data <- input[, 2:ncol(data)]

# DEGs Analysis with limma
m_sample <- data.frame(sample = rep(c("N", "T"), c(72,72)))
rownames(m_sample) <- colnames(data)
groups <- m_sample$sample
f <- factor(groups, levels = c("N","T"))
design <- model.matrix(~ 0 + f)
colnames(design) <- c("N","T")
rownames(design) <- colnames(data)
data.fit <- lmFit(data, design)
contrast.matrix <- makeContrasts(N-T,levels = design)
data.fit.con <- contrasts.fit(data.fit, contrast.matrix)
data.fit.eb <- eBayes(data.fit.con)
options(digits = 2)
top <- topTable(data.fit.eb, coef=1, number = nrow(data), adjust.method = "BH", lfc = 0)
top <- cbind(gene = rownames(top), top)
filtered <- subset(top, adj.P.Val < 0.0001)

data <- cbind(gene = rownames(data), data)

# Filter up-regulated and down-regulated genes
up <- subset(top, adj.P.Val < 0.01 & logFC > 1)
up_mRNA <- merge(mRNA_data, up, by = "gene")
down <- subset(top, adj.P.Val < 0.01 & logFC< -1)
down_mRNA <- merge(mRNA_data, down, by = "gene")

# GO Enrichment for up-regulated genes
target_gene_up <- up_mRNA[, 1]
display_number = c(50, 50, 50)
ego_BP <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_up,
                   pvalueCutoff = 0.05,
                   ont = "BP",
                   readable = T)
ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]

ego_MF <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_up,
                   pvalueCutoff = 0.05,
                   ont = "MF",
                   readable = T)
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[2], ]

ego_CC <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_up,
                   pvalueCutoff = 0.05,
                   ont = "CC",
                   readable = T)
ego_result_CC <- na.omit(as.data.frame(ego_CC)[1:display_number[3], ])

pdf("BP_up_mRNA.pdf", width = 8)
barplot(ego_BP, showCategory = 20, title = "Top 20 Biological Process for up-regulated mRNA")
dev.off()

#pdf("CC_up_mRNA.pdf", width = 8)
barplot(ego_CC, showCategory = 20, title = "Top 20 Cellular Component for up-regulated mRNA")
#dev.off()

#pdf("MF_up_mRNA.pdf", width = 8)
barplot(ego_MF, showCategory = 20, title = "Top 20 Molecular Function for up-regulated mRNA")
#dev.off()

ego_BP_plot <- ego_result_BP[1:10,]
ego_CC_plot <- ego_result_CC[1:10,]
ego_MF_plot <- ego_result_MF[1:10,]

go_enrich_df <- data.frame(ID = c(ego_BP_plot$ID, ego_CC_plot$ID, ego_MF_plot$ID),
                           Description = c(ego_BP_plot$Description, ego_CC_plot$Description, ego_MF_plot$Description),
                           GeneNumber = c(ego_BP_plot$Count, ego_CC_plot$Count, ego_MF_plot$Count),
                           PValue = c(ego_BP_plot$p.adjust, ego_CC_plot$p.adjust, ego_MF_plot$p.adjust),
                           Type = factor(c(rep("Biological Process", nrow(ego_BP_plot)), 
                                         rep("Cellular Component", nrow(ego_CC_plot)),
                                         rep("Molecular Function", nrow(ego_MF_plot))), 
                           levels = c("Biological Process", "Cellular Component", "Molecular Function")))

## numbers as data on x axis
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))

## shorten the names of GO terms
shorten_names <- function(x, n_word=4, n_char=40) {
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40)) {
    if (nchar(x) > 40) {
      x <- substr(x, 1, 40)
      x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)], collapse=" "), "...", sep="")
      return(x)
    }
  } 
  else
  {
    return(x)
  }
}

labels <- (sapply(levels(go_enrich_df$Description)[as.numeric(go_enrich_df$Description)], shorten_names))
names(labels) <- rev(1:nrow(go_enrich_df))

p <- ggplot(go_enrich_df, aes(x = number, y = GeneNumber, fill = Type)) +
            geom_bar(stat = "identity", width = 0.6) + 
            coord_flip() + 
            scale_fill_manual(values = c("#2f5688", "#66C3A5", "#CC0000")) + 
            theme_bw() + 
            scale_x_discrete(labels = labels) +
            xlab("GO terms") + 
            theme(axis.text = element_text(face = "bold", color = "gray50")) +
            labs(title = "The Most Enriched GO Terms") 
            #geom_text(aes(label = signif(PValue, digits = 2)), hjust = 0)


pdf("go_enrichment_of_up_mRNA.pdf")
p
dev.off()

# KEGG Enrichment Analysis
kk <- enrichKEGG(gene = target_gene_up,
                 organism = 'hsa',
                 pvalueCutoff = 0.05)

pdf("KEGG_dot_up_mRNA.pdf", width = 8)
dotplot(kk, showCategory = 10, title = "Top 10 KEGG Enriched Pathways")
dev.off()

pdf("KEGG_emapplot_up_mRNA.pdf")
emapplot(kk, showCategory = 10)
dev.off()

cnetplot(kk, showCategory = 5)


##################################################################################################################################################################################

# GO Enrichment for down-regulated genes
target_gene_down <- down_mRNA[, 1]
display_number = c(50, 50, 50)
ego_BP <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_down,
                   pvalueCutoff = 0.05,
                   ont = "BP",
                   readable = T)
ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]

ego_MF <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_down,
                   pvalueCutoff = 0.05,
                   ont = "MF",
                   readable = T)
ego_result_MF <- ego_MF@result[1:display_number[2], ]

ego_CC <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_down,
                   pvalueCutoff = 0.05,
                   ont = "CC",
                   readable = T)
ego_result_CC <- na.omit(as.data.frame(ego_CC)[1:display_number[3], ])

pdf("BP_down_mRNA.pdf", width = 12)
barplot(ego_BP, showCategory = 20, title = "Top 20 Biological Process for down-regulated mRNA")
dev.off()

barplot(ego_CC, showCategory = 20)
barplot(ego_MF, showCategory = 20)

ego_BP_plot <- ego_result_BP[1:10,]
ego_CC_plot <- ego_result_CC
ego_MF_plot <- ego_result_MF[1:10,]

go_enrich_df <- data.frame(ID = c(ego_BP_plot$ID, ego_CC_plot$ID, ego_MF_plot$ID),
                           Description = c(ego_BP_plot$Description, ego_CC_plot$Description, ego_MF_plot$Description),
                           GeneNumber = c(ego_BP_plot$Count, ego_CC_plot$Count, ego_MF_plot$Count),
                           Type = factor(c(rep("Biological Process", nrow(ego_BP_plot)), 
                                         rep("Cellular Component", nrow(ego_CC_plot)),
                                         rep("Molecular Function", nrow(ego_MF_plot))), 
                           levels = c("Biological Process", "Cellular Component", "Molecular Function")))

## numbers as data on x axis
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))

labels <- as.character(go_enrich_df$Description) 
names(labels) <- rev(1:nrow(go_enrich_df))

p <- ggplot(go_enrich_df, aes(x = number, y = GeneNumber, fill = Type)) +
            geom_bar(stat = "identity", width = 0.6) + 
            coord_flip() + 
            scale_fill_manual(values = c("#2f5688", "#66C3A5", "#CC0000")) + 
            theme_bw() + 
            scale_x_discrete(labels = labels) +
            xlab("GO terms") + 
            theme(axis.text = element_text(face = "bold", color = "gray50")) +
            labs(title = "The Most Enriched GO Terms")
            
pdf("go_enrichment_of_down_mRNA.pdf", width = 9)
p
dev.off()

# KEGG Enrichment Analysis
kk <- enrichKEGG(gene = target_gene_down,
                 organism = 'hsa',
                 pvalueCutoff = 0.05)

pdf("KEGG_dot_down_mRNA.pdf")
dotplot(kk, showCategory = 10, title = "Top 10 KEGG Enriched Pathways")
dev.off()

