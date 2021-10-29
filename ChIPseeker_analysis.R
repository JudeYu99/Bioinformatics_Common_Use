# @Author: Yu Zhu
# @Email: yzhu99@stu.suda.edu.cn
# @Timestamp for Creation: 2021-09-19 9:40:12
# @Last Modified by: Yu Zhu
# @Timestamp for Last Modification: 2021-10-21 20:33:10

########################################################################################################################
### Lab 1:
###       ChIP Seq Analysis with ChIPseeker package
########################################################################################################################

########################################################################################################################
#
# step0: Install and load the needed package.
#
########################################################################################################################

setwd("~/Epi")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChIPseeker")
BiocManager::install("ReactomePA")

library("ChIPseeker")
library("GenomicFeatures")
library("org.Hs.eg.db")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("ReactomePA")
library("clusterProfiler")

########################################################################################################################
#
# step1: Read peak files and prepare for further analysis.
#
########################################################################################################################

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

files_processed <- list(GSE14664 = "GSE14664.bed",
                        GSE19013 = "GSE19013.bed")

GSE14664 <- readPeakFile(files_processed[[1]])
GSE19013 <- readPeakFile(files_processed[[2]])



########################################################################################################################
#
# step2: ChIP peaks coverage plot
#
########################################################################################################################

pdf("GSE14664_peak_covplot.pdf", width = 4)
covplot(GSE14664)
dev.off()

pdf("GSE19013_peak_covplot.pdf", width = 4)
covplot(GSE19013)
dev.off()

pdf("covplot.pdf")
peak <- GenomicRanges::GRangesList(GSE14664 = readPeakFile(files_processed[[1]]), GSE19013 = readPeakFile(files_processed[[2]]))
covplot(peak, weightCol="V5") + facet_grid(chr ~ .id)
dev.off()


########################################################################################################################
#
# step3: Heatmap of ChIP binding to TSS regions
#
########################################################################################################################

promoter <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000) 

pdf("peakBinding_heatmap.pdf")
opar <- par(no.readonly = T)
par(mfrow=c(1, 2))

tagMatrix1 <- getTagMatrix(GSE14664, windows = promoter)
tagMatrix2 <- getTagMatrix(GSE19013, windows = promoter)

tagHeatmap(tagMatrix1, xlim = c(-3000, 3000), color = "blue", title = "GSE14664")
tagHeatmap(tagMatrix2, xlim = c(-3000, 3000), color = "purple", title = "GSE19013")

par(opar)
dev.off()

peakHeatmap(files_processed, TxDb=txdb, upstream=3000, downstream=3000, color=rainbow(length(files_processed)))



########################################################################################################################
#
# step4: Average Profile of ChIP peaks binding to TSS region
#
########################################################################################################################

tagMatrixList <- lapply(files_processed, getTagMatrix, windows = promoter)

pdf("plotAvgProf.pdf")
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
dev.off()

pdf("plotAvgProf_new.pdf")
##conf定义置信区间，facet决定从上到下还是从左往右
plotAvgProf(tagMatrixList, xlim =c(-3000, 3000), conf=0.95, resample=500, facet="row")
dev.off()



########################################################################################################################
#
# step5: Peak Annotation
#
########################################################################################################################

peakAnno1 <- annotatePeak(files_processed[[1]], 
                         tssRegion = c(-3000, 3000),
                         TxDb = txdb, 
                         annoDb = "org.Hs.eg.db")

peakAnno2 <- annotatePeak(files_processed[[2]], 
                          tssRegion = c(-3000, 3000),
                          TxDb = txdb, 
                          annoDb = "org.Hs.eg.db")



########################################################################################################################
#
# step6: Peak Annotation Visualization
#
########################################################################################################################

# Pie Plot
plotAnnoPie(peakAnno1)
plotAnnoPie(peakAnno2)


# Bar Plot
plotAnnoBar(peakAnno1)
plotAnnoBar(peakAnno2)


# Vennypie Plot
pdf("vennpie1.pdf", width = 11)
vennpie(peakAnno1)
dev.off()

pdf("vennpie2.pdf", width = 11)
vennpie(peakAnno2)
dev.off()


# plotDistToTSS Plot
plotDistToTSS(peakAnno1, title = "Distribution of transcription factor-binding loci relative to TSS")
plotDistToTSS(peakAnno2, title = "Distribution of transcription factor-binding loci relative to TSS")



########################################################################################################################
#
# step7: Functional enrichment analysis
#
########################################################################################################################

pathway1 <- enrichPathway(as.data.frame(peakAnno1)$geneId, organism = "human")
gene1 <- seq2gene(GSE14664, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb = txdb)
pathway2 <- enrichPathway(gene1, organism = "human")
pdf("GSE14664_KEGG.pdf", width = 10)
dotplot(pathway2, title = "Top 10 KEGG Enriched Pathways for GSE14664")
dev.off()

pathway3 <- enrichPathway(as.data.frame(peakAnno2)$geneId, organism = "human")
gene3 <- seq2gene(GSE19013, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb = txdb)
pathway4 <- enrichPathway(gene3, organism = "human")
pdf("GSE19013_KEGG.pdf", width = 10)
dotplot(pathway4, title = "Top 10 KEGG Enriched Pathways for GSE19013")
dev.off()



########################################################################################################################
#
# step8: ChIP peak annotation comparision
#
########################################################################################################################

peakAnnoList <- lapply(files_processed, 
                       annotatePeak, 
                       TxDb = txdb,
                       tssRegion = c(-3000, 3000), 
                       verbose = F)

pdf("plotAnnoBar_all.pdf", height = 3) 
plotAnnoBar(peakAnnoList) 
dev.off()

pdf("plotDistToTSS_all.pdf", height = 3) 
plotDistToTSS(peakAnnoList) 
dev.off()

genes <-lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) <- sub("_", "\n", names(genes)) 
compKEGG <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")

pdf("compKEGG.pdf", width = 10)
dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
dev.off()



########################################################################################################################
#
# step9: Overlap of peaks and annotated genes
#
########################################################################################################################

genes <- lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
pdf("VennPeaksGenes.pdf")
vennplot(genes)
dev.off()


