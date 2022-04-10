# Get PPIN (Gene name) input from String DB (Homo sapiens)

# Annotation file from ensembl. (http://uswest.ensembl.org/biomart/martview)
ENSP <- read.table("mart_export.txt", header = TRUE, sep = ",")
ENSP_filt <- ENSP
ENSP_filt[ENSP_filt == ""] <- NA
ENSP_filt <- na.omit(ENSP_filt)

# PPIN input with a threshold of combined score >= 400.
PPI <- read.table("PPI_400.txt", sep = " ", header = T)

temp1 <- data.frame(cbind(paste("9606.", ENSP_filt[,1], sep = "")), ENSP_filt[,2])
colnames(temp1) <- c("protein1","Gene.name")

temp2 <- data.frame(cbind(paste("9606.", ENSP_filt[,1], sep = "")), ENSP_filt[,2])
colnames(temp2) <- c("protein2","Gene.name")

temp3 <- merge(PPI, temp1, by = "protein1")
temp3 <- temp3[, -1]
colnames(temp3)[3] <- "protein1"

temp4 <- merge(temp3, temp2, by = "protein2")
temp4 <- temp4[, -1]
colnames(temp4)[3] <- "protein2"

# Output the PPIN with gene name annotation.
PPI_result <- cbind(temp4[, 2:3], temp4[, 1])
colnames(PPI_result)[3] <- "combined.score"
write.table(PPI_result, file = "PPI_input.csv", quote = F, row.names = F, sep = "\t")
