# Load required libraries
library(DESeq2)
library(readxl)

# Read counts data from Excel file
counts_data <- read_excel("C:/Users/lenovo/Downloads/merged.xlsx")

# Extract gene IDs and save them to a separate CSV file
gene_ids <- counts_data$Geneid
write.csv(gene_ids, file = "C:/Users/lenovo/Downloads/output.csv", row.names = FALSE)

# Convert counts_data to a data frame
counts_data <- as.data.frame(counts_data)

# Remove the first column assuming it contains gene IDs
counts_data <- counts_data[, -1]  

# Read metadata for DESeq
colData <- read.csv("C:/Users/lenovo/Downloads/metadata_for_DESeq.csv")
rownames(colData) <- colnames(counts_data)

# Check if column names match between counts_data and colData
all(colnames(counts_data) %in% rownames(colData))
all(colnames(counts_data) == rownames(colData))

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ condition)

# Estimate size factors
dds = estimateSizeFactors(dds)

sizeFactors(dds)
head(counts(dds, normalized = TRUE))

# Write a table of normalized counts
write.csv(counts(dds, normalized = TRUE), 
          file="C:/Users/lenovo/Downloads/DESeq_normalizedcounts.csv")

# Inspect sample relationships using PCA
ddsB = estimateDispersions(dds)
vsd = varianceStabilizingTransformation(ddsB)
p = plotPCA(vsd,intgroup=c("condition","LibraryLayout"))
p

# Estimate dispersion
dds = estimateDispersions(dds)
plotDispEsts(dds)

# Perform DESeq analysis
dds <- DESeq(dds)

# Perform Wald test for differential expression between conditions

res1 <- results(dds, contrast = c("condition", "noDA24", "noDA26"))
res2<- results(dds, contrast = c("condition", "noDA24", "noDA34"))
res3<- results(dds, contrast = c("condition", "noDA24", "noDA60"))
res4<- results(dds, contrast = c("condition", "noDA26", "noDA34"))
res5<- results(dds, contrast = c("condition", "noDA26", "noDA60"))
res6<- results(dds, contrast = c("condition", "noDA34", "noDA60"))
res7<- results(dds, contrast = c("condition", "DA26", "DA34"))
res8<- results(dds, contrast = c("condition", "noDA24", "DA26"))
res9<- results(dds, contrast = c("condition", "noDA24", "DA34"))
res10<- results(dds, contrast = c("condition", "DA26", "noDA26"))
res11<- results(dds, contrast = c("condition", "DA34", "noDA34"))
res12<- results(dds, contrast = c("condition", "DA34", "noDA60"))

# Write results to CSV files
write.csv(res1Sig, 
          file="C:/Users/lenovo/Downloads/res1Sig.csv")
write.csv(res2Sig, 
          file="C:/Users/lenovo/Downloads/res2Sig.csv")
write.csv(res3Sig, 
          file="C:/Users/lenovo/Downloads/res3Sig.csv")
write.csv(res4Sig, 
          file="C:/Users/lenovo/Downloads/res4Sig.csv")
write.csv(res5Sig, 
          file="C:/Users/lenovo/Downloads/res5Sig.csv")
write.csv(res6Sig, 
          file="C:/Users/lenovo/Downloads/res6Sig.csv")
write.csv(res7Sig, 
          file="C:/Users/lenovo/Downloads/res7Sig.csv")
write.csv(res8Sig, 
          file="C:/Users/lenovo/Downloads/res8Sig.csv")
write.csv(res9Sig, 
          file="C:/Users/lenovo/Downloads/res9Sig.csv")
write.csv(res10Sig, 
          file="C:/Users/lenovo/Downloads/res10Sig.csv")
write.csv(res12Sig, 
          file="C:/Users/lenovo/Downloads/res11Sig.csv")
write.csv(res12Sig, 
          file="C:/Users/lenovo/Downloads/res12Sig.csv")


# Plot differential expression against expression amount
plotMA(res1, ylim = c(-2, 2))
points(res$log2FoldChange[res1$padj > 0.01], res1$baseMean[res1$padj > 0.01], col = "blue")

# Threshold differential expression using padj < 0.01
res1Sig = res1[which(res1$padj < 0.01),]
res2Sig = res2[which(res2$padj < 0.01),]
res3Sig = res3[which(res3$padj < 0.01),]
res4Sig = res4[which(res4$padj < 0.01),]
res5Sig = res5[which(res5$padj < 0.01),]
res6Sig = res6[which(res6$padj < 0.01),]
res7Sig = res7[which(res7$padj < 0.01),]
res8Sig = res8[which(res8$padj < 0.01),]
res9Sig = res9[which(res9$padj < 0.01),]
res10Sig = res10[which(res10$padj < 0.01),]
res11Sig = res11[which(res11$padj < 0.01),]
res12Sig = res12[which(res12$padj < 0.01),]

# Tabulate the number of differentially expressed genes
table(res1$padj < 0.01)
table(res2$padj < 0.01)
table(res3$padj < 0.01)
table(res4$padj < 0.01)
table(res5$padj < 0.01)
table(res6$padj < 0.01)
table(res7$padj < 0.01)
table(res8$padj < 0.01)
table(res9$padj < 0.01)
table(res10$padj < 0.01)
table(res11$padj < 0.01)
table(res12$padj < 0.01)

# Generate histograms of p-values for each comparison
res1$pval <- as.numeric(as.character(res1$pval))

hist(res1$pvalue,breaks=100)
hist(res2$pval,breaks=100)
hist(res3$pval,breaks=100)
hist(res4$pval,breaks=100)
hist(res5$pval,breaks=100)
hist(res6$pval,breaks=100)
hist(res7$pval,breaks=100)
hist(res8$pval,breaks=100)
hist(res9$pval,breaks=100)
hist(res10$pval,breaks=100)
hist(res11$pval,breaks=100)
hist(res12$pval,breaks=100)

# Combine results with gene IDs and save to CSV
results_with_gene_ids <- data.frame(GeneID = gene_ids, Results = res1)
write.csv(results_with_gene_ids, file = "C:/Users/lenovo/Downloads/res1withgenes.csv", row.names = FALSE)

