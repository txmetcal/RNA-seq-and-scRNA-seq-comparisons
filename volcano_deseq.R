library(DESeq2)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library("EnhancedVolcano")
library("apeglm")
library(org.Hs.eg.db)
library(dplyr)



#setwd("~/Desktop")

# Load count matrix
counts <- read.csv("countstable.csv", row.names = 1)

# Load metadata
metadata <- read.csv("metadata.csv", row.names = 1)
metadata <- subset(metadata, select = c(Condition))


# Ensure row names of metadata match column names of counts
metadata <- metadata[match(colnames(counts), rownames(metadata)), ]
metadata <- data.frame(Condition = metadata, stringsAsFactors = FALSE)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ Condition
)

# remove low count genes
dds <- dds[rowSums(counts(dds)) > 10, ]



#set base condition
dds$Condition <- relevel(dds$Condition, ref = "Parental_ALDH_plus")
dds <- DESeq(dds)

#head(dds)



#create nonconverted v plus results
res_nonconv <- results(dds, contrast = c("Condition", "Non_converted", "Parental_ALDH_plus"))

##can do with other groups for example:
#res_conv <- results(dds, contrast = c("Condition", "Converted", "Parental_ALDH_plus"))

results=results(dds)

#all results
write_delim(results, "results_RNA_seq.csv")
head(results)



#map to symbol
symbol<- mapIds(org.Hs.eg.db,
                          keys = rownames(res_nonconv),
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")

#add new column
res_nonconv$gene <- symbol



##shrink LFC
resLFC <- lfcShrink(
  dds,
  coef = "Condition_Non_converted_vs_Parental_ALDH_plus",
  type = "apeglm"
)


##for nonconverted can repeat with each list
# The main function is named after the package
EnhancedVolcano(toptable = resLFC,              # We use the shrunken log2 fold change to reduce noise
                x = "log2FoldChange",           # contains the log2 fold changes
                y = "padj",                     # adjusted p-value
                lab = res_nonconv$gene
)
