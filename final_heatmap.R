
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(org.Hs.eg.db)
library(AnnotationDbi)


##set to working directory
#setwd("~/Desktop")


# Counts + metadata (RAW integer counts for DESeq2)


##run DESeq2
counts <- read.csv("countstable.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("metadata.csv", row.names = 1, check.names = FALSE)
metadata <- metadata[, "Condition", drop = FALSE]

stopifnot(all(colnames(counts) %in% rownames(metadata)))
metadata <- metadata[match(colnames(counts), rownames(metadata)), , drop = FALSE]
stopifnot(identical(rownames(metadata), colnames(counts)))

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ Condition
)

keep <- rowSums(counts(dds)) > 10
dds <- dds[keep, ]

dds$Condition <- relevel(dds$Condition, ref = "Parental_ALDH_plus")
dds <- DESeq(dds)

res_default <- as.data.frame(results(dds))
res_default$ENSG <- rownames(res_default)




# Normalized counts / heatmap annotation

norm_counts <- read.csv("GSF3454_normalized_countstable.csv", row.names = 1, check.names = FALSE)

# drop extra column if present
if (ncol(norm_counts) >= 13) {
  norm_counts <- norm_counts[, -13, drop = FALSE]
}

# column annotation must match count matrix used in heatmap
annot_info <- as.data.frame(colData(dds)[, "Condition", drop = FALSE])
rownames(annot_info) <- colnames(dds)


##for all of them 
rld <- rlog(dds, blind = FALSE)
vst_mat <- assay(vst(dds, blind = FALSE))

#Gene symbols (keep ENSG as rownames for matching dds/rld)
sym <- mapIds(
  org.Hs.eg.db,
  keys      = rownames(res_default),
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)
res_default$symbol <- ifelse(is.na(sym), rownames(res_default), sym)

# Filter significant results table
res_sig <- res_default[!is.na(res_default$padj), ]
sig_genes <- res_sig[(res_sig$padj < 0.05) & (abs(res_sig$log2FoldChange) >= 1), ]
dim(sig_genes)




# Z score heatmap

library(pheatmap)

mat <- assay(rld)   # or assay(vst(dds, blind = FALSE))

# z-scores 
Z <- t(scale(t(mat)))
Z[is.nan(Z)] <- 0   # genes with zero variance across samples

light_blue <- adjustcolor("blue", alpha.f = 0.4)
light_red  <- adjustcolor("red",  alpha.f = 0.4)
light_green<- adjustcolor("green",  alpha.f = 0.4)
light_purple<- adjustcolor("purple",  alpha.f = 0.4)
ann_colors <- list(
  Condition = c(
    Parental_ALDH_plus = "blue",
    Converted = "red",
    Non_converted= "green",
    ParentalALDH_minus= "purple"
  )
)

pheatmap(
  Z,
  cluster_rows    = TRUE,
  cluster_cols    = TRUE, ##FALSE if want samples grouped together 
  show_rownames   = FALSE,   # TRUE if few genes
  annotation_col  = annot_info,
  annotation_colors = ann_colors,
  main            = "Row z-scored expression (rlog)"
)
