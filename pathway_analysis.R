

library(DESeq2)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(pathview)
library(readr)
library(dplyr)
library(png)
library(grid)
library(pathview)
library(enrichplot)



setwd("~/Desktop")



##DESeq2

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

# Pre-filtering: Remove low count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

##pick base comparison
##can change depending on group of interest

dds$Condition <- relevel(dds$Condition, ref = "Parental_ALDH_plus")
dds <- DESeq(dds)
results <- results(dds)


##separate results into comparisons
res_nonconv <- results(dds, contrast = c("Condition", "Non_converted", "Parental_ALDH_plus"))
res_conv<- results(dds, contrast = c("Condition", "Converted", "Parental_ALDH_plus"))
res_minus<- results(dds, contrast = c("Condition", "ParentalALDH_minus", "Parental_ALDH_plus"))
result_list=list(as.data.frame(res_nonconv), as.data.frame(res_conv), as.data.frame(res_minus))






##save each result as data frame 
write.csv(distinct(as.data.frame(res_nonconv)), "res_nonconv.csv",row.names=TRUE)
write.csv(distinct(as.data.frame(res_conv)), "res_conv.csv",row.names=TRUE)
write.csv(distinct(as.data.frame(res_minus)), "res_minus.csv",row.names=TRUE)




##CAN start here if already have DESeq2 results files

##can also open instead

##import each DESeq2 results file for comparisons
res_nonconv<-read_csv("res_nonconv.csv")
res_conv<-read_csv("res_conv.csv")
res_minus<-read_csv("res_minus.csv")



result_names <- c("res_nonconv", "res_conv", "res_minus")

#storing results


go_results <- list()
plot_list <- list()    # store plots (optional)
gene_lists <-list()
for (i in seq_along(result_list)) {
  df <- result_list[[i]]
  name <- result_names[i]
  original_gene_list <- df$log2FoldChange
  
  # name the vector
  names(original_gene_list) <- rownames(df)
  
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  
  # sort the list in decreasing order
  gene_list = sort(gene_list, decreasing = TRUE)
  
  # use significance threshold
  sig_genes_df = subset(df, padj < 0.05)
  
  #filter on log2fold change
  genes <- sig_genes_df$log2FoldChange
  
  # Name the vector
  names(genes) <- rownames(sig_genes_df)
  
  # omit NA values
  genes <- na.omit(genes)
  
  # log2FoldChange > 2
  genes <- names(genes)[abs(genes) > 2]
  gse <- gseGO(geneList=gene_list, 
               ont ="ALL", 
               keyType = "ENSEMBL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Hs.eg.db, 
               pAdjustMethod = "none")
  
  gene_lists[[name]]<-gene_list
  go_results[[name]] <- gse
  require(DOSE)
  p<-dotplot(gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
  plot_list[[name]] <- p
  
}


#choose gene list of interest
#can replace with any of above groups
plot_list[["res_conv"]]

ego <- go_results[["res_conv"]]

ego <- pairwise_termsim(ego)

emapplot(ego)

cnetplot(go_results[["res_conv"]], foldChange=gene_lists[["res_conv"]])


# named vector
fc <- gene_lists[["res_conv"]]

# map ENSEMBL to SYMBOL
gene_map <- bitr(names(fc),
                 fromType = "ENSEMBL",
                 toType = "SYMBOL",
                 OrgDb = org.Hs.eg.db)

# match names
fc_mapped <- fc[gene_map$ENSEMBL]
names(fc_mapped) <- gene_map$SYMBOL





# KEGG plot with pathview

dme <- pathview(
  gene.data = fc_mapped,
  pathway.id = "hsa05200",  ##chose pathway id from KEGG
  species = "hsa",  ##species
  gene.idtype = "SYMBOL",
  out.suffix="conv"  ##outfile suffix
)

##outsuffix changes name from pathview to conv
knitr::include_graphics("hsa05200.conv.png")

img <- readPNG("hsa05200.conv.png")
grid::grid.raster(img)
