

library(dplyr)
library(readr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(msigdbr)
library(fgsea)
library(msigdbr)
library(dplyr)


##setwd
#setwd("~/Desktop")


##read in comparison files from DESeq2
res_nonconv <- read.csv("res_nonconv.csv", row.names = 1)
res_conv <-read.csv("res_conv.csv", row.names=1)
res_minus <- read.csv("res_minus.csv", row.names=1)
res_naive <- read.csv("mmc3_naive_v_NACT.csv", skip=1)

# Put datasets into a named list
res_list <- list(
  nonconv = res_nonconv,
  conv = res_conv,
  minus = res_minus 
)


# Get each collection separately
#halls <- msigdbr(species = "Homo sapiens", category = "H")
#c2    <- msigdbr(species = "Homo sapiens", category = "C2")
gene_sets_df<- msigdbr(species = "Homo sapiens", category = "H")

# Combine into one data frame
#gene_sets_df <- bind_rows(c5)

# names = gene set, values = gene symbols
gene_sets <- split(gene_sets_df$gene_symbol, gene_sets_df$gs_name)

# Store results
gsea_results_list <- list()

for (name in names(res_list)) {
  
  res <- res_list[[name]]
  
  # Map symbols
  symbol <- mapIds(org.Hs.eg.db,
                   keys = rownames(res),
                   column = "SYMBOL",
                   keytype = "ENSEMBL",
                   multiVals = "first")
  
  res$symbol <- symbol
  
  # Remove NA/duplicates
  gseaDat <- res %>%
    filter(!is.na(symbol), !duplicated(symbol))
  
  # ranking
  rankings <- sign(gseaDat$log2FoldChange) * (-log10(gseaDat$pvalue))
  names(rankings) <- gseaDat$symbol
  
  # Clean rankings
  max_ranking <- max(rankings[is.finite(rankings)])
  min_ranking <- min(rankings[is.finite(rankings)])
  
  rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
  rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
  rankings <- sort(rankings, decreasing = TRUE)
  
  # GSEA
  GSEAres <- fgsea(
    pathways = gene_sets,
    stats = rankings,
    scoreType = "std",
    minSize = 10,
    maxSize = 500,
    nproc = 1
  )
  
  # Get top pathways
  top_up <- GSEAres[ES > 0][head(order(pval), 10), pathway]
  top_down <- GSEAres[ES < 0][head(order(pval), 10), pathway]
  
  top_pathways <- c(top_up, rev(top_down))
  
  # Store results
  gsea_results_list[[name]] <- list(
    GSEAres = GSEAres,
    rankings = rankings,
    top_pathways = top_pathways
  )
  
  
}

##choose specific list for paired comparison gsea_results_list[["nonconv"]]

all_pathways <- unique(unlist(lapply(gsea_results_list, function(x) x$top_pathways)))


##using NES
get_NES <- function(gseaRes, pathways) {
  vals <- gseaRes$NES[match(pathways, gseaRes$pathway)]
  names(vals) <- pathways
  return(vals)
}

##make into matrix 
mat <- sapply(names(gsea_results_list), function(name) {
  get_NES(gsea_results_list[[name]]$GSEAres, all_pathways)
})

mat <- as.matrix(mat)
# Set row names
#rownames(mat) <- all_pathways


# Replace NA (missing pathways) with 0
mat[is.na(mat)] <- 0

# Order by strongest enrichment signal
mat <- mat[order(apply(abs(mat), 1, max), decreasing = TRUE), , drop = FALSE ]

# keep top 30 pathways
#mat_top <- mat[1:min(19, nrow(mat)), ]
#mat <- matrix(mat, ncol = 1)


mat_top<-mat[1:min(25, nrow(mat)), ]

# Pathways in columns 
conv_paths <- rownames(mat_top)[!is.na(mat_top[,"conv"])]
nonconv_paths <- rownames(mat_top)[!is.na(mat_top[,"nonconv"])]
minus_paths <- rownames(mat_top)[!is.na(mat_top[,"minus"])]


# Only keep pathways present in chosen comparision groups
common_paths <- intersect(nonconv_paths, conv_paths)
common_paths <- intersect(common_paths,minus_paths)

##can keep only specific group here if needed
mat_top_sub <- mat_top[common_paths,]   # keeps rows as matrix

# convert to 1-column matrix

library(pheatmap)

pheatmap(mat_top_sub,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = FALSE,  # keep your 3 conditions in order
         fontsize_row = 8,
         angle_col = 45,
         main = "Top Shared Pathways (NES)")







#--------------------------------------------------
#IF YOU WANT TO COMPARE TO NAIVE OR OTHER DATASETS
#----------------------------------



##use pathways from above and then work with new dataset to find shared pathways

# Remove NA + duplicates
res_naive <- subset(res_naive, celltype=="Epithelial_cells")
gseaDat <- res_naive %>%
  filter(!is.na(gene), !duplicated(gene))

# Ranking
rankings <- sign(gseaDat$avg_log2FC) * (-log10(gseaDat$p_val_adj))
names(rankings) <- gseaDat$gene

# Clean rankings
max_ranking <- max(rankings[is.finite(rankings)])
min_ranking <- min(rankings[is.finite(rankings)])

rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
rankings <- sort(rankings, decreasing = TRUE)

# Run GSEA
GSEAres <- fgsea(
  pathways = gene_sets,
  stats = rankings,
  scoreType = "std",
  minSize = 10,
  maxSize = 500,
  nproc = 1
)

# Get top pathways
top_up <- GSEAres[ES > 0][head(order(pval), 10), pathway]
top_down <- GSEAres[ES < 0][head(order(pval), 10), pathway]

top_pathways <- c(top_up, rev(top_down))

# Store results
gsea_results_list_naive <- list(
  GSEAres = GSEAres,
  rankings = rankings,
  top_pathways = top_pathways
)

all_pathways_naive <- unique(gsea_results_list_naive$top_pathways)


##naive NES
get_NES_naive <- function(GSEAres, pathways) {
  vals <- GSEAres$NES[match(pathways, GSEAres$pathway)]
  names(vals) <- pathways
  return(vals)
}

##same steps as above, make into own matrix
mat_naive <- get_NES_naive(GSEAres, all_pathways)

mat_naive <- as.matrix(mat_naive)

mat_naive <- mat_naive[order(abs(mat_naive[,1]), decreasing = TRUE), , drop = FALSE]

write.csv(mat, "pathways_naive.csv")



# Pathways in mat_naive
naive_paths <- rownames(mat_naive)


##here using NONCONV, replace NONCONV with set of interest

# Only keep pathways present in both
common_paths <- intersect(nonconv_paths, naive_paths)

##nonconv or conv or minus 
mat_top_sub <- mat_top[common_paths,"nonconv"]   # keeps rows as matrix
mat_naive_sub <- mat_naive[common_paths,]               # subset vector
mat_naive_sub <- as.matrix(mat_naive_sub)                # convert to 1-column matrix

# convert to 1-column matrix
new_mat <- cbind(new_mat <- cbind(mat_top_sub, mat_naive_sub))


##headings 
colnames(new_mat)<-c("Non-converted_v_ALDH+","Naive_v_NACT")
library(pheatmap)

pheatmap(new_mat,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = FALSE,  # keep your 3 conditions in order
         fontsize_row = 8,
         angle_col = 45,
         main = "Top Shared Pathways (NES)")




