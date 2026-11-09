

#install.packages("ggplot2")
#install.packages("ggVennDiagram")
#install.packages("VennDiagram")
#install.packages('ActivePathways')

library(ActivePathways)
library(ggVennDiagram)
library(VennDiagram)
library(limma)
library(magrittr)
library(ggplot2)
library(readr)
library(dplyr)
library(tibble)




setwd("~/Desktop")


##read in metcalfe paper signficant genes
counts_conv <-read_delim("Converted_vs_Parental_ALDH_plus.txt")
counts_non <-read_delim("Non_Converted_vs_Parental_ALDH_plus.txt")

#counts_conv <- read_csv("converted_only_genes_full.csv")


##read in comparison paper significant genes
counts_comp <- read_csv("mmc3_naive_v_NACT.csv",skip=1) ##sensitive v resistant
counts_comp <-subset(counts_comp,celltype=="Epithelial_cells")


##use genes as name
counts_comp <- counts_comp[,-1]



# thresholds
adjp_cutoff <- 0.05
pval_cutoff <-0.01
logfc_cutoff <- 2



##renaming for use later 
counts_conv <- counts_conv %>%
  dplyr::rename(`pval` = `Converted_vs_Parental_ALDH_plus pvalue`)


counts_conv <- counts_conv %>%
  dplyr::rename(`adj.P.val` = `Converted_vs_Parental_ALDH_plus adjusted_pvalue_BH`)

counts_conv <- counts_conv %>%
  dplyr::rename(`logFC` = `Converted_vs_Parental_ALDH_plus_log2fc`)



##if needed can make numeric

#counts_comp$FDR <- as.numeric(as.character(counts_comp$FDR))
#counts_comp$Log2FC <- as.numeric(as.character(counts_comp$avr_log2FC))

# Get genes unique to sig_conv
# Find genes unique to sig_conv
#sig_non <- subset(counts_non,adj.P.val < adjp_cutoff &   abs(logFC) > logfc_cutoff & pval< pval_cutoff)
#sig_conv <- subset(counts_conv, adj.P.val < adjp_cutoff & abs(logFC) > logfc_cutoff & pval< pval_cutoff)

#conv_only_df <- sig_conv[!(sig_conv$Gene_Name %in% sig_non$Gene_Name), ]

# Write full dataframe (with all columns) to CSV
#write.csv(conv_only_df,
#          "converted_only_genes_full.csv",
#          row.names = FALSE)



# significant genes in each dataset
sig_comp <- subset(counts_comp, p_val_adj < adjp_cutoff &  abs(avg_log2FC) > logfc_cutoff & p_val< pval_cutoff)
sig_conv <- subset(counts_conv, adj.P.val < adjp_cutoff & abs(logFC) > logfc_cutoff & pval< pval_cutoff)

venn_list_all <- list(
  Comp = sig_comp$gene,
  Conv = sig_conv$Gene_Name
)


##all genes  
##change columns for groups of interest, they are pulled from list above

venn.plot <- draw.pairwise.venn(
  area1 = length(venn_list_all[[1]]),
  area2 = length(venn_list_all[[2]]),
  cross.area = length(intersect(venn_list_all[[1]], venn_list_all[[2]])),
  category = c("Naive v NACT", "Converted v ALDH+"),
  main.cex = 1.5,
  fill = c("#F46474", "#4981bf"),
  alpha = c(0.5, 0.5),
  lwd = 0.2,
  cex = 1.5,
  cat.cex = 1.2,
  scaled=FALSE,
  cat.dist = c(0.06, 0.06),
  cat.just = list(c(0.5, 0.5), c(0.5, 0.5)),
  cat.pos = c(0, 0)
  
)


#subsetting only up genes

up_comp <- subset(sig_comp, p_val_adj < adjp_cutoff & avg_log2FC > logfc_cutoff & p_val< pval_cutoff)
up_conv <- subset(sig_conv, adj.P.val < adjp_cutoff & logFC > logfc_cutoff & pval< pval_cutoff)



##make list of genes

venn_list_up <- list(
  Patient_Tumors = up_comp$gene,
  Converted_Cells= up_conv$Gene_Name
)


#venn_list_up[1]


##change columns for groups
##venn diagram of genes above
venn.plot <- draw.pairwise.venn(
  area1 = length(venn_list_up[[1]]),
  area2 = length(venn_list_up[[2]]),
  cross.area = length(intersect(venn_list_up[[1]], venn_list_up[[2]])),
  category = c("Naive v NACT", "Converted v ALDH+"),
  main.cex = 1.5,
  fill = c("#F46474", "#4981bf"),
  alpha = c(0.5, 0.5),
  lwd = 0.2,
  cex = 1.5,
  cat.cex = 1.2,
  scaled=FALSE,
  cat.dist = c(0.06, 0.06),
  cat.just = list(c(0.5, 0.5), c(0.5, 0.5)),
  cat.pos = c(0, 0)
  
)

grid.newpage()
grid.draw(venn.plot)

grid.text("Overlap of ** Genes", y = 0.95, gp = gpar(fontfamily = "Times New Roman", fontsize = 16))


 

##do the same for the down comparison
down_comp <- subset(sig_comp, p_val_adj < adjp_cutoff & avg_log2FC < -logfc_cutoff & p_val< pval_cutoff)
down_conv <- subset(sig_conv, adj.P.val < adjp_cutoff & logFC <  -logfc_cutoff & pval< pval_cutoff)


#list of down genes

venn_list_down <- list(
  Patient = down_comp$gene,
  Converted_Down = down_conv$Gene_Name
)


##make down venn diagram

venn.plot <- draw.pairwise.venn(
  area1 = length(venn_list_down[[1]]),
  area2 = length(venn_list_down[[2]]),
  cross.area = length(intersect(venn_list_down[[1]], venn_list_down[[2]])),
  category = c("Naive v NACT", "Converted v ALDH+"),
  main.cex = 1.5,
  fill = c("#F46474", "#4981bf"),
  alpha = c(0.5, 0.5),
  lwd = 0.2,
  cex = 1.5,
  cat.cex = 1.2,
  scaled=FALSE,
  cat.dist = c(0.06, 0.06),
  cat.just = list(c(0.5, 0.5), c(0.5, 0.5)),
  cat.pos = c(0, 0)
  
)


# Use intersect() to see genes in both
genes_down_both <- intersect(down_comp$gene, down_conv$Gene_Name)
genes_up_both <- intersect(up_comp$gene,up_conv$Gene_Name)

#head(genes_up_both)

##save files as comparison 

write.csv(genes_up_both, "converted_shared_up_maia_NaiveNACT_pval.csv", row.names = TRUE) ##change for whatever file you wish to save
write.csv(genes_down_both, "converted_shared_down_maia_NaiveNACT_pval.csv", row.names = TRUE) ##change for whatever file you wish to save

