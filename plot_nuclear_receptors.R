#!/usr/bin/env Rscript
# Reproduce nuclear receptor dot plot from GSE293459 processed Seurat object.
# Usage: Rscript plot_nuclear_receptors.R path/to/GSE293459_patient_seurat.rds [output_dir]

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || !nzchar(args[[1]])) {
  stop('Usage: Rscript plot_nuclear_receptors.R <input.rds> [output_dir]', call. = FALSE)
}
input_file <- args[[1]]
output_dir <- if (length(args) >= 2L) args[[2]] else 'results'
if (!file.exists(input_file)) stop('Input file not found: ', input_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_file)
if (!inherits(obj, 'Seurat')) stop('Input RDS must contain a Seurat object.')
if (!'customtype' %in% colnames(obj[[]])) stop('Missing metadata column: customtype')
if (!'RNA' %in% Assays(obj)) stop('Missing RNA assay.')
if (!'seurat_clusters' %in% colnames(obj[[]])) stop('Missing original seurat_clusters metadata.')

epi <- subset(obj, cells = colnames(obj)[!is.na(obj$customtype) & obj$customtype == 'Tumor'])
if (ncol(epi) == 0L) stop('No cells with customtype == "Tumor".')
DefaultAssay(epi) <- 'RNA'

# Requires pre-existing normalized expression; do not change the original data.
if (!'data' %in% Layers(epi[['RNA']])) {
  stop('RNA assay has no normalized data layer. Normalize appropriately before running this script.')
}

features <- c('RARA', 'RXRA', 'VDR')
missing <- setdiff(features, rownames(epi[['RNA']]))
if (length(missing)) stop('Missing RNA genes: ', paste(missing, collapse = ', '))

# Preserve supplied clusters, sorting numeric labels when possible.
clusters <- as.character(epi$seurat_clusters)
cluster_levels <- unique(clusters)
if (all(grepl('^[0-9]+$', cluster_levels))) {
  cluster_levels <- as.character(sort(as.integer(cluster_levels)))
} else {
  cluster_levels <- sort(cluster_levels)
}
epi$original_cluster_for_plot <- factor(clusters, levels = cluster_levels)

# DotPlot color is average expression scaled per gene across groups;
# dot area corresponds to percentage of cells with detectable expression.
p <- DotPlot(
  object = epi,
  features = features,
  group.by = 'original_cluster_for_plot',
  assay = 'RNA',
  dot.scale = 8,
  scale = TRUE,
  col.min = -2.5,
  col.max = 2.5
) +
  scale_color_gradient2(low = 'lightgrey', mid = 'lightblue', high = 'darkblue', midpoint = 0) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12, face = 'italic', angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = 'bold', hjust = 0.5)
  ) +
  labs(
    title = 'Nuclear Receptor Expression',
    x = 'Nuclear receptor',
    y = 'Seurat cluster',
    color = 'Average scaled\nexpression',
    size = 'Percent expressed'
  )

print(p)
ggsave(file.path(output_dir, 'RARA_RXRA_VDR_DotPlot.pdf'), p, width = 5.5, height = 7, units = 'in')
ggsave(file.path(output_dir, 'RARA_RXRA_VDR_DotPlot.png'), p, width = 5.5, height = 7, units = 'in', dpi = 600)
write.csv(p$data, file.path(output_dir, 'RARA_RXRA_VDR_DotPlot_Data.csv'), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, 'sessionInfo.txt'))
message('Done. Outputs saved to: ', normalizePath(output_dir))
