library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(tibble)
library(clusterProfiler)
library(pheatmap)

sx <- readRDS('/storage/riosugimuralab/alexto/B1_objects/B1_Sx/0102_three_sample_macrophage.rds')

# visualize UMAP
options(repr.plot.width = 10, repr.plot.height = 3)
DimPlot(sx2, reduction = 'umap',group.by='layer_hao', pt.size=0.05) | 
DimPlot(sx2, group.by='orig.ident', pt.size=0.05)

# Visualize Heatmap
library(Seurat)

# choose the assay you used for plotting
assay_use <- DefaultAssay(sx2)

# pick the grouping used for the heatmap rows (change to your metadata column)
group_var <- "layer_hao"   # or "celltype", "anno", etc.

# genes you intend to plot
genes_use <- unique(unlist(gene_groups))

avg <- AverageExpression(
  sx2,
  assays   = assay_use,
  features = genes_use,
  group.by = group_var,
  slot     = "data"
)[[assay_use]]  # genes x groups

# scale per gene (row-wise), then transpose to match your code’s expectation
mat_scaled_t <- t(scale(t(avg)))  # groups x genes

library(Seurat)
library(pheatmap)
library(RColorBrewer)

## ---- User choices ----
assay_use <- DefaultAssay(sx2)

# Rows of heatmap (groups)
group_by <- "ident"  # or e.g. "celltype", "seurat_clusters", "orig.ident"

## ---- 1) Category order + gene groups ----
cat_order <- c("Mono_derived","Chemotaxis","Phagocytic_Mac",
               "MHC_I","MHC_II","Ag_processing","Co_stim")

gene_groups <- list(
  Mono_derived   = c("FCGR3A","NFKB1","CSF1R"),
  Chemotaxis     = c("CSF1","CCL18","CXCL9"),
  Phagocytic_Mac = c("C1QA","C1QB","C1QC"),
  MHC_I          = c("HLA-A","HLA-B","HLA-C"),
  MHC_II         = c("HLA-DRA","HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DMA","HLA-DMB","CD74"),
  Ag_processing  = c("TAP1","TAP2","TAPBP","PSMB9","B2M"),
  Co_stim        = c("CD86","CD40")
)
genes_use <- unique(unlist(gene_groups))

## ---- 2) Average expression from sx2 ----
avg_list <- tryCatch(
  AverageExpression(
    sx2, assays = assay_use, features = genes_use,
    group.by = group_by, layer = "data"
  ),
  error = function(e) {
    AverageExpression(
      sx2, assays = assay_use, features = genes_use,
      group.by = group_by, slot = "data"
    )
  }
)
avg <- avg_list[[assay_use]]  # genes x groups

# Optional: enforce a stable group order (prevents “mystery reordering”)
if (identical(group_by, "ident")) {
  lev <- levels(Idents(sx2))
  lev <- lev[lev %in% colnames(avg)]
  if (length(lev)) avg <- avg[, lev, drop = FALSE]
} else {
  v <- sx2[[group_by]][, 1]
  if (is.factor(v)) {
    lev <- levels(v)
    lev <- lev[lev %in% colnames(avg)]
    if (length(lev)) avg <- avg[, lev, drop = FALSE]
  }
}

## ---- 3) Scale per gene; matrix as groups x genes ----
avg_scaled <- t(scale(t(avg)))  # genes x groups
mat_scaled_t <- t(avg_scaled)   # groups x genes

## ---- 4) Gene order in x-axis blocks ----
present <- colnames(mat_scaled_t)
gene_order <- unlist(lapply(cat_order, function(ct) intersect(gene_groups[[ct]], present)))
gene_order <- c(gene_order, setdiff(present, gene_order))

mat_use <- mat_scaled_t[, gene_order, drop = FALSE]

## ---- 5) Reverse y-axis (rows) ONLY ----
mat_use <- mat_use[rev(rownames(mat_use)), , drop = FALSE]  # reverse y-axis [web:40][web:53]

## ---- 6) Build Category annotation for TOP + legend order ----
# Gene -> category mapping
gene_cat <- setNames(rep(NA_character_, length(present)), present)
for (ct in names(gene_groups)) gene_cat[intersect(gene_groups[[ct]], present)] <- ct
gene_cat[is.na(gene_cat)] <- "Unassigned"

# Factor levels control legend order; keep Unassigned last [web:48][web:57]
cat_levels <- c(cat_order, "Unassigned")
annotation_col <- data.frame(
  Category = factor(gene_cat[gene_order], ordered = TRUE),
  row.names = gene_order
)

stopifnot(identical(colnames(mat_use), rownames(annotation_col)))  # must match [web:57]

# Define colors in the SAME order as cat_levels to keep legend consistent [web:48]
pal <- colorRampPalette(brewer.pal(8, "Set2"))(length(cat_levels))
ann_colors <- list(Category = setNames(pal, cat_levels))

## ---- 7) Plot ----
pheatmap(mat_use,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
         border_color = NA,
         fontsize_row = 15,
         fontsize_col = 15,
         angle_col = 90,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         annotation_names_col = TRUE,
         drop_levels = FALSE)
