library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(tibble)
library(clusterProfiler)
library(pheatmap)

sx2 <- readRDS('/storage/riosugimuralab/alexto/B1_objects/B1_Sx/0102_three_sample_macrophage.rds')

hcc <- readRDS('HCCF1F4.RData')
hcc<-UpdateSeuratObject(hcc)

Idents(hcc) <- "CellType"
sharma_myeloid <- subset(hcc, idents = c("Myeloid"))

sx_carm <- subset(sx2, orig.ident == 'CAR_M')

Idents(sharma_myeloid) <- 'condition'
sharma_myeloidsub <- subset(sharma_myeloid, idents = c('Adj Normal','Tumor'))

sx_carm$condition <- 'CAR_M'

combined <- merge(sharma_myeloidsub, sx_carm)

colnames(hcc@meta.data)
unique(hcc$condition)

Idents(hcc) <- 'condition'
sharma_hcc <- subset(hcc, idents = c("Tumor"))

combined$condition <- factor(
  combined$condition,
  levels = c("CAR_M","Tumor","Adj Normal" )
)

sharma_myeloidnormal <- subset(sharma_myeloid, idents = c('Adj Normal'))
combined_only <- merge(sharma_myeloidnormal, sx_carm)


# 1) Make RNA the default for plotting
DefaultAssay(combined) <- "RNA"  # avoids default-search warnings

# 2) Check layer/cell consistency
Layers(combined[["RNA"]])
length(Cells(combined))
length(Cells(combined[["RNA"]]))

# 3) Join layers (Seurat v5) so features/cells are consistently available
combined <- JoinLayers(combined, assay = "RNA")

# 4) Ensure a data layer exists (optional but typical for expression plots)
combined <- NormalizeData(combined, assay = "RNA")

library(patchwork)
options(repr.plot.width=8, repr.plot.height=9)

genes_req <- c(
  "IL18BP",'CD274','MARCO','LRIG2',
  'IFITM1','IFITM2','IFIT2','IFR1'
)

missing_genes <- setdiff(genes_req, rownames(combined[["RNA"]]))  # genes that will be removed [web:79]
genes_req <- intersect(genes_req, rownames(combined[["RNA"]]))
missing_genes


DefaultAssay(combined) <- "RNA"
p_list <- VlnPlot(combined, features = genes_req, group.by = "condition",
                  pt.size = 0, combine = FALSE)

wrap_plots(p_list) + plot_layout(guides = "collect") &
  theme(legend.position = "right")
