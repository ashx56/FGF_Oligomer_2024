suppressPackageStartupMessages({ 
  library(Matrix)
  library(Seurat)
  library(tidyverse)
  library(monocle3)
  library(RColorBrewer)
  library(viridis)
  
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  # Set a seed to make umap and other non-deterministic steps consistent
  set.seed(seed = 42)
})

# Processed data files can be downloaded from Harvard Dataverse doi:10.7910/DVN/80RTTV
# Load raw counts, cell/gene metadata
counts <- readMM("FGF_Oligomer_scRNAseq_1_filtered_counts.mtx")
genes <- read.csv("FGF_Oligomer_scRNAseq_1_genes.tsv", header = F)
barcodes <- read.csv("FGF_Oligomer_scRNAseq_1_barcodes.tsv", header = F)

rownames(counts) <- genes[, 1]
colnames(counts) <- barcodes[, 1]

mdata <- read.csv("FGF_Oligomer_scRNAseq_1_metadata.csv", row.names = 1)

# Create Seurat object
seurat_FGF <- CreateSeuratObject(counts = counts,
                                 meta.data = mdata)

seurat_FGF

# QC
# Mitochondrial gene detection
# "MT-" genes sequenced from mitochondrial genome 
mito_genes <- rownames(seurat_FGF)[grep("^MT-", rownames(seurat_FGF))]
mito_genes

# Percentage of reads per cell mapping to mitochondrial genes
# <10% usually considered a good threshold
seurat_FGF <- PercentageFeatureSet(seurat_FGF, 
                                   "^MT-", 
                                   col.name = "percent_mito")
head(seurat_FGF@meta.data)

# Ribosomal gene detection
# "RP.." genes mapping to ribosomal proteins
ribo_genes <- rownames(seurat_FGF)[grep("^RP[SL]", rownames(seurat_FGF))]
head(ribo_genes)

# Percentage of reads per cell mapping to ribosomal genes
# <1% usually considered a good threshold
seurat_FGF <- PercentageFeatureSet(seurat_FGF, 
                                   "^RP[SL]", 
                                   col.name = "percent_ribo")
head(seurat_FGF@meta.data)

# Red blood cell contamination
# Red blood cell genes annotated "HB.." indicating RBC contamination
hb_genes <- rownames(seurat_FGF)[grep("^HB[^(P)]", rownames(seurat_FGF))]
hb_genes

# Percentage of reads per cell mapping to RBC genes
# <0.5% usually considered a good threshold
seurat_FGF <- PercentageFeatureSet(seurat_FGF, 
                                   "^HB[^(P)]", 
                                   col.name = "percent_hb")
head(seurat_FGF@meta.data)

# Plot QC metrics
Idents(seurat_FGF) <- "Day"
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
VlnPlot(seurat_FGF, 
                features = feats, 
                pt.size = 0, 
                ncol = 3) +
  NoLegend()

# Subset Seurat object to include cells within 95% percentile of QC metrics
vars_to_filter <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")

seurat_FGF[['pass_filter']] <- TRUE
for (var in vars_to_filter) {
  lower_quantile <- quantile(seurat_FGF@meta.data[[var]], 0.05, na.rm = TRUE)
  upper_quantile <- quantile(seurat_FGF@meta.data[[var]], 0.95, na.rm = TRUE)
  seurat_FGF[['pass_filter']] <- seurat_FGF[['pass_filter']] & 
    (seurat_FGF@meta.data[[var]] >= lower_quantile &
       seurat_FGF@meta.data[[var]] <= upper_quantile)
}
seurat_FGF <- subset(seurat_FGF, subset = pass_filter)
seurat_FGF

# Expression-based filtering
# Remove cells expressing <200 genes
selected_c <- WhichCells(seurat_FGF, expression = nFeature_RNA > 200)
# Remove genes expressed in <3 cells
selected_f <- rownames(seurat_FGF)[Matrix::rowSums(seurat_FGF@assays$RNA$counts) > 3]

# Filter genes to include protein coding only
df_gene <- readRDS("helper_files/df_gene.RDS")
genes_coding <- df_gene %>%
  dplyr::filter(gene_type == "protein_coding") %>%
  dplyr::select(gene_short_name)

seurat_FGF <- subset(seurat_FGF, 
                     features = intersect(selected_f, genes_coding$gene_short_name), 
                     cells = selected_c)

seurat_FGF

# Correlation between recovered UMIs per cell and detected genes
FeatureScatter(seurat_FGF, 
                       "nCount_RNA", 
                       "nFeature_RNA", 
                       pt.size = 0.5) +
  ggplot2::theme_classic()

# Cell-cycle scoring
seurat_FGF = NormalizeData(seurat_FGF, verbose = T)

# Percentage of counts mapping to G2M/S phase genes
seurat_FGF <- CellCycleScoring(object = seurat_FGF, 
                               g2m.features = cc.genes$g2m.genes,
                               s.features = cc.genes$s.genes)

VlnPlot(seurat_FGF, 
        features = c("S.Score", "G2M.Score"),
        ncol = 2, pt.size = 0)

# Convert to monocle3 CDS for further processing
# Extract gene metadata and format for monocle3 cds
gene_mdata <- as.data.frame(rownames(seurat_FGF))
colnames(gene_mdata) <- "gene_short_name"
rownames(gene_mdata) <- rownames(seurat_FGF)

# Create new monocle3 cds
cds_FGF <- new_cell_data_set(expression_data = seurat_FGF[["RNA"]]$counts,
                             cell_metadata = seurat_FGF@meta.data,
                             gene_metadata = gene_mdata)
cds_FGF

# Data processing
# Size factor estimation and normalization/scaling
cds_FGF <- estimate_size_factors(cds_FGF) 
cds_FGF <- preprocess_cds(cds_FGF, verbose = T)

# Batch correction
# Counts are adjusted for contaminant and cell cycle genes
cds_FGF <- align_cds(cds_FGF, residual_model_formula_str = "~percent_mito +
                     percent_ribo +
                     percent_hb +
                     S.Score +
                     G2M.Score", verbose = T)

# Dimensionality reduction and clustering
cds_FGF <- reduce_dimension(cds_FGF, umap.min_dist = 0.25, verbose = T)
cds_FGF <- cluster_cells(cds_FGF, k = 20, verbose = T)

# Visualize UMAP by day of harvest (Fig 5A)
plot_cells(cds_FGF, 
                   color_cells_by = "Day",
                   label_cell_groups = F,
                   cell_size = 0.5) +
  scale_color_viridis_d()

# Visualize expression of a few canonical cell type markers to determine cell identities
# iPSCs, Endothelial, Perivascular
plot_cells(cds_FGF,
           genes = c("SOX2", "PECAM1", "PDGFRB", "TTN"),
           cell_size = 0.75)

colData(cds_FGF)$cluster <- monocle3::clusters(cds_FGF)
plot_cells(cds_FGF)

# Match clusters to cell types, clusters might differ between runs
colData(cds_FGF)$cell_type <- dplyr::case_match(colData(cds_FGF)$cluster,
                                                c("1","5", "14") ~ "D14 Perivascular",
                                                c("2","8", "13") ~ "D14 Endothelial",
                                                c("12", "11") ~ "D5 Common Precursor",
                                                "10" ~ "D28",
                                                .default = "iPSC")

plot_cells(cds_FGF, 
           color_cells_by = "cell_type",
           cell_size = 0.25) 

# Cell type marker heatmap
# Create list of markers for each annotated cell type
manual_markers_by_cluster <- list(
  iPSC = c("POU5F1", "SOX2", "MYC", "NANOG", "POU2F1", "DNMT3B", "LIN28A", "PRDM14", "SALL4"),
  Endothelial = c("PECAM1", "CDH5", "VWF", "TEK", "PDGFB", "CLDN5", "KDR", "THSD7A", "NRP1", "FLT1", "EFNB2", "HSPG2", "ITGA9", "NOTCH4"),
  Perivascular = c("PDGFRB", "CSPG4", "PRRX1", "ACTA2", "COL1A1", "FN1", "COL1A2", "COL3A1", "VIM", "COL5A1", "COL4A1", "MYL9"))

all_manual_markers <- unique(unlist(manual_markers_by_cluster))

# Calculate average expression of all genes grouped by cell type
agg_expr <- aggregate_gene_expression(cds_FGF, 
                                      cell_group_df = data.frame(Cell_ID = colnames(cds_FGF),
                                                                 Cell_Type = colData(cds_FGF)$cell_type))
# Subset to visualize genes of interest
agg_expr <- agg_expr[all_manual_markers, -c(3)] # Remove D28

# Reorder columns 
column_order <- c("iPSC", "D5 Common Precursor", "D14 Endothelial", "D14 Perivascular")
agg_expr <- agg_expr[, column_order]

# Reorder rows 
row_order <- unlist(manual_markers_by_cluster)
agg_expr <- agg_expr[row_order,]

annotation_df <- data.frame(Cluster = factor(column_order, levels = column_order))
rownames(annotation_df) <- column_order

# Plot heatmap (Fig S23B)
# Expression values are scaled across rows aka genes
pheatmap::pheatmap(agg_expr,
                           annotation_col = annotation_df,
                           color = brewer.pal(n = 9, name = "Greens"),  
                           scale = "row",
                           cluster_rows = F,
                           cluster_cols = F,
                           show_rownames = TRUE,
                           show_colnames = TRUE)
