# Setup and load libraries
suppressPackageStartupMessages({ 
  library(Matrix)
  library(Seurat)
  library(tidyverse)
  library(monocle3)
  library(viridis)
  library(RColorBrewer)
  library(pheatmap)
  library(DESeq2)
  
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  # Set a seed to make umap and other non-deterministic steps consistent
  set.seed(seed = 23)
})

# Load raw counts, cell metadata
counts <- readMM("data/Exp2_raw_counts.mtx") # Raw counts from GEO
genes <- read.csv("data/genes.tsv", header = F)
barcodes <- read.csv("data/barcodes.tsv", header = F)

rownames(counts) <- genes[, 1]
colnames(counts) <- barcodes[, 1]

mdata <- read.csv("data/cell_metadata.csv", row.names = 1)

# Create Seurat object
seurat_FGF <- CreateSeuratObject(counts = counts,
                                 meta.data = mdata)

seurat_FGF

#QUALITY CONTROL
# Mitochondrial gene detection
# "MT-" genes sequenced from mitochondrial genome 
mito_genes <- rownames(seurat_FGF)[grep("^MT-", rownames(seurat_FGF))]
head(mito_genes)

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
                                   assay = "RNA",
                                   "^RP[SL]", 
                                   col.name = "percent_ribo")
head(seurat_FGF@meta.data)

# Red blood cell contamination
# Red blood cell genes annotated "HB.." indicating RBC contamination
hb_genes <- rownames(seurat_FGF)[grep("^HB[^(P)]", rownames(seurat_FGF))]
hb_genes

# Percentage of reads per cell mapping to RBC genes
# <0.5% usually considered a good threshold
# This should really be 0, check raw data again if %s are high
seurat_FGF <- PercentageFeatureSet(seurat_FGF, 
                                   "^HB[^(P)]", 
                                   col.name = "percent_hb")
head(seurat_FGF@meta.data)

# Plot QC metrics
Idents(seurat_FGF) <- "sample"
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
plot <- VlnPlot(seurat_FGF, 
                 features = feats, 
                 pt.size = 0, 
                 ncol = 3) +
  NoLegend()
plot

# Subset Seurat object based on calculated QC metrics
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

# Detection-based filtering 
# Remove cells expressing <200 genes
# Remove genes expressed in <3 cells
# Remove all non-protein coding genes
selected_c <- WhichCells(seurat_FGF, expression = nFeature_RNA > 200)
selected_f <- rownames(seurat_FGF)[Matrix::rowSums(seurat_FGF@assays$RNA@layers$counts) > 3]

df_gene <- readRDS("data/df_gene.RDS")
genes_coding <- df_gene %>%
  dplyr::filter(gene_type == "protein_coding") %>%
  dplyr::select(gene_short_name)

seurat_FGF <- subset(seurat_FGF, 
                     features = intersect(selected_f, genes_coding$gene_short_name), 
                     cells = selected_c)
seurat_FGF

# Revisualize QC metrics
Idents(seurat_FGF) <- "sample"
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
plot <- VlnPlot(seurat_FGF, 
                 features = feats, 
                 pt.size = 0, 
                 ncol = 3) +
  NoLegend()
plot

# Correlation between recovered UMIs per cell and detected genes
plot <- FeatureScatter(seurat_FGF, 
               "nCount_RNA", 
               "nFeature_RNA", 
               pt.size = 0.5) +
  ggplot2::theme_classic()
plot + NoLegend()

# Cell-cycle scoring
seurat_FGF = NormalizeData(seurat_FGF, verbose = T)

# Percentage of counts mapping to G2M/S phase genes
seurat_FGF <- CellCycleScoring(object = seurat_FGF, 
                             g2m.features = cc.genes$g2m.genes,
                             s.features = cc.genes$s.genes)
seurat_FGF@meta.data

plot <- VlnPlot(seurat_FGF, 
        features = c("S.Score", "G2M.Score"),
        ncol = 2, pt.size = 0)
plot

# Convert to monocle3 CDS for processing
gene_mdata <- data.frame(row.names = rownames(seurat_FGF),
                                                gene_short_name = rownames(seurat_FGF))

# Create new monocle3 cds
cds_FGF <- new_cell_data_set(expression_data = seurat_FGF[["RNA"]]@layers$counts,
                             cell_metadata = seurat_FGF@meta.data,
                             gene_metadata = gene_mdata)
cds_FGF

# DATA PROCESSING
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
cds_FGF <- reduce_dimension(cds_FGF, verbose = T, umap.min_dist = 0.75)
cds_FGF <- cluster_cells(cds_FGF, verbose = T)

plot_cells(cds_FGF, color_cells_by = "partition")

# Eliminate undesired cells for analysis
# Cluster 10 resembles D5 Perivascular cells (immature)
# Cluster 11 resembles D28 cells (unwanted diff)
colData(cds_FGF)$cluster <- clusters(cds_FGF)
cds_FGF <- cds_FGF[, !colData(cds_FGF)$cluster %in% c(10, 11)]

# Visualize expression of endothelial and perivascular markers to annotate clusters
plot_cells(cds_FGF, genes = c("CDH5", "PDGFRB"), cell_size = 0.75)

# Cell type annotation
colData(cds_FGF)$partition <- partitions(cds_FGF)
colData(cds_FGF)$cell_type <- dplyr::case_match(colData(cds_FGF)$partition,
                                                c("1") ~ "Endothelial",
                                                .default = "Perivascular")
plot_cells(cds_FGF, color_cells_by = "cell_type", cell_size = 1)

# Endothelial/Perivascular proportions grouped by D2-5 treatment (Fig 5B)
# Extract metadata
coldata <- as.data.frame(colData(cds_FGF))

original_df <- coldata %>%
  mutate(treatment = gsub("_[0-9]$","", sample)) %>% # Treatment without replicate number
  dplyr::count(treatment, sample, cell_type) %>% # Count number of cells of each type
  group_by(treatment, sample) %>%
  mutate(total = sum(n),
         percent = (n / total) * 100) %>% # Calculate percentages in each replicate
  ungroup() %>%
  dplyr::filter(treatment %in% c("C6_1.0", "FGF_1.0", "No_FGF", "C2_100.0", "mb7_10.0", "mb7-FGF_10.0")) %>% # Only include treatments shown in figure
  mutate(treatment = factor(treatment, levels = c("No_FGF", "FGF_1.0", "C6_1.0", "C2_100.0", "mb7_10.0", "mb7-FGF_10.0"))) %>% # Change order of treatments in plot
  rbind(data.frame(treatment = "mb7-FGF_10.0",  # Manually change percentages where count is 0 (to avoid division by zero)  
                   sample = "mb7-FGF_10.0_3",
                   cell_type = "Endothelial",
                   n = 0,
                   total = 2,
                   percent = 0))

# Create summary data frame with averages and SE
summary_df <- original_df %>%
  group_by(treatment, cell_type) %>%
  summarise(AVERAGE = mean(percent), # Mean and SEM from 3 replicates
            SE = ifelse(n() > 1, sd(percent) / sqrt(n()), 0), .groups = "drop") %>%
  dplyr::filter(treatment %in% c("C6_1.0", "FGF_1.0", "No_FGF", "C2_100.0", "mb7_10.0", "mb7-FGF_10.0")) %>%
  mutate(treatment = factor(treatment, levels = c("No_FGF", "FGF_1.0", "C6_1.0", "C2_100.0", "mb7_10.0", "mb7-FGF_10.0")))

# Create plot
plot <- ggplot() +
  geom_point(data = summary_df, 
             aes(x = treatment, y = AVERAGE, fill = cell_type),
             stat = "identity", 
             pch = 22,
             alpha = 0.5,
             stroke = 0,
             size = 8) + # Box for average
  geom_errorbar(data = summary_df, 
                aes(x = treatment, ymin = AVERAGE - SE, ymax = AVERAGE + SE, color = cell_type), 
                width = 0.1,
                size = 0.75) + # Error bars
  geom_point(data = original_df, 
             aes(x = treatment, y = percent, fill = cell_type), 
             size = 3, 
             pch = 21, 
             alpha = 0.5,
             position = position_jitter(width = 0.1),
             stroke = 0.25) + # Individual data points
  scale_fill_manual(values = c("darkred", "darkcyan")) +  # Manually set the fill colors
  scale_color_manual(values = c("darkred", "darkcyan")) +  # Manually set the color of error bars and average lines
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = c("No_FGF" = "No FGF", "C2_100.0" = "C2-58-2X", "mb7_10.0" = "mb7", "mb7-FGF_10.0" = "mb7 + FGF2", "C6_1.0" = "C6-79C", "FGF_1.0" = "FGF2")) + # Manually set x labels 
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  theme_classic()

plot

# Endogenously secreted FGFs (Supp Fig 36C)
plot_cells(cds_FGF, genes = c("FGF1", "FGF2", "FGF5", "FGF7", "FGF9", "FGF11", "FGF12", "FGF13", "FGF14", "FGF20", "FGF21", "FGF22", "FGF23"), 
           cell_size = 1,
           label_cell_groups = F) +
  scale_color_viridis(option = "inferno", direction = -1)

# Subclustering endothelial cells
# Subset only annotated endothelial cells and convert to Seurat object
cds_endo <- cds_FGF[,colData(cds_FGF)$cell_type %in% "Endothelial"]
seurat_endo <- as.Seurat(cds_endo, data = NULL)

# Processing Seurat object
seurat_endo <- NormalizeData(seurat_endo, 
                           normalization.method = "LogNormalize", 
                           scale.factor = 10000)
seurat_endo <- FindVariableFeatures(seurat_endo, 
                                    selection.method = "vst", 
                                    nfeatures = 2000)

# Scale and run PCA
seurat_endo <- ScaleData(seurat_endo, features = rownames(seurat_endo))
seurat_endo <- RunPCA(seurat_endo, features = VariableFeatures(object = seurat_endo), npcs = 100)

# Check number of PC components 
ElbowPlot(seurat_endo)

# Cluster and visualize
seurat_endo <- FindNeighbors(seurat_endo, dims = 1:100)
seurat_endo <- FindClusters(seurat_endo, resolution = 0.8)
seurat_endo <- RunUMAP(seurat_endo, dims = 1:100, min.dist = 0.1, seed.use = -3)

DimPlot(seurat_endo, reduction = 'umap')
#seurat_endo <- seurat_endo[, !(seurat_endo@meta.data$seurat_clusters %in% c(9, 11))]

umap_data <- Embeddings(seurat_endo, "umap")

# Define a rotation matrix for 45 degrees
angle <- -45
rotation_matrix <- matrix(c(cos(pi * angle / 180), -sin(pi * angle / 180), 
                            sin(pi * angle / 180), cos(pi * angle / 180)), 
                          nrow = 2)

# Apply the rotation matrix to the UMAP coordinates
rotated_umap <- t(rotation_matrix %*% t(umap_data))

# Update the Seurat object with the new rotated UMAP coordinates
seurat_endo[["umap"]] <- CreateDimReducObject(embeddings = rotated_umap, key = "umap_rotated_")

# Plot the rotated UMAP
DimPlot(seurat_endo, reduction = "umap")

# Process Artery-Vein bulkRNAseq data set (Ang et. al., Cell, 2022)
# Read raw counts
av_counts <- read.delim("data/artery_vein_counts.txt", 
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Filter duplicate gene rows
av_counts <- av_counts %>%
  mutate(gene = ifelse(gene_name == '---', gene_id, gene_name)) %>%
  dplyr::select(SRR19221939, SRR19221940, SRR19221941, SRR19222030, SRR19222034, SRR19222036, gene) %>%
  distinct(gene, .keep_all = T) 

rownames(av_counts) <- av_counts$gene
av_counts <- av_counts[,-c(7)] # Remove 'gene' column
av_counts

# Set up metadata column for DESeq2
col_meta <-
  data.frame(row.names = colnames(av_counts),
             sample = factor(
               c(rep("vein", 3),
                 rep("artery", 3))))

# Set up DESeq2 object 
dds <- DESeqDataSetFromMatrix(countData = av_counts,
                              colData = col_meta,
                              design= ~sample)

# Perform DE analysis
dds <- DESeq(dds)

# Calculate results --> artery vs vein
res <- results(dds, contrast = c("sample", "artery", "vein"))  
res$gene <- row.names(res)

# Only keep sig genes
res <- res %>%
  as.data.frame() %>%
  dplyr::filter(padj < 0.05) 
res

# Log transformation for PCA plotting
rld <- rlog(dds)
plotPCA(rld, "sample") +
  theme_minimal()

# Save DE genes for artery, vein separately
artery <- res$gene[res$log2FoldChange > 0]
vein <- res$gene[res$log2FoldChange < 0]
genes_list <- list(artery, vein)

# Calculate arterial/venous module score for each cell
seurat_endo <- AddModuleScore(seurat_endo,
                  features = genes_list,
                  name = "Endo_subtype_score")

# A-V specificity defined as artery module score - vein module score
seurat_endo[['AV_specificity']] <- seurat_endo[['Endo_subtype_score1']] - seurat_endo[['Endo_subtype_score2']]

# Plot A-V specificity
md <- seurat_endo[[]]
coords <- Embeddings(seurat_endo[["umap"]])
md <- cbind(md, coords)

midpoint <- mean(range(md$AV_specificity))

# Order by absolute difference from the midpoint
md <- md[order(abs(md$AV_specificity - midpoint)), ]

# Plot (Fig 6A, left)
ggplot(md, aes(x = umaprotated_1, y = umaprotated_2, fill = AV_specificity)) +
  geom_point(shape = 21, size = 2, alpha = 1, stroke = 0.01) + 
  #scale_fill_gradientn(colors = colorRampPalette(brewer.pal(9, "PiYG"))(100)) +
  scale_fill_viridis(option = 'B', direction = -1) +
  theme_classic() 

# Visualize Seurat cluster annotations
ggplot(md, aes(x = umaprotated_1, y = umaprotated_2, fill = seurat_clusters)) +
  geom_point(shape = 21, size = 2, alpha = 1, stroke = 0) + 
  theme_classic()

# Assign endothelial subtypes to each cluster
seurat_endo$Endo_subtype <- ifelse(seurat_endo@meta.data$seurat_clusters %in% c(0, 1, 4, 6, 10),
                                                'Venous', 'Arterial')

# Revisualize
md <- seurat_endo[[]]
coords <- Embeddings(seurat_endo[["umap"]])
md <- cbind(md, coords)

ggplot(md, aes(x = umaprotated_1, y = umaprotated_2, fill = Endo_subtype)) +
  geom_point(shape = 21, size = 2, alpha = 1, stroke = 0) + 
  theme_classic()

dat <- as.data.frame(seurat_endo@meta.data)
coords <- Embeddings(seurat_endo[["umap"]])
dat$umap1 = coords[,1]
dat$umap2 = coords[,2]

#Create a data frame that doesn't contain a "sample" column. This will allow us to facet the density layer without affecting the points
dat$treatment <- gsub("_[0-9]$","", dat$sample)
dat_bg <- dat[,-(which(colnames(dat) == "treatment"))]

# (Fig 6A, middle)
dat %>%
  dplyr::filter(treatment %in% c("No_FGF", "FGF_1.0", "C6_1.0")) %>%
  ggplot(aes(x=umap1, y=umap2)) +
  stat_density_2d(geom="raster", aes(fill=after_stat(ndensity)), contour=F) +  #ndensity calculates the normalized density for each sample--otherwise density would be affected by the number of cells for each sample, which is variable
  geom_point(data=dat_bg, shape=16, size=0.5, alpha=0.4, color="black") +
  scale_fill_gradientn(colours=brewer.pal(n = 9, name = "Greys"), name="Density") +
  #scale_fill_gradientn(colours=viridisLite::mako(1000), name="Density") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  facet_wrap(~treatment, ncol=8) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12, color="black"),
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=14, color="black"))

# Percentage of each treatment belonging to each endothelial subtype (Fig 6A, right)
seurat_endo[[]] %>%
  mutate(treatment = gsub("_[0-9]$","", sample)) %>%
  dplyr::filter(treatment %in% c("No_FGF", "FGF_1.0", "C6_1.0")) %>%
  dplyr::count(treatment, sample, Endo_subtype) %>%
  group_by(treatment, sample) %>%
  mutate(total = sum(n),
         percent = (n / total) * 100) %>%
    ungroup() %>%
  group_by(treatment, Endo_subtype) %>%
  summarise(AVERAGE = mean(percent),
            SE = ifelse(n() > 1, sd(percent) / sqrt(n()), 0), .groups = "drop") %>%
  ggplot(aes(x = "", y = AVERAGE, fill = Endo_subtype)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start=0) +
    theme_void() +
    theme(legend.position="bottom") +
    labs(title = paste("Treatment"),
         fill = "Category",
         y = "Average Percentage") +
  scale_color_brewer(palette = "Spectral") +
  facet_wrap(~treatment) 



