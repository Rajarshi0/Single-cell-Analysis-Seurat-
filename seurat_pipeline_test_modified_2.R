setwd('/home/rajarshi/')

#edits................

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(svglite)
library(ggridges)
library(pheatmap)
library(ggrepel)
library(reshape2)
library(grid)
library(cowplot)
library(tibble)


#edits................

# Data Loading

arranged_data <- read.csv("data/GSE246269_count_matrix.csv", row.names = NULL)

if (any(duplicated(rownames(arranged_data)))) {
  cat("Found duplicate Genes")
}

arranged_data <- arranged_data[!duplicated(arranged_data$gene_symbols), ]

head(arranged_data)
rownames(arranged_data) <- arranged_data$gene_symbols
arranged_data <- arranged_data[,-1]

metadata <- read.csv("data/metadata6.csv", row.names = 1)





 head(arranged_data)
 head(metadata)

metadata_rownames <- rownames(metadata)
arranged_data_colnames <- colnames(arranged_data)
matches <- metadata_rownames %in% arranged_data_colnames

if (!all(matches)) {
  stop("Mismatch between metadata row names and expression data column names.")
}

 #Ensure the column "Condition" exists
if (!"Condition" %in% colnames(metadata)) {
  stop("Column 'Condition' not found in the metadata.")
}

merged_data <- t(arranged_data)
merged_data <- data.frame(merged_data)
merged_data$SampleID <- rownames(merged_data)
merged_data <- merge(merged_data, metadata, by.x = "SampleID", by.y = "row.names")

melted_data <- melt(merged_data, id.vars = c("SampleID", "Condition"))


gene_ids <- rownames(arranged_data)
seurat <- CreateSeuratObject(counts = arranged_data, meta.data = metadata)
#seurat <- NormalizeData(seurat)
#edits................


#library(Seurat)

#counts <- Read10X(data.dir = "data/DS1/")
#seurat <- CreateSeuratObject(counts, project="DS1")

#library(Matrix)
#counts <- readMM("data/DS1/matrix.mtx.gz")
#barcodes <- read.table("data/DS1/barcodes.tsv.gz", stringsAsFactors=F)[,1]
#features <- read.csv("data/DS1/features.tsv.gz", stringsAsFactors=F, sep="\t", header=F)
#rownames(counts) <- make.unique(features[,2])
#colnames(counts) <- barcodes

#seurat <- CreateSeuratObject(counts, project="DS1")


seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT[-\\.]")

VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ggsave("data/plot_1.png", plot = last_plot(), width = 8, height = 6)

ggsave("data/plot_1.pdf", plot = last_plot(), width = 8, height = 6)

VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)

ggsave("data/plot_2.png", plot = last_plot(), width = 8, height = 6)

ggsave("data/plot_2.pdf", plot = last_plot(), width = 8, height = 6)

library(patchwork)
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave("data/plot_3.png", plot = last_plot(), width = 8, height = 6)

ggsave("data/plot_3.pdf", plot = last_plot(), width = 8, height = 6)

#seurat <- subset(seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5)

#Step 3 to 5 
seurat <- NormalizeData(seurat)


seurat <- FindVariableFeatures(seurat, nfeatures = 3000)



top_features <- head(VariableFeatures(seurat), 20)
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE)
plot1 + plot2

ggsave("data/plot_4.png", plot = last_plot(), width = 8, height = 6)

ggsave("data/plot_4.pdf", plot = last_plot(), width = 8, height = 6)

seurat <- ScaleData(seurat)

#seurat <- ScaleData(seurat, vars.to.regress = c("nFeature_RNA", "percent.mt"))


seurat <- RunPCA(seurat, npcs = 50)
ElbowPlot(seurat, ndims = ncol(Embeddings(seurat, "pca")))

ggsave("data/plot_5.png", plot = last_plot(), width = 8, height = 6)

ggsave("data/plot_5.pdf", plot = last_plot(), width = 8, height = 6)

png(file="data/plot_heatmap.png")
#pdf(file="data/plot_heatmap.pdf")
PCHeatmap(seurat, dims = 1:20, cells = 300, balanced = TRUE, ncol = 4)
dev.off()


pdf(file="data/plot_heatmap.pdf")
PCHeatmap(seurat, dims = 1:20, cells = 300, balanced = TRUE, ncol = 4)
dev.off()


#ggsave("data/plot_heatmap.png", plot = last_Heatmap(), width = 8, height = 6)

#ggsave("data/plot_heatmap.pdf", plot = last_Heatmap(), width = 8, height = 6)


#edit###############
viz_dim_loadings <- VizDimLoadings(seurat, dims = 1:2, reduction = "pca")
print(viz_dim_loadings)
ggsave("data/viz_dim_loadings.png", plot = last_plot(), width = 8, height = 6)

#DimPlot(pbmc, reduction = "pca") + NoLegend()


png(file="data/plot_heatmap2.png")
DimHeatmap(seurat, dims= 1:20, cells =300, balanced = TRUE)
dev.off()


pdf(file="data/plot_heatmap2.pdf")
DimHeatmap(seurat, dims= 1:20, cells =300, balanced = TRUE)
dev.off()

ggsave("data/plot_9_heatmap2.png", plot = last_plot(), width = 8, height = 6)


# Extract PCA feature loadings
loadings <- seurat[["pca"]]@feature.loadings

# Convert the loadings to a data frame
loadings_df <- as.data.frame(loadings)
loadings_df$gene <- rownames(loadings_df)

# Function to create a plot for a given PC
create_pc_plot <- function(pc_num) {
  pc <- paste0("PC_", pc_num)
  top_genes <- loadings_df %>% arrange(desc(abs(get(pc)))) %>% head(10)
  top_genes$PC <- pc
  
  # Melt the data frame for ggplot
  loadings_melted <- melt(top_genes, id.vars = c("gene", "PC"))
  
  # Plot the feature loadings
  plot <- ggplot(loadings_melted, aes(x = gene, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal(base_size = 15) +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 8)
    ) +
    labs(title = paste("Top 10 Genes for", pc), x = "Genes", y = "Loading Value")
  
  return(plot)
}

library(ggplot2)

for (i in 1:10) {
  plot <- create_pc_plot(i)
  ggsave(plot = plot, filename = paste0("data/PC_", i, "_top_10_genes.png"))
}

#for (i in 1:10) {
#  plot <- create_pc_plot(i)
#  save_plot(plot, paste0("data/PC_", i, "_top_10_genes"))
#}
#edit*************************


seurat <- RunTSNE(seurat, dims = 1:20, perplexity = 2)
seurat <- RunUMAP(seurat, dims = 1:20, perplexity = 2)


plot1 <- TSNEPlot(seurat)
plot2 <- UMAPPlot(seurat)
plot1 + plot2

ggsave("data/plot_6.png", plot = last_plot(), width = 8, height = 6)

ggsave("data/plot_6.pdf", plot = last_plot(), width = 8, height = 6)

#plot1 <- FeaturePlot(seurat, c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","TFAP2A"),
#                     ncol=3, reduction = "tsne")
#plot2 <- FeaturePlot(seurat, c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","TFAP2A"),
 #                    ncol=3, reduction = "umap")
#plot1 / plot2

#ggsave("data/plot_7.png", plot = last_plot(), width = 8, height = 6)

#ggsave("data/plot_7.pdf", plot = last_plot(), width = 8, height = 6)

seurat <- FindNeighbors(seurat, dims = 1:20)
seurat <- FindClusters(seurat, resolution = 1)

plot1 <- DimPlot(seurat, reduction = "tsne", label = TRUE)
plot2 <- DimPlot(seurat, reduction = "umap", label = TRUE)
plot1 + plot2

ggsave("data/plot_8.png", plot = last_plot(), width = 8, height = 6)

ggsave("data/plot_8.pdf", plot = last_plot(), width = 8, height = 6)   ##perfect done


seurat.pos.markers <- FindAllMarkers(seurat, only.pos = TRUE)
seurat.pos.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

seurat.neg.markers <- FindAllMarkers(seurat, only.neg = TRUE)
seurat.neg.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC < 1)






library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)  # For better text label placement

# Specify the directory where plots will be saved
plot_dir <- "data/"

# Ensure the directory exists
if (!dir.exists(plot_dir)) {
    dir.create(plot_dir)
}

# Detect clusters from Seurat object
seurat_clusters <- unique(seurat@meta.data$seurat_clusters)

# Loop through each cluster
for (cluster in seurat_clusters) {
    cat("Processing Cluster", cluster, "\n")
    
    # Extract markers for the current cluster
    markers_cluster <- FindMarkers(seurat, ident.1 = cluster)

    # Convert row names to a column named 'gene'
    markers_cluster <- markers_cluster %>%
        tibble::rownames_to_column(var = "gene")

    # Add cluster information
    markers_cluster$cluster <- cluster

    # Identify positive and negative markers
    markers_cluster <- markers_cluster %>%
        mutate(marker_type = ifelse(avg_log2FC > 1, "Positive", 
                                    ifelse(avg_log2FC < -1, "Negative", "None")))

    # Identify top 15-20 outermost markers based on fold change
    top_markers <- markers_cluster %>%
        arrange(desc(abs(avg_log2FC))) %>%
        slice_max(abs(avg_log2FC), n = 20) %>%
        mutate(outermost = TRUE) %>%
        select(gene, avg_log2FC, p_val_adj, marker_type, outermost)

    # Remove duplicates and join with original data
    markers_cluster <- markers_cluster %>%
        left_join(top_markers %>% select(gene, outermost), by = "gene") %>%
        mutate(outermost = ifelse(is.na(outermost), FALSE, TRUE))

    # Prepare data for plotting
    plot_data <- markers_cluster %>%
        mutate(label = ifelse(outermost, as.character(gene), NA))  # Only label outermost markers

    # Create the volcano plot
    volcano_plot <- ggplot(plot_data, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
        geom_point(aes(color = marker_type), alpha = 0.5) +
        scale_color_manual(values = c("Positive" = "red", "Negative" = "blue", "None" = "grey")) +
        geom_text_repel(aes(label = label), size = 3, box.padding = 0.5, na.rm = TRUE) +
        labs(title = paste("Cluster", cluster, "Volcano Plot"),
             x = "Log2 Fold Change",
             y = "-Log10 Adjusted p-value") +
        theme_minimal() +
        theme(legend.position = "right")

    # Save the volcano plot as PNG
    ggsave(filename = file.path(plot_dir, paste0("Cluster_", cluster, "_VolcanoPlot.png")), 
           plot = volcano_plot, 
           width = 10, 
           height = 6)

    # Save the volcano plot as PDF
    ggsave(filename = file.path(plot_dir, paste0("Cluster_", cluster, "_VolcanoPlot.pdf")), 
           plot = volcano_plot, 
           width = 10, 
           height = 6)
}



