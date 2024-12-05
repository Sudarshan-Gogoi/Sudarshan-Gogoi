# Here we have provided the complete instructions to apply the code for analysis of scRNA-seq datasets of 
# cancer obtained from 10x Genomics database by applying the methodology mentioned in our paper named 
# "An integrated computational framework utilizing single-cell genomics for precise classification and prediction of multiple cancer types" 
# by Sudarshan Gogoi, Soumen Bera and Amit Chakraborty
# RStudio Version: 2023.12.1 (402)
# Example file used for illustration: "V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5"
# The code provides analysis of the Breast cancer Dataset (b) as an example
# Load libraries
library(Seurat)
library(tidyverse)
library(data.table)
library(Matrix)
library(cluster)
library("Hmisc")
library(reshape2)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(cowplot)
#Step 1: Data preparation
# Load the dataset
Breast.sparse.m <- Read10X_h5(filename ="D:/My PhD information folder first Paper/scRNA seq data(cancer)/Filtered Data/Breast cancer/Data4/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5")
str(Breast.sparse.m)
cnts <-Breast.sparse.m
# Initialize the Seurat object with the filtered data collected from the 10x genomics database
Breastcancer.seurat.obj<- CreateSeuratObject(counts = cnts, project = "Breastcancer",min.cells = 3, min.features = 100)
str(Breastcancer.seurat.obj)
Breastcancer.seurat.obj
#Step 2: Quality control and filtering(Table 1)
View(Breastcancer.seurat.obj@meta.data)
# Calculation of mitochondrial genes percentage in each cell
Breastcancer.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(Breastcancer.seurat.obj, pattern = "^MT-")
View(Breastcancer.seurat.obj@meta.data)
# Violin plot showing the distribution of the data
dev.new(height=7, width=7)
VlnPlot(Breastcancer.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# Scatter plot illustrating the data distribution
dev.new(height=7, width=7)
FeatureScatter(Breastcancer.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
# Filtration to subset the dataset by setting the thresholds(as per requirement) from violin and scatter plots observation(See Table 1)
# Breastcancer.seurat.obj <- subset(Breastcancer.seurat.obj, subset = nFeature_RNA > Min features threshold & nFeature_RNA < Max features threshold & percent.mt < Mitochondrial genes% threshold)
# Breastcancer.seurat.obj
# Skip this step for filtered data within the range of the sample size(2500-5000 cells)
#Step 3: Data normalization
Breastcancer.seurat.obj<- NormalizeData(Breastcancer.seurat.obj,normalization.method = "LogNormalize", scale.factor = 10000)
str(Breastcancer.seurat.obj)
#Step4: Top 2000 high variable features identification
Breastcancer.seurat.obj <- FindVariableFeatures(Breastcancer.seurat.obj, selection.method = "vst", nfeatures = 2000)
# Identifying the 3 most highly variable genes
top3 <- head(VariableFeatures(Breastcancer.seurat.obj), 3)
# Plot variable features with and without labels
dev.new(height=7, width=7)
plot1 <- VariableFeaturePlot(Breastcancer.seurat.obj)
LabelPoints(plot = plot1, points = top3)
#Step5: Scaling data
all.genes <- rownames(Breastcancer.seurat.obj)
Breastcancer.seurat.obj <- ScaleData(Breastcancer.seurat.obj, features = all.genes)
str(Breastcancer.seurat.obj)
#Step6: Perform linear dimensionality reduction technique PCA to identify principal components
Breastcancer.seurat.obj<- RunPCA(Breastcancer.seurat.obj, features = VariableFeatures(object = Breastcancer.seurat.obj))
# Determine the ideal amount of principal components to keep for subsequent analyses by observing the elbow plot 
dev.new(height=7, width=7)
ElbowPlot(Breastcancer.seurat.obj)
#Step7: Clustering(See Table S1)
Breastcancer.seurat.obj <- FindNeighbors(Breastcancer.seurat.obj, dims = 1:10)
#Fitting a resolution to form 8 clusters in the dataset
Breastcancer.seurat.obj <- FindClusters(Breastcancer.seurat.obj, resolution = c(0.2))
View(Breastcancer.seurat.obj@meta.data)
#Dimensional plot showing 8 identified clusters in the dataset derived from PCA
dev.new(height=7, width=7)
DimPlot(Breastcancer.seurat.obj, group.by = "RNA_snn_res.0.2", label = TRUE)
# setting identity of the clusters
Idents(Breastcancer.seurat.obj)
Idents(Breastcancer.seurat.obj) <- "RNA_snn_res.0.2"
Idents(Breastcancer.seurat.obj)
#Step8: Perform non-linear dimensionality reduction technique UMAP
Breastcancer.seurat.obj <- RunUMAP(Breastcancer.seurat.obj, dims = 1:10)
#UMAP dimensional plot illustrating 8 identified clusters
dev.new(height=7, width=7)
DimPlot(Breastcancer.seurat.obj, reduction = "umap",label = TRUE)
# Clusters quality checking
# Get UMAP embeddings
umap_emb <- as.data.frame(Breastcancer.seurat.obj@reductions$umap@cell.embeddings)
# Get cluster information
cluster_info <- Breastcancer.seurat.obj$RNA_snn_res.0.2
# Check if cluster information is available
if (!is.null(cluster_info)) {
  # Convert cluster info to numeric
  cluster_info <- as.numeric(as.character(cluster_info))
  # Check if the cluster info is numeric
  if (!anyNA(cluster_info) && is.numeric(cluster_info)) {
    # Calculate silhouette scores
    silhouette_scores <- silhouette(cluster_info, dist(umap_emb))
    # View the silhouette scores
    print(silhouette_scores)
  } else {
    print("Cluster information is not numeric or contains missing values.")
  }
} else {
  print("Cluster information not found.")
}
# Calculate silhouette statistics
sil_summary <- summary(silhouette_scores)
# Extract the mean silhouette width
mean_silhouette <- sil_summary$avg.width
# View the mean silhouette score
print(paste("Mean Silhouette Score:", mean_silhouette))
#Step9: Cluster-wise Highly Expressed genes extraction
# Create an empty list to store Highly_Expressed genes for each cluster
Highly_Expressed_genes_list <- list()
# Get the cluster identities
cluster_ids <- levels(Idents(Breastcancer.seurat.obj))
# Loop over each cluster
for (cluster_id in cluster_ids) {
  # Subset the cells belonging to the current cluster
  cluster_cells <- WhichCells(Breastcancer.seurat.obj, ident = cluster_id)
  # Get the gene expression matrix for the current cluster
  cluster_expr_matrix <- GetAssayData(Breastcancer.seurat.obj, cells = cluster_cells, assay = "RNA")
  # Calculate the average expression of each gene across cells in the cluster
  avg_expr_per_gene <- rowMeans(cluster_expr_matrix)
  # Get the top 1000 highly expressed genes for this cluster
  Highly_Expressed_genes <- names(sort(avg_expr_per_gene, decreasing = TRUE)[1:1000])
  # Store Highly_Expressed genes for this cluster in the list with appropriate names
  Highly_Expressed_genes_list[[paste0("Highly_Expressed_cluster", cluster_id)]] <- Highly_Expressed_genes
}
#Step10: Gene correlation network construction
  # Extract the meta data for cells in clusters 0,1,2,3,4,5,6,7
  # Create an empty list to store metadata for each cluster
cluster_meta_list <- list()
# Loop through clusters from 0 to 7
for (i in 0:7) {
  # Extract metadata for the current cluster
  cluster_meta <- Breastcancer.seurat.obj@meta.data[Breastcancer.seurat.obj@meta.data$seurat_clusters == i, ]
  # Assign the metadata to a variable
  assign(paste0("cluster_", i, "_meta"), cluster_meta)
  # Store the metadata in the list
  cluster_meta_list[[paste0("cluster_", i, "_meta")]] <- cluster_meta
}
# Get the cell IDs from the Seurat object
cell_ids <- rownames(Breastcancer.seurat.obj@meta.data)
# Add the cell IDs to the Seurat object's metadata data.frame
Breastcancer.seurat.obj@meta.data$cell_id <- cell_ids
# Create an empty list to store filtered correlation matrices for each cluster
filtered_correlation_list <- list()
    # Loop through clusters from 0 to 7
    for (i in 0:7) {
      # Extract cell names for the current cluster
      cell_names <- rownames(get(paste0("cluster_", i, "_meta")))
      # Subset the Seurat object to include only the specified cells
      seurat_subset <- Breastcancer.seurat.obj[, Breastcancer.seurat.obj@meta.data$cell_id %in% cell_names]
      # Extract the Highly_Expressed genes for the current cluster
      Highly_Expressed_genes <- Highly_Expressed_genes_list[[paste0("Highly_Expressed_cluster", i)]]
      # Subsetting individual cluster's gene expression data by keeping only Highly_Expressed genes
      cell_data <- as.data.frame(x = GetAssayData(seurat_subset, assay = "RNA"))
      # Add the gene names column to the Seurat subset dataframe
      cell_data$genes <- rownames(cell_data)
      # Filter the Seurat subset dataframe by including only selected Highly_Expressed genes
      filtered_data <- cell_data[cell_data$genes %in% Highly_Expressed_genes, ]
      # Gene correlation matrix construction
      m <- filtered_data[, -c(1, ncol(filtered_data))]
      m1 <- t(m)
      expression_matrix <- as.matrix(m1)
      # Calculate the correlation matrix
      res <- rcorr(m1)
      # Extract correlation and p-value matrices
      n2 <- as.data.frame(res$r)
      n1 <- as.data.frame(res$P)
      # Set the significance level (alpha)
      alpha <- 0.05
      # Filter the correlation matrix based on p-values
      filtered_correlation <- n2
      filtered_correlation[n1 >= alpha] <- 0
      # Re-add gene names as row names and column names
      rownames(filtered_correlation) <- rownames(filtered_data)
      colnames(filtered_correlation) <- rownames(filtered_data)
      # Store the filtered correlation matrix in the list
      filtered_correlation_list[[paste0("filtered_correlation_", i)]] <- filtered_correlation
    }
    # Print the list to check its contents
    #print(filtered_correlation_list)
# Gene correlation network construction of each cluster in the form of adjacency list
# Create a list to store threshold-wise adjacency lists
thresholdwise_adj_lists <- list()
# Loop through thresholds
thresholds <- c(0.50, 0.53, 0.56, 0.59, 0.62, 0.65, 0.68)
for (threshold in thresholds) {
  # Create a list to store adjacency lists for this threshold
  threshold_adj_lists <- list()
  # Loop through filtered correlation matrices
  for (i in 0:7) {
    # Access the filtered correlation matrix for the current cluster
    filtered_correlation <- filtered_correlation_list[[paste0("filtered_correlation_", i)]]
    # Set diagonal to 0
    diag(filtered_correlation) <- 0
    # Replace upper triangle values
    filtered_correlation[upper.tri(filtered_correlation)] <- 42
    # Convert matrix to data frame
    filtered_correlation_df <- as.data.frame(filtered_correlation)
    filtered_correlation_df <- rownames_to_column(filtered_correlation_df, var = "variable1")
    # Melt the data frame
    my_cor_df <- reshape2::melt(filtered_correlation_df, id.vars = "variable1")
    # Filter out unnecessary values and self correlations
    my_cor_df <- filter(my_cor_df, value != 42) %>% filter(variable1 != variable)
    # Create adjacency list with correlations > threshold
    adjacency_list <- my_cor_df %>% filter(abs(value) > threshold)
    # Add adjacency list to the list
    threshold_adj_lists[[paste0("cluster_", i, "_adjacency")]] <- adjacency_list
  }
  # Add threshold-wise adjacency lists to the main list
  thresholdwise_adj_lists[[paste0("threshold_", threshold, "_adj_lists")]] <- threshold_adj_lists
}
# Step 11: Hausdorff distance matrix and norm calculation
# Create a list to store Hausdorff distance matrices for thresholds 0.50, 0.53, 0.56, 0.59, 0.62, 0.65, and 0.68
Hausdorff_distance_matrices <- list()
thresholds <- c(0.50, 0.53, 0.56, 0.59, 0.62, 0.65, 0.68)
for (threshold in thresholds) {
  # Extract data frames of adjacency lists for the current threshold
  A1 <- data.frame(thresholdwise_adj_lists[[paste0("threshold_", threshold, "_adj_lists")]]$cluster_0_adjacency, directed = FALSE)
  B1 <- data.frame(thresholdwise_adj_lists[[paste0("threshold_", threshold, "_adj_lists")]]$cluster_1_adjacency, directed = FALSE)
  C1 <- data.frame(thresholdwise_adj_lists[[paste0("threshold_", threshold, "_adj_lists")]]$cluster_2_adjacency, directed = FALSE)
  D1 <- data.frame(thresholdwise_adj_lists[[paste0("threshold_", threshold, "_adj_lists")]]$cluster_3_adjacency, directed = FALSE)
  E1 <- data.frame(thresholdwise_adj_lists[[paste0("threshold_", threshold, "_adj_lists")]]$cluster_4_adjacency, directed = FALSE)
  F1 <- data.frame(thresholdwise_adj_lists[[paste0("threshold_", threshold, "_adj_lists")]]$cluster_5_adjacency, directed = FALSE)
  G1 <- data.frame(thresholdwise_adj_lists[[paste0("threshold_", threshold, "_adj_lists")]]$cluster_6_adjacency, directed = FALSE)
  H1 <- data.frame(thresholdwise_adj_lists[[paste0("threshold_", threshold, "_adj_lists")]]$cluster_7_adjacency, directed = FALSE)
  # Extract the correlation values column from the data frames of the adjacency lists
  A <- as.matrix(A1$value)
  B <- as.matrix(B1$value)
  C <- as.matrix(C1$value)
  D <- as.matrix(D1$value)
  E <- as.matrix(E1$value)
  F <- as.matrix(F1$value)
  G <- as.matrix(G1$value)
  H <- as.matrix(H1$value)
  # List of data frames' correlation values columns
  dataframes <- list(A, B, C, D, E, F, G, H)
  # Function to compute relative Hausdorff distance
  d_H <- function(M, N) {
    d_M <- apply(M, 1, function(x) min(sqrt(colSums((x - t(N))^2))))
    d_N <- apply(N, 1, function(x) min(sqrt(colSums((x - t(M))^2))))
    max(c(max(d_M), max(d_N)))
  }
  # Create an empty matrix to store distances
  distance_matrix <- matrix(0, nrow = length(dataframes), ncol = length(dataframes))
  # Calculate Hausdorff distances and fill the matrix
  for (i in 1:length(dataframes)) {
    for (j in 1:length(dataframes)) {
      if (i != j) {
        haus_dist <- d_H(dataframes[[i]], dataframes[[j]])
        distance_matrix[i, j] <- haus_dist
      }
    }
  }
  # Assign the distance matrix to the list for the current threshold
  Hausdorff_distance_matrices[[paste0("threshold_", threshold, "_distance_matrix")]] <- distance_matrix
}
# Print the Hausdorff distance matrices
print(Hausdorff_distance_matrices)
# Function to calculate Frobenius norm of a matrix
frobenius_norm <- function(matrix) {
  sqrt(sum(matrix^2))
}
# Create a vector to store Frobenius norms
frobenius_norms <- numeric(length(Hausdorff_distance_matrices))
# Calculate Frobenius norms for each matrix in increasing threshold order
for (i in seq_along(Hausdorff_distance_matrices)) {
  matrix <- Hausdorff_distance_matrices[[i]]
  norm <- frobenius_norm(matrix)
  frobenius_norms[i] <- norm
}
# Print Frobenius norms in increasing threshold order
for (i in seq_along(Hausdorff_distance_matrices)) {
  threshold <- thresholds[i]
  norm <- frobenius_norms[i]
  cat("Threshold:", threshold, "- Frobenius Norm:", norm, "\n")
}
# Combine thresholds and Frobenius norms into a data frame
norms_df <- data.frame(threshold = thresholds, frobenius_norm = frobenius_norms)
# Save the combined norms data frame to a CSV file
write.csv(norms_df, file = "D:/My PhD information folder first Paper/Data Analysis in R Modified/Lung/Data3/frobenius_normsB(b).csv", row.names = FALSE)
#Step 12: Construction of stacked line plots showing the trends in norm values at different gene correlation thresholds in the cancer dataset
# Prepare data
norm_data <- norms_df %>% select(threshold, frobenius_norm)
data_long <- gather(norm_data, key = "Variable", value = "Value", -threshold)
# Plot
p <- ggplot(data_long, aes(x = threshold, y = Value, color = Variable, group = Variable)) +
  geom_line() +
  geom_point() +
  labs(
    x = "threshold",
    y = "Frobenius Norm",
    color = "Variables"
  ) +
  scale_color_manual(
    values = c("blue"),
    breaks = c("frobenius_norm"),
    name = "Variables"
  ) +
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0)) +  # Set Y-axis limits from 0 to 10
  theme_minimal() +
  theme(
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15, face = "bold"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    legend.spacing.y = unit(1, "cm")
  )
# Print the plot
dev.new(height=7, width=7)
print(p)
# Save the plot with 600 dpi
ggsave("plot.png", plot = p, height = 7, width = 9, dpi = 600)

# Additional Analysis
# (a)Construct the loaded dataset's gene correlation network of each cluster in the form of adjacency list at gene correlation threshold of 0.6
# Define the correlation threshold
threshold <- 0.6
# Initialize the list to store adjacency lists for different clusters
adjacency_lists <- list()
# Loop through filtered correlation matrices for the threshold
for (i in 0:7) {  # Adjust the loop to run from 0 to 7
  # Access the filtered correlation matrix for the current cluster
  filtered_correlation <- filtered_correlation_list[[paste0("filtered_correlation_", i)]]
  # Set diagonal to 0
  diag(filtered_correlation) <- 0
  # Replace upper triangle values
  filtered_correlation[upper.tri(filtered_correlation)] <- 42
  # Convert matrix to data frame
  filtered_correlation_df <- as.data.frame(filtered_correlation)
  filtered_correlation_df <- rownames_to_column(filtered_correlation_df, var = "variable1")
  # Melt the data frame
  my_cor_df <- reshape2::melt(filtered_correlation_df, id.vars = "variable1")
  # filter out the self correlations and unnecessary values
  my_cor_df <- filter(my_cor_df, value != 42) %>% filter(variable1 != variable)
  dim(my_cor_df)
  # create adjacency list with correlations > threshold
  adjacency_list <- my_cor_df %>% filter(abs(value) > threshold)
  # Add adjacency list to the main list
  adjacency_lists[[paste0("Network_Cluster", i, "_threshold_", gsub("\\.", "", as.character(threshold)))]] <- adjacency_list
}
# Now adjacency_lists contains the adjacency lists for gene correlation networks at threshold 0.6
# (b)Calculation of number of genes and edges present in the networks of different clusters
num_genes <- numeric(length(adjacency_lists))
num_edges <- numeric(length(adjacency_lists))
# Loop through adjacency lists
for (i in seq_along(adjacency_lists)) {
  # Access the adjacency list for the current cluster
  adjacency_list <- adjacency_lists[[i]]
  # Extract unique genes from both 'variable1' and 'variable' columns
  all_genes <- unique(c(adjacency_list$variable1, adjacency_list$variable))
  # Count the number of genes
  num_genes[i] <- length(all_genes)
  # Count the number of edges
  num_edges[i] <- nrow(adjacency_list)
}
# Create a data frame to store the results
Gene_Edge_Numbers <- data.frame(
  Network = names(adjacency_lists),
  NumGenes = num_genes,
  NumEdges = num_edges
)
# Print the Gene_Edge_Numbers
print(Gene_Edge_Numbers)
# (c) Shannon index calculation of the networks on the basis of the number of genes and edges present in the networks
# Calculate the total count of genes or edges (total abundance)
total_abundance_genes <- sum(Gene_Edge_Numbers$NumGenes)
total_abundance_edges <- sum(Gene_Edge_Numbers$NumEdges)
# Calculate the proportion of each network in the community
proportions_genes <- Gene_Edge_Numbers$NumGenes / total_abundance_genes
proportions_edges <- Gene_Edge_Numbers$NumEdges / total_abundance_edges
# Calculate the Shannon Index for genes
shannon_index_genes <- -sum(proportions_genes * log2(proportions_genes), na.rm = TRUE)
# Calculate the Shannon Index for edges
shannon_index_edges <- -sum(proportions_edges * log2(proportions_edges), na.rm = TRUE)
# Print the results
cat("Shannon Index for Number of Genes:", shannon_index_genes, "\n")
cat("Shannon Index for Number of Edges:", shannon_index_edges, "\n")
