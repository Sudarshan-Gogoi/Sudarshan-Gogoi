# Full R function with integration of gene correlation network construction, Hausdorff distances, and Cancer Classification 
# RStudio Version: 2023.12.1 (402)
# Copyright (C) 2007 Free Software Foundation
# Sudarshan Gogoi, Soumen Bera, and Amit Chakraborty
analyze10xgenomCancer <- function(
    file_path,
    min_features = 0, # Set default minimum features 
    max_features = 25000, # Set default maximum features
    mt_threshold = 100, # Default mitochondrial content threshold
    dims = 1:10, # Default dimensions for PCA and UMAP
    resolution = 0.2 # Default clustering resolution
) {  
  # Load necessary libraries
  necessary_packages <- c("Seurat", "ggplot2", "dplyr", "reshape2", "cluster", "Hmisc", "tibble", "tidyr", "tidyverse", "hdf5r")
  lapply(necessary_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  })  
  # Load the dataset
  Cancer.sparse.m <- Read10X_h5(filename = file_path)  
  # Initialize the Seurat object
  Cancer.seurat.obj <- CreateSeuratObject(counts = Cancer.sparse.m,
                                          project = "Cancer",
                                          min.cells = 3,
                                          min.features = 100)
  # Quality control and filtering
  Cancer.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(Cancer.seurat.obj, pattern = "^MT-")
  # Filtration
  Cancer.seurat.obj <- subset(Cancer.seurat.obj,
                              subset = nFeature_RNA > min_features &
                                nFeature_RNA < max_features &
                                percent.mt < mt_threshold)
  # Data normalization
  Cancer.seurat.obj <- NormalizeData(Cancer.seurat.obj,
                                     normalization.method = "LogNormalize",
                                     scale.factor = 10000)
  # Identification of top variable features
  Cancer.seurat.obj <- FindVariableFeatures(Cancer.seurat.obj,
                                            selection.method = "vst",
                                            nfeatures = 2000)
  # Scaling data
  all.genes <- rownames(Cancer.seurat.obj)
  Cancer.seurat.obj <- ScaleData(Cancer.seurat.obj, features = all.genes)
  # Perform PCA
  Cancer.seurat.obj <- RunPCA(Cancer.seurat.obj, features = VariableFeatures(object = Cancer.seurat.obj))
  # Clustering
  Cancer.seurat.obj <- FindNeighbors(Cancer.seurat.obj, dims = dims)
  Cancer.seurat.obj <- FindClusters(Cancer.seurat.obj, resolution = resolution)
  # Perform UMAP
  Cancer.seurat.obj <- RunUMAP(Cancer.seurat.obj, dims = dims)
  # Plot the UMAP
  umap_plot <- DimPlot(Cancer.seurat.obj, reduction = "umap",pt.size = 1, label = TRUE)
  dev.new(height=7, width=7)
  print(umap_plot)
  # Cluster quality checking
  umap_emb <- as.data.frame(Cancer.seurat.obj@reductions$umap@cell.embeddings)
  cluster_info_name <- paste0("RNA_snn_res.", as.character(resolution))
  cluster_info <- Cancer.seurat.obj@meta.data[[cluster_info_name]] # accessing cluster info from meta.data
  if (!is.null(cluster_info)) {
    cluster_info <- as.numeric(cluster_info)
    if (!anyNA(cluster_info) && is.numeric(cluster_info)) {
      silhouette_scores <- silhouette(cluster_info, dist(umap_emb))
      print(silhouette_scores)
      sil_summary <- summary(silhouette_scores)
      mean_silhouette <- sil_summary$avg.width
      print(paste("Mean Silhouette Score:", mean_silhouette))
    } else {
      print("Cluster information is not numeric or contains missing values.")
    }
  } else {
    print("Cluster information not found.")
  }
  # Subfunction to find marker genes for each cluster
  Find_Marker_Genes <- function(seurat_object) {
    Highly_Expressed_genes_list <- list()
    cluster_ids <- levels(Idents(seurat_object))
    for (cluster_id in cluster_ids) {
      cluster_cells <- WhichCells(seurat_object, ident = cluster_id)
      cluster_expr_matrix <- GetAssayData(seurat_object, cells = cluster_cells, assay = "RNA")
      avg_expr_per_gene <- rowMeans(cluster_expr_matrix)
      Highly_Expressed_genes <- names(sort(avg_expr_per_gene, decreasing = TRUE)[1:1000])
      Highly_Expressed_genes_list[[paste0("Highly_Expressed_cluster", cluster_id)]] <- Highly_Expressed_genes
    }
    return(Highly_Expressed_genes_list)
  }
  # Cluster-wise highly expressed genes extraction
  Highly_Expressed_genes_list <- Find_Marker_Genes(Cancer.seurat.obj)
  # Gene correlation network construction
  # Individual cluster's gene expression data extraction
  # Extract the meta data for cells in clusters 0,1,2,3,4,5,6,7
  # Create an empty list to store metadata for each cluster
  cluster_meta_list <- list()
  # Loop through clusters from 0 to 7
  for (i in 0:7) {
    # Extract metadata for the current cluster
    cluster_meta <- Cancer.seurat.obj@meta.data[Cancer.seurat.obj@meta.data$seurat_clusters == i, ]
    # Assign the metadata to a variable
    assign(paste0("cluster_", i, "_meta"), cluster_meta)
    # Store the metadata in the list (optional)
    cluster_meta_list[[paste0("cluster_", i, "_meta")]] <- cluster_meta
  }
  # Get the cell IDs from the Seurat object
  cell_ids <- rownames(Cancer.seurat.obj@meta.data)
  # Add the cell IDs to the Seurat object's metadata data.frame
  Cancer.seurat.obj@meta.data$cell_id <- cell_ids
  # Create an empty list to store filtered correlation matrices for each cluster
  filtered_correlation_list <- list()
  # Loop through clusters from 0 to 7
  for (i in 0:7) {
    # Extract cell names for the current cluster
    cell_names <- rownames(get(paste0("cluster_", i, "_meta")))
    # Subset the Seurat object to include only the specified cells
    seurat_subset <- Cancer.seurat.obj[, Cancer.seurat.obj@meta.data$cell_id %in% cell_names]
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
  print(filtered_correlation_list)
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
  # Create a list to store Hausdorff distance matrices for thresholds 0.50, 0.53, 0.56, 0.58, 0.59, 0.62, 0.65, and 0.68
  Hausdorff_distance_matrices <- list()
  # Add threshold 0.50 to the list of thresholds
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
}

# Example usage:
analyze10xgenomCancer("file_path\\V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5")
                   
