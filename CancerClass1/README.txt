# Here we have provided the complete instructions to apply the code for analysis of scRNA-seq datasets of cancer obtained from 10x Genomics database by applying the methodology mentioned in our paper named "An integrated computational framework utilizing single-cell genomics for precise classification and prediction of multiple cancer types" by Sudarshan Gogoi, Soumen Bera and Amit Chakraborty
# There are total 3 scripts written in R: 1st script (file name: script1.R) contains step 1 to step 12 and some additional analysis, 2nd script (file name:script2.R) contains step 13, and the 3rd script (script3.R) contains plotting script of the stacked line plots showing norm patterns in training datasets. All steps 1 to 13, described in the Methodology section of the paper, are brifed below:
# For filtered data within the range of the sample size(2500-5000 cells) as given in the example, omit the step 2
# script1.R (example file used for illustration: "V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5" which is Breast cancer Dataset (b) as mentioned in the paper)
# Load the required libraries Seurat, tidyverse, data.table, Matrix, Cluster, Hmisc, reshape2, dplyr, tibble, tidyr, ggplot2, cowplot
#Step 1: Data preparation
# Load the dataset
# Convert the dataset into a seurat object with the criteria of minimum 3 cells and minimum 100 features 
#Step 2: Quality control and filtering
# Calculate the mitochondrial genes percentage in each cell
# Plot 3 violin plots showing density of the cells w.r.t number of genes, number of RNA molecules and number of mitochondrial genes percentages in the cells
# Plot a scatter plot showing the linear relationships between number of genes and number of RNA molecules
# Use these plots to filter the dataset and maintain the sample size(2500-5000 cells)
# Use Table 1 from the paper to select the minimum and maximum features threshold to filter the datasets used in the paper                                       
#Step 3: Data normalization
#Step 4: Top 2000 high variable features identification
#Step 5: Scaling data
#Step 6: Perform linear dimensionality reduction technique PCA to identify principal components 
# Select the number of principal components to keep for downstream analysis of the dataset by observing the elbow plot
# Use Table S1 from the paper to select the number of PCs for the datasets used in the paper
#Step 7: Clustering
# Fit a resolution to form 8 clusters in the dataset
# Use Table S1 to select the resolution parameters for the datasets used in the paper
# Plot a Dimensional plot to show 8 identified clusters in the dataset derived from PCA
#Step 8: Perform non-linear dimensionality reduction technique UMAP for better visualisation of the clusters
# Plot UMAP dimensional plot illustrating 8 identified clusters
# Calculate Mean Silhouette Score to check the quality of the clusters produced    
#Step 9: Cluster-wise Marker genes extraction
# Extract individual cluster's gene expression data
# Determine the top 1000 highly expressed genes as Marker genes for each cluster
#Step 10: Gene correlation network construction
# Subset individual cluster's gene expression data by keeping only marker genes    
# Construct gene correlation matrix from each of the subsetted individual cluster's gene expression data
# Construct gene correlation network of each cluster at 7 different gene correlation thresholds (0.50,0.53,0.56,0.59,0.62,0.65,0.68) in the form of adjacency list    
#Step 11: Hausdorff distance matrix and norm calculation
# Convert the adjacency lists into data frames
# Extract the correlation values column from the data frames of the adjacency lists
# Calculate the relative hausdorff distances between the network data frames and construct Hausdorff distance matrices for the gene correlation thresholds (0.50,0.53,0.56,0.59,0.62,0.65,0.68)
# Calculate Frobenius norm of each Hausdorff distance matrix and thus note down the norm values at increasing values of gene correlation thresholds(0.50,0.53,0.56,0.59,0.62,0.65 and 0.68) and get the norm dataset
# Save the norm dataset in CSV format
#Step 12: Construction of stacked line plots showing the trends in norm values at different gene correlation thresholds in the cancer dataset
# Save the plot in png format
# Additional Analysis
# (a)Construct the loaded dataset's gene correlation network of each cluster in the form of adjacency list at gene correlation threshold of 0.6
# (b)Calculation of number of genes and edges present in the networks of different clusters
# (c) Shannon index calculation of the networks on the basis of the number of genes and edges present in the networks
-----------------------------------------------------------------------------------------------------------------------------
# script2.R (example file used for illustration:"CancersNormData.csv")
# "CancersNormData.csv" contains norm values of all the training and validation datasets at different correlation thresholds. 
#Step13: Interval wise and overall normalized similarity percentage calculation between two cancer datasets' norm values
# Load "CancersNormData.csv" in R studio
# Calculate the maximum RMSE value
# Note down the similarity percentage between two similar cancer type training datasets' norm values
# Check the similarity percentage of a validation or test dataset with the training datasets(norm values).
# Take the norm values column of a test(validation) dataset (Example: Take Norm.B.c. or Norm.B.d. or Norm.C.c. or Norm.L.c. or Norm.O.c.) and compare its similarity percentages with the training datasets
# In the script2.R we showed similarity percentages of Breast cancer dataset (c) (validation dataset) with the training datasets by selecting the column Norm.B.c. from the loaded example file "CancersNormData.csv"
# Replace Norm.B.c with Norm.B.d. or Norm.C.c. or Norm.L.c. or Norm.O.c. to get other validation datasets' similarity percenatges
# Save the results to a CSV file
-----------------------------------------------------------------------------------------------------------------------------
# script3.R (example file used for illustration:"CancersNormData.csv")
# "CancersNormData.csv" contains norm values of all the training and validation datasets at different correlation thresholds. 
# Load "CancersNormData.csv" in R studio
# Plot the stacked line plots showing norm patterns in training datasets
-----------------------------------------------------------------------------------------------------------------------------
