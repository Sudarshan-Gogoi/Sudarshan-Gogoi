
#Step13: Interval wise and overall normalized similarity percentage calculation between two cancer datasets' norm values

# First, checked the normalized similarity percentage between two similar cancer type training datasets' norm values
# Second, Checked the normalized similarity percentage of a validation or test dataset with the training datasets(norm values)

library(dplyr)
# Load CancersNormData.csv in R studio
data <- read.csv("D:/My PhD information folder first Paper/Data Analysis in R Modified/FinalNorms/FinalSript/CancersNormData.csv")
# Select all pairs of columns for comparison
column_pairs <- combn(names(data)[2:ncol(data)], 2, simplify = TRUE)
# Create an empty data frame to store results
result_df <- data.frame()
# Initialize vector to store maximum RMSE values
max_rmse_values <- numeric()
# Loop through each column pair
for (i in 1:ncol(column_pairs)) {
  col1 <- column_pairs[1, i]
  col2 <- column_pairs[2, i]
  # Define intervals based on Threshold column values
  threshold_intervals <- unique(sort(data$Threshold))
  intervals <- lapply(1:(length(threshold_intervals) - 1), function(i) c(threshold_intervals[i], threshold_intervals[i + 1]))
  # Initialize vectors to store RMSE and similarity percentages for intervals
  rmse_values <- numeric(length(intervals))
  similarity_percentages <- numeric(length(intervals))
  # Calculate RMSE and similarity percentage for each interval
  for (i in seq_along(intervals)) {
    start_threshold <- intervals[[i]][1]
    end_threshold <- intervals[[i]][2]
    # Extract data within the threshold interval for both datasets
    interval_data_plot1 <- data[data$Threshold >= start_threshold & data$Threshold <= end_threshold, col1]
    interval_data_plot2 <- data[data$Threshold >= start_threshold & data$Threshold <= end_threshold, col2]
    # Calculate RMSE
    rmse_value <- sqrt(mean((interval_data_plot1 - interval_data_plot2)^2))
    rmse_values[i] <- rmse_value
    # Calculate normalized similarity percentage (with check for division by zero)
    max_possible_rmse_value <- 6  # assuming max possible RMSE is 6 for each interval
    if (max_possible_rmse_value != 0) {
      similarity_percentage <- 100 * (1 - (rmse_value / max_possible_rmse_value))
    } else {
      similarity_percentage <- NA
    }
    similarity_percentages[i] <- similarity_percentage
  }
  
  # Calculate the overall similarity percentage (average of interval similarities)
  overall_similarity_percentage <- mean(similarity_percentages, na.rm = TRUE)
  # Find the maximum RMSE value for the current column pair
  max_rmse_value <- max(rmse_values, na.rm = TRUE)
  # Store the maximum RMSE value in the vector
  max_rmse_values <- c(max_rmse_values, max_rmse_value)
  # Create a data frame with results for the current column pair
  result <- data.frame(
    Column1 = col1,
    Column2 = col2,
    Interval_1 = similarity_percentages[1],
    Interval_2 = similarity_percentages[2],
    Interval_3 = similarity_percentages[3],
    Interval_4 = similarity_percentages[4],
    Interval_5 = similarity_percentages[5],
    Interval_6 = similarity_percentages[6],
    Overall_Similarity_Percentage = overall_similarity_percentage,
    Max_RMSE = max_rmse_value
  )
  # Bind the result to the main data frame
  result_df <- rbind(result_df, result)
}

# Filter the result_df to determine the interval wise and overall similarity percentages between two same type of training datasets
Similarity_results_1 <- result_df %>%
  filter((Column1 == "Norm.B.a." & Column2 == "Norm.B.b.") |
           (Column1 == "Norm.C.a." & Column2 == "Norm.C.b.") |
           (Column1 == "Norm.L.a." & Column2 == "Norm.L.b.") |
           (Column1 == "Norm.O.a." & Column2 == "Norm.O.b."))
# Round off the similarity percentage values to three decimal places
Similarity_results_1 <- Similarity_results_1 %>%
  mutate_at(vars(starts_with("Interval")), ~paste0(round(., 3), "%"))

# Print the rounded results
print(Similarity_results_1)

# Take the norm values column of a test(random) dataset and compare its similarity percentages with the training datasets
# For example, here we filter the result_df to check the similarity percentage of Breast cancer dataset 3(validation dataset) with the training datasets by selecting the column Norm.B.c. 
Similarity_results_2 <- result_df %>%
  filter((Column1 == "Norm.B.c." & Column2 %in% c("Norm.B.a.", "Norm.B.b.", "Norm.L.a.", "Norm.L.b.", "Norm.C.a.", "Norm.C.b.", "Norm.O.a.", "Norm.O.b.")) |   # Replace the Norm.B.c. with the required 
           (Column2 == "Norm.B.c." & Column1 %in% c("Norm.B.a.", "Norm.B.b.", "Norm.L.a.", "Norm.L.b.", "Norm.C.a.", "Norm.C.b.", "Norm.O.a.", "Norm.O.b."))) # test data Norm column like Norm.B.d., Norm.L.c. etc 
                                                                                                                                            
# Round off the similarity percentage values to three decimal places
Similarity_results_2$Interval_1 <- paste0(round(Similarity_results_2$Interval_1, 3), "%")
Similarity_results_2$Interval_2 <- paste0(round(Similarity_results_2$Interval_2, 3), "%")
Similarity_results_2$Interval_3 <- paste0(round(Similarity_results_2$Interval_3, 3), "%")
Similarity_results_2$Interval_4 <- paste0(round(Similarity_results_2$Interval_4, 3), "%")
Similarity_results_2$Interval_5 <- paste0(round(Similarity_results_2$Interval_5, 3), "%")
Similarity_results_2$Interval_6 <- paste0(round(Similarity_results_2$Interval_6, 3), "%")
Similarity_results_2$Overall_Similarity_Percentage <- paste0(round(Similarity_results_2$Overall_Similarity_Percentage, 3), "%")

# Print the rounded results
print(Similarity_results_2)
# Save the results to a CSV file
write.csv(Similarity_results_2, "D:/My Folder/Similarity_Percentages.csv", row.names = FALSE)

