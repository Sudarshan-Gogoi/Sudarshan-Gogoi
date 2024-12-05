# Stacked line plot of the training datasets for better visualisation of their norms patterns
#Load libraries
library(tidyr)
library(ggplot2)
library(cowplot)
# Load CancersNormData.csv in R studio (RStudio Version: 2023.12.1 (402))
data <- read.csv("D:/My PhD information folder first Paper/Data Analysis in R Modified/CancersNormData.csv")
# Stacked line plots of Breast cancer training datasets
data1 <- data[, c("Threshold", "Norm.B.a.", "Norm.B.b.")]
data_long <- gather(data1, key = "Variable", value = "Value", -Threshold)
filtered_data <- data_long %>%
  filter(Variable %in% c("Norm.B.a.", "Norm.B.b."))
p <- ggplot(filtered_data, aes(x = Threshold, y = Value, color = Variable, group = Variable)) +
  geom_area(fill = NA) +
  labs(
    x = "Threshold",
    y = "Value",
    color = "Variable"
  ) +
  scale_color_manual(
    values = c("Norm.B.a." = "blue", "Norm.B.b." = "red"),
    name = "Variables"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 20,face = "bold"),     # Increase axis title font size
    axis.text = element_text(size = 20,face = "bold"),      # Increase axis text font size
    legend.title = element_text(size = 20,face = "bold"),   # Increase legend title font size
    legend.text = element_text(size = 20,face = "bold"),    # Increase legend text font size
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    legend.spacing.y = unit(1, "cm")  # Add space between legend items in y-axis (adjust the value as needed)
  )
ggsave("stacked_line_plot_1_Breast.png", p, width = 8, height = 6, dpi = 1200)
# Stacked line plots of Colorectal cancer training datasets
data2 <- data[, c("Threshold", "Norm.C.a.", "Norm.C.b.")]
data_long <- gather(data2, key = "Variable", value = "Value", -Threshold)
filtered_data <- data_long %>%
  filter(Variable %in% c("Norm.C.a.", "Norm.C.b."))
p <- ggplot(filtered_data, aes(x = Threshold, y = Value, color = Variable, group = Variable)) +
  geom_area(fill = NA) +
  labs(
    x = "Threshold",
    y = "Value",
    color = "Variable"
  ) +
  scale_color_manual(
    values = c("Norm.C.a." = "blue", "Norm.C.b." = "red"),
    name = "Variables"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 20,face = "bold"),     
    axis.text = element_text(size = 20,face = "bold"),      
    legend.title = element_text(size = 20,face = "bold"),   
    legend.text = element_text(size = 20,face = "bold"),    
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    legend.spacing.y = unit(1, "cm")  # Add space between legend items in y-axis (adjust the value as needed)
  )
ggsave("stacked_line_plot_2_Colorectal.png", p, width = 8, height = 6, dpi = 1200)
# Stacked line plots of Lung cancer training datasets
data3 <- data[, c("Threshold", "Norm.L.a.", "Norm.L.b.")]
data_long <- gather(data3, key = "Variable", value = "Value", -Threshold)
filtered_data <- data_long %>%
  filter(Variable %in% c("Norm.L.a.", "Norm.L.b."))
p <- ggplot(filtered_data, aes(x = Threshold, y = Value, color = Variable, group = Variable)) +
  geom_area(fill = NA) +
  labs(
    x = "Threshold",
    y = "Value",
    color = "Variable"
  ) +
  scale_color_manual(
    values = c("Norm.L.a." = "blue", "Norm.L.b." = "red"),
    name = "Variables"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 20,face = "bold"),     
    axis.text = element_text(size = 20,face = "bold"),      
    legend.title = element_text(size = 20,face = "bold"),   
    legend.text = element_text(size = 20,face = "bold"),    
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    legend.spacing.y = unit(1, "cm")  # Add space between legend items in y-axis (adjust the value as needed)
  )
ggsave("stacked_line_plot_3_Lung.png", p, width = 8, height = 6, dpi = 1200)
# Stacked line plots of Ovarian cancer training datasets
data4 <- data[, c("Threshold", "Norm.O.a.", "Norm.O.b.")]
data_long <- gather(data4, key = "Variable", value = "Value", -Threshold)
filtered_data <- data_long %>%
  filter(Variable %in% c("Norm.O.a.", "Norm.O.b."))
p <- ggplot(filtered_data, aes(x = Threshold, y = Value, color = Variable, group = Variable)) +
  geom_area(fill = NA) +
  labs(
    x = "Threshold",
    y = "Value",
    color = "Variable"
  ) +
  scale_color_manual(
    values = c("Norm.O.a." = "blue", "Norm.O.b." = "red"),
    name = "Variables"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 20,face = "bold"),     
    axis.text = element_text(size = 20,face = "bold"),      
    legend.title = element_text(size = 20,face = "bold"),   
    legend.text = element_text(size = 20,face = "bold"),    
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    legend.spacing.y = unit(1, "cm")  # Add space between legend items in y-axis (adjust the value as needed)
  )
ggsave("stacked_line_plot_4_Ovarian.png", p, width = 8, height = 6, dpi = 1200)
