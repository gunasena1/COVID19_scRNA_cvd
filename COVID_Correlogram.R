install.packages("corrplot")
library(corrplot)
# Replace 'your_data.csv' with your actual data file path
data <- read.csv("T cells and oxLDL Recovered.csv")

numeric_data <- data[, 2:ncol(data)]
# Replace 'numeric_column_names' with the names of your numeric columns


correlation_matrix <- cor(numeric_data)

corrplot(correlation_matrix, method = "color")


corrplot(
  correlation_matrix,
  method = "color",  # Color-based representation
  type = "upper",   # Show only the upper triangle of the matrix
  tl.col = "black", # Color of variable labels
  tl.srt = 45,      # Angle of variable labels
  diag = FALSE      # Exclude diagonal values (correlation with itself)
)


# Save the correlogram as a PNG image
png("T cells correlogram_control.png", width = 800, height = 800)
corrplot(correlation_matrix, method = "color")
dev.off()


# Assuming you have already calculated the correlation_matrix as in your previous code

# Specify the file path where you want to save the CSV file
output_file <- "correlation_matrix_T cells_Recovered.csv"

# Write the correlation matrix to a CSV file
write.csv(correlation_matrix, file = output_file, row.names = FALSE)

# Print a message indicating where the file has been saved
cat("Correlation matrix has been saved to:", output_file, "\n")


#---p values table--------

num_vars <- ncol(correlation_matrix)
p_values <- matrix(NA, nrow = num_vars, ncol = num_vars)

for (i in 1:num_vars) {
  for (j in 1:num_vars) {
    if (i == j) {
      p_values[i, j] <- NA  # Skip diagonal (correlation with itself)
    } else {
      p_values[i, j] <- cor.test(numeric_data[, i], numeric_data[, j])$p.value
    }
  }
}


rownames(p_values) <- colnames(correlation_matrix)
colnames(p_values) <- colnames(correlation_matrix)
p_values_df <- as.data.frame(p_values)


# Specify the file path where you want to save the CSV file
csv_file_path <- "p_values_table_T cells_Recovered.csv"

# Export the p_values_df table as a CSV file
write.csv(p_values_df, file = csv_file_path, row.names = FALSE)

#------star marks---------------------

# Create the correlogram
corrplot(correlation_matrix, method = "color")

# Annotate the correlogram with stars for significant p-values
significance_level <- 0.05  # Set your desired significance level
p_values_significant <- p_values < significance_level

# Add star marks to the correlogram
for (i in 1:num_vars) {
  for (j in 1:num_vars) {
    if (i != j && p_values_significant[i, j]) {
      text(i, j, "*", cex = 1.5, col = "red")  # You can customize the star appearance
    }
  }
}


install.packages("ggcorrplot")
library(ggcorrplot)

# Create a correlogram using ggcorrplot
corr_plot <- ggcorrplot(correlation_matrix, method = "color")

# Create a ggplot2 object from the ggcorrplot object
gg_corr <- ggplot_gtable(ggplot_build(corr_plot))

# Annotate the ggplot2 object with stars for significant p-values
significance_level <- 0.05  # Set your desired significance level
p_values_significant <- p_values < significance_level
rows <- as.vector(row(p_values_significant))
cols <- as.vector(col(p_values_significant))
gg_corr <- gg_corr + 
  geom_text(aes(x = cols - 0.5, y = rev(rows) - 0.5, label = ifelse(p_values_significant, "*", "")), color = "red")

# Save the correlogram as a PNG file
ggsave("correlogram_with_stars.png", gg_corr, width = 6, height = 6)