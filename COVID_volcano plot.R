# Load required libraries
library(readr)
library(dplyr)
library(ggplot2)

# Load the CSV file
data <- read_csv("data.csv")  # Replace "your_data.csv" with your file's name and path

# Assuming the first column contains sample names and the rest contain cell frequency data
# Separate the data into two groups (group1 and group2)
group1 <- data[, 2:5]  # Assuming the first group has 4 samples (columns 2 to 5)
group2 <- data[, 6:9]  # Assuming the second group has 4 samples (columns 6 to 9)

# Calculate log-fold change (log2) between group1 and group2
log_fold_change <- log2(rowMeans(group1) / rowMeans(group2))

# Perform a statistical test to calculate p-values (e.g., t-test, Wilcoxon test)
# Here's an example using a t-test; you may need to adjust the test depending on your data
p_values <- apply(data[, 2:9], 1, function(x) t.test(x[1:4], x[5:8])$p.value)

# Create a data frame for the volcano plot
volcano_data <- data.frame(LogFC = log_fold_change, PValue = p_values)

# Set a significance threshold (adjust as needed)
significance_threshold <- 0.10

# Create the volcano plot
ggplot(volcano_data, aes(x = LogFC, y = -log10(PValue))) +
  geom_point(aes(color = ifelse(PValue < significance_threshold, "Significant", "Not Significant")), size = 1) +
  geom_hline(yintercept = -log10(significance_threshold), linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-log10(P-Value)", title = "Volcano Plot") +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
  theme(legend.position = "none")

#-------------------------------------

# Create a column for upregulated and downregulated genes
volcano_data$Regulation <- ifelse(volcano_data$LogFC > 0, "Upregulated", "Downregulated",)

# Create the volcano plot with different colors for upregulated and downregulated significant data
ggplot(subset(volcano_data, PValue < significance_threshold), aes(x = LogFC, y = -log10(PValue))) +
  geom_point(aes(color = Regulation), size = 1) +
  geom_hline(yintercept = -log10(significance_threshold), linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-log10(P-Value)", title = "Volcano Plot") +
  scale_color_manual(values = c("Upregulated" = "blue", "Downregulated" = "red")) +
  theme(legend.title = element_blank())


#-------------------------------
# Adjust the cutoff values as per your requirement
log2FC_cutoff <- 0.10
padj_cutoff <- 0.05


ggplot(volcano_data, aes(x = LogFC, y = -log10(PValue))) +
  geom_point(aes(color = ifelse(LogFC < 0, "red", "blue")), size = 2, alpha = 0.7) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), linetype = "dashed", color = "gray40") +
  #geom_text(data = significant_genes, aes(label = ensembl_gene_id), vjust = -0.5, hjust = 0, size = 3, color = "red") + # Adding gene labels for significant genes
  labs(title = "Volcano Plot", x = "Average Log2 Fold Change", y = "-log10(p-value)") +
  theme_minimal()
coord_cartesian(clip = 'off', xlim = c(-2, 500), ylim = c(-2, 500))  # Adjust the xlim and ylim values as needed
