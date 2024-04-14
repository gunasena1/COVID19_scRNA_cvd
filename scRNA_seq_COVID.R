obj<-readRDS("local (2).rds")

head(obj@meta.data)

# Save metadata as a comma-separated text file
write.csv(obj@meta.data, file = "metadata.csv", quote = FALSE, row.names = TRUE)




#CD16Mono_data <- subset(obj, subset = Cell.group == "CD16+ Monocyte")
NK_data <- subset(obj, subset = Cell.group == "NK cell")

#TotalMono_data <- merge(CD16Mono_data, y = CD14Mono_data, add.cell.ids = c("CD16+ Monocyte", "CD14+ Monocyte"))

control_data <- subset(NK_data, subset = disease_general == "Healthy/Control")

severe_data <- subset(NK_data, subset = disease_general == "COVID-19 Severe/Late stage/Vent")


merged_data <- merge(control_data, y = severe_data, add.cell.ids = c("Healthy/Control", "COVID-19 Severe/Late stage/Vent"))
head(merged_data@meta.data)
write.csv(merged_data@meta.data, file = "mergeddata_metadata.csv", quote = FALSE, row.names = TRUE)


VlnPlot(merged_data,features = c("nCounts_RNA","nFeaturess_RNA","percent_mito"),pt.size = 0.2)
FeatureScatter(merged_data,feature1="nCounts_RNA",feature2="nFeaturess_RNA")
FeatureScatter(merged_data,feature1="nCounts_RNA",feature2="percent_mito")


mt_cutH <- 10
obj_unfiltered <- merged_data
merged_data <- subset(merged_data,subset = percent_mito < mt_cutH)
merged_data


merged_data <- FindVariableFeatures(merged_data)
merged_data <- ScaleData(merged_data,vars.to.regress = c("percent_mito"))


set.seed(0)
merged_data <- RunPCA(merged_data, npcs = 30, verbose = FALSE)
ElbowPlot(merged_data,ndims = 30)

numPC <- 10
maxPC <- numPC
merged_data <- FindNeighbors(merged_data, reduction = "pca", dims = 1:maxPC)


merged_data <- FindClusters(merged_data, resolution = 0.5)

merged_data <- RunUMAP(merged_data, reduction = "pca", dims = 1:maxPC)

DimPlot(merged_data,reduction = "umap")
DimPlot(merged_data, group.by = "disease_general")




merged_data <- SetIdent(merged_data,value = "disease_general")

clust.markers <- FindAllMarkers(merged_data, 
                                only.pos = FALSE,
                                min.pct = 0.25,
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.1,
                                logfc.threshold = 0.25)

head(clust.markers)

# Specify the file path where you want to save the CSV file
csv_file_path <- "clust.markers.csv"

# Save the sorted data frame as a CSV file
write.csv(clust.markers, file = csv_file_path, row.names = FALSE)



topG <- clust.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 200, wt = avg_log2FC)

# Replace "path/to/topG.csv" with the desired file path and name
csv_file_path <- "NK_topG.csv"

# Write the topG data frame to a CSV file
write.csv(topG, file = csv_file_path, row.names = FALSE)


heatmap_plot <-DoHeatmap(merged_data, features = topG$gene) + NoLegend()

pdf("_NKfinalheatmap_plot.pdf", width = 8, height = 6)  # You can adjust the width and height as per your preference
print(heatmap_plot)
dev.off()
#-------------------------Remove control------------------------------------

# Filter out samples labeled "Healthy/Control" from topG
topG_filtered <- topG %>%
  filter(cluster != "Healthy/Control")

# Print the first few rows of the filtered data frame
head(topG_filtered)
# Specify the file path where you want to save the CSV file
csv_file_path <- "topG_filtered.csv"





# Save the sorted data frame as a CSV file
write.csv(topG_filtered, file = csv_file_path, row.names = FALSE)

#if you edit the above file use the below code to reupload
topG_filtered <- read.csv("ALVAC new GO pathway DEGs switched neg pos.csv")



gene_list <- topG_filtered$gene

enrich_result <- enrichGO(gene = gene_list, 
                          OrgDb = org.Mm.eg.db, 
                          keyType = "ENSEMBL",
                          ont = "BP",  # Choose the ontology you are interested in, e.g., "BP" for Biological Process, "MF" for Molecular Function, or "CC" for Cellular Component.
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05,
                          universe = NULL,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = TRUE)


# Sort the 'enrich_result' data frame by 'p.adjust' values
enrich_result_sorted <- enrich_result[order(enrich_result$p.adjust), ]

# Specify the file path where you want to save the CSV file
csv_file_path <- "enrich_result_sorted_ALVAC switch neg pos.csv"

# Save the sorted data frame as a CSV file
write.csv(enrich_result_sorted, file = csv_file_path, row.names = FALSE)

#if you edit the above file use the below code to reupload
enrich_result_sorted <- read.csv("GO_pathways.csv")
# Create the bar plot
ggplot(enrich_result_sorted, aes(x = reorder(Description, p.adjust), y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Enriched GO Terms",
       x = "GO Terms",
       y = "Count",
       fill = "Adjusted p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "right") +
  scale_y_continuous(labels = comma)



ggsave("NK_ALVAC switch neg pos__enriched_go_terms_plot.png", width = 20, height = 10, units = "in")



#-----------------CNET------------------------

# Filter enriched pathways based on p-value
filtered_enrich_result <- enrich_result[enrich_result$pvalue < 0.05, ]






cnetplot(enrich_result, showCategory = 20)



ggsave("NK_SEVERE_BP__CNET.png", width = 40, height = 40, units = "in")
#------------------------------------
goplot(enrich_result, showCategory = 20, label = FALSE)
ggsave("NK_SEVERE_BP__GOnetwork.png", width = 20, height = 10, units = "in")



#----------------------------------------------------------------
# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)


eg <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP")

# Convert gene IDs to gene symbols within the enrichment result
eg_symbols <- setReadable(eg, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")

# Display the first few rows of the enrichment results with gene symbols
head(eg_symbols)

# Create a dot plot visualization of the enrichment results with gene symbols
dotplot(eg_symbols, showCategory = 20)

# Create a cluster network plot
cnetplot(eg_symbols, showCategory = 20)








# Assuming the column 'p.adjust' contains the adjusted p-values for enrichment
sorted_eg_symbols <- eg_symbols[order(eg_symbols$p.adjust), ]

# Display the sorted enrichment results
head(sorted_eg_symbols)

# Create a dot plot visualization of the sorted enrichment results
dotplot(sorted_eg_symbols)

# Create a cluster network plot for the sorted results
cnetplot(sorted_eg_symbols)







#------------------- Volcano plot-------------------------------------------------------------
clust.markers_filtered <- clust.markers %>%
  filter(cluster != "Healthy/Control")

# Specify the file path where you want to save the CSV file
csv_file_path <- "clust.markers_filtered_NK.csv"

# Save the sorted data frame as a CSV file
write.csv(clust.markers_filtered, file = csv_file_path, row.names = FALSE)

#if you edit the above file use the below code to reupload
clust.markers_filtered <- read.csv("Differential_expression_analysis_table (9).csv")


# Adjust the cutoff values as per your requirement
log2FC_cutoff <- 0.25
padj_cutoff <- 0.05

# Subset the significant genes based on log2FC and p-value cutoffs
significant_genes <- clust.markers_filtered %>% filter(abs(avg_log2FC) >= log2FC_cutoff & p_val_adj < padj_cutoff)
# Create the volcano plot with gene labels
volcano_plot <- ggplot(clust.markers_filtered, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = ifelse(avg_log2FC < 0, "red", "blue")), size = 4, alpha = 0.7) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), linetype = "dashed", color = "gray40") +
  #geom_text(data = significant_genes, aes(label = ensembl_gene_id), vjust = -0.5, hjust = 0, size = 3, color = "red") + # Adding gene labels for significant genes
  labs(title = "Volcano Plot", x = "Average Log2 Fold Change", y = "-log10(p-value)") +
  theme_minimal()
coord_cartesian(clip = 'off', xlim = c(-2, 500), ylim = c(-2, 500))  # Adjust the xlim and ylim values as needed


# Set the desired width and height for the plot
plot_width <- 80  # Adjust the width as needed
plot_height <- 50  # Adjust the height as needed

# Display the plot with the specified size
print(volcano_plot + theme(
  plot.title = element_text(size = 14),
  axis.title = element_text(size = 12),
  axis.text = element_text(size = 10),
  axis.text.x = element_text(angle = 45, hjust = 1)
) + 
  coord_cartesian(clip = 'off') +
  theme(
    aspect.ratio = plot_height/plot_width
  ))


tiff("volcano_plot_NK.tiff", width = 40, height = 25, units = "in", res = 300)
print(volcano_plot)
dev.off()





#------------------ENsemble ID to ENTREZ ID-------------------------------


convert_ensembl_to_entrez <- function(topG_filtered) {
  # Create a biomaRt connection to Ensembl database
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Query for the mapping between Ensembl and Entrez IDs
  mapping <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                   filters = "ensembl_gene_id",
                   values = topG$gene,
                   mart = mart)
  
  # Filter out rows with missing Entrez IDs
  topG <- topG[match(mapping$ensembl_gene_id, topG$gene, nomatch = 0), ]
  
  # Replace the "gene" column with the mapped Entrez IDs
  topG$gene <- mapping$entrezgene_id[match(topG$gene, mapping$ensembl_gene_id)]
  
  # Return the modified topG dataframe
  return(topG)
}

# Replace Ensembl IDs in topG with Entrez IDs
topG_filtered <- convert_ensembl_to_entrez(topG)

#------------KEGG Pathway analysis--------------------------

# Assuming you have already loaded the 'clust.markers' dataframe with marker genes

# Extract the list of marker genes
marker_genes <- topG_filtered$gene

# Convert Ensembl IDs in marker_genes to Entrez IDs
gene_entrez <- bitr(marker_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Extract the Entrez IDs of the marker genes
entrez_ids <- gene_entrez$ENTREZID

# Perform the KEGG pathway analysis
kegg_enrichment <- enrichKEGG(gene = entrez_ids, organism = 'hsa')

# Print the top enriched pathways
print(kegg_enrichment)

# Alternatively, you can save the results to a CSV file
write.csv(kegg_enrichment, file = "kegg_enrichment_results.csv", row.names = FALSE)

# Assuming you have already performed the KEGG pathway analysis and stored the results in 'kegg_enrichment' object

# Plot bar plot of the top enriched pathways

names(kegg_enrichment)

# Assuming the "p.adjust" column exists in kegg_enrichment
barplot(kegg_enrichment$p.adjust, las = 2, col = "lightblue", main = "Top Enriched KEGG Pathways")

# Dot plot of the top enriched pathways
keggplot <-dotplot(kegg_enrichment, showCategory = 15)

tiff("kegg_enrichment_NK.tiff", width = 16, height = 12, units = "in", res = 300)
print(keggplot)
dev.off()




# pathway thing
browseKEGG(kegg_enrichment, 'hsa04662')

# Specify the pathway ID you want to visualize
pathway_id <- "hsa05417"  # Replace with the desired KEGG pathway ID

# Create pathway enrichment map for the specified pathway
pathview(gene.data = kegg_enrichment, pathway.id = pathway_id, species = "hsa")



cnetplot(kegg_enrichment)




#-----------------------------------------------















#-------------Gene expression on UMAP-------------------------------

# Replace 'ENSG00000110848' with the actual gene ID you want to plot
gene_of_interest <- "ENSG00000111796"

# Find the index of the gene in the gene list of the Seurat object
gene_index <- which(rownames(merged_data) == gene_of_interest)

# Ensure that the gene is found in the dataset
if (length(gene_index) == 0) {
  stop(paste("Gene", gene_of_interest, "not found in the dataset."))
}

# Plot the gene expression using FeaturePlot
FeaturePlot(merged_data, features = gene_of_interest, pt.size = 0.2)




#------------Ridge Plots
data <- read.csv("ridgeplots ALVAC.csv")
top_gene <- data$gene
  
# Create RidgePlot for the combined subset
ridge_plot_combined <- RidgePlot(merged_data, features = top_gene, group.by = "disease_general")

# Save the combined RidgePlot as a PDF file
pdf("RidgePlot_NK.pdf", width = 16, height = 12)
print(ridge_plot_combined)
dev.off()


tiff("RidgePlot_NK.tiff", width = 16, height = 12, units = "in", res = 300)
print(ridge_plot_combined)
dev.off()

