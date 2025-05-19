# This script performs WGCNA on vcar 

# Get a vector of genes recognised previously as expressed in species 
vcar_dds <- read_rds("~/2023-species_comparisons/DGE_pipelines/vcar_dds_filtered_052024.rds")

vcar_anno <- read.csv("~/2023-species_comparisons/tables/annotation/vcar_anno.csv")

vcar_all_genes <- rownames(vcar_dds)
vcar_all_genes <- gsub("gene-", "", vcar_all_genes)

# Upload metadata 
vcar_meta <- read.csv("~/2023-species_comparisons/tables/metadata/vcar_meta.csv")

# Check for outliers and general quality 
WGCNA::goodSamplesGenes(t(assay(vcar_dds)))
# $allOK
# [1] TRUE
# Hierarchical clustering 
htree <- hclust(dist(t(assay(vcar_dds))), method = "average")
plot(htree)

# Exclude sulfoxaflor 
# vcar_meta <- vcar_meta %>%
#   filter(CODE != "Sulfoxaflor")
# vcar_meta <- vcar_meta[which(vcar_meta$sample_name %in% colnames(vcar_dds)),]

# vcar_dds <- vcar_dds[,vcar_meta$sample_name]

# Perform variance stabilization
dds_norm <- vst(vcar_dds, blind = TRUE)

dds_norm <- assay(dds_norm) %>% as.data.frame()

#dds_norm <- dds_norm[,colnames(dds_norm) %in% vcar_meta$sample_name]

# Get normalized counts
vst_full <- dds_norm %>% 
  t()

colnames(vst_full) <- gsub("gene-", "", colnames(vst_full))
# Construct the network

# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 30, by = 1))
# Call the network topology analysis function
sft <- WGCNA::pickSoftThreshold(vst_full,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft.data <- sft$fitIndices
# Select maximum r-square and minimum mean-connectivity

# Visualisation to pick power
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

a1 / a2

# Convert matrix to numeric
vst_full[] <- sapply(vst_full, as.numeric)

soft_power <- 10
temp_cor <- cor
cor <- WGCNA::cor

# Memory estimate w.r.t blocksize
bwnet <- blockwiseModules(vst_full,
                          maxBlockSize = 16000, # Estimate for at least 32 RAM 
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

cor <- temp_cor

# Calculate module eigengenes (PC1 for each module)
module_eigengenes <- bwnet$MEs

# Get number of genes for each module
table(bwnet$colors)
length(table(bwnet$colors)) # 19
table(bwnet$colors) %>% sum() # 10480 - 6141
# black         blue        brown         cyan        green  greenyellow         grey 
# 130          563          507           26          285           27         6141 
# grey60    lightcyan   lightgreen      magenta midnightblue         pink       purple 
# 22           22           21           33           23           42           27 
# red       salmon          tan    turquoise       yellow 
# 199           27           27         2001          357 

# [1] BOTH: "MEred"  SULF: "MEblue"  "MEgreen"

as.data.frame(bwnet$colors)[(rownames(as.data.frame(bwnet$colors)) %in% (updown_list$vcar_sulf %>% filter(padj < 0.05) %>% pull(X))),] %>% table()
# black      blue     brown     green      grey      pink       red turquoise    yellow 
# 11       171        16       107       216         8       136        66         3 

as.data.frame(bwnet$colors)[(rownames(as.data.frame(bwnet$colors)) %in% (updown_list$vcar_clot %>% filter(padj < 0.05) %>% pull(X))),] %>% table()
# blue     brown     green      grey      pink       red turquoise 
# 2         1         6        28         1         9         5 


plotDendroAndColors(bwnet$dendrograms[[1]], bwnet$colors,
                    c("merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.08,
                    guideHang = 0.01)


# Module trait associations

# Subset only cont and clothianidin
vst_full <- vst_full[(grep("CLOT|CONT", rownames(vst_full))),]

# Create traits matrix - binarise categorical variables
traits <- as.integer(grepl("SULF|CLOT", rownames(vst_full))) # treatment
# 
# sub_traits <- str_extract(rownames(vst_full), "GC|_") # tissue
# 
# sub_traits.out <- binarizeCategoricalColumns(sub_traits,
#                                              includePairwise = FALSE,
#                                              includeLevelVsAll = TRUE,
#                                              minCount = 1)
# traits <- cbind(traits, sub_traits.out)
traits <- as.data.frame(traits)
rownames(traits) <- rownames(vst_full)

# Define numbers of genes and samples
nSamples <- nrow(vst_full)
nGenes <- ncol(vst_full)

# Subset only cont and sulfoxaflor
module_eigengenes <- module_eigengenes[(grep("CLOT|CONT", rownames(module_eigengenes))),]

module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# Visualize module-trait association as a heatmap
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')
colnames(heatmap.data)

CorLevelPlot::CorLevelPlot(heatmap.data,
                           x = names(heatmap.data)[20],
                           y = names(heatmap.data)[1:19],
                           col = c("#009E73", "#FFFFFF", "#D55E00"))

correlations <- sapply(colnames(heatmap.data)[-which(colnames(heatmap.data) == "traits")], function(col) {
  cor(heatmap.data[[col]], heatmap.data$traits, method="pearson")
})

# Create empty vectors for storing correlations and p-values
correlations <- numeric(length = length(heatmap.data) - 1)
pvalues <- numeric(length = length(heatmap.data) - 1)

for (i in colnames(heatmap.data)[-which(colnames(heatmap.data) == "traits")]) {
  test <- cor.test(heatmap.data[[i]], heatmap.data$traits, method="pearson")
  correlations[names(heatmap.data)[which(names(heatmap.data) == i)]] <- test$estimate
  pvalues[names(heatmap.data)[which(names(heatmap.data) == i)]] <- test$p.value
}

# Combine correlations and p-values into a data frame for better visualization
result <- data.frame(Correlation = correlations, Pvalue = pvalues)
p_values <- result %>% filter(Pvalue > 0) %>% pull(Pvalue)

# Adjust p-values using Benjamini-Hochberg procedure
adjusted_p_values <- p.adjust(p_values, method = "BH")
which(adjusted_p_values < 0.05)
names(heatmap.data)[which(adjusted_p_values < 0.05)]
# [1] "MEred"
adjusted_p_values[which(adjusted_p_values < 0.05)]
# [1] 0.04197147

# Sulfoxaflor 
vst_full <- dds_norm %>% 
  t()
colnames(vst_full) <- gsub("gene-", "", colnames(vst_full))
# Subset only cont and clothianidin
vst_full <- vst_full[(grep("SULF|CONT", rownames(vst_full))),]

# Create traits matrix - binarise categorical variables
traits <- as.integer(grepl("SULF|CLOT", rownames(vst_full))) # treatment

traits <- as.data.frame(traits)
rownames(traits) <- rownames(vst_full)

# Define numbers of genes and samples
nSamples <- nrow(vst_full)
nGenes <- ncol(vst_full)

# Subset only cont and sulfoxaflor
module_eigengenes <- bwnet$MEs
module_eigengenes <- module_eigengenes[(grep("SULF|CONT", rownames(module_eigengenes))),]

module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# Visualize module-trait association as a heatmap
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')
colnames(heatmap.data)

CorLevelPlot::CorLevelPlot(heatmap.data,
                           x = names(heatmap.data)[20],
                           y = names(heatmap.data)[1:19],
                           col = c("#009E73", "#FFFFFF", "#D55E00"))

correlations <- sapply(colnames(heatmap.data)[-which(colnames(heatmap.data) == "traits")], function(col) {
  cor(heatmap.data[[col]], heatmap.data$traits, method="pearson")
})

# Create empty vectors for storing correlations and p-values
correlations <- numeric(length = length(heatmap.data) - 1)
pvalues <- numeric(length = length(heatmap.data) - 1)

for (i in colnames(heatmap.data)[-which(colnames(heatmap.data) == "traits")]) {
  test <- cor.test(heatmap.data[[i]], heatmap.data$traits, method="pearson")
  correlations[names(heatmap.data)[which(names(heatmap.data) == i)]] <- test$estimate
  pvalues[names(heatmap.data)[which(names(heatmap.data) == i)]] <- test$p.value
}

# Combine correlations and p-values into a data frame for better visualization
result <- data.frame(Correlation = correlations, Pvalue = pvalues)
p_values <- result %>% filter(Pvalue > 0) %>% pull(Pvalue)

# Adjust p-values using Benjamini-Hochberg procedure
adjusted_p_values <- p.adjust(p_values, method = "BH")
which(adjusted_p_values < 0.05)
names(heatmap.data)[which(adjusted_p_values < 0.05)]
# [1] "MEred"   "MEblue"  "MEgreen"
adjusted_p_values[which(adjusted_p_values < 0.05)]
# [1] 8.313151e-05 3.296071e-03 3.296071e-03

module.gene.mapping <- as.data.frame(bwnet$colors)

vcar_red <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'red') %>% 
  rownames()
# 199 genes in the "MEred" module 
vcar_green <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'green') %>% 
  rownames()

keggSpecies(gsub("LOC", "", vcar_red), "vcd") %>% filter(pvalue < 0.05) %>% 
  ggplot(aes(x = -log10(pvalue), y = factor(Description, levels = (keggSpecies(gsub("LOC", "", vcar_red), "vcd") %>% filter(pvalue < 0.05) %>% arrange(pvalue) %>% pull(Description))))) + 
  geom_bar(stat = "identity", width = 0.5, alpha = 1) +
  labs(x="-log10q-value", y=" ") 

associated_annotation_red <- vcar_anno %>%
  dplyr::filter(Symbol %in% vcar_red) %>%
  dplyr::select(Symbol, description) %>% unique()

write.csv(associated_annotation_red, "~/2023-species_comparisons/2024_results/vcar_red_moduleWGCNA.csv")

# Plot the module

module_eigengenes <- bwnet$MEs

eigengene_df <- data.frame(Sample = rownames(module_eigengenes),
                           Eigengene = module_eigengenes$MEred) %>%
  dplyr::mutate(treatment = case_when(
    grepl("CONT", Sample) ~ "Control",
    grepl("CLOT", Sample) ~ "Clothianidin",
    grepl("SULF", Sample) ~ "Sulfoxaflor"
  )) 

# Plot the eigengene values
# MEred are downregulated 
group_means <- eigengene_df %>% 
  group_by(treatment) %>% 
  summarise(mean_eigengene = mean(Eigengene))
ggplot(eigengene_df, aes(x = Sample, y = Eigengene)) +
  geom_hline(data = group_means, aes(yintercept = mean_eigengene), colour = "#F7D240", linewidth = 1) +
  geom_point(aes(colour = treatment), size = 2) +
  xlab("Sample") +
  ylab("Eigengene Value") +
  theme(axis.text.x = element_blank()) +
  facet_grid(~factor(treatment), scales = "free", space = "free", drop = TRUE)


# Load required packages
library(igraph)

module_genes <- c(which(bwnet$colors == "red"))
# Calculate the adjacency matrix for the module
adjacency_matrix <- adjacency(vst_full[, module_genes], type = "signed", power = soft_power)
# Filter connections
# Plot node relathioship values
#hist(as.vector(adjacency_matrix), main="Distribution of relationship strength between nodes", xlab="Value", ylab="Frequency")
q95 <- quantile(as.vector(adjacency_matrix), 0.5)
threshold = q95 # set a threshold
adjacency_matrix_thresholded <- adjacency_matrix * (adjacency_matrix > threshold)
N <- 5 # top 4 connections
adjacency_matrix_topN <- apply(adjacency_matrix_thresholded, 1, function(row) {
  row[order(row, decreasing = TRUE)[-(1:N)]] <- 0
  return(row)
})
# Create a graph object
g <- graph_from_adjacency_matrix(adjacency_matrix_topN, mode = "undirected", weighted = TRUE, diag = FALSE)
# Plot the network
# Get edge weights
edge_weights <- E(g)$weight
# Apply min-max scaling to normalize edge weights
min_weight <- min(edge_weights)
max_weight <- max(edge_weights)
scaled_weights <- (edge_weights - min_weight) / (max_weight - min_weight)
# Select gene names
gene_names <- colnames(vst_full)[black_module_genes]
# Add the names as vertex attributes
V(g)$name <- gene_names

labels <- rep(NA, length(V(g)))
# Loop through the vertices and set the label if the gene name is in the selected list
for (i in 1:length(V(g))) {
  gene_name <- V(g)$name[i]
  if (gene_name %in% selected_genes) {
    labels[i] <- gene_name
  }
}
# Add the labels as a vertex attribute
V(g)$selected_label <- labels
# Find the vertices with zero degree (i.e., no connections)
isolated_vertices <- V(g)[degree(g) == 0]
# Delete these isolated vertices
g <- delete.vertices(g, isolated_vertices)
# Set layout
layout <- layout_nicely(g)

selected_genes <- updown_list$vcar_sulf %>% filter(padj < 0.05) %>% pull(X)
selected_genes2 <- updown_list$vcar_clot %>% filter(padj < 0.05) %>% pull(X)
selected_genes3 <- selected_genes2[selected_genes2 %in% selected_genes]

# node_colors <- ifelse(V(g)$selected_label %in% selected_genes2, "red", "blue")
node_colors <- ifelse(V(g)$selected_label %in% (selected_genes3), "#556B2F",
                      ifelse(V(g)$selected_label %in% selected_genes, "yellow",
                             ifelse(V(g)$selected_label %in% selected_genes2, "pink", "white")))

# Plot the graph with the selected labels for selected genes enriched for ATP-related genes
plot(g, layout = layout_nicely, #layout_with_lgl, 
     vertex.label = V(g)$selected_label, 
     vertex.size = 4, edge.width = scaled_weights * 3, 
     vertex.label.cex = 0.2, 
     vertex.color = node_colors)
# Plot the graph with vertex labels for all genes
plot(g, layout = layout, vertex.label = V(g)$name, vertex.size = 4, edge.width = scaled_weights * 4, vertex.label.cex = 0.4)



table(bwnet$colors[names((bwnet$colors)) %in% selected_genes2[(!(selected_genes2 %in% selected_genes))]])

vcar_deg <- read.csv("~/2023-species_comparisons/tables/deg_results/vcar_sulf_singleModel_31032024.csv") %>% 
  filter(padj < 0.05)
vcar_deg$X <- gsub("gene-", "", vcar_deg$X)
vcar_deg_clot <- read.csv("~/2023-species_comparisons/tables/deg_results/vcar_clot_17092023.csv") %>% 
  filter(padj < 0.05)
vcar_deg_clot$X <- gsub("gene-", "", vcar_deg_clot$X)

table(vcar_deg$X %in% wgcna_genes_clot_vcar) # 117 out of 721 / 213 genes in module (54.9% genes in the module de in sulfoxaflor)
table(vcar_deg_clot$X %in% wgcna_genes_clot_vcar)# 19 out of 59 / 8.22% deg under clot

selected_genes <- wgcna_genes_clot_vcar[wgcna_genes_clot_vcar %in% vcar_deg$X]
selected_genes2 <- wgcna_genes_clot_vcar[wgcna_genes_clot_vcar %in% vcar_deg_clot$X]

# 181 DEG out of the total 643 = 28%

vcar_anno %>% filter(Symbol %in% gene_names_mod[!(gene_names_mod %in% (vcar_deg$X[vcar_deg$X %in% gene_names_mod]))]) %>% dplyr::select(Symbol, description)
