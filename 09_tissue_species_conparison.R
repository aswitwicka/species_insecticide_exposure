# Uplad species
vcar_matrix <- read.csv("~/2023-species_comparisons/tables/star_counts/vcar_star_featureCounts.csv")
obic_matrix <- read.csv("~/2023-species_comparisons/tables/star_counts/FULL_SET_obic_star_featureCounts_notfilt.csv")
bter_matrix <- read.csv("~/2023-species_comparisons/tables/star_counts/bter_star_featureCounts.csv")
lser_matrix <- read_csv("~/2022-Lsericata_brainDGE/results_update/lser_FULL_STAR_counts_290324.csv")

colnames(vcar_matrix) <- gsub(".Aligned.out.sorted.bam", "", colnames(vcar_matrix))
colnames(obic_matrix) <- gsub(".Aligned.out.sorted.bam", "", colnames(obic_matrix))
colnames(lser_matrix) <- gsub(".Aligned.out.sorted.bam", "", colnames(lser_matrix))
colnames(lser_matrix) <- sub("_S.*", "", colnames(lser_matrix))
colnames(bter_matrix) <- gsub(".Aligned.out.sorted.bam", "", colnames(bter_matrix))

bter_meta <- read.csv("~/2023-species_comparisons/tables/metadata/bter_meta.csv")
vcar_meta <- read.csv("~/2023-species_comparisons/tables/metadata/vcar_meta.csv")
obic_meta <- read.csv("~/2023-species_comparisons/tables/metadata/obic_meta.csv")
lser_meta <-  read.csv("~/2023-species_comparisons/tables/metadata/lser_meta_update.csv")

bter_dds <- read_rds("~/2023-species_comparisons/DGE_pipelines/bter_dds_filtered_052024.rds")
vcar_dds <- read_rds("~/2023-species_comparisons/DGE_pipelines/vcar_dds_filtered_052024.rds")
lser_dds <- read_rds("~/2023-species_comparisons/DGE_pipelines/lser_dds_filtered_052024.rds")
obic_dds <- read_rds("~/2023-species_comparisons/DGE_pipelines/obic_dds_filtered_052024.rds")

# Subset samples in the raw matrices (only controls)
bter_g <- bter_matrix$X
vcar_g <- vcar_matrix$X
lser_g <- lser_matrix$...1

vcar_matrix_cont <- vcar_matrix[, grep("CONT", colnames(vcar_matrix))]
lser_matrix_cont <- lser_matrix[, grep("CONT", colnames(lser_matrix))]
bter_matrix_cont <- bter_matrix[, grep("CONT", colnames(bter_matrix))]

colnames(vcar_matrix_cont) <- gsub(".CONT", "-CONT", colnames(vcar_matrix_cont))
vcar_matrix_cont <- vcar_matrix_cont[, colnames(vcar_matrix_cont) %in% vcar_meta$SAMPLE_CODE]
colnames(bter_matrix_cont) <- gsub("\\.", "-", colnames(bter_matrix_cont))
bter_matrix_cont <- bter_matrix_cont[, colnames(bter_matrix_cont) %in% bter_meta$SAMPLE_NAME]
colnames(lser_matrix_cont) <- gsub("-", ".", colnames(lser_matrix_cont))
lser_matrix_cont <- lser_matrix_cont[, colnames(lser_matrix_cont) %in% colnames(lser_dds)]

lser_matrix_cont$X <- lser_g
vcar_matrix_cont$X <- vcar_g
bter_matrix_cont$X <- bter_g

obic_g <- obic_matrix$X
colnames(obic_matrix) <- gsub("\\.", "-", colnames(obic_matrix))
obic_matrix <- obic_matrix[, colnames(obic_matrix) %in% obic_meta$SAMPLE_SEQ]
sample_code_mapping <- setNames(obic_meta$SAMPLE_CODE, obic_meta$SAMPLE_SEQ)
colnames(obic_matrix) <- sample_code_mapping[colnames(obic_matrix)]

obic_matrix <- obic_matrix[, grep("cont", colnames(obic_matrix))]
obic_matrix$X <- obic_g

obic_matrix_cont <- obic_matrix

# Upload tissue ban files
# bamFiles <- list.files(pattern="*.bam",
#                        path = "~/2022-bter_tissues_chronic_exposure/results/2022-10-07-star/update/tmp/star_combined_run/")
# bamFiles <- bamFiles[grep("FCONT", bamFiles)]
# bamFiles <- paste("~/2022-bter_tissues_chronic_exposure/results/2022-10-07-star/update/tmp/star_combined_run", bamFiles, sep = "/")
# annotationFile <- "~/2023-species_comparisons/genomes/bter/bter.gtf"
# outDir <- ""
# 
# counts_bter <- featureCounts(
#   files = bamFiles,
#   annot.ext = annotationFile,
#   isGTFAnnotationFile = TRUE,
#   isPairedEnd = TRUE,
#   minMQS = 0,
# )
# 
# bter_tissues_matrix <- counts_bter$counts

muscle_matrix <- bter_tissues_matrix[, grep("LEG", colnames(bter_tissues_matrix))]
tubules_matrix <- bter_tissues_matrix[, grep("TUBULES", colnames(bter_tissues_matrix))]

colnames(muscle_matrix) <- gsub(".Aligned.out.sorted.bam", "", colnames(muscle_matrix))
colnames(tubules_matrix) <- gsub(".Aligned.out.sorted.bam", "", colnames(tubules_matrix))


# Subset genes that are conserved 1-1 orthologs: n = 2956
muscle_matrix_subset <- muscle_matrix[rownames(muscle_matrix) %in% ortholog_list$all_species_orthologs$Bter_gene, ]
tubules_matrix_subset <- tubules_matrix[rownames(tubules_matrix) %in% ortholog_list$all_species_orthologs$Bter_gene, ]
bter_matrix_subset <- bter_matrix_cont[bter_matrix_cont$X %in% ortholog_list$all_species_orthologs$Bter_gene, ]
vcar_matrix_cont$X <- gsub("gene-", "", vcar_matrix_cont$X) 
vcar_matrix_subset <- vcar_matrix_cont[vcar_matrix_cont$X %in% ortholog_list$all_species_orthologs$Vcar_gene, ]
lser_matrix_subset <- lser_matrix_cont[lser_matrix_cont$X %in% ortholog_list$all_species_orthologs$Lser_gene, ]
obic_matrix_subset <- obic_matrix_cont[obic_matrix_cont$X %in% ortholog_list$all_species_orthologs$Obic_gene, ]

# Change names to bter gene IDs
##Ensure 'X' column exists and not in rownames, and convert to data frames if necessary
muscle_matrix_subset <- as.data.frame(muscle_matrix_subset)  # Ensure it's a data frame
muscle_matrix_subset$X <- rownames(muscle_matrix_subset)  # Move rownames to a column 'X'
rownames(muscle_matrix_subset) <- NULL  # Remove rownames
tubules_matrix_subset <- as.data.frame(tubules_matrix_subset)  # Ensure it's a data frame
tubules_matrix_subset$X <- rownames(tubules_matrix_subset)  # Move rownames to a column 'X'
rownames(tubules_matrix_subset) <- NULL  # Remove rownames
bter_matrix_subset <- as.data.frame(bter_matrix_subset)  # Ensure it's a data frame
vcar_matrix_subset <- as.data.frame(vcar_matrix_subset)  # Ensure it's a data frame
lser_matrix_subset <- as.data.frame(lser_matrix_subset)  # Ensure it's a data frame
obic_matrix_subset <- as.data.frame(obic_matrix_subset) 

# Replace gene IDs
lser_gene_mapping <- setNames(ortholog_list$all_species_orthologs$Bter_gene, ortholog_list$all_species_orthologs$Lser_gene)
lser_matrix_subset$X <- lser_gene_mapping[lser_matrix_subset$X]  # Replace Lser gene IDs with Bter_gene IDs
vcar_gene_mapping <- setNames(ortholog_list$all_species_orthologs$Bter_gene, ortholog_list$all_species_orthologs$Vcar_gene)
vcar_matrix_subset$X <- vcar_gene_mapping[vcar_matrix_subset$X]  # Replace Vcar gene IDs with Bter_gene IDs
obic_gene_mapping <- setNames(ortholog_list$all_species_orthologs$Bter_gene, ortholog_list$all_species_orthologs$Obic_gene)
obic_matrix_subset$X <- obic_gene_mapping[obic_matrix_subset$X]  # Replace Obic gene IDs with Bter_gene IDs

### SVA batch correction
### Combine the data frames
### Run PCA 

# combined_matrix <- Reduce(function(x, y) merge(x, y, by = "X"), 
#                           list(muscle_matrix_scaled, tubules_matrix_scaled, bter_matrix_scaled, 
#                                vcar_matrix_scaled, lser_matrix_scaled, obic_matrix_scaled))

# OR:
combined_matrix <- Reduce(function(x, y) merge(x, y, by = "X"), 
                          list(muscle_matrix_subset, tubules_matrix_subset, bter_matrix_subset, 
                               vcar_matrix_subset, lser_matrix_subset, obic_matrix_subset))

scale_0_to_1 <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

# Apply the scaling function to each column (except the 'X' column)
combined_matrix[, -1] <- as.data.frame(lapply(combined_matrix[, -1], scale_0_to_1))



# Create metadata 
all_meta <- data.frame(sample = colnames(combined_matrix))
all_meta <- all_meta %>% 
  mutate(
    tissue = case_when(
      grepl("TUBULES", sample) ~ "tubules",
      grepl("LEG", sample) ~ "muscle",
      TRUE ~ "brain"
    ),
    species = case_when(
      grepl("FCONT", sample) ~ "bter",
      grepl("cont_", sample) ~ "obic",
      grepl("X", sample) ~ "lser",
      grepl("AU", sample) ~ "vcar",
      TRUE ~ "unknown"  # Use this to catch any unexpected cases
    )
  ) %>%
  mutate(batch = case_when(
    species == "bter" & tissue == "brain" ~ "batch_a",     # bter brain
    species == "obic"                     ~ "batch_b",     # obic species
    species == "lser"                     ~ "batch_c",     # lser species
    species == "vcar"                     ~ "batch_d",     # vcar species
    species == "bter" & (tissue == "tubules" | tissue == "muscle") ~ "batch_e",  # tubules and muscle (bter)
    TRUE ~ "unknown"  # Fallback in case there are any unhandled cases
  ))

all_meta <- all_meta %>% filter(sample != "X")

# Batch effect 
numeric_expression_data <- as.matrix(combined_matrix[, -1])

mod <- model.matrix(~ tissue, data = all_meta)

# Create a null model for the null hypothesis (no biological signal)
mod0 <- model.matrix(~ 1, data = all_meta)

# Estimate surrogate variables (SVs) using sva()
svobj <- sva(numeric_expression_data, mod, mod0)

batch_corrected <- limma::removeBatchEffect(numeric_expression_data, covariates = svobj$sv, design = mod0)

# batch_corrected <- ComBat(dat = batch_corrected, batch = all_meta$batch, mod = NULL)

data_for_pca <- t(batch_corrected)

# Perform PCA
pca_result <- prcomp(data_for_pca, scale. = TRUE)

# Summary of PCA (optional)
summary(pca_result)

# Create a data frame of PCA results (scores)
pca_scores <- as.data.frame(pca_result$x)

# Add sample names as a column to the PCA scores
pca_scores$sample <- rownames(pca_scores)

pca_scores <- pca_scores %>%
  mutate(
    tissue = case_when(
      grepl("TUBULES", sample) ~ "tubules",
      grepl("LEG", sample) ~ "muscle",
      TRUE ~ "brain"
    ),
    species = case_when(
      grepl("FCONT", sample) ~ "bter",
      grepl("cont_", sample) ~ "obic",
      grepl("X", sample) ~ "lser",
      grepl("AU", sample) ~ "vcar",
      TRUE ~ "unknown"  # Use this to catch any unexpected cases
    )
  )

# Plot the PCA using ggplot2 with color indicating species and shape indicating tissue
ggplot(pca_scores, aes(x = PC1, y = PC2, color = species, shape = tissue)) +
  geom_point(size = 3) +
  #geom_text(aes(label = sample), vjust = -1, size = 3) +
  labs(x = "PC1 [38.257%]", y = "PC2 [25.494%]") +
  scale_shape_manual(values = c(0,1,23))

