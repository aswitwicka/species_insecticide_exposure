# This script created dds object using DESeq2 library and count matrices. It also filters out genes with low expression levels. 

# Upload count matrices 
vcar_matrix <- read.csv("~/2023-species_comparisons/tables/star_counts/vcar_star_featureCounts.csv")
obic_matrix <- read.csv("~/2023-species_comparisons/tables/star_counts/obic_star_featureCounts.csv")
lser_matrix <- read.csv("~/2022-Lsericata_brainDGE/results_update/lser_NEW_STAR_counts_ALL_310324.csv")
bter_matrix <- read.csv("~/2023-species_comparisons/tables/star_counts/bter_star_featureCounts_052024.csv")

# Upload metadata 
bter_meta <- read.csv("~/2023-species_comparisons/tables/metadata/bter_meta.csv")
vcar_meta <- read.csv("~/2023-species_comparisons/tables/metadata/vcar_meta.csv")
obic_meta <- read.csv("~/2023-species_comparisons/tables/metadata/obic_meta.csv")
lser_meta <-  read.csv("~/2023-species_comparisons/tables/metadata/lser_meta_update.csv")

# Create DESeq2 objects 
# Bombus terrestris 
# Transform names
bter_meta$transformed_column <- gsub("-", ".", bter_meta$SAMPLE_NAME)
bter_meta$transformed_column <- paste0(bter_meta$transformed_column, ".Aligned.out.sorted.bam")
# Set the first column of bter_matrix as the row names
rownames(bter_matrix) <- bter_matrix[, 1]
# Then, remove the first column from the matrix
bter_matrix <- bter_matrix[, -1]
# Check
table(colnames(bter_matrix) %in% bter_meta$transformed_column)
# Find the matching indices
identical(colnames(bter_matrix), bter_meta$transformed_column)
match_indices <- match(colnames(bter_matrix), bter_meta$transformed_column)
identical(colnames(bter_matrix), bter_meta$transformed_column)

# Replace the colnames using the matched indices
colnames(bter_matrix) <- bter_meta$SAMPLE_CODE[match_indices]
identical(colnames(bter_matrix), bter_meta$SAMPLE_CODE)

bter_dds <- DESeq2::DESeqDataSetFromMatrix(bter_matrix,
                                           colData = bter_meta,
                                           design = ~  QUEENRIGHT + CODE # The factor of interest should be last 
                                           # update the model - preferably make it as simple as possible
)
bter_dds$CODE <- relevel(bter_dds$CODE, ref = "Control")
# Estimate the size factors
bter_dds <- DESeq2::estimateSizeFactors(bter_dds)

# Filter
bter_norm <- counts(bter_dds, normalized = TRUE, replaced = FALSE) %>%
  as.data.frame() # 12008
bter_norm <- filterByExpr(
  bter_norm, 
  group = bter_meta$CODE, 
  min.count = 10, min.total.count = 0, large.n = 100, min.prop = 0.5) 
bter_norm_filter <- bter_norm[bter_norm == TRUE] %>% names() # 10088
# Sulset
bter_dds <- bter_dds[bter_norm_filter, ] # 10088

# Check
colData(bter_dds)

saveRDS(bter_dds,
        file = "DGE_pipelines/bter_dds_filtered_052024.rds"
)

# Osmia bicornis
# Transform names
obic_meta$transformed_column <- gsub("-", ".", obic_meta$SAMPLE_SEQ)
obic_meta$transformed_column <- paste0(obic_meta$transformed_column, ".Aligned.out.sorted.bam")
# Set the first column of bter_matrix as the rownames
rownames(obic_matrix) <- obic_matrix[, 1]
# Then, remove the first column from the matrix
obic_matrix <- obic_matrix[, -1]
# Filter the data to just keep sulf and clot 
obic_matrix <- obic_matrix[, colnames(obic_matrix) %in% obic_meta$transformed_column]
# Check
table(colnames(obic_matrix) %in% obic_meta$transformed_column)
table(obic_meta$transformed_column %in% colnames(obic_matrix))
# Find the matching indices
identical(colnames(obic_matrix), obic_meta$transformed_column)
# Match 
obic_matrix <- obic_matrix[, (match(obic_meta$transformed_column, colnames(obic_matrix)))]

# Replace the colnames using the matched indices
match_indices <- match(colnames(obic_matrix), obic_meta$transformed_column)
colnames(obic_matrix) <- obic_meta$SAMPLE_CODE[match_indices]
identical(colnames(obic_matrix), obic_meta$SAMPLE_CODE)
# 
obic_dds <- DESeq2::DESeqDataSetFromMatrix(obic_matrix,
                                           colData = obic_meta,
                                           design = ~  CODE # The factor of interest should be last 
                                           # update the model - preferably make it as simple as possible
)
obic_dds$CODE <- relevel(obic_dds$CODE, ref = "Control")
# Estimate the size factors
obic_dds <- DESeq2::estimateSizeFactors(obic_dds)

identical(rownames(colData(obic_dds)), colData(obic_dds)$SAMPLE_CODE)

# Filter
obic_norm <- counts(obic_dds, normalized = TRUE, replaced = FALSE) %>%
  as.data.frame() # 12630
obic_norm <- filterByExpr(
  obic_norm,
  group = obic_meta$CODE,
  min.count = 10, min.total.count = 0, large.n = 100, min.prop = 0.5)
obic_norm_filter <- obic_norm[obic_norm == TRUE] %>% names() # 10205
# Subset
obic_dds <- obic_dds[obic_norm_filter, ] # 10163
# Check
colData(obic_dds)

# Save 
saveRDS(obic_dds,
        file = "DGE_pipelines/obic_dds_filtered_052024.rds"
)

# Lucilia sericata
lser_meta$CODE <- as.factor(lser_meta$CODE)
rownames(lser_matrix) <- lser_matrix[, 1]
lser_matrix <- lser_matrix[, -1]

# Transform
colnames(lser_matrix) <- gsub("^X|_.*$", "", colnames(lser_matrix))
lser_meta$transformed_column <- gsub("_", ".", lser_meta$SAMPLE_CODE)

table(lser_meta$transformed_column %in% colnames(lser_matrix))
table(colnames(lser_matrix) %in% lser_meta$transformed_column)

# Subset
lser_matrix <- lser_matrix[,which(colnames(lser_matrix) %in% lser_meta$transformed_column)]

# Sort 
identical(colnames(lser_matrix), lser_meta$transformed_column)
match_indices <- match(lser_meta$transformed_column, colnames(lser_matrix))
lser_matrix <- lser_matrix[, match_indices]
identical(colnames(lser_matrix), lser_meta$transformed_column)
#
lser_dds <- DESeq2::DESeqDataSetFromMatrix(as.matrix(lser_matrix),
                                                  colData = lser_meta,
                                                  design = ~  CODE # The factor of interest should be last 
                                                  # update the model - preferably make it as simple as possible
)
lser_dds$treatment <- relevel(lser_dds$CODE, ref = "Control")
# Estimate the size factors
lser_dds <- DESeq2::estimateSizeFactors(lser_dds)

lser_norm <- counts(lser_dds, normalized = TRUE, replaced = FALSE) %>%
  as.data.frame() # 16436
lser_norm <- filterByExpr(
  lser_norm, 
  group = lser_meta$CODE, 
  min.count = 5, min.total.count = 0, large.n = 100, min.prop = 0.35) 
lser_norm_filter <- lser_norm[lser_norm == TRUE] %>% names() # 11968
# Subset
lser_dds <- lser_dds[lser_norm_filter, ]
# Check 
colData(lser_dds)
# Save 
saveRDS(lser_dds,
        file = "DGE_pipelines/lser_dds_filtered_052024.rds"
)

# Vanessa cardui 
# Transform names
vcar_meta$transformed_column <- gsub("-", ".", vcar_meta$SAMPLE_CODE)
vcar_meta$transformed_column <- paste0(vcar_meta$transformed_column, ".Aligned.out.sorted.bam")
# Set the first column of bter_matrix as the rownames
rownames(vcar_matrix) <- vcar_matrix[, 1]
# Then, remove the first column from the matrix
vcar_matrix <- vcar_matrix[, -1]
# Filter the data to just keep sulf and clot 
vcar_matrix <- vcar_matrix[, colnames(vcar_matrix) %in% vcar_meta$transformed_column]

# Check
table(colnames(vcar_matrix) %in% vcar_meta$transformed_column)
table(vcar_meta$transformed_column %in% colnames(vcar_matrix))
# Find the matching indices
identical(colnames(vcar_matrix), vcar_meta$transformed_column)
# Match 
vcar_matrix <- vcar_matrix[, (match(vcar_meta$transformed_column, colnames(vcar_matrix)))]
match_indices <- match(colnames(vcar_matrix), vcar_meta$transformed_column)
colnames(vcar_matrix) <- vcar_meta$SAMPLE_CODE[match_indices]
identical(colnames(vcar_matrix), vcar_meta$SAMPLE_CODE)
# 
vcar_dds <- DESeq2::DESeqDataSetFromMatrix(vcar_matrix,
                                           colData = vcar_meta,
                                           design = ~  CODE
)
vcar_dds$CODE <- relevel(vcar_dds$CODE, ref = "Control")
# Estimate the size factors
vcar_dds <- DESeq2::estimateSizeFactors(vcar_dds)

# Filter
vcar_norm <- counts(vcar_dds, normalized = TRUE, replaced = FALSE) %>%
  as.data.frame() # 15000
vcar_norm <- filterByExpr(
  vcar_norm,
  group = vcar_meta$CODE,
  min.count = 10, min.total.count = 0, large.n = 100, min.prop = 0.35)
vcar_norm_filter <- vcar_norm[vcar_norm == TRUE] %>% names() # 14480
# Sulset
vcar_dds <- vcar_dds[vcar_norm_filter, ]

saveRDS(vcar_dds,
        file = "DGE_pipelines/vcar_dds_filtered_052024.rds"
)
