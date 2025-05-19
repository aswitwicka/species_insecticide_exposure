# Iterate 5 samples of Bombus terrestris (clot exposed): 

bter_meta <- read.csv("~/2023-species_comparisons/tables/metadata/bter_meta.csv")
bter_matrix <- read.csv("~/2023-species_comparisons/tables/star_counts/bter_star_featureCounts_052024.csv")

rownames(bter_matrix) <- bter_matrix[, 1]
bter_matrix <- bter_matrix[, -1]

# Remove the prefix "GC.AW.9791." and the suffix "_SXX.Aligned.out.sorted.bam"
colnames(bter_matrix) <- gsub("GC.AW.9791.", "", colnames(bter_matrix))
colnames(bter_matrix) <- gsub("_S[0-9]+\\.Aligned\\.out\\.sorted\\.bam", "", colnames(bter_matrix))
# Replace the dot after a number with an underscore
colnames(bter_matrix) <- gsub("(\\d)\\.", "\\1_", colnames(bter_matrix))

table(bter_meta$SAMPLE_CODE %in% colnames(bter_matrix))
identical(bter_meta$SAMPLE_CODE, colnames(bter_matrix))


resulting_counts <- vector("numeric", 100) 

set.seed(42)

for(i in 1:100) {
  # Splitting the data into CLOT and CONT groups
  clot_files <- grep("FCLOT", bter_meta$SAMPLE_NAME, value = TRUE)
  cont_files <- grep("CONT", bter_meta$SAMPLE_NAME, value = TRUE)
  # Extract numbers from the filenames
  get_number <- function(filename) {
    as.numeric(gsub("^.*-([0-9]+)-.*$", "\\1", filename))
  }
  # Randomly select 5 CLOT
  selected_clot <- sample(clot_files, 5)
  selected_numbers <- sapply(selected_clot, get_number)
  # Exclude CONT files with the selected numbers
  available_cont <- cont_files[!sapply(cont_files, get_number) %in% selected_numbers]
  # Randomly select 5 CONT
  selected_cont <- sample(available_cont, 5)
  # Combine 
  bter_subset <- c(selected_clot, selected_cont)
  # Subset data frame 
  bter_meta_Subset <- bter_meta %>% filter(SAMPLE_NAME %in% bter_subset)
  # Subset matrix 
  bter_matrix_Subset <- bter_matrix[, colnames(bter_matrix) %in% bter_meta_Subset$SAMPLE_CODE]
  
  print("Are the samples ordered:")
  print(identical(bter_meta_Subset$SAMPLE_CODE, colnames(bter_matrix_Subset)))
  
  bter_dds_Subset <- DESeq2::DESeqDataSetFromMatrix(bter_matrix_Subset,
                                                    colData = bter_meta_Subset,
                                                    design = ~  CODE 
  )
  bter_dds_Subset$CODE <- relevel(bter_dds_Subset$CODE, ref = "Control")
  # Estimate the size factors
  bter_dds_Subset <- DESeq2::estimateSizeFactors(bter_dds_Subset)
  
  # Filter
  bter_norm <- counts(bter_dds_Subset, normalized = TRUE, replaced = FALSE) %>%
    as.data.frame()
  bter_norm <- filterByExpr(
    bter_norm, 
    group = bter_meta$CODE, 
    min.count = 10, min.total.count = 0, large.n = 100, min.prop = 0.5) 
  bter_norm_filter <- bter_norm[bter_norm == TRUE] %>% names() 
  print("No. genes considered in the iteration:")
  print(length(bter_norm_filter))
  bter_dds_Subset <- bter_dds_Subset[bter_norm_filter, ]
  
  # Wald 
  bter_dds_Subset_wald <- DESeq2::DESeq(bter_dds_Subset, test = "Wald", betaPrior = FALSE)
  ## Explore the results 
  DESeq2::resultsNames(bter_dds_Subset_wald)
  # Clot
  bter_Subset_wald_clot <- results(bter_dds_Subset_wald,
                                   # independentFiltering = FALSE,
                                   pAdjustMethod = "BH",
                                   contrast = c(
                                     "CODE",
                                     "Clothianidin",
                                     "Control"
                                   )
  ) %>% as.data.frame()
  
  resulting_counts[i] <- bter_Subset_wald_clot %>% filter(padj < 0.05) %>% nrow()
  print("No DEGs in the iteration:")
  print(resulting_counts[i])
}

p1 <- ggplot() +
  geom_histogram(mapping=aes(x=resulting_counts), binwidth = 95, colour = "#8B8989") +
  geom_vline(xintercept=225, lty=2, color="#00CD66", linewidth=1) +
  theme(text=element_text(size=12)) +
  xlab("No DEGs") +
  ylab("Count")

# Now the same but the samples are paired according to the colony effect: 
resulting_counts_paired <- vector("numeric", 100) 

for(i in 1:100) {
  print("Iteration:")
  print(i)
  # Splitting the data into CLOT and CONT 
  clot_files <- grep("FCLOT", bter_meta$SAMPLE_NAME, value = TRUE)
  cont_files <- grep("CONT", bter_meta$SAMPLE_NAME, value = TRUE)
  # Extract numbers 
  get_number <- function(filename) {
    as.numeric(gsub("^.*-([0-9]+)-.*$", "\\1", filename))
  }
  # Randomly select 5 CLOT
  selected_clot <- sample(clot_files, 5)
  selected_numbers <- sapply(selected_clot, get_number)
  # Exclude already selected CONT 
  available_cont <- cont_files[sapply(cont_files, get_number) %in% selected_numbers]
  # Randomly select 5 available CONT
  selected_cont <- sample(available_cont, 5)
  # Combine 
  bter_subset <- c(selected_clot, selected_cont)
  # Subset data frame 
  bter_meta_Subset <- bter_meta %>% filter(SAMPLE_NAME %in% bter_subset)
  # Subset matrix 
  bter_matrix_Subset <- bter_matrix[, colnames(bter_matrix) %in% bter_meta_Subset$SAMPLE_CODE]
  
  print("Are the samples ordered:")
  print(identical(bter_meta_Subset$SAMPLE_CODE, colnames(bter_matrix_Subset)))
  
  bter_dds_Subset <- DESeq2::DESeqDataSetFromMatrix(bter_matrix_Subset,
                                                    colData = bter_meta_Subset,
                                                    design = ~  QUEENRIGHT + CODE 
  )
  bter_dds_Subset$CODE <- relevel(bter_dds_Subset$CODE, ref = "Control")
  # Estimate the size factors
  bter_dds_Subset <- DESeq2::estimateSizeFactors(bter_dds_Subset)
  
  # Filter
  bter_norm <- counts(bter_dds_Subset, normalized = TRUE, replaced = FALSE) %>%
    as.data.frame()
  bter_norm <- filterByExpr(
    bter_norm, 
    group = bter_meta$CODE, 
    min.count = 10, min.total.count = 0, large.n = 100, min.prop = 0.5) 
  bter_norm_filter <- bter_norm[bter_norm == TRUE] %>% names() 
  print("No. genes considered in the iteration:")
  print(length(bter_norm_filter))
  bter_dds_Subset <- bter_dds_Subset[bter_norm_filter, ]
  
  # Wald 
  bter_dds_Subset_wald <- DESeq2::DESeq(bter_dds_Subset, test = "Wald", betaPrior = FALSE)
  ## Explore the results 
  DESeq2::resultsNames(bter_dds_Subset_wald)
  # Clot
  bter_Subset_wald_clot <- results(bter_dds_Subset_wald,
                                   # independentFiltering = FALSE,
                                   pAdjustMethod = "BH",
                                   contrast = c(
                                     "CODE",
                                     "Clothianidin",
                                     "Control"
                                   )
  ) %>% as.data.frame()
  
  resulting_counts_paired[i] <- bter_Subset_wald_clot %>% filter(padj < 0.05) %>% nrow()
  print("No DEGs in the iteration:")
  print(resulting_counts_paired[i])
}

p2 <- ggplot() +
  geom_histogram(mapping=aes(x=resulting_counts_paired), binwidth = 50, colour = "#8B8989") +
  geom_vline(xintercept=543, lty=2, color="#EE5C42", linewidth=1) +
  theme(text=element_text(size=12)) +
  xlab("No DEGs") +
  ylab("Count")


p1 + p2

max(resulting_counts)
max(resulting_counts_paired)

mean(resulting_counts)
mean(resulting_counts_paired)

median(resulting_counts)
median(resulting_counts_paired)


# Calculate Empirical Cumulative Distribution Function (ECDFs)
ecdf_dist1 <- ecdf(resulting_counts)
ecdf_dist2 <- ecdf(resulting_counts_paired)

# Get the ECDF values at 2186
ecdf_value1 <- ecdf_dist1(2186)
ecdf_value2 <- ecdf_dist2(2186)

cat("ECDF value for distribution 1 at 2186:", ecdf_value1, "\n")
cat("ECDF value for distribution 2 at 2186:", ecdf_value2, "\n")

# Permutation test
set.seed(42) # For reproducibility
n_perm <- 10000
perm_diff <- numeric(n_perm)

for (i in 1:n_perm) {
  perm_sample <- sample(c(resulting_counts, resulting_counts_paired))
  perm_dist1 <- perm_sample[1:length(resulting_counts)]
  perm_dist2 <- perm_sample[(length(resulting_counts) + 1):length(perm_sample)]
  perm_diff[i] <- ecdf(perm_dist1)(2186) - ecdf(perm_dist2)(2186)
}

observed_diff <- ecdf_value1 - ecdf_value2
p_value <- mean(abs(perm_diff) >= abs(observed_diff))

cat("Observed difference in ECDF values at 2186:", observed_diff, "\n")
cat("Permutation test p-value:", p_value, "\n")

# Observed difference in ECDF values at 2186: -0.03 
# Permutation test p-value: 0.247 

#
summary(resulting_counts)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.00   39.25  225.50  527.96  737.00 3477.00 
summary(resulting_counts_paired)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 17.0   105.0   543.0   661.1  1214.5  2147.0
# Calculate skewness
library(e1071)
skewness(resulting_counts)
# [1] 1.990363 Highly positively skewed distribution
skewness(resulting_counts_paired)
# [1] 0.5320319 Moderately positively skewed

qqnorm(resulting_counts, main = "Q-Q Plot Distribution 1")
qqline(resulting_counts)
qqnorm(resulting_counts_paired, main = "Q-Q Plot Distribution 2")
qqline(resulting_counts_paired)

# Asymptotic two-sample Kolmogorov-Smirnov test
ks.test(resulting_counts, resulting_counts_paired)
# data:  resulting_counts and resulting_counts_paired
# D = 0.21, p-value = 0.02431
# alternative hypothesis: two-sided


proportion_above_2186_counts <- sum(resulting_counts > 2186) / length(resulting_counts)
proportion_above_2186_counts_paired <- sum(resulting_counts_paired > 2186) / length(resulting_counts_paired)

cat("Proportion above 2186 in resulting_counts:", proportion_above_2186_counts, "\n")
cat("Proportion above 2186 in resulting_counts_paired:", proportion_above_2186_counts_paired, "\n")


