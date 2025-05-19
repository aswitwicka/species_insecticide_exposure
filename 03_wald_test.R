# This script runs Wald test from DEseq2 library, the output are .csv files for each pairwise comparison between pesticide and control for each species. 
# The files contail adjusted p-values and log2fold change values

# Upload data 
bter_dds <- read_rds("~/2023-species_comparisons/DGE_pipelines/bter_dds_filtered_052024.rds")
vcar_dds <- read_rds("~/2023-species_comparisons/DGE_pipelines/vcar_dds_filtered_052024.rds")
lser_dds <- read_rds("~/2023-species_comparisons/DGE_pipelines/lser_dds_filtered_052024.rds")
obic_dds <- read_rds("~/2023-species_comparisons/DGE_pipelines/obic_dds_filtered_052024.rds")

bter_meta <- read.csv("~/2023-species_comparisons/tables/metadata/bter_meta.csv")
vcar_meta <- read.csv("~/2023-species_comparisons/tables/metadata/vcar_meta.csv")
obic_meta <- read.csv("~/2023-species_comparisons/tables/metadata/obic_meta.csv")
lser_meta <-  read.csv("~/2023-species_comparisons/tables/metadata/lser_meta_update.csv")

# Run DESeq2 Wald test 
# Bombus terrestris 
bter_dds_wald <- DESeq2::DESeq(bter_dds, test = "Wald", betaPrior = FALSE)
## Explore the results 
DESeq2::resultsNames(bter_dds_wald)
# Clot
bter_wald_clot <- results(bter_dds_wald,
                          # independentFiltering = FALSE,
                          pAdjustMethod = "BH",
                          contrast = c(
                            "CODE",
                            "Clothianidin",
                            "Control"
                          )
) %>% as.data.frame() # 1935
bter_wald_clot %>% filter(padj < 0.05) %>% nrow()
# Sulf
bter_wald_sulf <- results(bter_dds_wald,
                          # independentFiltering = FALSE,
                          pAdjustMethod = "BH",
                          contrast = c(
                            "CODE",
                            "Sulfoxaflor",
                            "Control"
                          )
)  %>% as.data.frame() # 1
bter_wald_sulf %>% filter(padj < 0.05) %>% nrow()
# Save
# write.csv(bter_wald_sulf,
#           file = "DGE_pipelines/bter_sulf_17092023.csv"
# )
# write.csv(bter_wald_clot,
#           file = "DGE_pipelines/bter_clot_17092023.csv"
# )
# Clot only to obtain more gene candidates (sulfoxaflor always returns just one gene - defensin)
bter_meta_clot <- bter_meta %>% filter(CODE == "Clothianidin" | CODE == "Control")
bter_meta_clot$CODE <- as.factor(bter_meta_clot$CODE)

bter_dds_clot <- bter_dds[, which(colData(bter_dds)$SAMPLE_NAME %in% bter_meta_clot$SAMPLE_NAME)]
bter_dds_clot$CODE <- droplevels(bter_dds_clot$CODE)
colData(bter_dds_clot)

bter_dds_wald <- DESeq2::DESeq(bter_dds_clot, test = "Wald", betaPrior = FALSE)
## Explore the results 
DESeq2::resultsNames(bter_dds_wald)
# Clot
bter_wald_clot <- results(bter_dds_wald,
                          # independentFiltering = FALSE,
                          pAdjustMethod = "BH",
                          contrast = c(
                            "CODE",
                            "Clothianidin",
                            "Control"
                          )
) %>% as.data.frame()
bter_wald_clot %>% filter(padj < 0.05) %>% nrow() # 2186
write.csv(bter_wald_clot,
          file = "~/2023-species_comparisons/2024_results/bter_clot_singleModel_052024.csv"
)

# Vanessa cardui
vcar_dds_wald <- DESeq2::DESeq(vcar_dds, test = "Wald", betaPrior = FALSE)
## Explore the results 
DESeq2::resultsNames(vcar_dds_wald)
# Clot
vcar_wald_clot <- results(vcar_dds_wald,
                          # independentFiltering = FALSE,
                          pAdjustMethod = "BH",
                          contrast = c(
                            "CODE",
                            "Clothianidin",
                            "Control"
                          )
) %>% as.data.frame() # 49 
vcar_wald_clot %>% filter(padj < 0.05) %>% nrow() 
# Sulf
vcar_wald_sulf <- results(vcar_dds_wald,
                          # independentFiltering = FALSE,
                          pAdjustMethod = "BH",
                          contrast = c(
                            "CODE",
                            "Sulfoxaflor",
                            "Control"
                          )
)  %>% as.data.frame() # 664
vcar_wald_sulf %>% filter(padj < 0.05) %>% nrow() 

table((vcar_wald_sulf %>% filter(padj < 0.05) %>% rownames()) %in% 
  (vcar_wald_clot %>% filter(padj < 0.05) %>% rownames()))
table((vcar_wald_clot %>% filter(padj < 0.05) %>% rownames()) %in% 
        (vcar_wald_sulf %>% filter(padj < 0.05) %>% rownames()))

# Save
# write.csv(vcar_wald_sulf,
#           file = "vcar_sulf_17092023.csv"
# )
# write.csv(vcar_wald_clot,
#           file = "vcar_clot_17092023.csv"
# )

# Seperately
# Clot
vcar_meta_clot <- vcar_meta %>% filter(CODE == "Clothianidin" | CODE == "Control")
vcar_meta_clot$CODE <- as.factor(vcar_meta_clot$CODE)
vcar_dds_clot <- vcar_dds[, which(colData(vcar_dds)$SAMPLE_CODE %in% vcar_meta_clot$SAMPLE_CODE)]
vcar_dds_clot$CODE <- droplevels(vcar_dds_clot$CODE)
colData(vcar_dds_clot)
#
vcar_dds_wald <- DESeq2::DESeq(vcar_dds_clot, test = "Wald", betaPrior = FALSE)
DESeq2::resultsNames(vcar_dds_wald)
vcar_wald_clot <- results(vcar_dds_wald,
                          # independentFiltering = FALSE,
                          pAdjustMethod = "BH",
                          contrast = c(
                            "CODE",
                            "Clothianidin",
                            "Control"
                          )
) %>% as.data.frame() # 52
vcar_wald_clot %>% filter(padj < 0.05) %>% nrow()

rownames(vcar_wald_clot) <- gsub("gene-", "", rownames(vcar_wald_clot))

write.csv(vcar_wald_clot, 
          "~/2023-species_comparisons/2024_results/deseq2_degs/vcar_clot_singleModel_052024.csv")

# Sulf
vcar_meta_sulf <- vcar_meta %>% filter(CODE == "Sulfoxaflor" | CODE == "Control")
vcar_meta_sulf$CODE <- as.factor(vcar_meta_sulf$CODE)
vcar_dds_sulf <- vcar_dds[, which(colData(vcar_dds)$SAMPLE_CODE %in% vcar_meta_sulf$SAMPLE_CODE)]
vcar_dds_sulf$CODE <- droplevels(vcar_dds_sulf$CODE)
colData(vcar_dds_sulf)
#
vcar_dds_wald <- DESeq2::DESeq(vcar_dds_sulf, test = "Wald", betaPrior = FALSE)
DESeq2::resultsNames(vcar_dds_wald)
vcar_wald_sulf <- results(vcar_dds_wald,
                          # independentFiltering = FALSE,
                          pAdjustMethod = "BH",
                          contrast = c(
                            "CODE",
                            "Sulfoxaflor",
                            "Control"
                          )
) %>% as.data.frame() # 734
vcar_wald_sulf %>% filter(padj < 0.05) %>% nrow()

table((vcar_wald_sulf %>% filter(padj < 0.05) %>% rownames()) %in% 
        (vcar_wald_clot %>% filter(padj < 0.05) %>% rownames()))
table((vcar_wald_clot %>% filter(padj < 0.05) %>% rownames()) %in% 
        (vcar_wald_sulf %>% filter(padj < 0.05) %>% rownames()))

rownames(vcar_wald_sulf) <- gsub("gene-", "", rownames(vcar_wald_sulf))

write.csv(vcar_wald_sulf, 
          "~/2023-species_comparisons/2024_results/deseq2_degs/vcar_sulf_singleModel_052024.csv")


# Lucilia sericata 
lser_dds_wald_update <- DESeq2::DESeq(lser_dds, test = "Wald", betaPrior = FALSE)
## Explore the results 
DESeq2::resultsNames(lser_dds_wald_update)

lser_wald_sulf_update <- results(lser_dds_wald_update,
                                 # independentFiltering = FALSE,
                                 pAdjustMethod = "BH",
                                 contrast = c(
                                   "CODE",
                                   "Sulfoxaflor",
                                   "Control"
                                 )
) %>% as.data.frame() # 1563
lser_wald_sulf_update %>% filter(padj < 0.05) %>% nrow()

lser_wald_clot_update <- results(lser_dds_wald_update,
                                 # independentFiltering = FALSE,
                                 pAdjustMethod = "BH",
                                 contrast = c(
                                   "CODE",
                                   "Clothianidin",
                                   "Control"
                                 )
) %>% as.data.frame() # 447
lser_wald_clot_update %>% filter(padj < 0.05) %>% nrow()

# Seperately
# Sulfoxaflor 
lser_meta_sulf <- lser_meta %>% filter(CODE == "Sulfoxaflor" | CODE == "Control")
lser_meta_sulf$CODE <- as.factor(lser_meta_sulf$CODE)

lser_dds_sulf <- lser_dds[, which(colData(lser_dds)$SAMPLE_CODE %in% lser_meta_sulf$SAMPLE_CODE)]
lser_dds_sulf$CODE <- droplevels(lser_dds_sulf$CODE)
colData(lser_dds_sulf)

lser_dds_wald <- DESeq2::DESeq(lser_dds_sulf, test = "Wald", betaPrior = FALSE)
## Explore the results 
DESeq2::resultsNames(lser_dds_wald)
lser_wald_sulf <- results(lser_dds_wald,
                          # independentFiltering = FALSE,
                          pAdjustMethod = "BH",
                          contrast = c(
                            "CODE",
                            "Sulfoxaflor",
                            "Control"
                          )
) %>% as.data.frame() # 1868
lser_wald_sulf %>% filter(padj < 0.05) %>% nrow()

write.csv(lser_wald_sulf, 
          "~/2023-species_comparisons/2024_results/deseq2_degs/lser_sulf_singleModel_052024.csv")

# Clothianidin 
lser_meta_clot <- lser_meta %>% filter(CODE == "Clothianidin" | CODE == "Control")
lser_meta_clot$CODE <- as.factor(lser_meta_clot$CODE)

lser_dds_clot <- lser_dds[, which(colData(lser_dds)$SAMPLE_CODE %in% lser_meta_clot$SAMPLE_CODE)]
lser_dds_clot$CODE <- droplevels(lser_dds_clot$CODE)
colData(lser_dds_clot)

lser_dds_wald <- DESeq2::DESeq(lser_dds_clot, test = "Wald", betaPrior = FALSE)
## Explore the results 
DESeq2::resultsNames(lser_dds_wald)
lser_wald_clot <- results(lser_dds_wald,
                          # independentFiltering = FALSE,
                          pAdjustMethod = "BH",
                          contrast = c(
                            "CODE",
                            "Clothianidin",
                            "Control"
                          )
) %>% as.data.frame() # 498
lser_wald_clot %>% filter(padj < 0.05) %>% nrow() 

table((lser_wald_clot %>% filter(padj < 0.05) %>% rownames()) %in%
        (lser_wald_sulf %>% filter(padj < 0.05) %>% rownames()))

write.csv(lser_wald_clot, 
          "~/2023-species_comparisons/2024_results/deseq2_degs/lser_clot_singleModel_052024.csv")

# Osmia bicornis
obic_dds_wald <- DESeq2::DESeq(obic_dds, test = "Wald", betaPrior = FALSE)
## Explore the results 
DESeq2::resultsNames(obic_dds_wald)
# Clot
obic_wald_clot <- results(obic_dds_wald,
                          # independentFiltering = FALSE,
                          pAdjustMethod = "BH",
                          contrast = c(
                            "CODE",
                            "Clothianidin",
                            "Control"
                          )
) %>% as.data.frame() # 235
obic_wald_clot %>% filter(padj < 0.05) %>% nrow()

# Sulfoxaflor
obic_wald_sulf <- results(obic_dds_wald,
                          # independentFiltering = FALSE,
                          pAdjustMethod = "BH",
                          contrast = c(
                            "CODE",
                            "Sulfoxaflor",
                            "Control"
                          )
)  %>% as.data.frame() # 4
obic_wald_sulf %>% filter(padj < 0.05) %>% nrow()

table((obic_wald_sulf %>% filter(padj < 0.05) %>% rownames()) %in%
        (obic_wald_clot %>% filter(padj < 0.05) %>% rownames()))

write.csv(obic_wald_sulf, 
          "~/2023-species_comparisons/2024_results/deseq2_degs/obic_sulf_052024.csv")

# Seperately:
# Clot
obic_meta_clot <- obic_meta %>% filter(CODE == "Control" | CODE == "Clothianidin")
obic_meta_clot$CODE <- as.factor(obic_meta_clot$CODE)
obic_dds_clot <- obic_dds[, which(colData(obic_dds)$SAMPLE_CODE %in% obic_meta_clot$SAMPLE_CODE)]
obic_dds_clot$CODE <- droplevels(obic_dds_clot$CODE)
colData(obic_dds_clot)

obic_dds_wald <- DESeq2::DESeq(obic_dds_clot, test = "Wald", betaPrior = FALSE)
DESeq2::resultsNames(obic_dds_wald)
obic_wald_clot <- results(obic_dds_wald,
                          # independentFiltering = FALSE,
                          pAdjustMethod = "BH",
                          contrast = c(
                            "CODE",
                            "Clothianidin",
                            "Control"
                          )
) %>% as.data.frame() # 758
obic_wald_clot %>% filter(padj < 0.05) %>% nrow()

write.csv(obic_wald_clot, 
          "~/2023-species_comparisons/2024_results/deseq2_degs/obic_clot_singleModel_052024.csv")
