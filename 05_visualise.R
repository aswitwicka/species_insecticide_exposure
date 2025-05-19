# Upload data:
bter_wald_clot <- read.csv("~/2023-species_comparisons/2024_results/deseq2_degs/bter_clot_singleModel_052024.csv")
vcar_wald_clot <- read.csv("~/2023-species_comparisons/2024_results/deseq2_degs/vcar_clot_singleModel_052024.csv")
lser_wald_clot <- read.csv("~/2023-species_comparisons/2024_results/deseq2_degs/lser_clot_singleModel_052024.csv")
obic_wald_clot <- read.csv("~/2023-species_comparisons/2024_results/deseq2_degs/obic_clot_singleModel_052024.csv")
bter_wald_sulf <- read.csv("~/2023-species_comparisons/2024_results/deseq2_degs/bter_sulf_17092023.csv")
vcar_wald_sulf <- read.csv("~/2023-species_comparisons/2024_results/deseq2_degs/vcar_sulf_singleModel_052024.csv")
lser_wald_sulf <- read.csv("~/2023-species_comparisons/2024_results/deseq2_degs/lser_sulf_singleModel_052024.csv")
obic_wald_sulf <- read.csv("~/2023-species_comparisons/2024_results/deseq2_degs/obic_sulf_052024.csv") 

bter_meta <- read.csv("~/2023-species_comparisons/tables/metadata/bter_meta.csv")
vcar_meta <- read.csv("~/2023-species_comparisons/tables/metadata/vcar_meta.csv")
obic_meta <- read.csv("~/2023-species_comparisons/tables/metadata/obic_meta.csv")
lser_meta <-  read.csv("~/2023-species_comparisons/tables/metadata/lser_meta_update.csv")

bter_dds <- read_rds("~/2023-species_comparisons/DGE_pipelines/bter_dds_filtered_052024.rds")
vcar_dds <- read_rds("~/2023-species_comparisons/DGE_pipelines/vcar_dds_filtered_052024.rds")
lser_dds <- read_rds("~/2023-species_comparisons/DGE_pipelines/lser_dds_filtered_052024.rds")
obic_dds <- read_rds("~/2023-species_comparisons/DGE_pipelines/obic_dds_filtered_052024.rds")

treatment_list <- list(bter_clot = bter_wald_clot, vcar_clot = vcar_wald_clot, 
                       lser_clot = lser_wald_clot, obic_clot = obic_wald_clot, 
                       bter_sulf = bter_wald_sulf, vcar_sulf = vcar_wald_sulf, 
                       lser_sulf = lser_wald_sulf, obic_sulf = obic_wald_sulf)

# UP- DOWN- REGULATED GENES TABLE 
# Get up and down -regulated genes
updown_table <- function(deg_table) {
  mutate_table <- deg_table %>%
    mutate(expression = case_when(
      log2FoldChange > 0.0 & padj <= 0.05 ~ "Up-regulated",
      log2FoldChange < 0.0 & padj <= 0.05 ~ "Down-regulated",
      TRUE ~ "Unchanged"
    ))
  return(mutate_table)
}
updown_list <- purrr::map(.x = treatment_list, .f = updown_table)

# Write up- down- regulated table 
tableFunction <- function(data_frame) {
  table_vector <- data_frame %>%
    pull(expression)
  table(factor(table_vector, levels = c("Down-regulated", "Unchanged", "Up-regulated")))
}
updown_table <- purrr::map(.x = updown_list, .f = tableFunction) %>%
  unlist() %>%
  matrix(ncol = 3, byrow = TRUE) %>%
  as.data.frame()
colnames(updown_table) <- c("Down-regulated", "Unchanged", "Up-regulated")
updown_table$Treatment <- c(
  "bter_clot", "vcar_clot", 
  "lser_clot", "obic_clot", 
  "bter_sulf", "vcar_sulf", 
  "lser_sulf", "obic_sulf"
)
updown_table %>% gt::gt(rowname_col = "Treatment")

# VOLCANO PLOTS
volcano_plot <- function(data_frame) {
  data_frame <- as.data.frame(data_frame)
  updown <- data.frame(
    logFC = data_frame$log2FoldChange,
    negLogPval = -log10(pmax(data_frame$padj, 1e-16)),
    padj = data_frame$padj,
    expression = data_frame$expression,
    gene = data_frame$X
  )
  plot <- ggplot(updown) +
    geom_point(aes(y = negLogPval, x = logFC, colour = expression), size = 2, pch = 20, alpha = 0.5) +
    theme_classic() +
    xlab(expression("log"[2] * "Fold Change")) +
    ylab(expression("-log"[10] * "Benjamini-Hochberg adjusted p-value")) +
    scale_colour_manual(values = c("#A97DE5", "#fecf79", "#A97DE5")) +
    theme(legend.position = "none") + xlim(-10, 10) + ylim(0, 16)
  
  return(list(updown, plot))
}
volcanos_species <- purrr::map(.x = updown_list, .f = volcano_plot)

# OVERLAPS
# Vanessa cardui 
# Clot: 10 | 42 | 692 :Sulf
((updown_list$vcar_sulf %>% filter(padj < 0.05) %>% pull(X)) %in% 
  (updown_list$vcar_clot %>% filter(padj < 0.05) %>% pull(X))) %>% table

# Lucilia sericata
#Clot: 73 | 425 | 1443 :Sulf
((updown_list$lser_sulf %>% filter(padj < 0.05) %>% pull(X)) %in% 
    (updown_list$lser_clot %>% filter(padj < 0.05) %>% pull(X))) %>% table

# Bombus terrestris 
#Clot: 2185 | 1 | 0 :Sulf
((updown_list$bter_clot %>% filter(padj < 0.05) %>% pull(X)) %in% 
    (updown_list$bter_sulf %>% filter(padj < 0.05) %>% pull(X))) %>% table

# Osmia bicornis
#Clot: 758 | 1 | 3 :Sulf
((updown_list$obic_clot %>% filter(padj < 0.05) %>% pull(X)) %in% 
    (updown_list$obic_sulf %>% filter(padj < 0.05) %>% pull(X))) %>% table

# REGRESSION MODELS 
# Vanessa cardui
folds_vcar_sulf <- updown_list$vcar_sulf %>% 
  dplyr::select(X, log2FoldChange)
folds_vcar_sulf$treatment <- "Sulfoxaflor"
folds_vcar_sulf$sulf_Deg <- folds_vcar_sulf$X %in% 
  (updown_list$vcar_sulf %>% filter(padj < 0.05) %>% pull(X))

folds_vcar_clot <- updown_list$vcar_clot %>% filter(X %in% folds_vcar_sulf$X) %>%
  dplyr::select(X, log2FoldChange)
folds_vcar_clot$treatment <- "Clothianidin"
vcar_clot_degs <- updown_list$vcar_clot %>% filter(padj < 0.05) %>% pull(X)

folds_vcar <- merge(folds_vcar_sulf, folds_vcar_clot, by = "X")
folds_vcar$clot_Deg <- folds_vcar$X %in% vcar_clot_degs

colnames(folds_vcar) <- c("X", "log2FoldChange.S", "treatment.S", "sulf_Deg", 
                          "log2FoldChange.C", "treatment.C", "clot_Deg")

wgcna_vcar <- read.csv("~/2023-species_comparisons/2024_results/vcar_red_moduleWGCNA.csv")

# Plot
## Both treatments DE
vcar_1 <- ggplot(data = (folds_vcar %>% filter(sulf_Deg == TRUE & clot_Deg == TRUE)), 
       aes(x = log2FoldChange.S, y = log2FoldChange.C)) + 
  geom_point(colour = "#078a61", alpha = 1, size = 1.5) + 
  geom_smooth(method = "lm", colour = "#078a61", fill = "#078a61") +
  ylab("log2fold Clothianidin") + xlab("log2fold Sulfxoaflor") +
  xlim(c(-2, 3)) + 
  ylim(c(-2, 3)) + 
  geom_point(data = (folds_vcar %>% filter(sulf_Deg == TRUE & clot_Deg == TRUE) %>% 
                       filter(X %in% vcar_green)), 
             aes(x = log2FoldChange.S, y = log2FoldChange.C),
             colour = "green", alpha = 1, size = 1.5) +
  geom_point(data = (folds_vcar %>% filter(sulf_Deg == TRUE & clot_Deg == TRUE) %>% 
                       filter(X %in% wgcna_vcar$Symbol)), 
             aes(x = log2FoldChange.S, y = log2FoldChange.C),
             colour = "red", alpha = 1, size = 1.5) 
# LM
model <- lm(data = (folds_vcar %>% filter(sulf_Deg == TRUE & clot_Deg == TRUE)), log2FoldChange.S~log2FoldChange.C)
slope <- coef(model)[2]
summary(model)
# Estimate Std. Error t value Pr(>|t|)
# 0.97720    0.02698  36.214   <2e-16 ***
# Pearson 
filtered_data <- folds_vcar %>% filter(sulf_Deg == TRUE & clot_Deg == TRUE)
pearson_test <- cor.test(filtered_data$log2FoldChange.S, filtered_data$log2FoldChange.C, method = "pearson")
pearson_test$estimate
pearson_test$p.value
# t = 36.214, df = 40, p-value < 2.2e-16 ([1] 3.384945e-32)
# 0.98509 

# The slope is close to 1 (0.97720), suggesting a nearly perfect proportional relationship between the treatments.
# Very high correlation (0.98509), indicating a strong linear relationship.
# Lower standard error (0.02698), suggesting more precise estimates.
# Model shows a nearly perfect linear relationship between the treatments with a slope close to 1.
  
## Only Sulfoxaflor
vcar_2 <- ggplot(data = (folds_vcar %>% filter(sulf_Deg == TRUE & clot_Deg == FALSE)), 
       aes(x = log2FoldChange.S, y = log2FoldChange.C)) + 
  geom_point(colour = "#FEE05A", alpha = 1, size = 1.5) + 
  geom_smooth(method = "lm", colour = "#FEE05A", fill = "#FEE05A") +
  ylab("log2fold Clothianidin") + xlab("log2fold Sulfxoaflor") +
  xlim(c(-2, 3)) + 
  ylim(c(-2, 3)) +
  geom_point(data = (folds_vcar %>% filter(sulf_Deg == TRUE & clot_Deg == FALSE) %>% 
                       filter(X %in% vcar_green)), 
             aes(x = log2FoldChange.S, y = log2FoldChange.C),
             colour = "green", alpha = 1, size = 1.5) +
  geom_point(data = (folds_vcar %>% filter(sulf_Deg == TRUE & clot_Deg == FALSE) %>% 
                       filter(X %in% wgcna_vcar$Symbol)), 
             aes(x = log2FoldChange.S, y = log2FoldChange.C),
             colour = "red", alpha = 1, size = 1.5) 

model <- lm(data = (folds_vcar %>% filter(sulf_Deg == TRUE & clot_Deg == FALSE)), log2FoldChange.S~log2FoldChange.C)
coef(model)[2]
summary(model)
# Estimate Std. Error t value Pr(>|t|)
# 1.505163   0.031249  48.167   <2e-16 ***
# Pearson 
filtered_data <- folds_vcar %>% filter(sulf_Deg == TRUE & clot_Deg == FALSE)
pearson_test <- cor.test(filtered_data$log2FoldChange.S, filtered_data$log2FoldChange.C, method = "pearson")
pearson_test$estimate
pearson_test$p.value
# t = 48.167, df = 690, p-value < 2.2e-16 ([1] 6.846784e-223)
# 0.8779324  

# The slope is significantly higher than 1 (1.505163), 
# indicating that for every unit increase in log2 fold change in treatment C, 
# the log2 fold change in treatment S increases by more than 1.5 units, suggesting a steeper relationship.
# High correlation but lower than Model 1 (0.8779324), indicating a less strong but still significant linear relationship.
# Slightly higher standard error (0.031249), indicating less precise estimates compared to Model 1.
# Model shows a significant linear relationship but with a steeper slope, 
# suggesting a different dynamic between the two treatments where the change in one treatment 
# leads to a more than proportional change in the other.

set.seed(2)
random_vcar <- folds_vcar %>% 
  filter(sulf_Deg == FALSE & clot_Deg == FALSE) %>% 
  slice_sample(n = 42)
vcar_3 <- ggplot(data = (random_vcar), 
       aes(x = log2FoldChange.S, y = log2FoldChange.C)) + 
  geom_point(colour = "#8B8989", alpha = 1, size = 1.5) + 
  geom_smooth(method = "lm", colour = "#8B8989", fill = "#8B8989") +
  ylab("log2fold Clothianidin") + xlab("log2fold Sulfxoaflor") +
  xlim(c(-2, 3)) + 
  ylim(c(-2, 3))
vcar_3

model <- lm(data = random_vcar, log2FoldChange.S~log2FoldChange.C)
coef(model)[2]
summary(model)
# Estimate Std. Error t value Pr(>|t|)
# 0.15310    0.18058   0.848    0.402
# Pearson 
pearson_test <- cor.test(random_vcar$log2FoldChange.S, random_vcar$log2FoldChange.C, method = "pearson")
pearson_test$estimate
pearson_test$p.value
# t = 0.84781, df = 40, p-value = 0.4016
# 0.1328621  

# Model demonstrates a very weak and statistically insignificant relationship, 
# as indicated by a low correlation coefficient (0.1328621) and high p-values for both the slope (0.402) 
# and the correlation test (0.4016). This suggests there is no meaningful linear relationship between 
# the two variables in this data set.

# Lucilia sericata
folds_lser_sulf <- updown_list$lser_sulf %>% 
  dplyr::select(X, log2FoldChange)
folds_lser_sulf$treatment <- "Sulfoxaflor"
folds_lser_sulf$sulf_Deg <- folds_lser_sulf$X %in% 
  (updown_list$lser_sulf %>% filter(padj < 0.05) %>% pull(X))

folds_lser_clot <- updown_list$lser_clot %>% filter(X %in% folds_lser_sulf$X) %>%
  dplyr::select(X, log2FoldChange)
folds_lser_clot$treatment <- "Clothianidin"
lser_clot_degs <- updown_list$lser_clot %>% filter(padj < 0.05) %>% pull(X)

folds_lser <- merge(folds_lser_sulf, folds_lser_clot, by = "X")
folds_lser$clot_Deg <- folds_lser$X %in% lser_clot_degs

colnames(folds_lser) <- c("X", "log2FoldChange.S", "treatment.S", "sulf_Deg", 
                          "log2FoldChange.C", "treatment.C", "clot_Deg")

# Plot
## Both treatments DE
lser_1 <- ggplot(data = (folds_lser %>% filter(sulf_Deg == TRUE & clot_Deg == TRUE)), 
       aes(x = log2FoldChange.S, y = log2FoldChange.C)) + 
  geom_point(colour = "#078a61", alpha = 1, size = 1.5) + 
  geom_smooth(method = "lm", colour = "#078a61", fill = "#078a61") +
  ylab("log2fold Clothianidin") + xlab("log2fold Sulfxoaflor") +
  scale_x_continuous(limits = c(-4, 8.5), breaks = seq(-4, 8, 2)) +
  scale_y_continuous(limits = c(-4, 8.5), breaks = seq(-4, 8, 2)) +
  geom_point(data = (folds_lser %>% filter(sulf_Deg == TRUE & clot_Deg == TRUE) %>% 
                       filter(X %in% lser_turquoise)), 
             aes(x = log2FoldChange.S, y = log2FoldChange.C),
             colour = "darkslategray1", alpha = 1, size = 1.5) +
  geom_point(data = (folds_lser %>% filter(sulf_Deg == TRUE & clot_Deg == TRUE) %>% 
                       filter(X %in% lser_green)), 
             aes(x = log2FoldChange.S, y = log2FoldChange.C),
             colour = "#C1FFC1", alpha = 1, size = 1.5) 

model <- lm(data = (folds_lser %>% filter(sulf_Deg == TRUE & clot_Deg == TRUE)), log2FoldChange.S~log2FoldChange.C)
coef(model)[2]
summary(model)
# Estimate Std. Error t value Pr(>|t|)
# 0.99043    0.01537  64.451  < 2e-16 ***
# Pearson 
filtered_data <- folds_lser %>% filter(sulf_Deg == TRUE & clot_Deg == TRUE)
pearson_test <- cor.test(filtered_data$log2FoldChange.S, filtered_data$log2FoldChange.C, method = "pearson")
pearson_test$estimate
pearson_test$p.value
summary(pearson_test)
# t = 64.451, df = 423, p-value < 2.2e-16 ([1] 7.409663e-221)
# 0.9526695  

# Shows a nearly perfect linear relationship (slope ~ 1).
# Very strong correlation (0.9526695), highly significant (p < 2.2e-16).
  
## Only Sulfoxaflor
lser_2 <- ggplot(data = (folds_lser %>% filter(sulf_Deg == TRUE & clot_Deg == FALSE)), 
       aes(x = log2FoldChange.S, y = log2FoldChange.C)) + 
  geom_point(colour = "#FEE05A", alpha = 1, size = 1.5) + 
  geom_smooth(method = "lm", colour = "#FEE05A", fill = "#FEE05A") +
  ylab("log2fold Clothianidin") + xlab("log2fold Sulfxoaflor") +
  scale_x_continuous(limits = c(-4, 8.5), breaks = seq(-4, 8, 2)) +
  scale_y_continuous(limits = c(-4, 8.5), breaks = seq(-4, 8, 2)) +
  geom_point(data = (folds_lser %>% filter(sulf_Deg == TRUE & clot_Deg == FALSE) %>% 
                       filter(X %in% lser_turquoise)), 
             aes(x = log2FoldChange.S, y = log2FoldChange.C),
             colour = "darkslategray1", alpha = 1, size = 1.5) +
  geom_point(data = (folds_lser %>% filter(sulf_Deg == TRUE & clot_Deg == FALSE) %>% 
                       filter(X %in% lser_green)), 
             aes(x = log2FoldChange.S, y = log2FoldChange.C),
             colour = "#C1FFC1", alpha = 1, size = 1.5) 

model <- lm(data = (folds_lser %>% filter(sulf_Deg == TRUE & clot_Deg == FALSE)), log2FoldChange.S~log2FoldChange.C)
coef(model)[2]
summary(model)
# Estimate Std. Error t value Pr(>|t|)
# 1.31394    0.01505  87.296  < 2e-16 ***
# Pearson 
filtered_data <- folds_lser %>% filter(sulf_Deg == TRUE & clot_Deg == FALSE)
pearson_test <- cor.test(filtered_data$log2FoldChange.S, filtered_data$log2FoldChange.C, method = "pearson")
pearson_test$estimate
pearson_test$p.value
# t = 87.296, df = 1441, p-value < 2.2e-16 ([1] 0)
# 0.9170473

# The slope is significantly higher than 1 (1.31394), indicating a steeper relationship.
# Strong correlation (0.9170473), highly significant (p < 2.2e-16).

## Only Clothianidin
lser_3 <- ggplot(data = (folds_lser %>% filter(sulf_Deg == FALSE & clot_Deg == TRUE)), 
       aes(x = log2FoldChange.S, y = log2FoldChange.C)) + 
  geom_point(colour = "mediumpurple1", alpha = 1, size = 1.5) + 
  geom_smooth(method = "lm", colour = "mediumpurple1", fill = "mediumpurple1") +
  ylab("log2fold Clothianidin") + xlab("log2fold Sulfxoaflor") +
  scale_x_continuous(limits = c(-4, 8.5), breaks = seq(-4, 8, 2)) +
  scale_y_continuous(limits = c(-4, 8.5), breaks = seq(-4, 8, 2)) +
  geom_point(data = (folds_lser %>% filter(sulf_Deg == FALSE & clot_Deg == TRUE) %>% 
                       filter(X %in% lser_turquoise)), 
             aes(x = log2FoldChange.S, y = log2FoldChange.C),
             colour = "darkslategray1", alpha = 1, size = 1.5) +
  geom_point(data = (folds_lser %>% filter(sulf_Deg == FALSE & clot_Deg == TRUE) %>% 
                       filter(X %in% lser_green)), 
             aes(x = log2FoldChange.S, y = log2FoldChange.C),
             colour = "#C1FFC1", alpha = 1, size = 1.5)

model <- lm(data = (folds_lser %>% filter(sulf_Deg == FALSE & clot_Deg == TRUE)), log2FoldChange.S~log2FoldChange.C)
coef(model)[2]
summary(model)
# Estimate Std. Error t value Pr(>|t|)
# 0.30422    0.04036   7.538  1.2e-10 ***
# Pearson 
filtered_data <- folds_lser %>% filter(sulf_Deg == FALSE & clot_Deg == TRUE)
pearson_test <- cor.test(filtered_data$log2FoldChange.S, filtered_data$log2FoldChange.C, method = "pearson")
pearson_test$estimate
pearson_test$p.value
# t = 7.5385, df = 71, p-value = 1.197e-10 ([1] 1.196952e-10)
# 0.6667589 

# The slope is much lower (0.30422), indicating a weaker relationship.
# Moderate correlation (0.6667589), significant but less strong compared to Models 1 and 2.

random_lser <- folds_lser %>% 
  filter(sulf_Deg == FALSE & clot_Deg == FALSE) %>% 
  slice_sample(n = 425)
lser_4 <- ggplot(data = (random_lser), 
       aes(x = log2FoldChange.S, y = log2FoldChange.C)) + 
  geom_point(colour = "#8B8989", alpha = 1, size = 1.5) + 
  geom_smooth(method = "lm", colour = "#8B8989", fill = "#8B8989") +
  ylab("log2fold Clothianidin") + xlab("log2fold Sulfxoaflor") +
  scale_x_continuous(limits = c(-4, 8.5), breaks = seq(-4, 8, 2)) +
  scale_y_continuous(limits = c(-4, 8.5), breaks = seq(-4, 8, 2))
lser_4

model <- lm(data = random_lser, log2FoldChange.S~log2FoldChange.C)
coef(model)[2]
summary(model)
# Estimate Std. Error t value Pr(>|t|)
# 0.551851   0.039163  14.091   <2e-16 ***
# Pearson 
pearson_test <- cor.test(random_lser$log2FoldChange.S, random_lser$log2FoldChange.C, method = "pearson")
pearson_test$estimate
pearson_test$p.value
pearson_test
# t = 14.091, df = 423, p-value < 2.2e-16 (3.041825e-37)
# R = 0.5652051 

# The slope is 0.5652051, indicating a moderate relationship. (The correlation coefficient)
# Moderate correlation (0.551851), significant but the weakest among the models.

# Models 1 and 2 show strong and significant linear relationships, with Model 1 having a slope closer to 1, indicating near-perfect proportionality.
# Model 2 shows a stronger effect (slope > 1) but still a strong correlation.
# Model 3 shows a weaker relationship with a lower slope and moderate correlation.
# Model 4 (Random Data) shows the weakest relationship, indicating less meaningful correlation in this dataset.

# Combine plots
(lser_1 + lser_2 + lser_4) /
(vcar_1 + vcar_2 + vcar_3 )

# LOG2FOLDS
vcar_logs <- data.frame(
  logfolds = updown_list$vcar_sulf %>% filter(padj < 0.05) %>% pull(log2FoldChange),
  logfolds_abs = updown_list$vcar_sulf %>% filter(padj < 0.05) %>% pull(log2FoldChange) %>% abs(),
  species = "Vanessa cardui") %>%
  arrange(desc(logfolds_abs)) %>%
  mutate(top_five = row_number() <= n() * 0.1)
bter_logs <- data.frame(
  logfolds = updown_list$bter_clot %>% filter(padj < 0.05) %>% pull(log2FoldChange),
  logfolds_abs = updown_list$bter_clot %>% filter(padj < 0.05) %>% pull(log2FoldChange) %>% abs(),
  species = "Bombus terrestris") %>%
  arrange(desc(logfolds_abs)) %>%
  mutate(top_five = row_number() <= n() * 0.1)
obic_logs <- data.frame(
  logfolds = updown_list$obic_clot %>% filter(padj < 0.05) %>% pull(log2FoldChange),
  logfolds_abs = updown_list$obic_clot %>% filter(padj < 0.05) %>% pull(log2FoldChange) %>% abs(),
  species = "Osmia bicornis") %>%
  arrange(desc(logfolds_abs)) %>%
  mutate(top_five = row_number() <= n() * 0.1)
lser_logs <- data.frame(
  logfolds = updown_list$lser_sulf %>% filter(padj < 0.05) %>% pull(log2FoldChange),
  logfolds_abs = updown_list$lser_sulf %>% filter(padj < 0.05) %>% pull(log2FoldChange) %>% abs(),
  species = "Lucilia sericata") %>%
  arrange(desc(logfolds_abs)) %>%
  mutate(top_five = row_number() <= n() * 0.1)

lser_logs %>%
  ggplot(aes(x = logfolds, fill = top_five)) + 
  geom_histogram(binwidth = 0.2, aes(y = ..count.. + 1)) + 
  scale_y_continuous(trans = 'log10') + 
  labs(y = "Count (log scale)", x = "Logfolds")

vcar_logs %>%
  ggplot(aes(x = logfolds, fill = top_five)) + 
  geom_histogram(binwidth = 0.1, aes(y = ..count.. + 1)) + 
  scale_y_continuous(trans = 'log10') + 
  labs(y = "Count (log scale)", x = "Logfolds")

obic_logs %>%
  ggplot(aes(x = logfolds, fill = top_five)) + 
  geom_histogram(binwidth = 0.2, aes(y = ..count.. + 1)) + 
  scale_y_continuous(trans = 'log10') + 
  labs(y = "Count (log scale)", x = "Logfolds")

bter_logs %>%
  ggplot(aes(x = logfolds, fill = top_five)) + 
  geom_histogram(binwidth = 0.1, aes(y = ..count.. + 1)) + 
  scale_y_continuous(trans = 'log10') + 
  labs(y = "Count (log scale)", x = "Logfolds")

