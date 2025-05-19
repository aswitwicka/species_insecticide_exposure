# Check the 1-1 orthologs and orthogroup overlaps 

# Get the DEGs 
model_results <- list.files(pattern="*.csv", 
                            path = "~/2023-species_comparisons/2024_results/deseq2_degs/")
files_full <- paste(
  "~/2023-species_comparisons/2024_results/deseq2_degs/",
  model_results, sep = ""
)
for (file in 1:length(model_results)) {
  assign(model_results[file], 
         read.csv(files_full[file])
  )
}

# DEGs detected by DESeq2 in independent testing 
dge_res_list <- list(
  bter_clot = `bter_clot_singleModel_052024.csv`,
  obic_clot = `obic_clot_singleModel_052024.csv`,
  bter_sulf = `bter_sulf_17092023.csv`,
  obic_sulf = `obic_sulf_052024.csv`,
  lser_clot = `lser_clot_singleModel_052024.csv`,
  lser_sulf = `lser_sulf_singleModel_052024.csv`,
  vcar_clot = `vcar_clot_singleModel_052024.csv`,
  vcar_sulf = `vcar_sulf_singleModel_052024.csv`
)
deg_only <- function(deg_table) {
  mutate_table <- deg_table %>%
    mutate(expression = case_when(
      padj <= 0.05 & log2FoldChange > 0 ~ "Up-regulated",
      padj <= 0.05 & log2FoldChange < 0 ~ "Down-regulated",
      TRUE ~ "Unchanged"
    )) %>% filter( expression != "Unchanged")
  mutate_table$X <- gsub("gene-", "", mutate_table$X)
  return(mutate_table)
}
dge_res_list <- purrr::map(.x = dge_res_list, .f = deg_only)

# Get ortholog lists
orthofinder_results <- list.files(pattern="*.csv", 
                                  path = "~/2023-species_comparisons/tables/one_one_orthologs/")
files_full <- paste(
  "~/2023-species_comparisons/tables/one_one_orthologs/",
  orthofinder_results, sep = ""
)
for (file in 1:length(orthofinder_results)) {
  assign(orthofinder_results[file], 
         read.csv(files_full[file])
  )
}

ortholog_list <- list(
  all_species_orthologs = `prot2gene_all_species_orthologs.csv`,
  vcar_bter_orthologs = `vcar_bter_1to1_ort_prot2gene.csv`,
  obic_bter_orthologs = `obic_bter_1to1_ort_prot2gene.csv`,
  lser_bter_orthologs = `lser_bter_1to1_ort_prot2gene.csv`,
  vcar_obic_orthologs = `vcar_obic_1to1_ort_prot2gene.csv`,
  lser_obic_orthologs = `lser_obic_1to1_ort_prot2gene.csv`,
  vcar_lser_orthologs = `vcar_lser_1to1_ort_prot2gene.csv`
)

# Check how many 1-1 orthologs present in all species are differentially expressed 

de_ort_bter_clot <- dge_res_list$bter_clot[dge_res_list$bter_clot$X %in% 
                                             ortholog_list$all_species_orthologs$Bter_gene,] # 623 / 2186 = 28.5%
de_ort_vcar_clot <- dge_res_list$vcar_clot[dge_res_list$vcar_clot$X %in% 
                                             ortholog_list$all_species_orthologs$Vcar_gene,] # 13 / 52 = 25%
de_ort_lser_clot <- dge_res_list$lser_clot[dge_res_list$lser_clot$X %in% 
                                             ortholog_list$all_species_orthologs$Lser_gene,] # 49 / 498 = 9%
de_ort_obic_clot <- dge_res_list$obic_clot[dge_res_list$obic_clot$X %in% 
                                             ortholog_list$all_species_orthologs$Obic_gene,] # 220 / 759 = 29% 

# Sulfoxaflor
de_ort_vcar_sulf <- dge_res_list$vcar_sulf[dge_res_list$vcar_sulf$X %in% 
                                             ortholog_list$all_species_orthologs$Vcar_gene,] # 181 / 734 = 24.7%
de_ort_lser_sulf <- dge_res_list$lser_sulf[dge_res_list$lser_sulf$X %in% 
                                             ortholog_list$all_species_orthologs$Lser_gene,] # 426 / 1868 = 23%
de_ort_obic_sulf <- dge_res_list$obic_sulf[dge_res_list$obic_sulf$X %in% 
                                             ortholog_list$all_species_orthologs$Obic_gene,] # 0

# Compare orthologs function: 
deg_ort_overlap <- function(species_gene1, species_gene2, deg_vector_species1, deg_vector_species2, ortholog_df){
  comparison_genes1 <- ortholog_df %>%
    filter({{ species_gene1 }} %in% deg_vector_species1) %>%
    pull({{ species_gene2 }})
  comparison_genes2 <- ortholog_df %>%
    filter({{ species_gene2 }} %in% deg_vector_species2) %>%
    pull({{ species_gene2 }})
  output_table1 <- table(comparison_genes1 %in% comparison_genes2)
  output_table2 <- table(comparison_genes2 %in% comparison_genes1)
  output_df <- comparison_genes1[comparison_genes1 %in% comparison_genes2]
  return(list(output_df, output_table1, output_table2))
}

# Enrichment function:
keggSpeciesEnrich <- function(genes_of_interest, kegg_organism){
  enrichment_result <- enrichKEGG(
    gene         = genes_of_interest,
    organism     = kegg_organism,
    pvalueCutoff = 0,
    pAdjustMethod= "BH",
    qvalueCutoff = 0,
    maxGSSize    = 200000
  )

  # Check if the result is NULL
  if (is.null(enrichment_result)) {
    warning("No enriched pathways found. Returning NULL.")
    return(NULL)
  }

  return(enrichment_result@result)
}

ort_bter_obic <- deg_ort_overlap(species_gene1 = Obic_gene, # 136
                                 species_gene2 = Bter_gene, 
                                 deg_vector_species1 = dge_res_list$obic_clot$X, 
                                 deg_vector_species2 = dge_res_list$bter_clot$X,
                                 ortholog_df = ortholog_list$obic_bter_orthologs)
keggSpeciesEnrich(gsub("LOC", "", ort_bter_obic[[1]]), "bter") %>% filter(pvalue < 0.05)

ort_bter_vcar <- deg_ort_overlap(species_gene1 = Vcar_gene, # 75
                                 species_gene2 = Bter_gene, 
                                 deg_vector_species1 = dge_res_list$vcar_clot$X, 
                                 deg_vector_species2 = dge_res_list$bter_clot$X,
                                 ortholog_df = ortholog_list$vcar_bter_orthologs)
keggSpeciesEnrich(gsub("LOC", "", ort_bter_vcar[[1]]), "bter") %>% filter(pvalue < 0.05)

ort_vcar_obic <- deg_ort_overlap(species_gene1 = Vcar_gene, # 23
                                 species_gene2 = Obic_gene, 
                                 deg_vector_species1 = dge_res_list$vcar_clot$X, 
                                 deg_vector_species2 = dge_res_list$obic_clot$X,
                                 ortholog_df = ortholog_list$vcar_obic_orthologs)

ort_lser_obic <- deg_ort_overlap(species_gene1 = Obic_gene, # 42
                                 species_gene2 = Lser_gene, 
                                 deg_vector_species1 = dge_res_list$obic_clot$X, 
                                 deg_vector_species2 = dge_res_list$lser_clot$X,
                                 ortholog_df = ortholog_list$lser_obic_orthologs)
keggSpeciesEnrich(gsub("LOC", "", ort_lser_obic[[1]]), "lsq") %>% filter(pvalue < 0.05)

ort_lser_bter <- deg_ort_overlap(species_gene1 = Bter_gene, # 151
                                 species_gene2 = Lser_gene, 
                                 deg_vector_species1 = dge_res_list$bter_clot$X, 
                                 deg_vector_species2 = dge_res_list$lser_clot$X,
                                 ortholog_df = ortholog_list$lser_bter_orthologs)
keggSpeciesEnrich(gsub("LOC", "", ort_lser_bter[[1]]), "lsq") %>% filter(pvalue < 0.05)

ort_lser_vcar <- deg_ort_overlap(species_gene1 = Vcar_gene, # 54
                                 species_gene2 = Lser_gene, 
                                 deg_vector_species1 = dge_res_list$vcar_sulf$X, 
                                 deg_vector_species2 = dge_res_list$lser_sulf$X,
                                 ortholog_df = ortholog_list$vcar_lser_orthologs)
keggSpeciesEnrich(gsub("LOC", "", ort_lser_vcar[[1]]), "lsq") %>% filter(pvalue < 0.05)



# Parameters for the hypergeometric test
# N <- length(unique(ortholog_list$lser_bter_orthologs$Lser_gene))  # Population size: total number of 1-1 orthologs
# K <- as.numeric(unname(ort_lser_bter[[2]]["TRUE"])) + as.numeric(unname(ort_lser_bter[[2]]["FALSE"])) # Number of successes in population: DE orthologs in Species A
# n <- as.numeric(unname(ort_lser_bter[[3]]["TRUE"])) + as.numeric(unname(ort_lser_bter[[3]]["FALSE"])) # Sample size: DE orthologs in Species B
# x <- as.numeric(unname(ort_lser_bter[[2]]["TRUE"])) # Observed successes: overlapping DE orthologs

hyperTest(N = length(unique(ortholog_list$vcar_lser_orthologs$Lser_gene)),
          K = as.numeric(unname(ort_lser_vcar[[2]]["FALSE"])) + as.numeric(unname(ort_lser_vcar[[2]]["TRUE"])),
          n = as.numeric(unname(ort_lser_vcar[[3]]["FALSE"])) + as.numeric(unname(ort_lser_vcar[[3]]["TRUE"])),
          x = as.numeric(unname(ort_lser_vcar[[2]]["TRUE"])))

hyperTest(N = length(unique(ortholog_list$obic_bter_orthologs$Bter_gene)), #
          K = as.numeric(unname(ort_bter_obic[[2]]["FALSE"])) + as.numeric(unname(ort_bter_obic[[2]]["TRUE"])),
          n = as.numeric(unname(ort_bter_obic[[3]]["FALSE"])) + as.numeric(unname(ort_bter_obic[[3]]["TRUE"])),
          x = as.numeric(unname(ort_bter_obic[[2]]["TRUE"])))

hyperTest(N = length(unique(ortholog_list$lser_bter_orthologs$Lser_gene)),
          K = as.numeric(unname(ort_lser_bter[[2]]["FALSE"])) + as.numeric(unname(ort_lser_bter[[2]]["TRUE"])),
          n = as.numeric(unname(ort_lser_bter[[3]]["FALSE"])) + as.numeric(unname(ort_lser_bter[[3]]["TRUE"])),
          x = as.numeric(unname(ort_lser_bter[[2]]["TRUE"])))

hyperTest(N = length(unique(ortholog_list$vcar_bter_orthologs$Bter_gene)),
          K = as.numeric(unname(ort_bter_vcar[[2]]["FALSE"])) + as.numeric(unname(ort_bter_vcar[[2]]["TRUE"])),
          n = as.numeric(unname(ort_bter_vcar[[3]]["FALSE"])) + as.numeric(unname(ort_bter_vcar[[3]]["TRUE"])),
          x = as.numeric(unname(ort_bter_vcar[[2]]["TRUE"])))

hyperTest(N = length(unique(ortholog_list$lser_obic_orthologs$Lser_gene)),
          K = as.numeric(unname(ort_lser_obic[[2]]["FALSE"])) + as.numeric(unname(ort_lser_obic[[2]]["TRUE"])),
          n = as.numeric(unname(ort_lser_obic[[3]]["FALSE"])) + as.numeric(unname(ort_lser_obic[[3]]["TRUE"])),
          x = as.numeric(unname(ort_lser_obic[[2]]["TRUE"])))

hyperTest(N = length(unique(ortholog_list$vcar_obic_orthologs$Obic_gene)),
          K = as.numeric(unname(ort_vcar_obic[[2]]["FALSE"])) + as.numeric(unname(ort_vcar_obic[[2]]["TRUE"])),
          n = as.numeric(unname(ort_vcar_obic[[3]]["FALSE"])) + as.numeric(unname(ort_vcar_obic[[3]]["TRUE"])),
          x = as.numeric(unname(ort_vcar_obic[[2]]["TRUE"])))


# Specify the function: 
hyperTest <- function(N, K, n, x){
# Calculating the p-value
# Note: 'lower.tail = FALSE' calculates the upper tail probability,
# which is the p-value for at least 'x' successes.
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
set.seed(42)
data <- rhyper(x - 1, # gene overlap - 1
       K, # list 2
       N - K, # universe - list 2
       n
) %>%
  data.frame() %>%
  pivot_longer(everything(), names_to = "Hypergeometric", values_to = "Overlap")
plot <- data %>%
  ggplot(aes(x=Overlap, fill = Hypergeometric)) +
  geom_histogram(binwidth = 2, colour = "white") +
  geom_vline(xintercept = x, lty = 2, color="grey", linewidth = 1) +
  theme(text = element_text(size = 10), legend.position = "none") +
  ylab("Count")
return(list(p_value, plot, data))
}

# ORTHOGROUPS 
# Tables come from orthogroups_search.R script in the project code directory 

gr_obic_vcar <- read.csv("~/2023-species_comparisons/tables/orthogroup_tables/obic_vcar_orthogroup_gene.csv") 
gr_obic_lser <- read.csv("~/2023-species_comparisons/tables/orthogroup_tables/obic_lser_orthogroup_gene.csv") 
gr_obic_bter <- read.csv("~/2023-species_comparisons/tables/orthogroup_tables/bter_obic_orthogroup_gene.csv") 
gr_lser_vcar <- read.csv("~/2023-species_comparisons/tables/orthogroup_tables/vcar_lser_orthogroup_gene.csv") 
gr_bter_vcar <- read.csv("~/2023-species_comparisons/tables/orthogroup_tables/bter_vcar_orthogroup_gene.csv") 
gr_bter_lser <- read.csv("~/2023-species_comparisons/tables/orthogroup_tables/bter_lser_orthogroup_gene.csv") 

# Plot for bter and obic 

bter_logs <- updown_list$bter_clot %>% 
  filter(X %in% ortholog_list$obic_bter_orthologs$Bter_gene)
bter_logs$Bter_gene <- bter_logs$X
bter_logs <- merge(bter_logs, ortholog_list$obic_bter_orthologs, by = "Bter_gene")

obic_logs <- updown_list$obic_clot %>% 
  filter(X %in% ortholog_list$obic_bter_orthologs$Obic_gene)
obic_logs$Obic_gene <- obic_logs$X
obic_logs <- merge(obic_logs, ortholog_list$obic_bter_orthologs, by = "Obic_gene")

bees_logs <- merge(bter_logs, obic_logs, by = "Obic_gene")

min_value <- min(c(bees_logs$log2FoldChange.x, bees_logs$log2FoldChange.y, nonbees_logs$log2FoldChange.x, nonbees_logs$log2FoldChange.y))
max_value <- max(c(bees_logs$log2FoldChange.x, bees_logs$log2FoldChange.y, nonbees_logs$log2FoldChange.x, nonbees_logs$log2FoldChange.y))

bter_obic_orth <- ggplot() +
  # Bter significant 
  geom_point(data = (bees_logs %>% filter(padj.x < 0.05)), 
             aes(x = log2FoldChange.x, y = log2FoldChange.y), alpha = 0.7, colour = "#EDCB5C", size = 1
             ) +
  # Obic significant 
  geom_point(data = (bees_logs %>% filter(padj.y < 0.05)), 
               aes(x = log2FoldChange.x, y = log2FoldChange.y), alpha = 0.7, colour = "#5CAFDC", size = 1
             ) +
  geom_point(data = (bees_logs %>% filter(padj.y < 0.05) %>% filter(padj.x < 0.05) %>% 
                       filter(Bter_gene.x %in% ort_bter_obic[[1]])),
             aes(x = log2FoldChange.x, y = log2FoldChange.y), colour = "#6E6E6E", alpha = 1, size = 1.2) +
  geom_text(data = (bees_logs %>% filter(padj.y < 0.05) %>% filter(padj.x < 0.05) %>% 
                      filter(Bter_gene.x %in% ort_bter_obic[[1]])),
            aes(x = log2FoldChange.x, y = log2FoldChange.y, label = Obic_gene), 
            colour = "red", hjust = 0, vjust = 0, size = 1) +
    ylab("log2fold solitary bee") + xlab("log2fold bumblebee") +
  xlim(c(min_value, max_value)) + 
  ylim(c(min_value, max_value)) 

# Genes with the same direction 
same_direction_bees <- bees_logs %>% filter(padj.x < 0.05) %>% filter(padj.y < 0.05) %>%
  filter((log2FoldChange.x > 0 & log2FoldChange.y > 0) | (log2FoldChange.x < 0 & log2FoldChange.y < 0))

bter_anno$Bter_gene.x <- bter_anno$Symbol

same_direction_bees <- merge(same_direction_bees, bter_anno, by = "Bter_gene.x")


# Plot for vcar and lser 

vcar_logs <- updown_list$vcar_sulf %>% 
    filter(X %in% ortholog_list$vcar_lser_orthologs$Vcar_gene)
  vcar_logs$Vcar_gene <- vcar_logs$X
  vcar_logs <- merge(vcar_logs, ortholog_list$vcar_lser_orthologs, by = "Vcar_gene")
  
  lser_logs <- updown_list$lser_sulf %>% 
    filter(X %in% ortholog_list$vcar_lser_orthologs$Lser_gene)
  lser_logs$Lser_gene <- lser_logs$X
  lser_logs <- merge(lser_logs, ortholog_list$vcar_lser_orthologs, by = "Lser_gene")
  
nonbees_logs <- merge(vcar_logs, lser_logs, by = "Lser_gene")

  lser_vcar_orth <- ggplot() +
    # Bter significant 
    geom_point(data = (nonbees_logs %>% filter(padj.x < 0.05)), 
               aes(x = log2FoldChange.y, y = log2FoldChange.x), alpha = 0.7, colour = "#EDCB5C", size = 1
    ) +
    # Obic significant 
    geom_point(data = (nonbees_logs %>% filter(padj.y < 0.05)), 
               aes(x = log2FoldChange.y, y = log2FoldChange.x), alpha = 0.7, colour = "#5CAFDC", size = 1
    ) +
    geom_point(data = (nonbees_logs %>% filter(padj.y < 0.05) %>% filter(padj.x < 0.05) %>% 
                         filter(Lser_gene %in% ort_lser_vcar[[1]])),
               aes(x = log2FoldChange.y, y = log2FoldChange.x), colour = "#6E6E6E", alpha = 1, size = 1.2) +
    geom_text(data = (nonbees_logs %>% filter(padj.y < 0.05) %>% filter(padj.x < 0.05) %>% 
                        filter(Lser_gene %in% ort_lser_vcar[[1]])),
              aes(x = log2FoldChange.y, y = log2FoldChange.x, label = Vcar_gene.y), 
              colour = "red", hjust = 0, vjust = 0, size = 1) +
    ylab("log2fold butterfly") + xlab("log2fold fly") +
    xlim(c(min_value, max_value)) + 
    ylim(c(min_value, max_value)) 
  

# Genes with the same direction 
same_direction_nonbees <- nonbees_logs %>% filter(padj.x < 0.05) %>% filter(padj.y < 0.05) %>%
  filter((log2FoldChange.x > 0 & log2FoldChange.y > 0) | (log2FoldChange.x < 0 & log2FoldChange.y < 0))

vcar_anno$Vcar_gene.x <- vcar_anno$Symbol
same_direction_nonbees <- merge(same_direction_nonbees, vcar_anno, by = "Vcar_gene.x")
same_direction_nonbees %>% dplyr::select(Vcar_gene.x, description)


bter_obic_orth + lser_vcar_orth

# Tafazzin (LOC124537509)
taf_vcar <- single_gene_plot(single_gene = "gene-LOC124537509", vcar_dds, "vcar")
taf_bter <- single_gene_plot(single_gene = "LOC100651201", bter_dds, "bter")
taf_obic <- single_gene_plot(single_gene = "LOC114872286", obic_dds, "obic")
taf_lser <- single_gene_plot(single_gene = "LOC119615633", lser_dds, "lser")
  
taf_vcar + taf_lser + taf_obic + taf_bter
