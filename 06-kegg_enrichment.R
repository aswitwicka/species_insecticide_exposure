library(enrichplot)
library(clusterProfiler)
library(pathview)

bter_pathways <- keggList("pathway", "bter")
obb_pathways <- keggList("pathway", "obb")
vcd_pathways <- keggList("pathway", "vcd")
lsq_pathways <- keggList("pathway", "lsq")

# Specify a new function 
keggSpecies <- function(genes_of_interest, kegg_organism){
  enrichment_result <- gseKEGG(
    geneList = genes_of_interest,
    organism = kegg_organism,
    keyType = "kegg",
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 500,
    eps = 1e-10,
    pvalueCutoff = 1,
    pAdjustMethod = "BH"
  )
  
  # Check if the result is NULL
  if (is.null(enrichment_result)) {
    warning("No enriched pathways found. Returning NULL.")
    return(NULL)
  }
  
  return(enrichment_result@result)
}

downGSE <- function(updown_list, species, treatment){
  bter_trial <- updown_list %>% filter(log2FoldChange < 0) %>% arrange(pvalue) %>% pull(pvalue)
  bter_trial <- -log10(bter_trial)
  names(bter_trial) <- gsub("gene-", "", (updown_list %>% filter(log2FoldChange < 0) %>% arrange(pvalue) %>% pull(X)))
  names(bter_trial) <- gsub("LOC", "", names(bter_trial))
  geneList_bter_clot <- sort(bter_trial, decreasing = TRUE)
  bter_clot_gse <- keggSpecies(geneList_bter_clot, species) %>% filter(pvalue < 0.05)
  bter_clot_gse$species <- species
  bter_clot_gse$treatment <- treatment
  return(bter_clot_gse)
}

upGSE <- function(updown_list, species, treatment){
  bter_trial <- updown_list %>% filter(log2FoldChange > 0) %>% arrange(pvalue) %>% pull(pvalue)
  bter_trial <- -log10(bter_trial)
  names(bter_trial) <- gsub("gene-", "", (updown_list %>% filter(log2FoldChange > 0) %>% arrange(pvalue) %>% pull(X)))
  names(bter_trial) <- gsub("LOC", "", names(bter_trial))
  geneList_bter_clot <- sort(bter_trial, decreasing = TRUE)
  bter_clot_gse <- keggSpecies(geneList_bter_clot, species) %>% filter(pvalue < 0.05)
  bter_clot_gse$species <- species
  bter_clot_gse$treatment <- treatment
  return(bter_clot_gse)
}

allGSE <- function(updown_list, species, treatment){
  bter_trial <- updown_list %>% arrange(pvalue) %>% pull(pvalue)
  bter_trial <- -log10(bter_trial)
  names(bter_trial) <- gsub("gene-", "", (updown_list %>% arrange(pvalue) %>% pull(X)))
  names(bter_trial) <- gsub("LOC", "", names(bter_trial))
  geneList_bter_clot <- sort(bter_trial, decreasing = TRUE)
  bter_clot_gse <- keggSpecies(geneList_bter_clot, species) %>% filter(pvalue < 0.05)
  bter_clot_gse$species <- species
  bter_clot_gse$treatment <- treatment
  return(bter_clot_gse)
}


bter_clot_gse <- allGSE(updown_list = updown_list$bter_clot, species = "bter", treatment = "clot")
bter_clot_gse_down <- upGSE(updown_list = updown_list$bter_clot, species = "bter", treatment = "clot")
bter_clot_gse_up <- downGSE(updown_list = updown_list$bter_clot, species = "bter", treatment = "clot")

bter_sulf_gse <- allGSE(updown_list = updown_list$bter_sulf, species = "bter", treatment = "sulf")
bter_sulf_gse_down <- upGSE(updown_list = updown_list$bter_sulf, species = "bter", treatment = "sulf")
bter_sulf_gse_up <- downGSE(updown_list = updown_list$bter_clot, species = "bter", treatment = "sulf")

obic_clot_gse <- allGSE(updown_list = updown_list$obic_clot, species = "obb", treatment = "clot")
obic_clot_gse_down <- upGSE(updown_list = updown_list$obic_clot, species = "obb", treatment = "clot")
obic_clot_gse_up <- downGSE(updown_list = updown_list$obic_clot, species = "obb", treatment = "clot")

obic_sulf_gse <- allGSE(updown_list = updown_list$obic_sulf, species = "obb", treatment = "sulf")
obic_sulf_gse_down <- upGSE(updown_list = updown_list$obic_sulf, species = "obb", treatment = "sulf")
obic_sulf_gse_up <- downGSE(updown_list = updown_list$obic_sulf, species = "obb", treatment = "sulf")

vcar_sulf_gse <- allGSE(updown_list = updown_list$vcar_sulf, species = "vcd", treatment = "sulf")
vcar_sulf_gse_down <- upGSE(updown_list = updown_list$vcar_sulf, species = "vcd", treatment = "sulf")
vcar_sulf_gse_up <- downGSE(updown_list = updown_list$vcar_sulf, species = "vcd", treatment = "sulf")

lser_sulf_gse <- allGSE(updown_list = updown_list$lser_sulf, species = "lsq", treatment = "sulf")
lser_sulf_gse_down <- upGSE(updown_list = updown_list$lser_sulf, species = "lsq", treatment = "sulf")
lser_sulf_gse_up <- downGSE(updown_list = updown_list$lser_sulf, species = "lsq", treatment = "sulf")

vcar_clot_gse <- allGSE(updown_list = updown_list$vcar_clot, species = "vcd", treatment = "clot")
vcar_clot_gse_down <- upGSE(updown_list = updown_list$vcar_clot, species = "vcd", treatment = "clot")
vcar_clot_gse_up <- downGSE(updown_list = updown_list$vcar_clot, species = "vcd", treatment = "clot")

lser_clot_gse <- allGSE(updown_list = updown_list$lser_clot, species = "lsq", treatment = "clot")
lser_clot_gse_down <- upGSE(updown_list = updown_list$lser_clot, species = "lsq", treatment = "clot")
lser_clot_gse_up <- downGSE(updown_list = updown_list$lser_clot, species = "lsq", treatment = "clot")

all_full_kegg <- rbind(bter_clot_gse, obic_clot_gse,
                       lser_sulf_gse, lser_clot_gse, vcar_sulf_gse, vcar_clot_gse,
                       bter_clot_gse_down, bter_clot_gse_up,
                       obic_clot_gse_up, obic_clot_gse_down,
                       lser_clot_gse_up, lser_clot_gse_down, vcar_clot_gse_up, vcar_clot_gse_down,
                       lser_sulf_gse_up, vcar_sulf_gse_down, vcar_sulf_gse_down, vcar_sulf_gse_up)

all_full_kegg <- rbind(bter_clot_gse, obic_clot_gse, # bter_sulf_gse, obic_sulf_gse, 
                       lser_sulf_gse, lser_clot_gse, vcar_sulf_gse, vcar_clot_gse)

all_full_kegg$Description <- sub(" -.*", "", all_full_kegg$Description)

table(table(all_full_kegg$Description))
table(all_full_kegg$Description)[table(all_full_kegg$Description) == 4] # Carbon metabolism Phototransduction 
table(all_full_kegg$Description)[table(all_full_kegg$Description) == 5] # Fatty acid degradation Hippo signaling pathway 
table(all_full_kegg$treatment)

all_full_kegg$species_treatment <- paste(all_full_kegg$species, all_full_kegg$treatment, sep = "_")

# Define the specific order of species_treatment
species_treatment_order <- c("bter_clot", "obb_clot", "vcd_sulf", "vcd_clot", "lsq_sulf", "lsq_clot")

# Calculate the count of unique species_treatment for each description and sort within groups
description_order <- all_full_kegg %>%
  group_by(Description) %>%
  summarize(count = n_distinct(species_treatment), 
            species_treatment_order = min(match(species_treatment, species_treatment_order))) %>%
  arrange(desc(count), species_treatment_order) %>%
  pull(Description) %>% rev()

description_order <- c("Hedgehog signaling pathway", 
  "Glyoxylate and dicarboxylate metabolism", 
  "Glycine, serine and threonine metabolism",
  "Glutathione metabolism", 
  "SNARE interactions in vesicular transport",
  "Mitophagy", 
  "Fatty acid metabolism",
  "Ether lipid metabolism",
  "Arachidonic acid metabolism",
  "Pentose and glucuronate interconversions",
  "Notch signaling pathway",
  "Motor proteins",
  "Lysosome",
  "Inositol phosphate metabolism",
  "Histidine metabolism",
  "Cysteine and methionine metabolism", 
  "Polycomb repressive complex",
  "Lipoic acid metabolism",
  "Fanconi anemia pathway",
  "Drug metabolism",
  "Biosynthesis of cofactors",
  "ABC transporters",
  "Ribosome biogenesis in eukaryotes",
  "Ribosome",
  "Nucleotide metabolism",
  "DNA replication",
  "Apoptosis",
  "Aminoacyl-tRNA biosynthesis", 
  "beta-Alanine metabolism",
  "Valine, leucine and isoleucine degradation",
  "Tryptophan metabolism",
  "Starch and sucrose metabolism",
  "MAPK signaling pathway",
  "Hippo signaling pathway",
  "Glycolysis / Gluconeogenesis",
  "Glycerolipid metabolism",
  "Fructose and mannose metabolism",
  "Biosynthesis of nucleotide sugars",
  "Biosynthesis of amino acids",
  "Arginine and proline metabolism",
  "Amino sugar and nucleotide sugar metabolism",
  "Alanine, aspartate and glutamate metabolism",
  "Tyrosine metabolism",
  "Pentose phosphate pathway",
  "Fatty acid elongation", 
  "Efferocytosis",
  "ECM-receptor interaction",
  "Lysine degradation",
  "Ascorbate and aldarate metabolism",
  "Pantothenate and CoA biosynthesis",
  "Oxidative phosphorylation",
  "Nicotinate and nicotinamide metabolism", 
  "Propanoate metabolism",
  "Phototransduction",
  "Phagosome",
  "Insect hormone biosynthesis",
  "Arginine biosynthesis",
  "Toll and Imd signaling pathway",
  "One carbon pool by folate",
  "Pyruvate metabolism",
  "Pyrimidine metabolism",
  "Carbon metabolism",
  "Fatty acid degradation")

all_full_kegg %>%
  ggplot(aes(y = factor(species_treatment, levels = c("bter_clot", "obb_clot", 
                                                      "vcd_sulf", "vcd_clot", "lsq_sulf", "lsq_clot")), 
             x = factor(Description, levels = rev(description_order)), fill = -log10(pvalue + 0.000000001))) + 
  #geom_bar(stat = "identity", width = 0.5, alpha = 1) +
  geom_tile(colour = "#E7E7E7") + 
  # facet_grid(~factor(species, levels = c("bter", "obb", "vcd", "lsq")), drop = TRUE,
  #            scales = "free_x", space = "free") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    panel.spacing = unit(2.5, "mm")
  ) +
  scale_fill_gradientn(colors = rev(c("#000b0d", "#003842", "#05c3e6")))


# labs(x="log10p-value", y=" ") +
# facet_grid(~factor(species, levels = c("bter", "obb", "vcd", "lsq")) + 
#              factor(treatment, levels = c( "clot", "sulf")),
#            scales = "free_y", space = "free", drop = TRUE) +
# theme_minimal()

ggsave("~/2023-species_comparisons/paper_figures/GSEkegg_heatmap_reduced.pdf", units = "cm", width = 15, height = 20)
