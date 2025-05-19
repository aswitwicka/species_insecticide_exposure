# This script creates a summerised transcript per gene per sample table using .bam files created in STAR alignments and saves the output as a .csv file 

# Osmia bicornis 
bamFiles <- list.files(pattern="*.bam", 
                            path = "~/2023-osmia_bicornis/results/2023-04-star") # List of BAM files
bamFiles <- paste("~/2023-osmia_bicornis/results/2023-04-star", bamFiles, sep = "/")
annotationFile <- "~/2023-osmia_bicornis/input/obic/obic_genomic.gtf" # Annotation file in GTF format
outDir <- "" # Output directory

obic_counts <- featureCounts(
  files = bamFiles[-5],
  annot.ext = annotationFile,
  isGTFAnnotationFile = TRUE,
  isPairedEnd = TRUE,  # Set to TRUE if your data is paired-end
  minMQS = 0,          # Minimum mapping quality score
)
write.csv(obic_counts$counts, "~/2023-species_comparisons/tables/star_counts/FULL_SET_obic_star_featureCounts_notfilt.csv")

# Lucilia sericata 
bamFiles <- list.files(pattern="*.bam", 
                       path = "~/2022-Lsericata_brainDGE/results_update/star_alignment/tmp/star_combined_run") # List of BAM files
bamFiles <- paste("~/2022-Lsericata_brainDGE/results_update/star_alignment/tmp/star_combined_run", bamFiles, sep = "/")
annotationFile <- "~/2023-species_comparisons/genomes/lser/lser_genomic.gtf"              # Annotation file in GTF format
outDir <- ""                            # Output directory

lser_star_counts <- featureCounts(
  files = bamFiles,
  annot.ext = annotationFile,
  isGTFAnnotationFile = TRUE,
  isPairedEnd = TRUE,  # Set to TRUE if your data is paired-end
  minMQS = 0,          # Minimum mapping quality score
)

write.csv(lser_star_counts$counts, "~/2022-Lsericata_brainDGE/results_update/lser_NEW_STAR_counts_ALL_310324.csv")

# Bombus terrestris 
bamFiles <- list.files(pattern="*.bam", 
                       path = "~/2021-bter_concentration_exposure_experiment/output/2022-04-11-star/tmp") 
bamFiles <- paste("~/2021-bter_concentration_exposure_experiment/output/2022-04-11-star/tmp", bamFiles, sep = "/")
annotationFile <- "~/2023-species_comparisons/genomes/bter/bter.gtf"            
outDir <- ""                         

# Subset samples 
bamFiles <- bamFiles[grepl("FCLOT|FCONT|FSULF", bamFiles)]

counts_bter <- featureCounts(
  files = bamFiles,
  annot.ext = annotationFile,
  isGTFAnnotationFile = TRUE,
  isPairedEnd = TRUE,  
  minMQS = 0,          
)

write.csv(counts_bter$counts, "~/2023-species_comparisons/tables/star_counts/bter_star_featureCounts.csv")

# Vanessa cardui  
bamFiles <- list.files(pattern="*.bam", 
                       path = "~/2023-butterfly_data/results/2023-04-star/tmp/star_combined_run/") 
bamFiles <- paste("~/2023-butterfly_data/results/2023-04-star/tmp/star_combined_run", bamFiles, sep = "/")
annotationFile <- "~/2023-species_comparisons/genomes/vcar/vcar_genommic.gtf"         
outDir <- ""                      

vcar_counts <- featureCounts(
  files = bamFiles,
  annot.ext = annotationFile,
  isGTFAnnotationFile = TRUE,
  isPairedEnd = TRUE,  
  minMQS = 0,         
)

write.csv(vcar_counts$counts, "~/2023-species_comparisons/tables/star_counts/vcar_star_featureCounts.csv")

