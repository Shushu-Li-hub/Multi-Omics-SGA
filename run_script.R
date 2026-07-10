
save_dir <- 'H:/LSS/save_result'
############ Data process ###################
library(readxl)
# Data preprocessed
metabolomic_dir <- 'H:/LSS/data/metabolomic'
metagenomic_dir <- 'H:/LSS/data/metagenomic'
metadata_dir <- 'H:/LSS/data/metadata'
# Read data
# metabolomic
metabolomic_df <- read_xlsx(paste0(metabolomic_dir, '/combine.intensity.xlsx'))
metabolomic_postive <- metabolomic_df[grepl('pos', metabolomic_df$ID), ]
metabolomic_negative <- metabolomic_df[grepl('neg', metabolomic_df$ID), ]
# metagenomic
metagenomic_abun <- read.table(paste0(metagenomic_dir, '/taxonomy_abund.txt'), header = TRUE)
metagenomic_sample_filter <- read_excel(paste0(metagenomic_dir, '/sample_map.xlsx'))
# Count used for diversity
metagenomic_count <- read_xlsx(paste0(metagenomic_dir, '/1_Species_count_format.xlsx'))
# Unigenes
Unigenes_abund <- read.table(paste0(metagenomic_dir, '/4_Unigenes_abund.txt'), header = TRUE)
# clinical data
metadata <- read_xlsx(paste0(metadata_dir, '/meta_raw_filter.xlsx'))


# Get the abundance 
colnames(metabolomic_postive)
metabolomic_abun <- metabolomic_df
filter_colnames <- grepl("QC|AGA|SGA", colnames(metabolomic_abun))
metabolomic_abun[, filter_colnames] <- as.data.frame(lapply(metabolomic_abun[, filter_colnames], function(x) (x / sum(x)) * 100))
colSums(metabolomic_abun[, filter_colnames])

colnames(metabolomic_postive)
metabolomic_pos_abun <- metabolomic_postive
filter_colnames <- grepl("QC|SGA|AGA", colnames(metabolomic_pos_abun))
metabolomic_pos_abun[, filter_colnames] <- as.data.frame(lapply(metabolomic_pos_abun[, filter_colnames], function(x) (x / sum(x)) * 100))
colSums(metabolomic_pos_abun[, filter_colnames])
metabolomic_neg_abun <- metabolomic_negative
filter_colnames <- grepl("QC|SGA|AGA", colnames(metabolomic_neg_abun))
metabolomic_neg_abun[, filter_colnames] <- as.data.frame(lapply(metabolomic_neg_abun[, filter_colnames], function(x) (x / sum(x)) * 100))
colSums(metabolomic_neg_abun[, filter_colnames])

colnames(metagenomic_count)
metagenomic_abun <- metagenomic_count
filter_colnames <- grepl("AGA|SGA", colnames(metagenomic_count))
metagenomic_abun[, filter_colnames] <- as.data.frame(lapply(metagenomic_abun[, filter_colnames], function(x) (x / sum(x)) * 100))
colSums(metagenomic_abun[, filter_colnames])

# Get the average relative abundance
metabolomic_abun[['Average_Mean']] <- rowMeans(metabolomic_abun[, grepl('AGA|SGA', colnames(metabolomic_abun))])
metabolomic_pos_abun[['Average_Mean']] <- rowMeans(metabolomic_pos_abun[, grepl('AGA|SGA', colnames(metabolomic_pos_abun))])
metabolomic_neg_abun[['Average_Mean']] <- rowMeans(metabolomic_neg_abun[, grepl('AGA|SGA', colnames(metabolomic_neg_abun))])
metagenomic_abun[['Average_Mean']] <- rowMeans(metagenomic_abun[, grepl('AGA|SGA', colnames(metagenomic_abun))])
sum(metabolomic_pos_abun[['Average_Mean']])
Unigenes_abund[['Average_Mean']] <- rowMeans(Unigenes_abund[, grepl('AGA|SGA', colnames(Unigenes_abund))])
sum(Unigenes_abund[['Average_Mean']])

# Filter: species or KOs with an average relative abundance above 10^−7(0.00001%) were considered in the analyses.
metabolomic_pos_filter <- metabolomic_pos_abun[metabolomic_pos_abun$Average_Mean > 0.00001, ]
metabolomic_neg_filter <- metabolomic_neg_abun[metabolomic_neg_abun$Average_Mean > 0.00001, ]
metagenomic_filter <- metagenomic_abun[metagenomic_abun$Average_Mean > 0.00001, ]
print(paste0("metabolomic_pos_filter: ", as.character(nrow(metabolomic_pos_abun)), " to ",  as.character(nrow(metabolomic_pos_filter))))
print(paste0("metabolomic_neg_filter: ", as.character(nrow(metabolomic_neg_abun)), " to ",  as.character(nrow(metabolomic_neg_filter))))
print(paste0("metagenomic_filter: ", as.character(nrow(metagenomic_abun)), " to ",  as.character(nrow(metagenomic_filter))))

# metagenomic filtered count used for diversity
metagenomic_filter_count <- metagenomic_count[metagenomic_count$Species %in% metagenomic_filter$Species, ]
# Get unigenes_kegg data
library(readr)
unigenes_kegg_df <- read_delim(paste0(metagenomic_dir, '/Unigenes_KEGG.txt'), delim = "\t", col_names = TRUE, col_types = cols(), trim_ws = TRUE)

# script run
directory <- 'H:/LSS/project'
source(paste0(directory, '/figure1.R'))  # microbiome 
source(paste0(directory, '/figure2.R')) # KOs 
source(paste0(directory, '/figure3.R')) # metaolites
source(paste0(directory, '/figure4.R')) # sankey
source(paste0(directory, '/figure5.R')) # heatmap

