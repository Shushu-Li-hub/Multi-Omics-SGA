
# save dir
setwd(save_dir)
# Figure 4 Sankey. Species-Functions-Metabiolites
############################
#### Taxa Diff Analysis  ###
############################
log2fc <- 0.585
test_method <- 'rfit.p_value'

# Load script
utils_direcotry <- 'H:/LSS/project/utils'
source(paste0(utils_direcotry, '/new_calc_microbiome_diff.R'))
# Get Kos Data
colnames(metagenomic_abun)
# Select samples
sample_dir <- 'H:/LSS/data'
sample_origin <- read_excel(paste0(sample_dir, '/sample_origin.xlsx'))
colnames(sample_origin)
microbiome_sample <- na.omit(sample_origin[['P20230730001_metagenomic']])
SGA_samples <- microbiome_sample[grepl('SGA', microbiome_sample)]
AGA_samples <- microbiome_sample[grepl('AGA', microbiome_sample)]
# Filtered data
AGA_SGA_group <- rbind(data.frame(sample = AGA_samples, group = rep('AGA', length(AGA_samples))), 
                       data.frame(sample = SGA_samples, group = rep('SGA', length(SGA_samples))))
AGA_SGA_group$group <- factor(AGA_SGA_group$group, levels = c('AGA', 'SGA'))
Taxa_diff <- new_calc_microbiome_diff(metagenomic_abun, AGA_SGA_group, metadata_covariates)
Taxa_diff[['Species']] <- metagenomic_abun$Species
# relative abundance > 0.00001 Taxa_diff
Taxa_diff <- Taxa_diff[(Taxa_diff$mean_AGA > 0.00001) & (Taxa_diff$mean_SGA > 0.00001) , ]
Taxa_diff <- Taxa_diff[!grepl('unclassified|uncultured', Taxa_diff$Species), ]
# FDR 20260620
Taxa_diff$rfit.fdr_values <- p.adjust(Taxa_diff[[test_method]], method = "fdr")  # FDR correction
# Difference taxa
Taxa_filter2 <- Taxa_diff[Taxa_diff[[test_method]] < 0.05 & abs(Taxa_diff$log_fold_change) >= log2fc, ]
nrow(Taxa_filter2)
# lefse + diff taxa
lefse_marker <- lefse_diff_select$Species
explore_taxa <- unique(c(lefse_marker, Taxa_filter2$Species))


############################
######## Sankey   ##########
############################
######### Species-Functions ##################
utils_direcotry <- 'H:/LSS/project/utils'
source(paste0(utils_direcotry, '/calculate_spearman.R'))
source(paste0(utils_direcotry, '/calculate_sankey_data.R'))
database <- list()
# Note: Here Need filter species and select
Taxa_filter_select <- Taxa_diff[Taxa_diff$Species %in% explore_taxa, ]
microbiome_data <- Taxa_filter_select
rownames(microbiome_data) <- Taxa_filter_select$Species
# Unigenes diff
Unigenes_Difference_all <- Unigenes_Difference_all[order(Unigenes_Difference_all[[test_method]]), ]
Unigenes_Difference_select <-  Unigenes_Difference_all
functions_data <- Unigenes_Difference_select
colnames(microbiome_data)
colnames(functions_data)
microbiome_sample <- na.omit(sample_origin[!is.na(sample_origin$P20230730001_metagenomic), ][['P20230730001_metagenomic']])
functions_sample <- microbiome_sample
result <- calculate_spearman(microbiome_data, functions_data, microbiome_sample, functions_sample)
database$taxa_function$spearman_pvalue <- result[["pvalue"]]
database$taxa_function$spearman_p_adjusted_value <- result[["adjusted_p_values"]]
database$taxa_function$spearman_correlation <- result[["correlation"]]
spearman_correlation <- result[["correlation"]]
spearman_correlation_abs <- abs(spearman_correlation)
Filter_matrix <- spearman_correlation_abs > 0.5 & result[["pvalue"]] < 0.05
spearman_filter_correlation <- as.matrix(spearman_correlation) * Filter_matrix
database$taxa_function$spearman_filter_correlation <- spearman_filter_correlation


######### Species-Metabiolites ##################
microbiome_data <- Taxa_filter_select
metabolites_data <- annot_kegg_diff             
rownames(microbiome_data) <- Taxa_filter_select$Species
rownames(metabolites_data) <- metabolites_data$KEGG_Metabolite
colnames(microbiome_data)
colnames(metabolites_data)
colnames(sample_origin)
microbiome_sample <- na.omit(sample_origin[!is.na(sample_origin$P20230730001_metagenomic), ][['P20230730001_metagenomic']])
metabolites_sample <- na.omit(sample_origin[!is.na(sample_origin$P20230730001_metagenomic), ][['P20230730002_metabolism']])
result <- calculate_spearman(microbiome_data, metabolites_data, microbiome_sample, metabolites_sample)
database$taxa_metabiolites$spearman_pvalue <- result[["pvalue"]]
database$taxa_metabiolites$spearman_p_adjusted_value <- result[["adjusted_p_values"]]
database$taxa_metabiolites$spearman_correlation <- result[["correlation"]]
spearman_correlation <- result[["correlation"]]
spearman_correlation_abs <- abs(spearman_correlation)
Filter_matrix <- spearman_correlation_abs > 0.5 & result[["pvalue"]] < 0.05
spearman_filter_correlation <- as.matrix(spearman_correlation) * Filter_matrix
database$taxa_metabiolites$spearman_filter_correlation <- spearman_filter_correlation


######### Functions-Metabiolites ##################
functions_data <- Unigenes_Difference_select
metabolites_data <- annot_kegg_diff
rownames(metabolites_data) <- metabolites_data$KEGG_Metabolite
colnames(functions_data)
colnames(metabolites_data)
metabolites_sample <- na.omit(sample_origin[!is.na(sample_origin$P20230730001_metagenomic), ][['P20230730002_metabolism']])
microbiome_sample <- na.omit(sample_origin[!is.na(sample_origin$P20230730001_metagenomic), ][['P20230730001_metagenomic']])
functions_sample <- microbiome_sample
result <- calculate_spearman(functions_data, metabolites_data, functions_sample, metabolites_sample)
database$functions_metabiolites$spearman_pvalue <- result[["pvalue"]]
database$functions_metabiolites$spearman_p_adjusted_value <- result[["adjusted_p_values"]]
database$functions_metabiolites$spearman_correlation <- result[["correlation"]]
spearman_correlation <- result[["correlation"]]
spearman_correlation_abs <- abs(spearman_correlation)
Filter_matrix <- spearman_correlation_abs > 0.5 & result[["pvalue"]] < 0.05
spearman_filter_correlation <- as.matrix(spearman_correlation) * Filter_matrix
database$functions_metabiolites$spearman_filter_correlation <- spearman_filter_correlation


########## Microbiome Sankey: ############
# sankey data
library(ggplot2)
library(ggalluvial)
sankey_data_A_B <- database$taxa_function$spearman_correlation
relations_A_B <- data.frame(
  source = rep(rownames(sankey_data_A_B), each = ncol(sankey_data_A_B)),
  target = rep(colnames(sankey_data_A_B), times = nrow(sankey_data_A_B)),
  value = as.vector(t(sankey_data_A_B))
)
relations_A_B <- relations_A_B[relations_A_B$value != 0, ]
colnames(relations_A_B) <- c("Microbial_species", "Functions", "Correlation_A_B")
# sankey data
sankey_data_B_C <- database$functions_metabiolites$spearman_correlation
relations_B_C <- data.frame(
  source = rep(rownames(sankey_data_B_C), each = ncol(sankey_data_B_C)),
  target = rep(colnames(sankey_data_B_C), times = nrow(sankey_data_B_C)),
  value = as.vector(t(sankey_data_B_C))
)
relations_B_C <- relations_B_C[relations_B_C$value != 0, ]
colnames(relations_B_C) <- c("Functions", "Metabolites", "Correlation_B_C")
# combind
colnames(relations_A_B)
merge_relations <- merge(relations_A_B, relations_B_C, by = "Functions", all = FALSE)
merge_relations[['Count']] <- rep(1, nrow(merge_relations))
output_merge_relations <- merge_relations

# pvalue and correlation
A_B_spearman_pvalue <- database$taxa_function$spearman_pvalue
B_C_spearman_pvalue <- database$functions_metabiolites$spearman_pvalue
A_C_spearman_pvalue <- database$taxa_metabiolites$spearman_pvalue
A_C_spearman_correlation <- database$taxa_metabiolites$spearman_correlation

colnames(output_merge_relations)

colnames(A_B_spearman_pvalue)
rownames(A_B_spearman_pvalue)
output_merge_relations[["pvalue_A_B"]] <- unlist(lapply(1:nrow(output_merge_relations), 
                                                        function(x) A_B_spearman_pvalue[output_merge_relations$Microbial_species[x], output_merge_relations$Functions[x]]))
colnames(B_C_spearman_pvalue)
rownames(B_C_spearman_pvalue)
output_merge_relations[["pvalue_B_C"]] <- unlist(lapply(1:nrow(output_merge_relations), 
                                                        function(x) B_C_spearman_pvalue[output_merge_relations$Functions[x], output_merge_relations$Metabolites[x]]))
colnames(A_C_spearman_pvalue)
rownames(A_C_spearman_pvalue)
output_merge_relations[["pvalue_A_C"]] <- unlist(lapply(1:nrow(output_merge_relations), 
                                                        function(x) A_C_spearman_pvalue[output_merge_relations$Microbial_species[x], output_merge_relations$Metabolites[x]]))

colnames(A_C_spearman_correlation)
rownames(A_C_spearman_correlation)
output_merge_relations[["Correlation_A_C"]] <- unlist(lapply(1:nrow(output_merge_relations), 
                                                             function(x) A_C_spearman_correlation[output_merge_relations$Microbial_species[x], output_merge_relations$Metabolites[x]]))

# FDR
A_B_spearman_fdr_pvalue <- database$taxa_function$spearman_p_adjusted_value
B_C_spearman_fdr_pvalue <- database$functions_metabiolites$spearman_p_adjusted_value
A_C_spearman_fdr_pvalue <- database$taxa_metabiolites$spearman_p_adjusted_value

colnames(A_B_spearman_fdr_pvalue)
rownames(A_B_spearman_fdr_pvalue)
output_merge_relations[["fdr_pvalue_A_B"]] <- unlist(lapply(1:nrow(output_merge_relations), 
                                                        function(x) A_B_spearman_fdr_pvalue[output_merge_relations$Microbial_species[x], output_merge_relations$Functions[x]]))
colnames(B_C_spearman_fdr_pvalue)
rownames(B_C_spearman_fdr_pvalue)
output_merge_relations[["fdr_pvalue_B_C"]] <- unlist(lapply(1:nrow(output_merge_relations), 
                                                        function(x) B_C_spearman_fdr_pvalue[output_merge_relations$Functions[x], output_merge_relations$Metabolites[x]]))
colnames(A_C_spearman_fdr_pvalue)
rownames(A_C_spearman_fdr_pvalue)
output_merge_relations[["fdr_pvalue_A_C"]] <- unlist(lapply(1:nrow(output_merge_relations), 
                                                        function(x) A_C_spearman_fdr_pvalue[output_merge_relations$Microbial_species[x], output_merge_relations$Metabolites[x]]))


# filter
# value
filter_output_merge_relations <- output_merge_relations[abs(output_merge_relations$Correlation_A_B) >= 0.5 &
                                                          abs(output_merge_relations$Correlation_B_C) >= 0.5 &
                                                          abs(output_merge_relations$Correlation_A_C) >= 0.5 &
                                                          abs(output_merge_relations$pvalue_A_B) < 0.05 &
                                                          abs(output_merge_relations$pvalue_B_C) < 0.05 &
                                                          abs(output_merge_relations$pvalue_A_C) < 0.05 ,  ]


# write.csv(filter_output_merge_relations, 'filter_output_merge_relations.csv')
# Select speceis lesfe and taxa sum.
# + -
# View focus species 
relations_species <- unique(filter_output_merge_relations$Microbial_species)
relations_species
Taxa_filter_select[['Sum']] <- rowSums(Taxa_filter_select[, c(AGA_samples, SGA_samples)])
Taxa_filter_select <- Taxa_filter_select[order(Taxa_filter_select$Sum, decreasing = TRUE), ]
Select_speceis <- Taxa_filter_select$Species[1:50]
merge_relations <- filter_output_merge_relations[filter_output_merge_relations$Microbial_species %in% Select_speceis, ]
merge_relations <- merge_relations[merge_relations$Functions %in% Unigenes_Difference$id, ]
unique(merge_relations$Microbial_species)
unique(merge_relations$Metabolites)
merge_relations$Metabolites <- factor(merge_relations$Metabolites, levels = c("5'-Methylthioadenosine", "Taurocholate","Gibberellin A12 aldehyde", "beta-Alanine"))


# colnames(merge_relations)
# library(dplyr) 
# relation_data <- group_by(merge_relations , `Microbial species`, Functions, Metabolites) %>% summarise(., count = n())
library(ggsci)
colors <- pal_lancet("lanonc")(length(unique(merge_relations$Metabolites)))
colors <- c("#ED0000FF", colors[!colors %in% "#ED0000FF"])
colnames(unigenes_kegg_df)
kos2genename <- unigenes_kegg_df[, c("GeneName", "KOEntry")]
kos2genename <- kos2genename[!duplicated(kos2genename$KOEntry), ]
merge_relations[['Genename']] <- replace_names(merge_relations$Functions, kos2genename$KOEntry, kos2genename$GeneName)
merge_relations$Functions <- paste0(merge_relations$Functions, " (", merge_relations$Genename, ")")
merge_relations$Microbial_species <- gsub("s__", "", merge_relations$Microbial_species)
species_count <- as.data.frame(table(merge_relations$Microbial_species))
species_count <- species_count[species_count$Freq == 1, ]
merge_relations[merge_relations$Microbial_species %in% species_count$Var1, ][['Count']] <- 2
# merge_relations$Count <- 1

# Create the Sankey diagram
sankey <- ggplot(data = merge_relations,
                 aes(axis1 = Microbial_species, axis2 = Functions, axis3 = Metabolites,
                     y = Count)) +
  scale_x_discrete(limits = c("Microbial_species", "Functions", "Metabolites"), position = "top") +
  geom_alluvium(aes(fill = Metabolites), curve_type = "quintic", alpha = 0.8) +
  geom_stratum(width = .42) + # must 0.4
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  theme(legend.position = "none" , 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        #axis.text = element_blank(),
        axis.text.y =  element_blank(),
        axis.text.x = element_text(size = 24, color = "black"),
        axis.title = element_blank()) +
  scale_fill_manual(values = colors)  +
  ggtitle(" ")
print(sankey)
png('Sankey_plot.png', width = 4000, height = 3000, res=300)  # Adjust width, height, and resolution as needed
print(sankey)
dev.off()
tiff('Sankey_plot.tif', width = 13.3, height = 10, units = "in", res = 300)
print(sankey)
dev.off()


library(openxlsx)
wb <- createWorkbook()
options(warn = -1)
addWorksheet(wb, 'spearman_pvalue')
addWorksheet(wb, 'spearman_p_adjusted_value')
addWorksheet(wb, 'spearman_correlation')
writeData(wb, 'spearman_pvalue', database$taxa_metabiolites$spearman_pvalue, startCol = 1, startRow = 1, rowNames = TRUE, colNames = TRUE)  
writeData(wb, 'spearman_p_adjusted_value', database$taxa_metabiolites$spearman_p_adjusted_value, startCol = 1, startRow = 1, rowNames = TRUE, colNames = TRUE) 
writeData(wb, 'spearman_correlation', database$taxa_metabiolites$spearman_correlation, startCol = 1, startRow = 1, rowNames = TRUE, colNames = TRUE)  
options(warn = 1)
saveWorkbook(wb, "taxa_metabiolites_spearman.xlsx", overwrite = TRUE)

wb <- createWorkbook()
options(warn = -1)
addWorksheet(wb, 'spearman_pvalue')
addWorksheet(wb, 'spearman_p_adjusted_value')
addWorksheet(wb, 'spearman_correlation')
writeData(wb, 'spearman_pvalue', database$functions_metabiolites$spearman_pvalue, startCol = 1, startRow = 1, rowNames = TRUE, colNames = TRUE)  
writeData(wb, 'spearman_p_adjusted_value', database$functions_metabiolites$spearman_p_adjusted_value, startCol = 1, startRow = 1, rowNames = TRUE, colNames = TRUE) 
writeData(wb, 'spearman_correlation', database$functions_metabiolites$spearman_correlation, startCol = 1, startRow = 1, rowNames = TRUE, colNames = TRUE)  
options(warn = 1)
saveWorkbook(wb, "functions_metabiolites_spearman.xlsx", overwrite = TRUE)


wb <- createWorkbook()
options(warn = -1)
addWorksheet(wb, 'spearman_pvalue')
addWorksheet(wb, 'spearman_p_adjusted_value')
addWorksheet(wb, 'spearman_correlation')
writeData(wb, 'spearman_pvalue', database$taxa_function$spearman_pvalue, startCol = 1, startRow = 1, rowNames = TRUE, colNames = TRUE)  
writeData(wb, 'spearman_p_adjusted_value', database$taxa_function$spearman_p_adjusted_value, startCol = 1, startRow = 1, rowNames = TRUE, colNames = TRUE) 
writeData(wb, 'spearman_correlation', database$taxa_function$spearman_correlation, startCol = 1, startRow = 1, rowNames = TRUE, colNames = TRUE)  
options(warn = 1)
saveWorkbook(wb, "taxa_function_spearman.xlsx", overwrite = TRUE)
write.xlsx(annot_kegg_diff, "Compound_kegg_name.xlsx", overwrite = TRUE)


# Note: Here only select p_value < 0.05 ;
# library(officer)
# doc <- read_docx()
# # Add the table to the document
# doc <- body_add_table(doc, value = table_data)
# # Save the document
# print(doc, target = "example.docx")
# Species
Table_S2 <- data.frame(
  Species = gsub("s__", "", Taxa_filter2$Species),
  `P value` = round(Taxa_filter2[[test_method]] , 4),
  `FDR-adjusted P value` = round(Taxa_filter2[['rfit.fdr_values']], 4) ,
  `Fold change (SGA vs AGA)` = round(Taxa_filter2$fold_change, 4)
)


# Calculate different of genus 
metagenomic_count_directory <- 'H:/LSS/data/metagenomic'
metagenomic_count_taxa <- read_xlsx(paste0(metagenomic_count_directory, '/Species_count_format_taxa.xlsx'))
SGA_samples <- microbiome_sample[grepl('SGA', microbiome_sample)]
AGA_samples <- microbiome_sample[grepl('AGA', microbiome_sample)]
library(stringr)
taxo_name <- data.frame(str_split(metagenomic_count_taxa[["Taxonomy"]], ";", simplify = TRUE))
colnames(taxo_name) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
metagenomic_abun_cbind <- cbind(metagenomic_abun, taxo_name)
colSums(metagenomic_abun_cbind[, c(AGA_sample, SGA_sample)])
temp_data <- metagenomic_abun_cbind[, c("Genus", AGA_sample, SGA_sample)]
metagenomic_abun2 <- aggregate(. ~ Genus, data = temp_data, FUN = sum)
colSums(metagenomic_abun2[, c(AGA_samples, SGA_samples)])
AGA_SGA_group <- rbind(data.frame(sample = AGA_samples, group = rep('AGA', length(AGA_samples))), 
                       data.frame(sample = SGA_samples, group = rep('SGA', length(SGA_samples))))
AGA_SGA_group$group <- factor(AGA_SGA_group$group, levels = c('AGA', 'SGA'))
Genus_Taxa_diff <- new_calc_microbiome_diff(metagenomic_abun2, AGA_SGA_group, metadata_covariates)
Genus_Taxa_diff[['Genus']] <- metagenomic_abun2$Genus

# relative abundance > 0.00001 Genus_Taxa_diff
Genus_Taxa_diff <- Genus_Taxa_diff[(Genus_Taxa_diff$mean_AGA > 0.00001) & (Genus_Taxa_diff$mean_SGA > 0.00001) , ]
Genus_Taxa_diff <- Genus_Taxa_diff[!grepl('unclassified|uncultured', Genus_Taxa_diff$Genus), ]
# FDR 20260620
Genus_Taxa_diff$rfit.fdr_values <- p.adjust(Genus_Taxa_diff[[test_method]], method = "fdr")  # FDR correction
# Difference taxa
Genus_Taxa_diff <- Genus_Taxa_diff[Genus_Taxa_diff[[test_method]] < 0.05 & abs(Genus_Taxa_diff$log_fold_change) >= log2fc, ]
# Genus
Table_S4 <- data.frame(
  Genus = gsub("g__", "", Genus_Taxa_diff$Genus),
  `P value` = round(Genus_Taxa_diff[[test_method]] , 4),
  `FDR-adjusted P value` = round(Genus_Taxa_diff[['rfit.fdr_values']], 4) ,
  `Fold change (SGA vs AGA)` = round(Genus_Taxa_diff$fold_change, 4)
)

# Calculate differenct of Phylum
metagenomic_count_directory <- 'H:/LSS/data/metagenomic'
metagenomic_count_taxa <- read_xlsx(paste0(metagenomic_count_directory, '/Species_count_format_taxa.xlsx'))
SGA_samples <- microbiome_sample[grepl('SGA', microbiome_sample)]
AGA_samples <- microbiome_sample[grepl('AGA', microbiome_sample)]
library(stringr)
taxo_name <- data.frame(str_split(metagenomic_count_taxa[["Taxonomy"]], ";", simplify = TRUE))
colnames(taxo_name) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
metagenomic_abun_cbind <- cbind(metagenomic_abun, taxo_name)
colSums(metagenomic_abun_cbind[, c(AGA_sample, SGA_sample)])
temp_data <- metagenomic_abun_cbind[, c("Phylum", AGA_sample, SGA_sample)]
metagenomic_abun3 <- aggregate(. ~ Phylum, data = temp_data, FUN = sum)
colSums(metagenomic_abun3[, c(AGA_samples, SGA_samples)])
AGA_SGA_group <- rbind(data.frame(sample = AGA_samples, group = rep('AGA', length(AGA_samples))), 
                       data.frame(sample = SGA_samples, group = rep('SGA', length(SGA_samples))))
AGA_SGA_group$group <- factor(AGA_SGA_group$group, levels = c('AGA', 'SGA'))
Phylum_Taxa_diff <- new_calc_microbiome_diff(metagenomic_abun3, AGA_SGA_group, metadata_covariates)
Phylum_Taxa_diff[['Phylum']] <- metagenomic_abun3$Phylum

# relative abundance > 0.00001 Phylum_Taxa_diff
Phylum_Taxa_diff <- Phylum_Taxa_diff[(Phylum_Taxa_diff$mean_AGA > 0.00001) & (Phylum_Taxa_diff$mean_SGA > 0.00001) , ]
Phylum_Taxa_diff <- Phylum_Taxa_diff[!grepl('unclassified|uncultured', Phylum_Taxa_diff$Phylum), ]
# FDR 20260620
Phylum_Taxa_diff$rfit.fdr_values <- p.adjust(Phylum_Taxa_diff[[test_method]], method = "fdr")  # FDR correction
# Difference taxa
Phylum_Taxa_diff <- Phylum_Taxa_diff[Phylum_Taxa_diff[[test_method]] < 0.05 & abs(Phylum_Taxa_diff$log_fold_change) >= log2fc, ]

# phylum 
Table_S3 <- data.frame(
  Phylum = gsub("p__", "", Phylum_Taxa_diff$Phylum),
  `P value` = round(Phylum_Taxa_diff[[test_method]] , 4),
  `FDR-adjusted P value` = round(Phylum_Taxa_diff$rfit.fdr_values, 4) ,
  `Fold change (SGA vs AGA)` = round(Phylum_Taxa_diff$fold_change, 4)
)


# Unigenes Kos
Filter_Unigenes_Difference <- Unigenes_filter[(Unigenes_filter[[test_method]] < 0.05) & (abs(Unigenes_filter$log_fold_change)>log2fc), ]
Table_S5 <- data.frame(
  Taxon = gsub("p__", "", Filter_Unigenes_Difference$id),
  `P value` = round(Filter_Unigenes_Difference[[test_method]] , 4),
  `FDR-adjusted P value` = round(Filter_Unigenes_Difference$rfit.fdr_values, 4) ,
  `Fold change (SGA vs AGA)` = round(Filter_Unigenes_Difference$fold_change, 4)
)


# Spearman correlation
merge_relations <- merge_relations[order(merge_relations$Functions), ]
merge_relations <- merge_relations[order(merge_relations$Metabolites), ]
Table_S6 <- data.frame(
  Species = merge_relations$Microbial_species,
  Fuction = merge_relations$Functions,
  Metabolite = merge_relations$Metabolites,
  `Species-function correlation` = round(merge_relations$Correlation_A_B, 4),
  `Species-function p_value` = round(merge_relations$pvalue_A_B, 4),
  `Species-function fdr p_value` = round(merge_relations$fdr_pvalue_A_B, 4),
  `Function-metabolite correlation` = round(merge_relations$Correlation_B_C, 4),
  `Function-metabolite p_value` = round(merge_relations$pvalue_B_C, 4),
  `Function-metabolite fdr p_value` = round(merge_relations$fdr_pvalue_B_C, 4),
  `Species-metabolite correlation` = round(merge_relations$Correlation_A_C, 4),
  `Species-metabolite p_value` = round(merge_relations$pvalue_A_C, 4),
  `Species-metabolite fdr p_value` = round(merge_relations$fdr_pvalue_A_C, 4)
)
unique(Table_S6$Fuction)


library(openxlsx)
wb <- createWorkbook()
## Add worksheets
addWorksheet(wb, "Table 1")
addWorksheet(wb, "Table 2")
addWorksheet(wb, "Table S1")
addWorksheet(wb, "Table S2")
addWorksheet(wb, "Table S3")
addWorksheet(wb, "Table S4")
addWorksheet(wb, "Table S5")
addWorksheet(wb, "Table S6")
addWorksheet(wb, "Table S7")
addWorksheet(wb, "Table S8")
addWorksheet(wb, "Table S9")
addWorksheet(wb, "Table S10")

# Species
Table_S1 <- data.frame(
  Species = gsub("s__", "", lefse_diff_select$Species),
  `P.adj` = round(lefse_diff_select$P.adj , 4),
  `LDA` = round(lefse_diff_select$LDA, 4) 
)

column_name1 <- "Table S1. Comparison of relative abundance of gut microbiota between term SGA and term AGA groups based on LEfSe analysis (P value < 0.05 and LDA score > 3.0) at the species level"
column_name2 <- colnames(Table_S1)
Table_S1 <- rbind(c(column_name1, rep("", 3)), column_name2, Table_S1)
column_name1 <- "Table S2. Comparison of relative abundance of gut microbiota between term SGA and term AGA groups at the species level"
column_name2 <- colnames(Table_S2)
Table_S2 <- rbind(c(column_name1, rep("", 3)), column_name2, Table_S2)
column_name1 <- "Table S3. Comparison of relative abundance of gut microbiota between term SGA and term AGA groups at the phylum level"
column_name2 <- colnames(Table_S3)
Table_S3 <- rbind(c(column_name1, rep("", 3)), column_name2, Table_S3)
column_name1 <- "Table S4. Comparison of relative abundance of gut microbiota between term SGA and term AGA groups at the genus level"
column_name2 <- colnames(Table_S4)
Table_S4 <- rbind(c(column_name1, rep("", 3)), column_name2, Table_S4)
column_name1 <- "Table S5. Comparison of relative abundance of KOs between term SGA and term AGA group"
column_name2 <- colnames(Table_S5)
Table_S5 <- rbind(c(column_name1, rep("", 2)), column_name2, Table_S5)
column_name1 <-  "Table S6. Interrelationship between gut microbial species, functions and metabolites"
# Remove unused levels from factors
Table_S6$Metabolite <- as.character(Table_S6$Metabolite)
column_name2 <- colnames(Table_S6)
Table_S6 <- rbind(c(column_name1, rep("", 5)), column_name2, Table_S6)

writeData(wb, "Table 1", table1.1, startCol = 1, startRow = 1, rowNames = FALSE,  colNames = FALSE)
writeData(wb, "Table 2", table2.1, startCol = 1, startRow = 1, rowNames = FALSE,  colNames = FALSE)
writeData(wb, "Table S1", Table_S1, startCol = 1, startRow = 1, rowNames = FALSE,  colNames = FALSE)
writeData(wb, "Table S2", Table_S2, startCol = 1, startRow = 1, rowNames = FALSE,  colNames = FALSE)
writeData(wb, "Table S3", Table_S3, startCol = 1, startRow = 1, rowNames = FALSE,  colNames = FALSE)
writeData(wb, "Table S4", Table_S4, startCol = 1, startRow = 1, rowNames = FALSE,  colNames = FALSE)
writeData(wb, "Table S5", Table_S5, startCol = 1, startRow = 1, rowNames = FALSE,  colNames = FALSE)
writeData(wb, "Table S6", Table_S6, startCol = 1, startRow = 1, rowNames = FALSE,  colNames = FALSE)









