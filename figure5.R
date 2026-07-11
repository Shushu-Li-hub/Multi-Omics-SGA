# save dir
setwd(save_dir)

# Color
pheamap_gradient_color2 <- rev(brewer.rdylbu(n = 100))
# Figure 5.2 Spearman: Metabolites and clinical indicators(Heatmap)
############################
######## Heatmap  ##########
############################
######### Metabiolites-Phenotype #############
Heatmap_database <- list()
colnames(annot_kegg_diff)
metabolites_data <- annot_kegg_diff[!duplicated(annot_kegg_diff$KEGG_Metabolite), ]
rownames(metabolites_data) <- metabolites_data$KEGG_Metabolite
phenotyop_directory <- 'H:/LSS/data/metadata'
# Note here use raw data.
phenotype_data <- read.xlsx(paste0(phenotyop_directory, '/meta_raw_zscore.xlsx'))
phenotype_data <- data.frame(t(phenotype_data))
colnames(phenotype_data) <- phenotype_data['SampleID',]
#  "WLZ" "WAZ""LAZ""BAZ""HAZ" 'gesell'
######born######
born_names <- rownames(phenotype_data)[5:12]
born_names
rename_born_names <- c("Neonatal_weight","Head_circumference", "Height" ,
                       "WLZ", "WAZ", "LAZ", "BAZ", "HAZ")
rownames(phenotype_data[born_names, ])
phenotype_data_born <- as.data.frame(lapply(phenotype_data[born_names, ], as.numeric))
rownames(phenotype_data_born) <- rename_born_names
######3 Month######
# add bfa-3
m3_names <- rownames(phenotype_data)[13:21]
m3_names
rename_m3_names <- c("3Mo_Weight", "3Mo_Height", "3Mo_Head_circumference",
                     "3Mo_WLZ", "3Mo_WAZ", "3Mo_LAZ", 
                     "3Mo_BAZ", "3Mo_HAZ", "3Mo_GV")
rownames(phenotype_data[m3_names, ])
phenotype_data_3mo <- as.data.frame(lapply(phenotype_data[m3_names, ], as.numeric))
rownames(phenotype_data_3mo) <- rename_m3_names
######6 Month######
m6_names <- rownames(phenotype_data)[22:30]
m6_names
rename_m6_names <- c("6Mo_Weight", "6Mo_Height", "6Mo_Head_circumference",
                     "6Mo_WLZ", "6Mo_WAZ", "6Mo_LAZ", 
                     "6Mo_BAZ", "6Mo_HAZ", "6Mo_GV")
rownames(phenotype_data[m6_names, ])
phenotype_data_6mo <- as.data.frame(lapply(phenotype_data[m6_names, ], as.numeric))
rownames(phenotype_data_6mo) <- rename_m6_names
######gesell Neural development######
# combine_pheno <- rbind(phenotype_data_born, phenotype_data_3mo, phenotype_data_6mo, phenotype_data_gesell)
# select 6mo indicator
combine_pheno <- rbind(phenotype_data_6mo)
rownames(combine_pheno) <- c("Weight","Height","Head_circumference",
                             "WLZ", "WAZ" , "LAZ",
                             "BAZ","HAZ","GV" )

phenotype_data <- combine_pheno
colnames(metabolites_data)
colnames(phenotype_data)
sample_dir <- 'H:/LSS/data'
sample_origin <- read_excel(paste0(sample_dir, '/sample_origin.xlsx'))
colnames(sample_origin)
phenotype_sample <- na.omit(sample_origin[!is.na(sample_origin$Clinical_features) & !is.na(sample_origin$P20230730002_metabolism) , ][['Clinical_features']])
metabolism_sample <- na.omit(sample_origin[!is.na(sample_origin$Clinical_features) & !is.na(sample_origin$P20230730002_metabolism), ][['P20230730002_metabolism']])
# 
result <- calculate_spearman(metabolites_data, phenotype_data, metabolism_sample, phenotype_sample)
database$metabiolites_phenotype$spearman_pvalue <- result[["pvalue"]]
database$metabiolites_phenotype$spearman_p_adjusted_value <- result[["adjusted_p_values"]]
database$metabiolites_phenotype$spearman_correlation <- result[["correlation"]]

################### Heatmap ####################
sort_heatmap <- function(spearman, pvalue){
  # spearman <- spearman1
  # pvalue <- pvalue2
  df <- spearman
  df <- df[order(df[[1]], decreasing = TRUE), ]
  ordered <- rownames(df)
  t.df <- as.data.frame(t(df))
  df <- t.df[order(t.df[[1]], decreasing = TRUE), ]
  ordered.t <- rownames(df)
  # sort spearman
  spearman <- spearman[, ordered.t]
  spearman.t <- t(spearman)
  spearman <- spearman.t[, ordered]
  # sort pvalue
  pvalue <- pvalue[, ordered.t]
  pvalue.t <- t(pvalue)
  pvalue <- pvalue.t[, ordered]
  result <- list()
  result[['spearman']] <- spearman
  result[['pvalue']] <- pvalue
  return(result)
}
# heatmap  and fill string '*'  if p<0.05
df_metabolism_phenotype_spearman <- database$metabiolites_phenotype$spearman_correlation
df_metabolism_phenotype_p_value <- database$metabiolites_phenotype$spearman_pvalue
df_metabolism_phenotype_fdr_p_value <- database$metabiolites_phenotype$spearman_p_adjusted_value
df_metabolism_phenotype_spearman <- t(as.data.frame(df_metabolism_phenotype_spearman))
df_metabolism_phenotype_p_value <- t(as.data.frame(df_metabolism_phenotype_p_value))
df_metabolism_phenotype_fdr_p_value <- t(as.data.frame(df_metabolism_phenotype_fdr_p_value))
# Check
rownames(df_metabolism_phenotype_spearman) == rownames(df_metabolism_phenotype_p_value)
colnames(df_metabolism_phenotype_spearman) == colnames(df_metabolism_phenotype_p_value)
# Rename colname
# deluplicated 
df_metabolism_phenotype_spearman <- df_metabolism_phenotype_spearman[!duplicated(rownames(df_metabolism_phenotype_spearman)), ]
df_metabolism_phenotype_p_value <- df_metabolism_phenotype_p_value[!duplicated(rownames(df_metabolism_phenotype_p_value)), ]
df_metabolism_phenotype_fdr_p_value <- df_metabolism_phenotype_fdr_p_value[!duplicated(rownames(df_metabolism_phenotype_fdr_p_value)), ]
# Check
rownames(df_metabolism_phenotype_spearman) == rownames(df_metabolism_phenotype_p_value)
library('stringr')
rownames(df_metabolism_phenotype_spearman) <- str_sub(rownames(df_metabolism_phenotype_spearman), start = 1, end = 40)
# plot
# sum(display_text == 0)
display_text <- df_metabolism_phenotype_p_value
df_metabolism_phenotype_spearman_abs <- abs(df_metabolism_phenotype_spearman)
Filter_matrix <- df_metabolism_phenotype_spearman_abs > 0.2 & df_metabolism_phenotype_p_value < 0.05
# Modified FALSE to NA
# Replace FALSE with NA
Filter_matrix <- as.data.frame(Filter_matrix)
Filter_matrix <- data.frame(lapply(Filter_matrix, function(x) ifelse(x == FALSE, NA, x)))
display_text <- display_text * Filter_matrix

# Apply the conditions using ifelse
display_text <- ifelse(display_text < 0.001 , "***", 
                       ifelse(display_text < 0.01, "**", 
                              ifelse(display_text < 0.05, "*", '')))
display_text[is.na(display_text)] <- ""

Heatmap_database$metabolism_phenotype_spearman <- df_metabolism_phenotype_spearman
Heatmap_database$metabolism_phenotype_p_value <- df_metabolism_phenotype_p_value
Heatmap_database$metabolism_phenotype_fdr_p_value <- df_metabolism_phenotype_fdr_p_value
library(pheatmap)
p <-  pheatmap(df_metabolism_phenotype_spearman,
               cellheight = 20, 
               cellwidth = 20, 
               border_color= "grey",  
               show_colnames = TRUE, 
               show_rownames=TRUE,  
               fontsize=9, 
               color = pheamap_gradient_color2,
               annotation_legend=TRUE, 
               cluster_rows = FALSE,
               cluster_cols = TRUE, 
               number_color = "white",
               display_numbers = display_text, 
               fontsize_number = 15, 
               angle_col = "45"
)

p
png('metabolism_phenotype_heatmap.png', width = 2800, height = 1600, res=300)  # Adjust width, height, and resolution as needed
print(p)
dev.off()
tiff('metabolism_phenotype_heatmap.tif', width = 10, height = 5, units = "in", res = 300)
print(p)
dev.off()

# Figure 5.1 Spearman: Gut bacterial species and clinical indicators(Heatmap)
############################
######## Heatmap  ##########
############################
###
# select_species
filter_select_species <- unique(c(paste0("s__", merge_relations$Microbial_species), lefse_diff_select$Species))
Taxa_filter_select <- Taxa_diff[Taxa_diff$Species %in% filter_select_species, ]
microbiome_data <- Taxa_filter_select
rownames(microbiome_data) <- Taxa_filter_select$Species
phenotype_data <- combine_pheno
colnames(microbiome_data)
colnames(phenotype_data)
colnames(sample_origin)
microbiome_sample <- na.omit(sample_origin[!is.na(sample_origin$P20230730001_metagenomic), ][['P20230730001_metagenomic']])
phenotype_sample <- na.omit(sample_origin[!is.na(sample_origin$P20230730001_metagenomic), ][['Clinical_features']])

result <- calculate_spearman(microbiome_data, phenotype_data, microbiome_sample, phenotype_sample)
database$microbiome_phenotype$spearman_pvalue <- result[["pvalue"]]
database$microbiome_phenotype$spearman_p_adjusted_value <- result[["adjusted_p_values"]]
database$microbiome_phenotype$spearman_correlation <- result[["correlation"]]

# heatmap  and fill string '*'  if p<0.05
df_metabolism_phenotype_spearman <- database$microbiome_phenotype$spearman_correlation
df_metabolism_phenotype_p_value <- database$microbiome_phenotype$spearman_pvalue
df_metabolism_phenotype_fdr_p_value <- database$microbiome_phenotype$spearman_p_adjusted_value
df_metabolism_phenotype_spearman <- t(as.data.frame(df_metabolism_phenotype_spearman))
df_metabolism_phenotype_p_value <- t(as.data.frame(df_metabolism_phenotype_p_value))
df_metabolism_phenotype_fdr_p_value <- t(as.data.frame(df_metabolism_phenotype_fdr_p_value))
# Check
rownames(df_metabolism_phenotype_spearman) == rownames(df_metabolism_phenotype_p_value)
colnames(df_metabolism_phenotype_spearman) == colnames(df_metabolism_phenotype_p_value)
# Rename colname
# deluplicated 
df_metabolism_phenotype_spearman <- df_metabolism_phenotype_spearman[!duplicated(rownames(df_metabolism_phenotype_spearman)), ]
df_metabolism_phenotype_p_value <- df_metabolism_phenotype_p_value[!duplicated(rownames(df_metabolism_phenotype_p_value)), ]
df_metabolism_phenotype_fdr_p_value <- df_metabolism_phenotype_fdr_p_value[!duplicated(rownames(df_metabolism_phenotype_fdr_p_value)), ]

# Check
rownames(df_metabolism_phenotype_spearman) == rownames(df_metabolism_phenotype_p_value)
library('stringr')
rownames(df_metabolism_phenotype_spearman) <- str_sub(rownames(df_metabolism_phenotype_spearman), start = 1, end = 40)
# plot
# sum(display_text == 0)
display_text <- df_metabolism_phenotype_p_value
df_metabolism_phenotype_spearman_abs <- abs(df_metabolism_phenotype_spearman)
Filter_matrix <- df_metabolism_phenotype_spearman_abs > 0.3 & df_metabolism_phenotype_p_value < 0.05
# Modified FALSE to NA
# Replace FALSE with NA
Filter_matrix <- as.data.frame(Filter_matrix)
Filter_matrix <- data.frame(lapply(Filter_matrix, function(x) ifelse(x == FALSE, NA, x)))
display_text <- display_text * Filter_matrix

# Apply the conditions using ifelse
display_text <- ifelse(display_text < 0.001 , "***", 
                       ifelse(display_text < 0.01, "**", 
                              ifelse(display_text < 0.05, "*", '')))
display_text[is.na(display_text)] <- ""

# 
Heatmap_database$microbiome_phenotype_spearman <- df_metabolism_phenotype_spearman
Heatmap_database$microbiome_phenotype_p_value <- df_metabolism_phenotype_p_value
Heatmap_database$microbiome_phenotype_fdr_p_value <- df_metabolism_phenotype_fdr_p_value
library(pheatmap)
p <-  pheatmap(df_metabolism_phenotype_spearman,
               cellheight = 20, 
               cellwidth = 20, 
               border_color= "grey", 
               show_colnames = TRUE, 
               show_rownames=TRUE,  
               fontsize=9,
               color = pheamap_gradient_color2,
               annotation_legend=TRUE, 
               cluster_rows = FALSE,
               cluster_cols = TRUE,
               number_color = "white",
               display_numbers = display_text, 
               fontsize_number = 15, 
               angle_col = "45"
)

p
png('microbiome_phenotype_heatmap.png', width = 4000, height = 2000, res=300)  # Adjust width, height, and resolution as needed
print(p)
dev.off()
tiff('microbiome_phenotype_heatmap.tif', width = 13, height = 6, units = "in", res = 300)
print(p)
dev.off()

##################################
######## Micro_Metabolites #######
##################################
# Note: Here Need filter species by yourself
metabolites_data <- annot_kegg_diff[!duplicated(annot_kegg_diff$KEGG_Metabolite), ]
rownames(metabolites_data) <- metabolites_data$KEGG_Metabolite
Taxa_filter_select <- Taxa_diff[Taxa_diff$Species %in% filter_select_species, ]
microbiome_data <- Taxa_filter_select
rownames(microbiome_data) <- Taxa_filter_select$Species

colnames(microbiome_data)
colnames(metabolites_data)
colnames(sample_origin)
microbiome_sample <- na.omit(sample_origin[!is.na(sample_origin$P20230730001_metagenomic), ][['P20230730001_metagenomic']])
metabolism_sample <- na.omit(sample_origin[!is.na(sample_origin$P20230730001_metagenomic), ][['P20230730002_metabolism']])

result <- calculate_spearman(microbiome_data, metabolites_data, microbiome_sample, metabolism_sample)
database$microbiome_metabolites$spearman_pvalue <- result[["pvalue"]]
database$microbiome_metabolites$spearman_p_adjusted_value <- result[["adjusted_p_values"]]
database$microbiome_metabolites$spearman_correlation <- result[["correlation"]]

# heatmap  and fill string '*'  if p<0.05
df_microbiome_metabolites_spearman <- database$microbiome_metabolites$spearman_correlation
df_microbiome_metabolites_p_value <- database$microbiome_metabolites$spearman_pvalue
df_microbiome_metabolites_fdr_p_value  <- database$microbiome_metabolites$spearman_p_adjusted_value
df_microbiome_metabolites_spearman <- t(as.data.frame(df_microbiome_metabolites_spearman))
df_microbiome_metabolites_p_value <- t(as.data.frame(df_microbiome_metabolites_p_value))
df_microbiome_metabolites_fdr_p_value <- t(as.data.frame(df_microbiome_metabolites_fdr_p_value))
# Check
rownames(df_microbiome_metabolites_spearman) == rownames(df_microbiome_metabolites_p_value)
colnames(df_microbiome_metabolites_spearman) == colnames(df_microbiome_metabolites_p_value)
# Rename colname
# deluplicated 
df_microbiome_metabolites_spearman <- df_microbiome_metabolites_spearman[!duplicated(rownames(df_microbiome_metabolites_spearman)), ]
df_microbiome_metabolites_p_value <- df_microbiome_metabolites_p_value[!duplicated(rownames(df_microbiome_metabolites_p_value)), ]
df_microbiome_metabolites_fdr_p_value <- df_microbiome_metabolites_fdr_p_value[!duplicated(rownames(df_microbiome_metabolites_fdr_p_value)), ]
# Check
rownames(df_microbiome_metabolites_spearman) == rownames(df_microbiome_metabolites_p_value)
library('stringr')
rownames(df_microbiome_metabolites_spearman) <- str_sub(rownames(df_microbiome_metabolites_spearman), start = 1, end = 40)
# plot
# Check
# sum(display_text == 0)
display_text <- df_microbiome_metabolites_p_value
df_microbiome_metabolites_spearman_abs <- abs(df_microbiome_metabolites_spearman)
Filter_matrix <- df_microbiome_metabolites_spearman_abs > 0.3 & df_microbiome_metabolites_p_value < 0.05
# Modified FALSE to NA
# Replace FALSE with NA
Filter_matrix <- as.data.frame(Filter_matrix)
Filter_matrix <- data.frame(lapply(Filter_matrix, function(x) ifelse(x == FALSE, NA, x)))
display_text <- display_text * Filter_matrix

# Apply the conditions using ifelse
display_text <- ifelse(display_text < 0.001 , "***", 
                       ifelse(display_text < 0.01, "**", 
                              ifelse(display_text < 0.05, "*", '')))
display_text[is.na(display_text)] <- ""

# 
Heatmap_database$microbiome_metabolites_spearman <- df_microbiome_metabolites_spearman
Heatmap_database$microbiome_metabolites_p_value <- df_microbiome_metabolites_p_value
Heatmap_database$microbiome_metabolites_fdr_p_value <- df_microbiome_metabolites_fdr_p_value
library(pheatmap)
p <-  pheatmap(df_microbiome_metabolites_spearman,
               cellheight = 20, 
               cellwidth = 20, 
               border_color= "grey",  
               show_colnames = TRUE, 
               show_rownames=TRUE,  
               fontsize=9, 
               color = pheamap_gradient_color2,
               annotation_legend=TRUE,
               cluster_rows = FALSE,
               cluster_cols = TRUE, 
               number_color = "white",
               display_numbers = display_text, 
               fontsize_number = 15, 
               angle_col = "45"
)

p
png('microbiome_metabolites_heatmap.png', width = 4200, height = 2400, res=300)  # Adjust width, height, and resolution as needed
print(p)
dev.off()
tiff('mmicrobiome_metabolites_heatmap.tif', width = 14, height = 8, units = "in", res = 300)
print(p)
dev.off()


# Microbiome - Phenotype
microbiome_phenotype_spearman <- Heatmap_database$microbiome_phenotype_spearman
microbiome_phenotype_p_value <- Heatmap_database$microbiome_phenotype_p_value
microbiome_phenotype_fdr_p_value <- Heatmap_database$microbiome_phenotype_fdr_p_value
output_micro_pheno_list <- data.frame(
  source = rep(rownames(microbiome_phenotype_spearman), each = ncol(microbiome_phenotype_spearman)),
  target = rep(colnames(microbiome_phenotype_spearman), times = nrow(microbiome_phenotype_spearman)),
  spearman = as.vector(t(microbiome_phenotype_spearman)),
  p_value = as.vector(t(microbiome_phenotype_p_value)),
  fdr_p_value = as.vector(t(microbiome_phenotype_fdr_p_value))
)

output_micro_pheno_list[['Relation']] <- ifelse(as.numeric(output_micro_pheno_list$spearman) > 0, "Postive", "Negative")
output_micro_pheno_list <- output_micro_pheno_list[abs(output_micro_pheno_list$p_value) < 0.05 & abs(output_micro_pheno_list$spearman) > 0.3, ]
# Order
output_micro_pheno_list <- output_micro_pheno_list[order(abs(output_micro_pheno_list$p_value)), ]
output_micro_pheno_list <- output_micro_pheno_list[order(output_micro_pheno_list$Relation), ]

Table_S7 <- data.frame(
  Species = output_micro_pheno_list$target,
  `Clinical indicators` = output_micro_pheno_list$source,
  `R value` = round(output_micro_pheno_list$spearman, 4),
  `P value` = round(output_micro_pheno_list$p_value, 4),
  `FDR P Value`  = round(output_micro_pheno_list$fdr_p_value, 4),
  Relation = output_micro_pheno_list$Relation
)

column_name1 <- "Table S7. Spearman correlation between species and clinical indicators"
column_name2 <- colnames(Table_S7)
Table_S7 <- rbind(c(column_name1, rep("", 5)), column_name2, Table_S7)

# Metabolites - Phenotype
metabolism_phenotype_spearman <- Heatmap_database$metabolism_phenotype_spearman
metabolism_phenotype_p_value <- Heatmap_database$metabolism_phenotype_p_value
metabolism_phenotype_fdr_p_value <- Heatmap_database$metabolism_phenotype_fdr_p_value
output_meta_pheno_list <- data.frame(
  source = rep(rownames(metabolism_phenotype_spearman), each = ncol(metabolism_phenotype_spearman)),
  target = rep(colnames(metabolism_phenotype_spearman), times = nrow(metabolism_phenotype_spearman)),
  spearman = as.vector(t(metabolism_phenotype_spearman)),
  p_value = as.vector(t(metabolism_phenotype_p_value)),
  fdr_p_value = as.vector(t(metabolism_phenotype_fdr_p_value))
)

output_meta_pheno_list[['Relation']] <- ifelse(as.numeric(output_meta_pheno_list$spearman) > 0, "Postive", "Negative")
output_meta_pheno_list <- output_meta_pheno_list[(abs(output_meta_pheno_list$p_value) < 0.05 & abs(output_meta_pheno_list$spearman) > 0.2), ]
# Order
output_meta_pheno_list <- output_meta_pheno_list[order(abs(output_meta_pheno_list$p_value)), ]
output_meta_pheno_list <- output_meta_pheno_list[order(output_meta_pheno_list$Relation), ]

Table_S8 <- data.frame(
  Metabolite = output_meta_pheno_list$target,
  `Clinical indicators` = output_meta_pheno_list$source,
  `R value` = round(output_meta_pheno_list$spearman, 4),
  `P value` = round(output_meta_pheno_list$p_value, 4),
  `FDR P Value`  = round(output_meta_pheno_list$fdr_p_value, 4),
  Relation = output_meta_pheno_list$Relation
)
column_name1 <- "Table S8. Spearman correlation between metabolites and clinical indicators"
column_name2 <- colnames(Table_S8)
Table_S8 <- rbind(c(column_name1, rep("", 5)), column_name2, Table_S8)

# Metabolites - Phenotype
microbiome_metabolites_spearman <- Heatmap_database$microbiome_metabolites_spearman
microbiome_metabolites_p_value <- Heatmap_database$microbiome_metabolites_p_value
microbiome_metabolites_fdr_p_value <- Heatmap_database$microbiome_metabolites_fdr_p_value
output_meta_pheno_list <- data.frame(
  source = rep(rownames(microbiome_metabolites_spearman), each = ncol(microbiome_metabolites_spearman)),
  target = rep(colnames(microbiome_metabolites_spearman), times = nrow(microbiome_metabolites_spearman)),
  spearman = as.vector(t(microbiome_metabolites_spearman)),
  p_value = as.vector(t(microbiome_metabolites_p_value)),
  fdr_p_value = as.vector(t(microbiome_metabolites_fdr_p_value))
  
)

output_meta_pheno_list[['Relation']] <- ifelse(as.numeric(output_meta_pheno_list$spearman) > 0, "Postive", "Negative")
output_meta_pheno_list <- output_meta_pheno_list[(abs(output_meta_pheno_list$p_value) < 0.05 & abs(output_meta_pheno_list$spearman) > 0.3), ]
# Order
output_meta_pheno_list <- output_meta_pheno_list[order(abs(output_meta_pheno_list$p_value)), ]
output_meta_pheno_list <- output_meta_pheno_list[order(output_meta_pheno_list$Relation), ]

Table_S9 <- data.frame(
  Metabolite = output_meta_pheno_list$target,
  `Clinical indicators` = output_meta_pheno_list$source,
  `R value` = round(output_meta_pheno_list$spearman, 4),
  `P value` = round(output_meta_pheno_list$p_value, 4),
  `FDR P value` = round(output_meta_pheno_list$fdr_p_value, 4),
  Relation = output_meta_pheno_list$Relation
)
column_name1 <- "Table S9. Spearman correlation between species and metabolites"
column_name2 <- colnames(Table_S9)
Table_S9 <- rbind(c(column_name1, rep("", 5)), column_name2, Table_S9)



writeData(wb, "Table S7", Table_S7, startCol = 1, startRow = 1, rowNames = FALSE,  colNames = FALSE)
writeData(wb, "Table S8", Table_S8, startCol = 1, startRow = 1, rowNames = FALSE,  colNames = FALSE)
writeData(wb, "Table S9", Table_S9, startCol = 1, startRow = 1, rowNames = FALSE,  colNames = FALSE)

#  Metabolite Formula m/z RT(min) MSI Level Annotation Method HMDB ID KEGG ID
metabolomic_df_filter <- metabolomic_df[metabolomic_df$ID %in% annot_kegg_diff$`Primary ID`, ]
metabolomic_df_filter$MS2kegg %in% annot_kegg_diff$MS2kegg
metabolomic_df_filter_merge <- merge(metabolomic_df_filter, annot_kegg_diff, by.x = 'ID', by.y = 'Primary ID')
Table_S10 <- data.frame(
  Metabolite = metabolomic_df_filter_merge$KEGG_Metabolite,
  Formula = metabolomic_df_filter_merge$MS2MetaboliteFormula,
  MZ = metabolomic_df_filter_merge$MZ,
  RT = metabolomic_df_filter_merge$RT,
  KEGG_ID = metabolomic_df_filter_merge$MS2kegg.x, 
  MS2_Score = metabolomic_df_filter_merge$MS2MetaboliteScore, 
  Annotation_Method = 'MS²'
)
column_name1 <- "Table S10. Metabolite annotations and supporting evidence (retention time and MS²)."
column_name2 <- colnames(Table_S10)
Table_S10 <- rbind(c(column_name1, rep("", 6)), column_name2, Table_S10)

writeData(wb, "Table S10", Table_S10, startCol = 1, startRow = 1, rowNames = FALSE,  colNames = FALSE)
saveWorkbook(wb, "Supplement.xlsx", overwrite = TRUE)












