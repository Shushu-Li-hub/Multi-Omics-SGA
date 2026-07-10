# save dir
setwd(save_dir)

# Color
FC_up_colour <- "#FF82AB" 
FC_down_colour <- "#00BFFF"
FC_up_colour <- "red" 
FC_down_colour <- "blue"

# Figure 3.2 Volcano
############################
######## Volcano  ##########
############################
log2fc <- 0.585
test_method <- 'wilcox.test.p_value'

utils_direcotry <- 'H:/LSS/project/utils'
source(paste0(utils_direcotry, '/new_calc_metabolism_diff.R'))

# Get Kos Data
colnames(metabolomic_abun)
sample_dir <- 'H:/LSS/data'
sample_origin <- read_excel(paste0(sample_dir, '/sample_origin.xlsx'))
metabolism_sample <- na.omit(sample_origin[['P20230730002_metabolism']])
AGA_samples <- metabolism_sample[grepl('AGA', metabolism_sample)]
SGA_samples <- metabolism_sample[grepl('SGA', metabolism_sample)]
# Filtered data 
AGA_SGA_metabolomi_Data <- metabolomic_abun[, c(AGA_samples, SGA_samples)]
rownames(AGA_SGA_metabolomi_Data) <- metabolomic_abun[["ID"]]
AGA_SGA_group <- rbind(data.frame(sample = AGA_samples, group = rep('AGA', length(AGA_samples))), 
                       data.frame(sample = SGA_samples, group = rep('SGA', length(SGA_samples))))
AGA_SGA_group$group <- factor(AGA_SGA_group$group, levels = c('AGA', 'SGA'))
Metabolomi_diff <- new_calc_metabolism_diff(AGA_SGA_metabolomi_Data, AGA_SGA_group)
Metabolomi_diff[['ID']] <- rownames(AGA_SGA_metabolomi_Data)
# relative abundance > 0.00001
# & (Unigenes_diff$ratio_AGA > 30) & (Unigenes_diff$ratio_SGA > 30)
Metabolomi_diff <- Metabolomi_diff[(Metabolomi_diff$mean_AGA > 0.00001) & (Metabolomi_diff$mean_SGA > 0.00001) , ]
nrow(Metabolomi_diff)

# merge vip value of simca
simca_directory <- 'H:/LSS/simca_process'
postive_vip <- read_xlsx(paste0(simca_directory, '/simca_postive_result/postive_vip.xlsx'))
negative_vip <- read_xlsx(paste0(simca_directory, '/simca_negative_result/negative_vip.xlsx'))
rbind_vip <- rbind(postive_vip, negative_vip)
colnames(Metabolomi_diff)
colnames(rbind_vip)
merge_diff <- merge(Metabolomi_diff, rbind_vip, by.x ='ID', by.y = 'Primary ID')

########## plot Volcano ##############
########## Postive    ################
origin_df <- merge_diff[grepl('pos', merge_diff$ID), ]
# Function: get_color
get_color <- function(value, log_fold_change, vip){
  if(value <= 0.05 & log_fold_change >= log2fc & vip >= 1){
    return('Sig_Up')
  }
  if(value <= 0.05 & log_fold_change <= -log2fc & vip >= 1){
    return('Sig_Down')
  }
  if(value > 0.05 | abs(log_fold_change) < log2fc | vip < 1){
    return('No_Diff')
  }
}

group <- unlist(lapply(1:nrow(origin_df), function(x) get_color(origin_df[[test_method]][x], 
                                                                origin_df$log_fold_change[x], 
                                                                origin_df$M1.VIPpred[x])))
point_df <- data.frame(id = origin_df$id, 
                       x.axis = origin_df$log_fold_change, 
                       y.axis = -log10(origin_df[[test_method]]), 
                       groups = group)

Sig_Down = Sig_down_color  # down
Sig_Up = Sig_up_color # UP
No_Diff = No_diff_color
Sig_Down_count <- origin_df[origin_df[[test_method]] <= 0.05 & origin_df$log_fold_change <= -log2fc & origin_df$M1.VIPpred >= 1, ]
Sig_Up_count <- origin_df[origin_df[[test_method]] <= 0.05 & origin_df$log_fold_change >= log2fc & origin_df$M1.VIPpred >= 1, ]
No_Diff_count <- origin_df[origin_df[[test_method]] > 0.05 | abs(origin_df$log_fold_change) < log2fc | origin_df$M1.VIPpred < 1, ]
nrow(Sig_Down_count) + nrow(Sig_Up_count) + nrow(No_Diff_count)
nrow(origin_df)
point_df$groups <- factor(point_df$groups, levels = c("Sig_Down","No_Diff","Sig_Up"))
# X . log2(Fold change)
# Y . -log10(P value)
library(ggplot2)
# Create the volcano plot
volcano_plot <- ggplot(point_df, aes(x = x.axis, y = y.axis, colour = groups)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Sig_Down" = Sig_Down, 'No_Diff' = No_Diff, "Sig_Up" = Sig_Up), 
                     labels = c(paste0("Sig_Down(", nrow(Sig_Down_count), ")"), 
                                paste0("No_Diff(", nrow(No_Diff_count), ")"), 
                                paste0("Sig_Up(", nrow(Sig_Up_count), ")"))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#7F7F7FFF", size = 1) +
  geom_vline(xintercept = -log2fc, linetype = "dashed", color = "#7F7F7FFF", size = 1) +
  geom_vline(xintercept = log2fc, linetype = "dashed", color = "#7F7F7FFF", size = 1) +
  theme_minimal() +
  labs(
    title = "",
    x = "Log2 (Fold Change)",
    y = "-log10 (P value)"
  )  +
  # scale_x_continuous(limits = c(-6, 6))+
  theme(
    legend.title = element_blank(), 
    legend.key = element_rect(fill = "white", color = "black", size = 0.5),
    legend.position = c(0.998, 0.998), legend.justification = c(1, 1),  # Adjust legend position
    legend.background = element_rect(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks = element_line(color = "black"), 
    axis.text = element_text(face = 'bold', size = 12, color = "black"),
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# Print the plot
print(volcano_plot)
png('Postive_volcano_plot.png', width = 2000, height = 2000, res=300)  # Adjust width, height, and resolution as needed
print(volcano_plot)
dev.off()
tiff('Postive_volcano_plot.tif', width = 6, height = 6, units = "in", res = 300)
print(volcano_plot)
dev.off()

########## Negative   ################
origin_df <- merge_diff[grepl('neg', merge_diff$ID), ]
# Function: get_color
get_color <- function(value, log_fold_change, vip){
  if(value <= 0.05 & log_fold_change >= log2fc & vip >= 1){
    return('Sig_Up')
  }
  if(value <= 0.05 & log_fold_change <= -log2fc & vip >= 1){
    return('Sig_Down')
  }
  if(value > 0.05 | abs(log_fold_change) < log2fc | vip < 1){
    return('No_Diff')
  }
}

group <- unlist(lapply(1:nrow(origin_df), function(x) get_color(origin_df[[test_method]][x], 
                                                                origin_df$log_fold_change[x], 
                                                                origin_df$M1.VIPpred[x])))
point_df <- data.frame(id = origin_df$id, 
                       x.axis = origin_df$log_fold_change, 
                       y.axis = -log10(origin_df[[test_method]]), 
                       groups = group)
point_df$groups <- factor(point_df$groups, levels = c("Sig_Down","No_Diff","Sig_Up"))

Sig_Down = Sig_down_color  # down
Sig_Up = Sig_up_color # UP
No_Diff = No_diff_color
Sig_Down_count <- origin_df[origin_df[[test_method]] <= 0.05 & origin_df$log_fold_change <= -log2fc & origin_df$M1.VIPpred >= 1, ]
Sig_Up_count <- origin_df[origin_df[[test_method]] <= 0.05 & origin_df$log_fold_change >= log2fc & origin_df$M1.VIPpred >= 1, ]
No_Diff_count <- origin_df[origin_df[[test_method]] > 0.05 | abs(origin_df$log_fold_change) < log2fc | origin_df$M1.VIPpred < 1, ]
nrow(Sig_Down_count) + nrow(Sig_Up_count) + nrow(No_Diff_count)
nrow(origin_df)

# X . log2(Fold change)
# Y . -log10(P value)
library(ggplot2)
# Create the volcano plot
volcano_plot <- ggplot(point_df, aes(x = x.axis, y = y.axis, colour = groups)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Sig_Down" = Sig_Down, 'No_Diff' = No_Diff, "Sig_Up" = Sig_Up), 
                     labels = c(paste0("Sig_Down(", nrow(Sig_Down_count), ")"), 
                                paste0("No_Diff(", nrow(No_Diff_count), ")"), 
                                paste0("Sig_Up(", nrow(Sig_Up_count), ")"))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#7F7F7FFF", size = 1) +
  geom_vline(xintercept = -log2fc, linetype = "dashed", color = "#7F7F7FFF", size = 1) +
  geom_vline(xintercept = log2fc, linetype = "dashed", color = "#7F7F7FFF", size = 1) +
  theme_minimal() +
  labs(
    title = "",
    x = "Log2 (Fold Change)",
    y = "-log10 (P value)"
  )  +
  scale_x_continuous(limits = c(-6, 6))+
  theme(
    legend.title = element_blank(), 
    legend.key = element_rect(fill = "white", color = "black", size = 0.5),
    legend.position = c(0.998, 0.998), legend.justification = c(1, 1),  # Adjust legend position
    legend.background = element_rect(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks = element_line(color = "black"), 
    axis.text = element_text(face = 'bold', size = 12, color = "black"),
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# Print the plot
print(volcano_plot)
png('Negative_volcano_plot.png',  width = 2000, height = 2000, res=300)  # Adjust width, height, and resolution as needed
print(volcano_plot)
dev.off()
tiff('Negative_volcano_plot.tif', width = 6, height = 6, units = "in", res = 300)
print(volcano_plot)
dev.off()


# Figure 3.3 Metabolites annotated as known metabolites
############################
######## Histogram  ########
############################
# Note: Data sum by KEGGID
colnames(metabolomic_abun)
# Filtered data
KEGG_abun <- metabolomic_abun[!is.na(metabolomic_abun$MS2kegg), ]
KEGG_abun <- KEGG_abun[!KEGG_abun$MS2kegg %in% c("null", "-", "NA", "[M-H]-"), ]
KEGG_abun <- KEGG_abun[, c('MS2kegg', AGA_samples, SGA_samples)]
KEGG_abun_summarized <- aggregate(. ~ MS2kegg, data = KEGG_abun, FUN = sum)
colnames(Unigenes_kos_summarized)
AGA_SGA_metabolomi_Data <- KEGG_abun_summarized[, c(AGA_samples, SGA_samples)]
rownames(AGA_SGA_metabolomi_Data) <- KEGG_abun_summarized[["MS2kegg"]]
AGA_SGA_group <- rbind(data.frame(sample = AGA_samples, group = rep('AGA', length(AGA_samples))), 
                       data.frame(sample = SGA_samples, group = rep('SGA', length(SGA_samples))))
AGA_SGA_group$group <- factor(AGA_SGA_group$group, levels = c('AGA', 'SGA'))
KEGG_diff <- new_calc_metabolism_diff(AGA_SGA_metabolomi_Data, AGA_SGA_group)

colnames(KEGG_diff)
rbind_vip <- rbind(postive_vip, negative_vip)
merge_rbind_vip <- merge(rbind_vip, metabolomic_abun[, c("ID", "MS2kegg")], by.x = "Primary ID", by.y = 'ID')
merge_rbind_vip <- merge_rbind_vip[order(merge_rbind_vip$M1.VIPpred, decreasing = TRUE), ]
merge_rbind_vip <- merge_rbind_vip[!duplicated(merge_rbind_vip$MS2kegg), ]
KEGG_diff <- merge(KEGG_diff, merge_rbind_vip, by.x = 'id', by.y = 'MS2kegg')
annot_diff <- KEGG_diff[KEGG_diff[[test_method]] < 0.05 & abs(KEGG_diff$log_fold_change) >= log2fc & KEGG_diff$M1.VIPpred >= 1, ]
length(unique(KEGG_diff$id))
kegg_directory <- 'H:/LSS/kegg'
library(readr)
cpd2name_data <- read_delim(paste0(kegg_directory, '/kegg_cpd2name.txt'), delim = "\t", col_names = FALSE, col_types = cols(), trim_ws = TRUE)
library(stringr)
utils_direcotry <- 'H:/LSS/project/utils'
source(paste0(utils_direcotry, '/define_function.R'))
select_name <- function(x){
  nchar_length <- nchar(strsplit(x, ';')[[1]])
  select_names <- strsplit(x, ';')[[1]][nchar_length <= 30]
  return(select_names[1])
  if (length(select_names) == 0){
    return(strsplit(x, ';')[[1]][1])
  }
  
}
cpd2name_data[['name']] <- lapply(cpd2name_data$X2, function(x) select_name(x))
annot_diff[['KEGG_Metabolite']] <- replace_names(annot_diff$id, cpd2name_data$X1, cpd2name_data$name)
annot_diff[['MS2Metabolite']] <- replace_names(annot_diff$`Primary ID`, metabolomic_abun$ID, metabolomic_abun$MS2Metabolite)

annot_diff <- annot_diff[order(annot_diff$fold_change, decreasing = TRUE), ]
annot_kegg_diff <- annot_diff
annot_kegg_diff[['MS2kegg']] <- annot_kegg_diff[['id']]

########## Plot bar ############
library(ggplot2)
annot_diff
Histogram_data <- data.frame(
  name =  annot_diff$KEGG_Metabolite,
  value = annot_diff$log_fold_change,  # Random normal values
  group = ifelse(annot_diff$log_fold_change > 0, 'UP', 'Down')
)
Histogram_data <- Histogram_data[!duplicated(Histogram_data$name), ]
Histogram_data$name <- factor(Histogram_data$name, levels = Histogram_data$name)


# View the first few rows of the dataset
head(Histogram_data)
p <- ggplot(Histogram_data, aes(x = name, y = value, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5, color = "black") +
  labs(
    title = " ",
    x = " ",
    y = "SGA vs AGA log2(FC)"
  ) +
  scale_fill_manual(values = c("UP" = FC_up_colour, "Down" = FC_down_colour)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.y = element_text(face = "bold"), 
        axis.ticks = element_line(size = 0.5), 
        panel.background = element_rect(fill = NA),
        legend.position = "none"
  )
print(p)
png('Histogram_plot.png', width = 2400, height = 1600, res=300)  # Adjust width, height, and resolution as needed
print(p)
dev.off()
tiff('Histogram_plot.tif', width = 9, height = 6, units = "in", res = 300)
print(p)
dev.off()

# Figure 2.4 Enrichment
############################
#### Bubble diagram  #######
############################
metabolites_set <- annot_kegg_diff
metabolites_set <- unique(metabolites_set$MS2kegg)
length(metabolites_set)

# N pathway Metabolism total
result_combine_directroy <- 'H:/LSS/kegg/hsa_pathway_kgml'
kegg.levels.directory <- 'H:/LSS/kegg/hsa_pathway'
library(readxl)
result_combine <- read_xlsx(paste0(result_combine_directroy, '/hsa_pathway_nodes.xlsx'))
kegg.pathway.levels <- read_xlsx(paste0(kegg.levels.directory, '/hsa.kegg.pathway.levels.xlsx'))
######## Note: level1: 'Metabolism'  level2: not include "Global and overview maps" for reduce false negatives/positive and reduce hyper test false positive(i think 4117-2797 compounds few appear in other pathways)
# Metabolism.pathway <- kegg.pathway.levels[kegg.pathway.levels$Level1 == 'Metabolism' & kegg.pathway.levels$Level2!= "Global and overview maps", ][['hsa.kegg.pathway.database']]
Metabolism.pathway <- kegg.pathway.levels[kegg.pathway.levels$Level1 == 'Metabolism', ][['hsa.kegg.pathway.database']]
######## NOTE :Filter result_combine by levels ########
result_combine <- result_combine[result_combine$name %in% Metabolism.pathway, ]
result_combine <- result_combine[!is.na(result_combine$cpd), ]
all.cpd.hsa <- unique(unlist(strsplit(unlist(result_combine$cpd), ';')))
N <- length(all.cpd.hsa)
# filter cpd in result_combine
metabolites_set_filter <- metabolites_set[metabolites_set %in% all.cpd.hsa]
# n diff number
n <- length(metabolites_set_filter)

phyper_result <- data.frame()
for (variable in result_combine$name) {
  # M specific pathway number 
  # variable <- 'hsa00513'
  specific_pathway_compound <- unlist(strsplit(result_combine[result_combine$name == variable, ][['cpd']], ';'))
  M <- length(specific_pathway_compound)
  # m Overlapping number 
  inter_compound <- intersect(specific_pathway_compound, metabolites_set_filter)
  m <- length(inter_compound)
  p.value <- stats::phyper(q = m - 1, m = M, n = (N - M), 
                           k = n, lower.tail = FALSE)
  # erichment ratio
  prop_reference <- length(metabolites_set_filter) / length(all.cpd.hsa)
  hits <- length(inter_compound)
  expected_hits <- prop_reference * length(specific_pathway_compound)
  new_data <- c(name = variable, 
                compounds.dem.length = m,
                compounds.dem = paste0(inter_compound, collapse = ';'), 
                compounds.Total.length = M, 
                compounds.Total = paste0(specific_pathway_compound, collapse = ';'),
                p.value = p.value, 
                enrichment_ratio = hits/expected_hits)
  if (nrow(phyper_result) == 0){
    phyper_result <- as.data.frame(t(new_data))
  } else {
    phyper_result <- rbind(phyper_result, t(new_data))
  }
}
# merge
enrichment_phyper_result <- merge(result_combine, phyper_result, by = "name", all = TRUE, sort = FALSE)
# Filter hit
enrichment_filter_hit <-  enrichment_phyper_result[enrichment_phyper_result$compounds.dem.length >= 1, ]
# adjust by FDR
enrichment_filter_hit[['p.value.FDR']] <- p.adjust(enrichment_filter_hit[['p.value']], method = "fdr")  # FDR correction
# merge 'Descirption'
merge_df <- merge(enrichment_filter_hit, kegg.pathway.levels, by.x = 'name', by.y = 'hsa.kegg.pathway.database')
colnames(merge_df)
numeirc_vairables <- c('compounds.dem.length', 'compounds.Total.length', 'p.value', 'enrichment_ratio', 'p.value.FDR')
merge_df[numeirc_vairables] <- lapply(merge_df[numeirc_vairables], as.numeric)

####Top 25 Overview#####
# Filter p.value < 0.05
merge_df <- merge_df[order(merge_df$p.value), ]
if (nrow(merge_df) >25){
  merge_df <- merge_df[1:25, ]
} 
# merge_df <- merge_df[merge_df$p.value < 0.05, ]

# Run plot
library(ggplot2)
#########PlotEnrichDotPlot############
# mSetObj=NA, enrichType = "ora", imgName, format="png", dpi=72, width=NA
enrich_df <- merge_df
# Function: GetMyHeatCols
GetMyHeatCols <- function(len){
  if(len > 50){
    ht.col <- c(substr(heat.colors(50), 0, 7), rep("#FFFFFF", len-50));
  }else{
    # reduce to hex by remove the last character so HTML understand
    ht.col <- substr(heat.colors(len), 0, 7);
  }
}
my.cols <- GetMyHeatCols(nrow(enrich_df))

if(nrow(enrich_df) > 25){
  enrich_df <- enrich_df[1:25,]
  my.cols <- my.cols[1:25]
}

df <- data.frame(Name = factor(enrich_df$Level3, levels = rev(enrich_df$Level3)),
                 rawp = enrich_df$p.value,
                 logp = -log10(enrich_df$p.value),
                 folds = enrich_df$enrichment_ratio)
maxp <- max(df$rawp)

p <- ggplot(df,  aes(x = logp, y = Name)) + 
  geom_vline(xintercept = -log10(0.05), color = 'gray', size = 1, linetype = 2) + 
  geom_point(aes(size = folds, color = rawp)) + scale_size_continuous(range = c(2, 8)) +
  theme_bw(base_size = 14.5) +
  scale_colour_gradient(limits=c(0, maxp), low=my.cols[1], high = my.cols[length(my.cols)]) +
  ylab(NULL) + xlab("-log10 (p-value)") + 
  ggtitle("Overview of Enriched Metabolite Sets") +
  theme(plot.title = element_text(face = "bold", hjust = 0.3),
        panel.border = element_rect(color = "black", fill = NA, size = 1.25),
        legend.text=element_text(face = 'bold', size=14),
        legend.title=element_text(face = 'bold', size=15),
        axis.title = element_text(color = 'black', face = 'bold'),
        axis.text.x = element_text(color = 'black', face = 'bold'),
        axis.text.y = element_text(color = 'black', face = 'bold'),
        axis.ticks = element_line(color = 'black', size = 1.25)
  )

p$labels$colour <- "P-value"
p$labels$size <- "Enrichment Ratio"
p
png('metabolites.kegg.enrichemnt.dot.png', width = 3100, height = 2500, res=300)  # Adjust width, height, and resolution as needed
print(p)
dev.off()
tiff('metabolites.kegg.enrichemnt.dot.tif', width = 10, height = 8, units = "in", res = 300)
print(p)
dev.off()

