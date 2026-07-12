# save dir
setwd(save_dir)


# Color
library(pals)
# Color of simca
# c("#3A5FCD", "#CD2626", "#FFFFFF")
# 58,95,205 205,38,38
Sig_down_color =  "#4876FF"# down
Sig_up_color = "#FF4040"
Sig_down_color =  "navy"# down
Sig_up_color = "firebrick3"
No_diff_color = "#7F7F7FFF"
# pheamap_gradient_color = coolwarm(n=100)

# Figure 2.1 KEGG annotation and classifcation of microbial functional genes.
############################
##### KEGG annotation ######
############################
log2fc <- 0.585
test_method <- 'rfit.p_value'

# KEGGLevel2
# Remove specified pathway
library(ggsci)
unigenes_kegg_kos <- unigenes_kegg_df[!is.na(unigenes_kegg_df$KOEntry), ]
kegg_annotation_df <- unigenes_kegg_kos[unigenes_kegg_kos$KEGGLevel2 != 'Global and overview maps', ]
colnames(kegg_annotation_df)
kegg_annotation_df <- kegg_annotation_df[, c('Query', 'KOEntry', 'KEGGLevel1', 'KEGGLevel2')]
kegg_annotation_count <- as.data.frame(table(kegg_annotation_df$KEGGLevel2))
colnames(kegg_annotation_count) <- c('KEGGLevel2', 'Count')
kegg_annotation_count[['Mean_Proportions']] <- kegg_annotation_count[['Count']] / sum(kegg_annotation_count[['Count']])
# KEGGLevel2 to
colnames(kegg_annotation_df)
kegg_level2_level1 <- kegg_annotation_df[, c("KEGGLevel1", "KEGGLevel2")]
kegg_level2_level1 <- kegg_level2_level1[!duplicated(kegg_level2_level1$KEGGLevel2), ]
merge_kegg_annotation <- merge(kegg_annotation_count, kegg_level2_level1, by = 'KEGGLevel2', all = TRUE)
colnames(merge_kegg_annotation)
merge_kegg_annotation <- merge_kegg_annotation[order(merge_kegg_annotation$Mean_Proportions, decreasing = TRUE), ]
unique(merge_kegg_annotation$KEGGLevel1)
specified_order <- c("Metabolism", "Environmental Information Processing", 
                     "Genetic Information Processing", "Cellular Processes", "Human Diseases")
merge_kegg_annotation <- merge_kegg_annotation[order(match(merge_kegg_annotation$KEGGLevel1, specified_order)), ]
merge_kegg_annotation$KEGGLevel1 <- factor(merge_kegg_annotation$KEGGLevel1, levels = specified_order)
merge_kegg_annotation$KEGGLevel2 <- factor(merge_kegg_annotation$KEGGLevel2, levels = merge_kegg_annotation$KEGGLevel2)
fill_colors <- pal_lancet("lanonc")(length(specified_order))

# Plot
library(viridis)
library(ggsci)
p <- ggplot(merge_kegg_annotation, aes(x = KEGGLevel2, y = Mean_Proportions, fill = KEGGLevel1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = fill_colors) + 
  # scale_fill_brewer(palette = "Set1") + 
  # scale_fill_viridis(discrete = TRUE) +
  labs(title = " ",
       x = " ",
       y = expression(bold("Mean Proportions")),
       fill = " ") +
  theme_minimal() +
  theme(axis.text.x = element_text(face = "bold",angle = 90, hjust = 1),
        axis.text.y = element_text(face = "bold", color = "black"), 
        axis.ticks = element_line(size = 0.5), 
        legend.title = element_text(face = "bold"),
        legend.text  = element_text(face = "bold"),
        # panel.background = element_rect(fill = panel_colors),
        panel.spacing = unit(0, "cm"), 
        # strip.text = element_blank(),
        # strip.text = element_text(size = 12, color = "black", face = "bold"),
        # strip.background = element_rect(fill = 'lightgray', color = "black"),
        panel.background = element_rect(fill = NA), 
        legend.position = c(0.8, 0.8)
        # legend.position = "none"
        # panel.grid.major = element_line(color = "gray")
        # strip.text.x = element_text(angle = 90, hjust = 1, color = "black", size = 12, face = "bold")
  )

# Print the plot
print(p)
png('Kegg_annotation_classification.png', width = 2400, height = 1800, res=300)  # Adjust width, height, and resolution as needed
print(p)
dev.off()
tiff('Kegg_annotation_classification.tif', width = 9, height = 6, units = "in", res = 300)
print(p)
dev.off()


# Figure 2.2 Volcano
############################
######## Volcano  ##########
############################
utils_direcotry <- 'H:/LSS/project/utils'
source(paste0(utils_direcotry, '/new_calc_kos_diff.R'))
# Get Kos Data
unigenes_kegg_kos <- unigenes_kegg_df[!is.na(unigenes_kegg_df$KOEntry), ]

# unigenes_kegg_kos <- unigenes_kegg_kos[unigenes_kegg_kos$Identity >= 100, ]
colnames(unigenes_kegg_kos)
unigenes_kegg_info <- unigenes_kegg_kos[, c('Query', 'Identity', 'KOEntry')]
unigenes_kegg_info <- unigenes_kegg_info[!duplicated(unigenes_kegg_info$Query), ]

# Filter the unigenes in kos
Unigenes_kos_abun <- merge(Unigenes_abund, unigenes_kegg_info, by.x = 'Unigene_ID', by.y = 'Query')
# identity>95% 
Unigenes_kos_abun <- Unigenes_kos_abun[Unigenes_kos_abun$Identity >= 95, ]
# Get sum of KOs.
colnames(Unigenes_kos_abun)
Unigenes_kos_abun <- Unigenes_kos_abun[, c('KOEntry', colnames(Unigenes_kos_abun)[grepl('AGA|SGA', colnames(Unigenes_kos_abun))])]
library(dplyr)
# Aggregate counts by KO terms and perform differential expression analysis 
Unigenes_kos_summarized <- aggregate(. ~ KOEntry, data = Unigenes_kos_abun, FUN = sum)
colnames(Unigenes_kos_summarized)

# Select samples
colnames(sample_origin)
sample_dir <- 'H:/LSS/data'
sample_origin <- read_excel(paste0(sample_dir, '/sample_origin.xlsx'))
microbiome_sample <- na.omit(sample_origin[['P20230730001_metagenomic']])
SGA_samples <- microbiome_sample[grepl('SGA', microbiome_sample)]
AGA_samples <- microbiome_sample[grepl('AGA', microbiome_sample)]

# Filtered data
AGA_SGA_Unigenes_Data <- Unigenes_kos_summarized[, c(AGA_samples, SGA_samples)]
rownames(AGA_SGA_Unigenes_Data) <- Unigenes_kos_summarized$KOEntry
AGA_SGA_group <- rbind(data.frame(sample = AGA_samples, group = rep('AGA', length(AGA_samples))), 
                       data.frame(sample = SGA_samples, group = rep('SGA', length(SGA_samples))))
AGA_SGA_group$group <- factor(AGA_SGA_group$group, levels = c('AGA', 'SGA'))
Unigenes_diff <- new_calc_kos_diff(AGA_SGA_Unigenes_Data, AGA_SGA_group, metadata_covariates)
Unigenes_diff[['KOEntry']] <- rownames(AGA_SGA_Unigenes_Data)
# relative abundance > 0.00001
# & (Unigenes_diff$ratio_AGA > 30) & (Unigenes_diff$ratio_SGA > 30)
# Filter: Only species or KOs with an average relative abundance above 10^−7(0.00001%) were considered in the analyses.
Unigenes_filter <- Unigenes_diff[(Unigenes_diff$mean_AGA > 0.00001) & (Unigenes_diff$mean_SGA > 0.00001) , ]
# FDR
Unigenes_filter$rfit.fdr_values <- p.adjust(Unigenes_filter[[test_method]], method = "fdr")  # FDR correction


########## plot Volcano ##############
origin_df <- Unigenes_filter
# Function: get_color
get_color <- function(value, log_fold_change){
  if(value <= 0.05 & log_fold_change >= log2fc){
    return('Sig_Up')
  }
  if(value <= 0.05 & log_fold_change <= -log2fc){
    return('Sig_Down')
  }
  if(value > 0.05 | abs(log_fold_change) < log2fc){
    return('No_Diff')
  }
}

group <- unlist(lapply(1:nrow(origin_df), function(x) get_color(origin_df[[test_method]][x],
                                                                origin_df$log_fold_change[x])))
point_df <- data.frame(id = origin_df$id,
                       x.axis = origin_df$log_fold_change,
                       y.axis = -log10(origin_df[[test_method]]),
                       groups = group)

Sig_Down = Sig_down_color  # down
Sig_Up = Sig_up_color # UP
No_Diff = No_diff_color
Sig_Down_count <- origin_df[origin_df[[test_method]] < 0.05 & origin_df$log_fold_change <= -log2fc, ]
Sig_Up_count <- origin_df[origin_df[[test_method]] < 0.05 & origin_df$log_fold_change >= log2fc, ]
No_Diff_count <- origin_df[origin_df[[test_method]] > 0.05 | abs(origin_df$log_fold_change) < log2fc, ]
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
    y = "-log10 (P.value)"
  )  +
  scale_x_continuous(limits = c(-8, 8))+
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
png('Unigenes_volcano_plot.png', width = 2000, height = 2000, res=300)  # Adjust width, height, and resolution as needed
print(volcano_plot)
dev.off()
tiff('Unigenes_volcano_plot.tif', width = 6, height = 6, units = "in", res = 300)
print(volcano_plot)
dev.off()


# Figure 2.3 Heatmap
############################
######## Heatmap  ##########
############################
Unigenes_Difference <- Unigenes_filter[(Unigenes_filter[[test_method]] < 0.05)&(abs(Unigenes_filter$log_fold_change)>log2fc), ]
nrow(Unigenes_Difference)
Unigenes_Difference_all <- Unigenes_Difference


# The top 50 most abundance
Unigenes_Difference <- Unigenes_Difference[order(Unigenes_Difference$Sum, decreasing = TRUE), ]
if (nrow(Unigenes_Difference) >= 50){
  Unigenes_Difference <- Unigenes_Difference[1:50, ]
} 


library(pheatmap)
df <- Unigenes_Difference
# SORT FC
df_1 <- df[order(df$log_fold_change), ]
df_2 <- df_1[, c(AGA_samples, SGA_samples)]
rownames(df_2) <- df_1$KOEntry
# group
dfSample_1 <- data.frame(Group=c(rep('AGA_samples', length(AGA_samples)),
                                 rep('SGA_samples', length(SGA_samples))),
                         row.names = colnames(df_2))
p <- pheatmap(df_2,
              # annotation_row=dfGene, # （可选）指定行分组文件
              # annotation_col=dfSample_1, # （可选）指定列分组文件
              cellheight = 10, 
              cellwidth = 10, 
              show_colnames = TRUE, # 是否显示列名
              show_rownames=TRUE,  # 是否显示行名
              fontsize=9, # 字体大小
              color = colorRampPalette(c("blue", "white", "red"))(50), # pheamap_gradient_color
              annotation_legend=TRUE, # 是否显示图例
              border_color= "grey",  # 边框颜色 NA表示没有
              # scale normalized
              scale="row",  # 指定归一化的方式。"row"按行归一化，"column"按列
              cluster_rows = TRUE, # 是否对行聚类
              cluster_cols = FALSE
)
p
png('Kos_heatmap.png', width = 1400, height = 2800, res=300)  # Adjust width, height, and resolution as needed
print(p)
dev.off()
tiff('Kos_heatmap.tif', width = 5, height = 10, units = "in", res = 300)
print(p)
dev.off()


# Figure 2.4 Enrichment
############################
#### Bubble diagram  #######
############################
Unigenes_Difference <- Unigenes_filter[(Unigenes_filter[[test_method]] < 0.05)&(abs(Unigenes_filter$log_fold_change)>log2fc), ]
Kos_set <- unique(Unigenes_Difference$KOEntry)
length(Kos_set)

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
result_combine <- result_combine[!is.na(result_combine$ko), ]
all.ko.hsa <- unique(unlist(strsplit(unlist(result_combine$ko), ';')))
N <- length(all.ko.hsa)
# filter cpd in result_combine
new.Kos_set <- Kos_set[Kos_set %in% all.ko.hsa]
# n diff number
n <- length(new.Kos_set)

phyper_result <- data.frame()
for (variable in result_combine$name) {
  # M specific pathway number 
  # variable <- 'hsa00513'
  specific_pathway_compound <- unlist(strsplit(result_combine[result_combine$name == variable, ][['ko']], ';'))
  M <- length(specific_pathway_compound)
  # m Overlapping number 
  inter_compound <- intersect(specific_pathway_compound, new.Kos_set)
  m <- length(inter_compound)
  p.value <- stats::phyper(q = m - 1, m = M, n = (N - M), 
                           k = n, lower.tail = FALSE)
  # erichment ratio
  prop_reference <- sum(all.ko.hsa %in% new.Kos_set) / length(all.ko.hsa)
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
  ggtitle(" Overview of Enriched Kos Sets (top25)") +
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
png('kos.kegg.enrichemnt.dot.png', width = 3100, height = 2500, res=300)  # Adjust width, height, and resolution as needed
print(p)
dev.off()
tiff('kos.kegg.enrichemnt.dot.tif', width = 10, height = 8, units = "in", res = 300)
print(p)
dev.off()











