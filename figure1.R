# Color
library(ggsci)
pal_jco("default")(10)
AGA_colour_point_alpha_diversity <- "#27408B"
SGA_colour_point_alpha_diversity <- "#FF3030"
AGA_colour_fill_alpha_diversity <- "#7EC0EE"
SGA_colour_fill_alpha_diversity <- "#EEA2AD"
AGA_colour_Veen <- "#1E90FF"
SGA_colour_Veen <- "#FF6A6A"
AGA_colour_PCA <- "#4876FF"
SGA_colour_PCA  <- "#FF4040"
AGA_colour_LDA  <- "#00008B"
SGA_colour_LDA  <- "#8B0000"
AGA_colour_LEFSE  <- "blue"
SGA_colour_LEFSE  <- "red"

# save dir
setwd(save_dir)

# Figure 1.1 Alpha-diversity
############################
######Alpha diversity#######
############################
library(vegan)
colnames(metagenomic_filter_count)
filter_names <- colnames(metagenomic_filter_count)[grepl('AGA|SGA', colnames(metagenomic_filter_count))]
otu_table <- metagenomic_filter_count[, filter_names]
rownames(otu_table) <- paste0("Otu", 1:nrow(otu_table))
otu_table <- t(otu_table)

# Combine OTU table and group information into a data frame
# Filter data by samples
sample_dir <- 'H:/LSS/data'
sample_origin <- read_excel(paste0(sample_dir, '/sample_origin.xlsx'))
colnames(sample_origin)
metabolism_sample <- na.omit(sample_origin[['P20230730001_metagenomic']])
otu_table <- otu_table[metabolism_sample, ]
group <- gsub("AGA\\d+", "AGA", rownames(otu_table))
group <- gsub("SGA\\d+", "SGA", group)
otu_df <- data.frame(Sample = rownames(otu_table), otu_table, Group = group)

# Calculate Shannon Simpson Chao1 ACE
shannon_values <- diversity(otu_table, index = 'shannon')
simpson_values <- diversity(otu_table, index = 'simpson')
# Calculate Chao1 for each sample
chao1_values <- apply(otu_table, 1, function(x) estimateR(x)["S.chao1"])
# Calculate ACE for each sample
ace_values <- apply(otu_table, 1, function(x) estimateR(x)["S.ACE"])
# Shannon Diversity up↑
# Richness
richness <- apply(t(otu_table), 2, function(x) sum(x > 0))
######## phylogenetic tree,################
# Faith's PD

# Combine calculated values into a data frame
alpha_diversity <- data.frame(Sample = rownames(otu_table), 
                              Shannon = shannon_values,
                              Simpson = simpson_values,
                              Chao1 = chao1_values,
                              ACE = ace_values,
                              Richness = richness,
                              Group = group)

library(tidyr)
# Reshape the dataframe to long format
df_long <- pivot_longer(alpha_diversity, cols = c('Shannon', 'Simpson', 'Chao1', 'ACE', 'Richness'),
                        names_to = "variable", values_to = "value")

# Chao1 Shannon ACE Simpson
plot_alpha_diversity <- function(df_long_chao1, my_comparisons, ylab_title){
  p <-  ggplot(df_long_chao1, aes(Group, value, color = 'Group', fill = Group)) +
    stat_boxplot(geom = "errorbar", size=1, color = 'black') +
    geom_boxplot(size=1, width=0.7, outlier.shape = NA, color = 'black') +
    scale_fill_manual(values = c("AGA" = AGA_colour_fill_alpha_diversity, "SGA" = SGA_colour_fill_alpha_diversity))  +
    theme_bw() +
    geom_jitter(size=5, aes(color = Group), show.legend = FALSE, width = 0.25)+ 
    scale_color_manual(values =  c("AGA" = AGA_colour_point_alpha_diversity, "SGA" = SGA_colour_point_alpha_diversity)) +
    theme(axis.text.x = element_text(face = "bold", color = "black"),
          axis.text.y = element_text(face = "bold", color = "black"),
          panel.border = element_rect(color = "black", size = 1.5),
          axis.ticks = element_line(size = 1, color = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.line = element_line(size = 0.5, colour = "black"),
          axis.ticks.y = element_line(color = "black"), 
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title.y = element_text(face = "bold"))+
    # ggsci::scale_color_igv() + # 配色
    xlab('') +
    ylab(ylab_title) +
    stat_compare_means(comparisons = my_comparisons,
                       method = "wilcox.test",
                       paired = F,
                       size = 4,
                       label = "p.signif",
                       # label = "p.format",
                       tip.length = 0, 
                       bracket.size = 1) # 显著性
  return(p)
}

library(ggplot2)
library(ggpubr)
my_comparisons <- list(c("AGA", "SGA"))
df_long_chao1 <- df_long[df_long$variable == 'Chao1', ]
df_long_Shannon <- df_long[df_long$variable == 'Shannon', ]
df_long_ACE <- df_long[df_long$variable == 'ACE', ]
df_long_Simpson <- df_long[df_long$variable == 'Simpson', ]
df_long_Richness <- df_long[df_long$variable == 'Richness', ]
p_chao1 <- plot_alpha_diversity(df_long_chao1, my_comparisons, 'Chao1')
p_Shannon <- plot_alpha_diversity(df_long_Shannon, my_comparisons, 'Shannon')
p_ACE <- plot_alpha_diversity(df_long_ACE, my_comparisons, 'ACE')
p_Simpson <- plot_alpha_diversity(df_long_Simpson, my_comparisons, 'Simpson')
p_Richness <- plot_alpha_diversity(df_long_Richness, my_comparisons, 'Richness')
library(patchwork)
combined_plot <- p_chao1 + p_Shannon + p_ACE + p_Simpson  + plot_layout(nrow = 1)
# combined_plot <- p_chao1 + p_Shannon + p_ACE + p_Simpson + p_Richness + plot_layout(nrow = 1)
combined_plot
png('alpha_diversity.png', width = 4000, height = 1000, res= 300)  # Adjust width, height, and resolution as needed
print(combined_plot)
dev.off()
library('tiff')
tiff('alpha_diversity.tif', width = 16, height = 4, units = "in", res = 300)
print(combined_plot)
dev.off()


# Figure 1.2 NMDS Beta diversity
############################
####NMDS Beta diversity#####
############################
library(vegan)
library(ggplot2)
library(readxl)
microbiome_sample <- na.omit(sample_origin[['P20230730001_metagenomic']])
SGA_sample <- microbiome_sample[grepl('SGA', microbiome_sample)]
AGA_sample <- microbiome_sample[grepl('AGA', microbiome_sample)]


# use abundance file or count file have same result
colnames(metagenomic_filter)
dune <- metagenomic_filter[, c(AGA_sample, SGA_sample)]
dune <- as.data.frame(t(dune))
# sample info
Group <- gsub("AGA\\d+", "AGA", rownames(dune))
Group <- gsub("SGA\\d+", "SGA", Group)
dune.env <- data.frame(SampleID = rownames(dune), 
                       Groups = Group)
#cal bray_curtis distance
distance <- vegdist(dune, method = 'bray')
distance
nmds <- metaMDS(distance, k = 2, trymax = 100)
summary(nmds)
# get stress
stress <- nmds$stress
# get plot data
df <- as.data.frame(nmds$points)
# cbind
df <- cbind(df, dune.env)
# include covariables
dune.env <- merge(dune.env, metadata_covariates, by.x = 'SampleID', by.y = 'SampleID')
dune.env <- dune.env[, -1]

set.seed(123)
permanova_result <- vegan::adonis2(distance ~ Groups + ., data = dune.env,  permutations = 999, method = distance)
PVALUE <- permanova_result$`Pr(>F)`[1]
PVALUE
R2_text <- round(permanova_result$R2[1], 4)
R2_text
stress_text <- paste("Stress  =", round(stress, 4))
stress_text

# Plot PCA results with two groups colored in red and blue
p <- ggplot(df, aes(x = MDS1, y= MDS2, color = Groups)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = Groups), level = 0.95, type = "norm", geom = "polygon", alpha = 0.25, show.legend = FALSE) +
  scale_color_manual(values = c("SGA" = SGA_colour_PCA, "AGA"= AGA_colour_PCA)) +
  labs(title = paste0("NMDS Analysis (", stress_text, ")"),
       x = paste0("NMDS1"),
       y = paste0("NMDS2"),
       color = "Groups") +
  theme_minimal() +
  theme(
    # panel.background = element_rect(fill = "white"),  # Set background color to white
    panel.grid = element_blank(),  # Remove grid lines
    legend.key = element_rect(fill = "white", color = "black", size = 0.5),
    legend.position = c(0.998, 0.998), legend.justification = c(1, 1),  # Adjust legend position
    legend.background = element_rect(color = "black"),
    # plot.background = element_rect(fill = "white"),  # Set plot background color
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Set panel border color and size
    # axis.text = element_text(face = "bold"),   # Adjust the axis text
    axis.title = element_text(face = "bold"),
    axis.ticks = element_line(size = 0.5),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15)  # Center the title
  ) +
  annotate("text", x = -Inf, y = Inf, label = paste0("PERMANOVA P=", PVALUE, '\nAdonis R2=', R2_text), hjust = -0.05, vjust = 1.2, size = 6)
p

png('NMDS_PERMANOVA_analysis.png', width = 2000, height = 2000, res=300)  # Adjust width, height, and resolution as needed
print(p)
dev.off()

tiff('NMDS_PERMANOVA_analysis.tif', width = 6, height = 6, units = "in", res = 300)
print(p)
dev.off()



# Figure 1.4 Veen
############################
########## Veen ############
############################
library(VennDiagram)
library(venn)
library(ggsci)
# select data
veen_set <- list()
metagenomic_filter_count_AGA <- metagenomic_abun[, grepl('AGA', colnames(metagenomic_abun))]
metagenomic_filter_count_AGA[['Rowsum']] <- rowSums(metagenomic_filter_count_AGA)
metagenomic_filter_count_AGA[['Species']] <- metagenomic_abun$Species
metagenomic_filter_count_SGA <- metagenomic_abun[, grepl('SGA', colnames(metagenomic_abun))]
metagenomic_filter_count_SGA[['Rowsum']] <- rowSums(metagenomic_filter_count_SGA)
metagenomic_filter_count_SGA[['Species']] <- metagenomic_abun$Species
# dat <- veen_set
# group species
veen_set <- list()
veen_set[['AGA']] <- metagenomic_filter_count_AGA[['Species']][metagenomic_filter_count_AGA$Rowsum != 0]
veen_set[['SGA']] <- metagenomic_filter_count_SGA[['Species']][metagenomic_filter_count_SGA$Rowsum != 0]

groups <- c('AGA', 'SGA')
# intersect
library(VennDiagram)
library(venn)
library(ggsci)
# select data
veen_set_names <- names(veen_set)
filter_names <- veen_set_names[grepl(paste0(groups, collapse = '|'), veen_set_names)]
dat <- veen_set[filter_names]
# select colors
colors <- pal_lancet("lanonc", alpha = 0.5)(length(dat))
# Veen plot
veen_count <- as.data.frame(extractInfo(dat))
veen_count[['combination']] <- paste0(veen_count$AGA,  veen_count$SGA)
# AGA\ b65\ SGA
veen_color <-  list(
  AGA = AGA_colour_Veen,  
  SGA = SGA_colour_Veen)
venn(2,  sncs = 0.000001, borders = FALSE,  box = FALSE, par = FALSE, lty = "blank")
# A- AGA
area <- getZones("1-") 
transparent_color <- adjustcolor(veen_color[['AGA']], alpha.f = 0.6)
polygon(area[[1]], col = transparent_color, lty = "blank" )
# B- Model
area <- getZones("-1") 
transparent_color <- adjustcolor(veen_color[['SGA']], alpha.f = 0.6)
polygon(area[[1]], col = transparent_color, lty = "blank")
# intersect counts
veen_count$counts <- c(0,  758, 189, 3572)  # 
for (variable in veen_count$combination[2:nrow(veen_count)]) {
  centroid <- getCentroid(getZones(variable))[[1]]
  count_value <- veen_count[veen_count$combination == variable, ][['counts']]
  text(centroid[1], centroid[2], labels = count_value, cex = 1, font = 2)
}
dat_labels <- list(
  AGA = c(100, 700), 
  SGA = c(900, 700)
)
for (variable in names(dat)) {
  # variable = names(dat)[1]
  labels_point <- dat_labels[[variable]]
  text(labels_point[1], labels_point[2], labels = variable, cex = 1, font = 2)
}

p <- recordPlot()
png('Veen.png', width = 2000, height = 2000, res=300)  # Adjust width, height, and resolution as needed
print(p)
dev.off()
tiff('Veen.tif', width = 6, height = 6, units = "in", res = 300)
print(p)
dev.off()


# Figure 1.6 LDA / 1.5 cladogram
############################
########## LDA  ###########
############################
set.seed(123)
library(microeco)
class_directory <- 'H:/LSS/project'
source(paste0(class_directory, '/private.R'))
LDA_SCORE <- 3
head(dataset$otu_table)
head(dataset$tax_table)
head(dataset$sample_table)
# Read data
library(readxl)
otu_tax_table <- read_xlsx(paste0(metagenomic_dir, "/Species_count_format_taxa.xlsx"))
microbiome_sample <- na.omit(sample_origin[['P20230730001_metagenomic']])
AGA_sample <- microbiome_sample[grepl('AGA', microbiome_sample)]
SGA_sample <- microbiome_sample[grepl('SGA', microbiome_sample)]
otu_tax_table[['Sum']] <- rowSums(otu_tax_table[,  c(AGA_sample, SGA_sample)])
otu_tax_table <- otu_tax_table[otu_tax_table$Sum != 0, ]
otu_tax_table[['OTU']] <- paste0('OTU_', 1:nrow(otu_tax_table))
otu_table <- otu_tax_table[, c(AGA_sample, SGA_sample)]
otu_table <- as.data.frame(otu_table)
rownames(otu_table) <- paste0('OTU_', 1:nrow(otu_table))
library(stringr)
tax_table <- str_split(otu_tax_table$Taxonomy, ";", simplify = TRUE)
tax_table <- as.data.frame(tax_table)
rownames(tax_table) <- paste0('OTU_', 1:nrow(tax_table))
colnames(tax_table) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")[1:length(tax_table)]
GROUP <- data.frame(SampleID = c(AGA_sample, SGA_sample), 
                    Group = c(rep("AGA", length(AGA_sample)), rep("SGA", length(SGA_sample))))
rownames(GROUP) <- c(AGA_sample, SGA_sample)
# Create microbiome dataset object
df <- microtable$new(sample_table = GROUP,
                     otu_table = otu_table,
                     tax_table = tax_table,
                     auto_tidy = FALSE)

# Initialize LEfSe analysis
# taxa_level = "Genus" "Species"
# Only save species levels
lefse <- trans_diff$new(dataset = df,
                        method = "lefse",
                        group = "Group",
                        alpha = 0.05,
                        taxa_level = "all",
                        p_adjust_method = "none"
)

refse_diff <- lefse[["res_diff"]]
lefse_diff_select <- refse_diff[refse_diff$LDA > LDA_SCORE, ]
lefse_diff_select[lefse_diff_select$Group == 'AGA', ][['LDA']] <- -lefse_diff_select[lefse_diff_select$Group == 'AGA', ][['LDA']]
lefse_diff_select <- lefse_diff_select[grepl('s__', lefse_diff_select$Taxa), ]
lefse_diff_select[['Species']] <- str_split(lefse_diff_select$Taxa, "\\|", simplify = TRUE)[, 7]
lefse_diff_select <- lefse_diff_select[!grepl('unclassified', lefse_diff_select$Species),  ]
lefse_diff_select <- lefse_diff_select[!grepl('uncultured', lefse_diff_select$Species), ]
lefse_diff_select$Group <- factor(lefse_diff_select$Group, levels = c('AGA', 'SGA'))
lefse_diff_select <- lefse_diff_select[order(lefse_diff_select$LDA), ]


# Plot LDA bar
library(ggplot2)
Histogram_data <- data.frame(
  name = gsub("s__", "",lefse_diff_select$Species),
  value = lefse_diff_select$LDA,  # Random normal values
  group = lefse_diff_select$Group
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
    y = "LDA Score"
  ) +
  scale_fill_manual(name = "Groups", values = c("AGA" = AGA_colour_LDA, "SGA" = SGA_colour_LDA)) +
  theme_minimal() +
  theme(axis.text.x = element_text(face = "bold", angle = 90, hjust = 1),
        axis.text.y = element_text(face = "bold"), 
        axis.title.y = element_text(face = "bold"), 
        axis.ticks = element_line(size = 0.5), 
        panel.background = element_rect(fill = NA),
        # strip.text = element_blank(),
        # strip.text = element_text(size = 12, color = "black", face = "bold"),
        # strip.background = element_rect(fill = 'lightgray', color = "black"),
        legend.position = "none"
        # legend.position = c(0.05, 0.8), 
        # panel.grid.major = element_line(color = "gray")
        # strip.text.x = element_text(angle = 90, hjust = 1, color = "black", size = 12, face = "bold")
  )
print(p)
png('LDA_plot.png', width = 1800, height = 1200, res=200)  # Adjust width, height, and resolution as needed
print(p)
dev.off()
tiff('LDA_plot.tif', width = 9, height = 6, units = "in", res = 300)
print(p)
dev.off()

############ cladogram #############
self <- lefse
color = c(AGA_colour_LEFSE, SGA_colour_LEFSE)
group_order = NULL
use_taxa_num = 200
filter_taxa = NULL
clade_label_level = 4
select_show_labels = NULL
only_select_show = FALSE
sep = "|"
branch_size = 0.2
alpha = 0.2
clade_label_size = 1.5
clade_label_size_add = 5
clade_label_size_log = exp(1)
node_size_scale = 1
node_size_offset = 1
annotation_shape = 22
annotation_shape_size = 5


# plot_diff_cladogram 
# developed based on microbiome Marker
abund_table <- self$abund_table
marker_table <- self$res_diff %>% dropallfactors
# Filter LDA > 3 p< 0.05
marker_table <- marker_table[marker_table$LDA > LDA_SCORE, ]
use_feature_num = nrow(marker_table)
library(stringr)
str_split(marker_table$Taxa, "\\|")
library(dplyr)
# Vector of strings to split
strings_to_split <- unlist(marker_table$Taxa)
# Function to generate prefixes for a single string
generate_prefixes <- function(input_string) {
  split_vector <- unlist(strsplit(input_string, "\\|"))
  prefixes <- c()
  for (i in 1:length(split_vector)) {
    prefixes[i] <- paste(split_vector[1:i], collapse = "|")
  }
  return(prefixes)
}
# Initialize an empty list to store all prefixes
all_prefixes <- list()
# Process each string
for (str in strings_to_split) {
  all_prefixes[[str]] <- generate_prefixes(str)
}
# Convert the list to a dataframe
df <- bind_rows(lapply(all_prefixes, function(x) data.frame(prefixes = x)), .id = "original_string")
abund_table_select <- unique(df$prefixes)
abund_table <-abund_table[rownames(abund_table) %in% abund_table_select, ]
use_taxa_num = nrow(abund_table)
method <- self$method

library(magrittr)
if(! method %in% c("lefse", "rf")){
  stop("This function currently can only be used for method = 'lefse' or 'rf' !")
}
# only Species 
if(self$taxa_level != "all"){
  stop("This function is available only when taxa_level = 'all' !")
}
if(!is.null(use_feature_num)){
  if(use_feature_num > nrow(marker_table)){
    message("Input use_feature_num ", use_feature_num, " larger than available features number ", nrow(marker_table), " ! ", 
            "Use ", nrow(marker_table), " instead of it ...")
  }else{
    message("Select ", use_feature_num, " significant features ...")
    marker_table %<>% .[1:use_feature_num, ]
  }
}
if(only_select_show == T){
  marker_table %<>% .[.$Taxa %in% select_show_labels, ]
}
# color legend order settings
if(is.null(group_order)){
  if(! is.null(self$group_order)){
    color_groups <- self$group_order
  }else{
    color_groups <- marker_table$Group %>% as.character %>% as.factor %>% levels
  }
}else{
  color_groups <- group_order
}

# filter redundant groups
if(! all(color_groups %in% unique(marker_table$Group))){
  tmp_message <- color_groups[! color_groups %in% unique(marker_table$Group)]
  message("Part of groups in group_order, ", paste(tmp_message, collapse = " "), ", not found in the feature table ...")
  color_groups %<>% .[. %in% unique(marker_table$Group)]
}
if(! all(unique(marker_table$Group) %in% color_groups)){
  tmp_message <- unique(marker_table$Group)[! unique(marker_table$Group) %in% color_groups]
  message("Part of groups in the feature table, ", paste(tmp_message, collapse = " "), ", not found in the group_order ...")
  marker_table %<>% .[.$Group %in% color_groups, ]
}
# get the color palette
if(length(color) < length(unique(marker_table$Group))){
  stop("Please provide enough color palette! There are ", length(unique(marker_table$Group)), 
       " groups, but only ", length(color), " colors provideed in color parameter!")
}else{
  color <- color[1:length(unique(marker_table$Group))]
}
library(magrittr)
tree <- private$plot_backgroud_tree(abund_table = abund_table, use_taxa_num = use_taxa_num, filter_taxa = filter_taxa, sep = sep)
tree
# generate annotation
annotation <- private$generate_cladogram_annotation(marker_table, tree = tree, color = color, color_groups = color_groups, sep = sep)
# check again filtered groups
if(!is.null(use_taxa_num)){
  if(! all(color_groups %in% unique(annotation$enrich_group))){
    tmp_message <- color_groups[! color_groups %in% unique(annotation$enrich_group)]
    message("Biomarkers in group(s), ", paste(tmp_message, collapse = " "), ", not found in the background tree!", 
            " Please try to enlarge parameter use_taxa_num from ", use_taxa_num, " to a larger one ...")
    color_groups %<>% .[. %in% unique(annotation$enrich_group)]
  }
}

annotation_info <- dplyr::left_join(annotation, tree$data, by = c("node" = "label")) %>%
  dplyr::mutate(label = .data$node, id = .data$node.y, level = as.numeric(.data$node_class))
hilight_para <- dplyr::transmute(
  annotation_info,
  node = .data$id,
  fill = .data$color,
  alpha = alpha,
  extend = private$get_offset(.data$level)
)
hilights_g <- purrr::pmap(hilight_para, ggtree::geom_hilight)
tree <- purrr::reduce(hilights_g, `+`, .init = tree)

# hilight legend
hilights_df <- dplyr::distinct(annotation_info, .data$enrich_group, .data$color)
hilights_df$x <- 0
hilights_df$y <- 1
# resort the table used for the legend color and text
hilights_df %<>% `row.names<-`(.$enrich_group) %>% .[color_groups, ]
# make sure the right order in legend
hilights_df$enrich_group %<>% factor(., levels = color_groups)

# add legend
tree <- tree + 
  geom_rect(aes(xmin = x, xmax = x, ymax = y, ymin = y, fill = enrich_group), data = hilights_df, inherit.aes = FALSE) +
  guides(fill = guide_legend(title = NULL, order = 1, override.aes = list(fill = hilights_df$color)))


# set nodes color and size
nodes_colors <- rep("white", nrow(tree$data))
nodes_colors[annotation_info$id] <- annotation_info$color
node_size <- node_size_scale*log(tree$data$abd) + node_size_offset
tree$data$node_size <- node_size
tree <- tree + ggtree::geom_point2(aes(size = I(node_size)), fill = nodes_colors, shape = 21)

## add clade labels
clade_label <- dplyr::transmute(
  annotation_info,
  node = .data$id,
  offset = private$get_offset(.data$level)-0.4,
  offset.text = 0,
  angle = purrr::map_dbl(.data$id, private$get_angle, tree = tree),
  label = .data$label,
  fontsize = clade_label_size + log(.data$level + clade_label_size_add, base = clade_label_size_log),
  barsize = 0,
  extend = 0.2,
  hjust = 0.5,
  level = .data$level
) %>% dplyr::arrange(desc(.data$level))

clade_label$offset.text <- unlist(lapply(seq_len(nrow(clade_label)), function(x){
  if(clade_label$angle[x] < 180){ 0.2 }else{ 0 }}))
clade_label$angle <- unlist(lapply(clade_label$angle, function(x){
  if(x < 180){ x - 90 }else{ x + 90 }}))
clade_label_new <- clade_label

# add letters label to replace long taxonomic label
# set ind all TRUE  select_show_labels <- NULL
ind <- ! clade_label$label %in% select_show_labels

ind_num <- sum(ind)

if(ind_num > 0){
  if(ind_num < 27){
    use_letters <- letters
  }else{
    if(ind_num < 326){
      use_letters <- apply(combn(letters, 2), 2, function(x){paste0(x, collapse = "")})
    }else{
      stop("Too much features to be labelled with letters, consider to use use_feature_num parameter to reduce the number!")
    }
  }
  clade_label_new$label_legend <- clade_label_new$label_show <- clade_label_new$label_raw <- clade_label_new$label
  clade_label_new$label_show[ind] <- use_letters[1:ind_num]
  clade_label_new$label_legend[ind] <- paste0(clade_label_new$label_show[ind], ": ", clade_label_new$label[ind])
  clade_label_new$label <- clade_label_new$label_show
  # delete redundant columns to avoid warnings
  clade_label <- clade_label_new %>% .[, which(! colnames(.) %in% c("label_raw", "label_show", "label_legend", "level"))]
}

clade_label_g <- purrr::pmap(clade_label, ggtree::geom_cladelabel)
tree <- purrr::reduce(clade_label_g, `+`, .init = tree)


# if letters are used, add guide labels
if(ind_num > 0){
  guide_label <- clade_label_new[ind, ] %>%
    dplyr::mutate(color = annotation_info$color[match(.data$label_raw, annotation_info$label)])
  tree <- tree + 
    geom_point(data = guide_label, inherit.aes = FALSE, aes_(x = 0, y = 0, shape = ~label_legend), size = 0, stroke = 0) +
    scale_shape_manual(values = rep(annotation_shape, nrow(guide_label))) +
    guides(shape = guide_legend(override.aes = list(
      size = annotation_shape_size, shape = annotation_shape, fill = guide_label$color), ncol = 3))
} #  Set legend to have ncol = 3 columns by guide_legend
tree <- tree + theme(legend.position = "right", legend.title = element_blank()) 

png('Cladogram_plot.png', width = 4000, height = 2000, res=200)  # Adjust width, height, and resolution as needed
print(tree)
dev.off()
tiff('Cladogram_plot.tif', width = 24, height = 12, units = "in", res = 300)
print(tree)
dev.off()



# Figure 1.3 Plot Circos of abundance 
# Phylum
############# Genus #############
metagenomic_count_taxa <- read_xlsx(paste0(metagenomic_dir, '/Species_count_format_taxa.xlsx'))
microbiome_sample <- na.omit(sample_origin[['P20230730001_metagenomic']])
AGA_sample <- microbiome_sample[grepl('AGA', microbiome_sample)]
SGA_sample <- microbiome_sample[grepl('SGA', microbiome_sample)]
# Calculate abundance 
metagenomic_abundance <- metagenomic_count_taxa
metagenomic_abundance[, c(AGA_sample, SGA_sample)] <- as.data.frame(lapply(metagenomic_abundance[, c(AGA_sample, SGA_sample)], function(x) (x / sum(x)) * 100))
# Filter 0.00001
metagenomic_abundance[['Average_Mean']] <- rowMeans(metagenomic_abundance[, grepl('AGA|SGA', colnames(metagenomic_abundance))])
metagenomic_abundance <- metagenomic_abundance[metagenomic_abundance$Average_Mean > 0.00001, ]
library(stringr)
taxo_name <- data.frame(str_split(metagenomic_abundance[["Taxonomy"]], ";", simplify = TRUE))
colnames(taxo_name) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")[1:length(taxo_name)]
#rownames(taxo_name) <- rownames(otu_table_phy)
# Set species replace genus convenient for use
metagenomic_Phylum <- cbind(metagenomic_abundance[, c(AGA_sample, SGA_sample)], Phylum = taxo_name[['Phylum']])
# calculate sum by Genus 
metagenomic_Phylum_summarized <- aggregate(. ~ Phylum, data = metagenomic_Phylum, FUN = sum)
colSums(metagenomic_Phylum_summarized[, AGA_sample]) == colSums(metagenomic_Phylum[, AGA_sample])

Phylum_data <- metagenomic_Phylum_summarized
Phylum_data$row_sum <- rowSums(Phylum_data[, -1])
Phylum_data <- Phylum_data[order(Phylum_data$row_sum, decreasing = TRUE), ]
Phylum_data_classified <- Phylum_data[!grepl('unclassified', Phylum_data$Phylum), ]
Phylum_data_top_10 <- Phylum_data_classified$Phylum[1:5]
sum_Others <- colSums(Phylum_data[!Phylum_data$Phylum %in% Phylum_data_top_10, ][, c(AGA_sample, SGA_sample)])
data_vector <- rbind(Phylum_data_classified[1:5,c(AGA_sample, SGA_sample)], sum_Others)
colSums(data_vector)
data_vector <- as.data.frame(data_vector)
data_vector[['Phylum']] <- c(Phylum_data_top_10, "Others")
data_vector[['SGA_sum']] <- rowSums(data_vector[, grepl("SGA", colnames(data_vector))])
data_vector[['AGA_sum']] <- rowSums(data_vector[, grepl("AGA", colnames(data_vector))])
# Plot by circos
Phylum_abundance_df <- data.frame(
  Samples = data_vector[['Phylum']], 
  SGA = data_vector[['SGA_sum']], 
  AGA = data_vector[['AGA_sum']]
)


############# Genus #############
metagenomic_count_taxa <- read_xlsx(paste0(metagenomic_dir, '/Species_count_format_taxa.xlsx'))
microbiome_sample <- na.omit(sample_origin[['P20230730001_metagenomic']])
AGA_sample <- microbiome_sample[grepl('AGA', microbiome_sample)]
SGA_sample <- microbiome_sample[grepl('SGA', microbiome_sample)]
# Calculate abundance 
metagenomic_abundance <- metagenomic_count_taxa
metagenomic_abundance[, c(AGA_sample, SGA_sample)] <- as.data.frame(lapply(metagenomic_abundance[, c(AGA_sample, SGA_sample)], function(x) (x / sum(x)) * 100))
# Filter 0.00001
metagenomic_abundance[['Average_Mean']] <- rowMeans(metagenomic_abundance[, grepl('AGA|SGA', colnames(metagenomic_abundance))])
metagenomic_abundance <- metagenomic_abundance[metagenomic_abundance$Average_Mean > 0.00001, ]
library(stringr)
taxo_name <- data.frame(str_split(metagenomic_abundance[["Taxonomy"]], ";", simplify = TRUE))
colnames(taxo_name) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")[1:length(taxo_name)]
#rownames(taxo_name) <- rownames(otu_table_phy)
# Set species replace genus convenient for use
metagenomic_genus <- cbind(metagenomic_abundance[, c(AGA_sample, SGA_sample)], Genus = taxo_name[['Genus']])
# calculate sum by Genus 
metagenomic_genus_summarized <- aggregate(. ~ Genus, data = metagenomic_genus, FUN = sum)
colSums(metagenomic_genus_summarized[, AGA_sample]) == colSums(metagenomic_genus[, AGA_sample])

Genus_data <- metagenomic_genus_summarized
Genus_data$row_sum <- rowSums(Genus_data[, -1])
Genus_data <- Genus_data[order(Genus_data$row_sum, decreasing = TRUE), ]
Genus_data_classified <- Genus_data[!grepl('unclassified', Genus_data$Genus), ]
Genus_data_top_10 <- Genus_data_classified$Genus[1:10]
sum_Others <- colSums(Genus_data[!Genus_data$Genus %in% Genus_data_top_10, ][, c(AGA_sample, SGA_sample)])
taxa <- c(Genus_data_top_10, 'Others')
data_vector <- rbind(Genus_data_classified[1:10,c(AGA_sample, SGA_sample)], sum_Others)
colSums(data_vector)
data_vector <- as.data.frame(data_vector)
data_vector[['Genus']] <- c(Genus_data_top_10, "Others")
data_vector[['SGA_sum']] <- rowSums(data_vector[, grepl("SGA", colnames(data_vector))])
data_vector[['AGA_sum']] <- rowSums(data_vector[, grepl("AGA", colnames(data_vector))])
# Plot by circos
Genus_abundance_df <- data.frame(
  Samples = data_vector[['Genus']], 
  SGA = data_vector[['SGA_sum']], 
  AGA = data_vector[['AGA_sum']]
)

library(openxlsx)
write.xlsx(Genus_abundance_df, 'Genus_abundance_df.xlsx')
write.xlsx(Phylum_abundance_df, 'Phylum_abundance_df.xlsx')
# circos plot



# Histogram plot
############################
#### Taxa Diff Analysis  ###
############################
log2fc <- 0.585
test_method <- 'rfit.p_value'
# Histogram illustrating the differential species significant 
# fold change ≥ 1.5 or ≤ 0.67; P < 0.05
# Load script
utils_direcotry <- 'H:/LSS/project/utils'
source(paste0(utils_direcotry, '/new_calc_microbiome_diff.R'))
metagenomic_Taxonomy_count <- read_xlsx(paste0(metagenomic_dir, '/Species_count_format_taxa.xlsx'))
sample_dir <- 'H:/LSS/data'
sample_origin <- read_excel(paste0(sample_dir, '/sample_origin.xlsx'))
microbiome_sample <- na.omit(sample_origin[['P20230730001_metagenomic']])
AGA_samples <- microbiome_sample[grepl('AGA', microbiome_sample)]
SGA_samples <- microbiome_sample[grepl('SGA', microbiome_sample)]
filter_colnames <- c('Taxonomy', c(AGA_samples, SGA_samples))
metagenomic_Taxonomy_count <- metagenomic_Taxonomy_count[, filter_colnames]
AGA_SGA_columns <- colnames(metagenomic_Taxonomy_count)[grepl('AGA|SGA', colnames(metagenomic_Taxonomy_count))]
metagenomic_Taxonomy_count[, AGA_SGA_columns] <- as.data.frame(lapply(metagenomic_Taxonomy_count[, AGA_SGA_columns], function(x) (x / sum(x)) * 100))
colSums(metagenomic_Taxonomy_count[, AGA_SGA_columns])
AGA_SGA_group <- rbind(data.frame(sample = AGA_samples, group = rep('AGA', length(AGA_samples))), 
                       data.frame(sample = SGA_samples, group = rep('SGA', length(SGA_samples))))
AGA_SGA_group$group <- factor(AGA_SGA_group$group, levels = c('AGA', 'SGA'))

Taxa_diff <- new_calc_microbiome_diff(metagenomic_Taxonomy_count, AGA_SGA_group, metadata_covariates)
Taxa_diff[['Taxonomy']] <- metagenomic_Taxonomy_count$Taxonomy
Taxa_diff[['Species']] <- as.data.frame(str_split(metagenomic_Taxonomy_count$Taxonomy, ';', simplify = TRUE))[['V7']]
# relative abundance > 0.00001 prevalence > 0.1 top
Taxa_diff <- Taxa_diff[(Taxa_diff$mean_AGA > 0.00001) & (Taxa_diff$mean_SGA > 0.00001) , ]
Taxa_diff <- Taxa_diff[!grepl('unclassified|uncultured', Taxa_diff$Species), ]


# Difference taxa
library(stringr)
taxonomy_names <-  as.data.frame(str_split(Taxa_diff$Taxonomy, ';', simplify = TRUE))
# SuperKingdom	Phylum	Class	Order	Family	Genus	Species
Taxa_diff$Phylum <- taxonomy_names$V2
# rename Phylum
Phylum_abundance_df$Samples[1:4]
Taxa_diff$Phylum <- sapply(Taxa_diff$Phylum, function(x) ifelse(x %in% Phylum_abundance_df$Samples[1:4], x,  'Others'))
Taxa_diff$Phylum <- gsub('p__', '', Taxa_diff$Phylum)
Taxa_diff$Species <- taxonomy_names$V7
Taxa_diff$Species <- gsub('s__', '', Taxa_diff$Species)
Taxa_filter2 <- Taxa_diff[Taxa_diff[[test_method]] < 0.05 & abs(Taxa_diff$log_fold_change) >= log2fc, ]
unique(Taxa_filter2$Phylum)
Taxa_filter2$Phylum <- factor(Taxa_filter2$Phylum, levels = c('Actinobacteria', 'Bacteroidetes', 'Firmicutes', 'Proteobacteria', 'Others'))
Taxa_filter2 <- Taxa_filter2[order(Taxa_filter2$log_fold_change),]
define_order <- c('Actinobacteria', 'Bacteroidetes', 'Firmicutes', 'Proteobacteria', 'Others')
sorted_indices <- order(match(Taxa_filter2$Phylum, define_order))
Taxa_filter2 <- Taxa_filter2[sorted_indices, ]
Taxa_filter2$Species <- factor(Taxa_filter2$Species, levels = Taxa_filter2$Species)
table(Taxa_filter2$Phylum)
# Plot hisgram
colnames(Taxa_filter2)
p <- ggplot(Taxa_filter2, aes(x = Species, y = log_fold_change, fill = Phylum)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6, color = NA) +
  # set the color to NA
  labs(
    title = " ",
    x = " ",
    y = "SGA vs AGA log2(FC)"
  ) +
  scale_fill_manual(name = "", values = c("Actinobacteria" = "#A52A2A", 
                                          "Bacteroidetes" = "#FF6EB4", 
                                          "Firmicutes" = "#63B8FF", 
                                          "Proteobacteria" = "#0000FF", 
                                          "Others" = "#00008B")) +
  annotate("text", x = 7, y = -12, label = "Species", color = "black", size = 12, angle = 0, fontface = "bold") +
  annotate("text", x = 13, y = 12, label = "Actinobacteria", color = "#A52A2A", size = 12, angle = 0, fontface = "bold") +
  annotate("text", x = 30, y = 12, label = "Bacteroidetes", color = "#FF6EB4", size = 12, angle = 0, fontface = "bold") +
  annotate("text", x = 130, y = 12, label = "Firmicutes", color = "#63B8FF", size = 12, angle = 0, fontface = "bold") +
  annotate("text", x = 220, y = 12, label = "Proteobacteria", color = "#0000FF", size = 12, angle = 0, fontface = "bold") +
  annotate("text", x = 238, y = 12, label = "Others", color = "#00008B", size = 12, angle = 0, fontface = "bold") +
  theme_minimal() +
  theme(# panel.grid.major.y = element_blank(),
    # panel.grid.major = element_blank(),
    panel.grid.major.y = element_line(color = 'grey80', linetype = "solid"), 
    axis.text.x = element_text(face = "bold", angle = 90, hjust = 1, size = 12),
    axis.text.y = element_text(face = "bold", size = 12), 
    axis.title.y = element_text(face = "bold", size = 30), 
    axis.ticks = element_line(size = 0.5), 
    panel.background = element_rect(fill = NA),
    legend.position = "none"
  )

print(p)

png('Histogram_illustrating_the diferential_species.png', width = 10000, height = 2000, res=200)  # Adjust width, height, and resolution as needed
print(p)
dev.off()
tiff('Histogram_illustrating_the diferential_species.tif', width = 50, height = 10, units = "in", res = 300)
print(p)
dev.off()













