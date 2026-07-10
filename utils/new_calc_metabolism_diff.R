# Function: Calculate different between two group
new_calc_metabolism_diff <- function(gene_expression, group) {
  # gene_expression <- AGA_SGA_metabolomi_Data
  # group <- AGA_SGA_group
  sample_rownames <- rownames(gene_expression)
  level_names <- levels(group[['group']])
  level_sample1 <- group[group$group %in% level_names[1], ][['sample']]
  level_sample2 <- group[group$group %in% level_names[2], ][['sample']]
  gene_expression <- gene_expression[, group[['sample']]]
  data_frame <- gene_expression

  # Perform t-wilcox.test for each gene
  wilcox_p_values <- apply(gene_expression, 1, function(expression) {
    # correct = FALSE
    wilcox_result <- wilcox.test(expression[level_sample1], 
                                 expression[level_sample2])
    return(wilcox_result$p.value)
  })
  
  # Calculate the mean expression for each group
  mean_Group1 <- rowMeans(gene_expression[, level_sample1])
  mean_Group2 <- rowMeans(gene_expression[, level_sample2])
  # Function to calculate the percentage of zeros in a row
  percentage_zeros <- function(row) {
    return(sum(row != 0) / length(row) * 100)
  }
  ratio_Group1<- apply(gene_expression[, level_sample1], 1, percentage_zeros)
  ratio_Group2<-  apply(gene_expression[, level_sample2], 1, percentage_zeros)
  # Calculate the fold change for each gene
  fold_change <- mean_Group2 / mean_Group1
  # Calculate the log-fold change for each gene
  log_fold_change <- log2(fold_change)

  data_frame[paste0('mean_', level_names[1])] <- mean_Group1
  data_frame[paste0('mean_', level_names[2])] <- mean_Group2
  data_frame[paste0('ratio_', level_names[1])] <- ratio_Group1
  data_frame[paste0('ratio_', level_names[2])] <- ratio_Group2
  data_frame$fold_change <- fold_change
  data_frame$log_fold_change <- log_fold_change
  data_frame$wilcox.test.p_value <- wilcox_p_values

  sum_abudance <- rowSums(gene_expression)
  data_frame$Sum <- sum_abudance
  data_frame$id <- sample_rownames
  result <- data_frame[c('id', 'fold_change', 'log_fold_change', 
                          'wilcox.test.p_value',
                          'Sum', paste0('mean_', level_names[1]), paste0('mean_', level_names[2]), paste0('ratio_', level_names[1]), paste0('ratio_', level_names[2]),
                          colnames(gene_expression))]
  return(result)
}




