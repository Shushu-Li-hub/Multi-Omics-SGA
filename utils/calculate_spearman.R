# Function:
calculate_spearman <- function(microbiome_data, metabolism_data, microbiome_sample, metabolism_sample){
  # Create an empty data frame with NULL values
  # microbiome_data=microbiome_data
  # metabolism_data= functions_data
  # microbiome_sample=microbiome_sample
  # metabolism_sample=functions_sample
  row_names = rownames(microbiome_data)
  column_names = rownames(metabolism_data)
  spearman_table <- data.frame(matrix(NA, nrow = length(row_names), ncol = length(column_names)))
  rownames(spearman_table) <- row_names
  colnames(spearman_table) <- column_names
  
  # Function to calculate Spearman  p-value
  calculate_spearman_pvalue <- function(row1, row2) {
    cor_test_result <- cor.test(row1, row2, method = "spearman")
    return(p_value = cor_test_result$p.value)
  }
  # Function to calculate Spearman correlation
  calculate_spearman_correlation <- function(row1, row2) {
    cor_test_result <- cor.test(row1, row2, method = "spearman")
    return(correlation  = cor_test_result$estimate)
  }
  # Calculate Spearman pvalue between corresponding rows
  spearman_p_value <- apply(microbiome_data[microbiome_sample], 1, function(row1) {
    apply(metabolism_data[metabolism_sample], 1, function(row2) {
      calculate_spearman_pvalue(row1, row2)
    })
  })
  # Calculate Spearman correlation between corresponding rows
  spearman_correlation <- apply(microbiome_data[microbiome_sample], 1, function(row1) {
    apply(metabolism_data[metabolism_sample], 1, function(row2) {
      calculate_spearman_correlation(row1, row2)
    })
  })
  
  # Apply Benjamini-Hochberg (BH) correction
  spearman_adjusted_p_values <- matrix(
    p.adjust(as.vector(spearman_p_value), method = "BH"),
    nrow = nrow(spearman_p_value),
    ncol = ncol(spearman_p_value)
  )
  rownames(spearman_adjusted_p_values) <- rownames(spearman_p_value)
  
  spear_p_result <- as.data.frame(t(spearman_p_value))
  colnames(spear_p_result) <- rownames(metabolism_data)
  rownames(spear_p_result) <- rownames(microbiome_data)
  spear_p_adjust_result <- as.data.frame(t(spearman_adjusted_p_values))
  colnames(spear_p_adjust_result) <- rownames(metabolism_data)
  rownames(spear_p_adjust_result) <- rownames(microbiome_data)
  spear_corr_result <- as.data.frame(t(spearman_correlation))
  colnames(spear_corr_result) <- rownames(metabolism_data)
  rownames(spear_corr_result) <- rownames(microbiome_data)
  return(list('pvalue' = spear_p_result, 
              'adjusted_p_values' = spear_p_adjust_result, 
              'correlation' = spear_corr_result))
}


