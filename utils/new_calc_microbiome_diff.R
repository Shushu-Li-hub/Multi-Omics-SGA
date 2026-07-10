# Function: Calculate different between two group
new_calc_microbiome_diff <- function(gene_expression, group, covariables_df) {
  # gene_expression <- metagenomic_Taxonomy_count
  # group <- AGA_SGA_group
  # covariables_df <- metadata_covariates
  level_names <- levels(group[['group']])
  level_sample1 <- group[group$group %in% level_names[1], ][['sample']]
  level_sample2 <- group[group$group %in% level_names[2], ][['sample']]
  gene_expression <- gene_expression[, group[['sample']]]
  data_frame <- gene_expression
  num_samples <- ncol(gene_expression)
  
  # Performed multivariable linear regression models, with CLR-transformed taxon abundances
  library(compositions)
  taxa_matrix <- t(gene_expression)
  taxa_clr <- as.data.frame(clr(as.matrix(taxa_matrix + 1e-8)))
  # combine with metadata
  sample_id <- covariables_df$SampleID
  covariables_df <- covariables_df[, -1]
  covariables_df <- cbind(data.frame(SGA = gsub("[0-9]+$", "", sample_id)), covariables_df)
  rownames(covariables_df) = sample_id
  covariables_df <- covariables_df[rownames(taxa_clr), , drop = FALSE]
  analysis_df <- cbind(
    covariables_df,
    taxa_clr
  )
  
  library(future.apply)
  lmfun <- function(taxon) {
    # taxon<- 'V1'
    covars <- colnames(covariables_df)
    form <- as.formula(
      paste(taxon, "~ ", paste(covars, collapse = " + "))
    )
    fit <- lm(form, data = analysis_df)
    return(summary(fit)$coefficients["SGASGA", "Pr(>|t|)"])
  }
  future::plan(future::multisession, workers = parallel::detectCores() - 2)
  lm_adjust_p_values <- future_apply(
    X = matrix(colnames(taxa_clr), nrow = length(colnames(taxa_clr))),
    MARGIN = 1,
    FUN = lmfun
  )
  # Performed rank-based regression models, with CLR-transformed taxon abundances
  library(Rfit)
  rfitfun <- function(taxon) {
    # taxon<- 'V602'
    covars <- colnames(covariables_df)
    form <- as.formula(
      paste(taxon, "~ ", paste(covars, collapse = " + "))
    )
    tryCatch(
      {
        fit <- rfit(form, data = analysis_df)
        return(summary(fit)$coef["SGASGA", "p.value"])
      },
      error = function(e) {
        return(1)
      }
    )
  }
  future::plan(future::multisession, workers = parallel::detectCores() - 2)
  rfit_adjust_p_values <- future_apply(
    X = matrix(colnames(taxa_clr), nrow = length(colnames(taxa_clr))),
    MARGIN = 1,
    FUN = rfitfun
  )
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
  data_frame$rfit.p_value <- rfit_adjust_p_values
  data_frame$lm.p_value <- lm_adjust_p_values
  sum_abudance <- rowSums(gene_expression)
  data_frame$Sum <- sum_abudance
  data_frame$id <- rownames(gene_expression)
  result <- data_frame[c('id', 'fold_change', 'log_fold_change', 'lm.p_value',
                         'rfit.p_value',
                         'Sum', paste0('mean_', level_names[1]), paste0('mean_', level_names[2]), paste0('ratio_', level_names[1]), paste0('ratio_', level_names[2]),
                         colnames(gene_expression))]
  return(result)
}




