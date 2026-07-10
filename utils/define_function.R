# Replace values to metabolism name
replace_names <- function(names, original_names, modified_names) {
  # Loop through each original name and its corresponding modified name
  query_data <- cbind(original_names, modified_names)
  colnames(query_data) <- c('original_names', 'modified_names')
  query_data <- as.data.frame(query_data)
  # Split names by ', ' delimiter note： space char
  split_names <- strsplit(names, "; ")
  # Replace original names with modified names
  names <- sapply(split_names, function(x) {
    for (variable in x) {
      data_filtered <- subset(query_data, original_names == variable, select = c('modified_names'))
      if(length(unlist(data_filtered)) == 1){
        x[x == variable] <- data_filtered[[1]]
      } else if(length(unlist(data_filtered)) > 1){
        x[x == variable] <- paste0(unique(data_filtered$modified_names), collapse = '; ')
      }
    }
    return(paste0(x, collapse = "; "))
  })
  return(names)
}