# Function: calculate sankey data
# Create a new data frame for translation
calculate_sankey_data <- function(sankey_data, group_names){
  # group_names : rownames, colnames
  # group_names <- c("microbiome", "metabolism")
  translated_df <- data.frame(
    source = rep(rownames(sankey_data), each = ncol(sankey_data)),
    target = rep(colnames(sankey_data), times = nrow(sankey_data)),
    value = as.vector(t(sankey_data))
  )
  # Remove rows where the value is equal to 0
  translated_df <- subset(translated_df, value != 0)
  variable_name <- as.data.frame(unique(c(translated_df$source, translated_df$target)))
  # variable_name$id <- rownames(variable_name)
  colnames(variable_name) <- c('name')
  # Variable name to numbers and create json data
  name_translation <- structure(1:nrow(variable_name)-1, names = variable_name$name)
  json_data <- list() 
  json_data$nodes  <- variable_name
  # set node color two group
  nametogroup <- structure(c(rep(group_names[1], nrow(sankey_data)), 
                             rep(group_names[2], ncol(sankey_data))), 
                           names = c(rownames(sankey_data), colnames(sankey_data)))
  json_data$nodes$type <- nametogroup[variable_name$name]
  # Variable name to numbers
  # translated_df$source <- unlist(name_translation[translated_df$source])
  # translated_df$target <- unlist(name_translation[translated_df$target])
  translated_df$source <- unlist(translated_df$source)
  translated_df$target <- unlist(translated_df$target)
  json_data$links <- translated_df
  # set link color
  json_data$links$type <- ifelse(json_data[["links"]][["value"]] > 0, 'Postive', 'Negative')
  # set all value to 1
  json_data$links$value <- 1
  return(json_data)
}
