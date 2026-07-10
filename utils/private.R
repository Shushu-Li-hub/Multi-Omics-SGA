
private = list(
  check_taxa_level_all = function(taxa_level){
    if(taxa_level == "all"){
      stop("The taxa_level parameter cannot be 'all' for this method! Please provide a taxonomic level, such as 'Genus' !")
    }
  },
  # group test for lefse or rf
  test_mark = function(dataframe, group, min_num_nonpara = 1, method = NULL){
    d1 <- as.data.frame(t(dataframe))
    taxaname <- colnames(d1)[1]
    d1$Group <- group
    colnames(d1)[1] <- "Value"
    formu <- reformulate("Group", "Value")
    if(any(table(as.character(group))) < min_num_nonpara){
      list(p_value = NA, med = NA)
    }else{
      if(! is.null(method)){
        method <- match.arg(method, c("wilcox.test", "kruskal.test"))
        if(method == "wilcox.test"){
          res1 <- wilcox.test(formula = formu, data = d1)
        }else{
          res1 <- kruskal.test(formula = formu, data = d1)
        }
      }else{
        if(length(unique(as.character(group))) == 2){
          res1 <- wilcox.test(formula = formu, data = d1)
        }else{
          res1 <- kruskal.test(formula = formu, data = d1)
        }
      }
      if(is.nan(res1$p.value)){
        res1$p.value <- 1
      }
      med <- tapply(d1[,1], group, median) %>% as.data.frame
      colnames(med) <- taxaname
      list(p_value = res1$p.value, med = med)
    }
  },
  check_taxa_number = function(sel_taxa, p_adjust_method){
    if(sum(sel_taxa) == 0){
      if(p_adjust_method == "none"){
        stop('No significant feature found!')
      }else{
        stop('No significant feature found! To disable p value adjustment, please use p_adjust_method = "none"!')
      }
    }
    if(sum(sel_taxa) == 1){
      if(p_adjust_method == "none"){
        stop('Only one significant feature found! Stop running subsequent process!')
      }else{
        stop('Only one significant feature found! Stop running subsequent process! To disable p value adjustment, please use p_adjust_method = "none"!')
      }
    }
  },
  generate_microtable_unrel = function(use_dataset, taxa_level, filter_thres, filter_features){
    use_dataset$tidy_dataset()
    suppressMessages(use_dataset$cal_abund(rel = FALSE))
    use_feature_table <- use_dataset$taxa_abund[[taxa_level]]
    if(filter_thres > 0){
      use_feature_table %<>% .[! rownames(.) %in% names(filter_features), ]
    }
    message("Available feature number: ", nrow(use_feature_table))
    newdata <- microtable$new(otu_table = use_feature_table, sample_table = use_dataset$sample_table)
    newdata$tidy_dataset()
    newdata
  },
  # plot the background tree according to raw abundance table
  plot_backgroud_tree = function(abund_table, use_taxa_num = NULL, filter_taxa = NULL, sep = "|"){
    # filter the taxa with unidentified classification or with space, in case of the unexpected error in the following operations
    abund_table %<>% {.[!grepl("\\|.__\\|", rownames(.)), ]} %>%
      {.[!grepl("\\s", rownames(.)), ]} %>%
      # also filter uncleared classification to make it in line with the lefse above
      {.[!grepl("Incertae_sedis|unculture", rownames(.), ignore.case = TRUE), ]}
    if(nrow(abund_table) <= 2){
      stop("After filtering out non-standard taxonomy information, the abundance table only has ", nrow(abund_table), " feature(s)! ", 
           "Is there an issue with the taxonomy table? ", 
           "Please first use the tidy_taxonomy function to process the taxonomy information table before constructing the microtable object.")
    }
    if(!is.null(use_taxa_num)){
      if(use_taxa_num < nrow(abund_table)){
        message("Select ", use_taxa_num, " most abundant taxa as the background cladogram ...")
        abund_table %<>% .[names(sort(apply(., 1, mean), decreasing = TRUE)[1:use_taxa_num]), ]
      }else{
        message("Provided use_taxa_num: ", use_taxa_num, " >= ", " total effective taxa number. Skip the selection ...")
      }
    }
    if(!is.null(filter_taxa)){
      abund_table %<>% .[apply(., 1, mean) > (self$lefse_norm * filter_taxa), ]
    }
    abund_table %<>% .[sort(rownames(.)), ]
    tree_table <- data.frame(taxa = row.names(abund_table), abd = rowMeans(abund_table), stringsAsFactors = FALSE) %>%
      dplyr::mutate(taxa = paste("r__Root", .data$taxa, sep = sep), abd = .data$abd/max(.data$abd)*100)
    taxa_split <- strsplit(tree_table$taxa, split = sep, fixed = TRUE)
    nodes <- purrr::map_chr(taxa_split, utils::tail, n = 1)
    # check whether some nodes duplicated from bad classification
    if(any(duplicated(nodes))){
      del <- nodes %>% .[duplicated(.)] %>% unique
      for(i in del){
        tree_table %<>% .[!grepl(paste0("\\|", i, "($|\\|)"), .$taxa), ]
      }
      taxa_split <- strsplit(tree_table$taxa, split = sep, fixed = TRUE)
      nodes <- purrr::map_chr(taxa_split, utils::tail, n = 1)
    }
    # add root node
    nodes %<>% c("r__Root", .)
    # levels used for extending clade label
    label_levels <- purrr::map_chr(nodes, ~ gsub("__.*$", "", .x)) %>%
      factor(levels = rev(unlist(lapply(taxa_split, function(x) gsub("(.)__.*", "\\1", x))) %>% .[!duplicated(.)]))
    
    # root must be a parent node
    nodes_parent <- purrr::map_chr(taxa_split, ~ .x[length(.x) - 1]) %>% c("root", .)
    
    # add index for nodes
    is_tip <- !nodes %in% nodes_parent
    index <- vector("integer", length(is_tip))
    index[is_tip] <- 1:sum(is_tip)
    index[!is_tip] <- (sum(is_tip) + 1):length(is_tip)
    
    edges <- cbind(parent = index[match(nodes_parent, nodes)], child = index)
    edges <- edges[!is.na(edges[, 1]), ]
    if(! inherits(edges, "matrix")){
      stop("The parent and child nodes are not correctly recognized! ",
           "Please try to use tidy_taxonomy function to tidy taxonomic table in your microtable object! Then rerun cal_abund function and trans_diff class!")
    }
    # donot label tips
    node_label <- nodes[!is_tip]
    phylo <- structure(list(
      edge = edges, 
      node.label = node_label, 
      tip.label = nodes[is_tip], 
      edge.length = rep(1, nrow(edges)), 
      Nnode = length(node_label)
    ), class = "phylo")
    mapping <- data.frame(
      node = index, 
      abd = c(100, tree_table$abd),
      node_label = nodes, 
      stringsAsFactors = FALSE)
    mapping$node_class <- label_levels
    tree <- tidytree::treedata(phylo = phylo, data = tibble::as_tibble(mapping))
    tree <- ggtree::ggtree(tree, size = 0.2, layout = 'circular')
    tree
  },
  # generate the cladogram annotation table
  generate_cladogram_annotation = function(marker_table, tree, color, color_groups, sep = "|") {
    use_marker_table <- marker_table
    feature <- use_marker_table$Taxa
    label <- strsplit(feature, split = sep, fixed = TRUE) %>% 
      purrr::map_chr(utils::tail, n =1)
    plot_color <- use_marker_table$Group %>% 
      as.character
    for(i in seq_along(color_groups)){
      plot_color[plot_color == color_groups[i]] <- color[i]
    }
    annotation <- data.frame(
      node = label,
      color = plot_color,
      enrich_group = use_marker_table$Group,
      stringsAsFactors = FALSE
    )
    # filter the feature with bad classification
    annotation %<>% .[label %in% tree$data$label, ]
    annotation
  },
  get_angle = function(tree, node){
    if (length(node) != 1) {
      stop("The length of `node` must be 1")
    }
    tree_data <- tree$data
    sp <- tidytree::offspring(tree_data, node)$node
    sp2 <- c(sp, node)
    sp.df <- tree_data[match(sp2, tree_data$node),]
    mean(range(sp.df$angle))
  },
  get_offset = function(x){(x*0.2+0.2)^2},
  # dependent functions in utility
  calculate_metastat = function(inputdata, g, pflag = FALSE, threshold = NULL, B = NULL){
    trans_data <- load_frequency_matrix(input = inputdata)
    res <- detect_differentially_abundant_features(jobj = trans_data, g = g, pflag = pflag, threshold = threshold, B = B)
    res
  }
)
