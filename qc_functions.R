# TODO: Something like this for UMAP by cluster:
# seurat@reductions$umap@cell.embeddings %>% as_tibble(rownames = "barcode") %>% left_join(y = seurat[[]] %>% as_tibble(rownames = "barcode")) %>% ggplot(aes(x = umap_1, y = umap_2, colour = demux_assignment == "patient_348")) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() + NoLegend()
# So that we can see the actual distribution of each patient



# Function to import list of barcode:patient assignments ----------------------------------------------------------

import_demux_meta <- function(sample, doublet_stringency, barcode_base_directory = barcode_base_dir){
  # Identify barcode files for all samples
  message("Finding barcode files for sample (library): ", sample)
  message("Filtering barcodes using doublet stringency: ", doublet_stringency)
  
  target_files <- list.files(path = file.path(barcode_base_dir, sample),
                             pattern = paste0("dbl_stringency_", doublet_stringency, ".barcodes.csv"),
                             full.names = TRUE)
  
  if(length(target_files) == 0) {
    stop("No matching files found for the given sample name and doublet stringency")
  }
  
  # Read barcode files
  demux_meta <-
    lapply(target_files, function(x){
      patient <- basename(x)
      # patient <- str_match(string = patient, pattern = "^([^.]+)")[,2]
      patient <- str_extract(string = patient, pattern = "^[^.]+")
      message("Reading barcode file for: ", patient)
      read_csv(file = x, col_names = "barcode", show_col_types = FALSE) %>%
        mutate(demux_assignment = patient)
    }) %>%
    # Condense to a single tibble
    reduce(bind_rows)
  
  message("Processing complete")
  
  demux_meta %>%
    group_by(demux_assignment) %>%
    summarise(n_cells = n()) %>%
    print()
  
  return(demux_meta)
}
# Note that unassigned are not present in this list, but will just be the remaining cells in the Seurat object




# Function to add metadata to a Seurat object ---------------------------------------------------------------------

add_metadata <- function(object, metadata){
  object[[]] <-
    object[[]] %>%
    as_tibble(rownames = "barcode") %>%
    left_join(y = metadata, by = "barcode") %>%
    column_to_rownames(var = "barcode")
  
  return(object)
}




# Function to create QC metrics then plot them --------------------------------------------------------------------

perform_qc <- function(object, n_dims = 30, clustering_resolution = 0.8, grouping_variable = NULL, dim_reduction = FALSE, plot_out_dir){
  
  ##########
  ### Checks
  ##########
  
  assertthat::assert_that(class(object) == "Seurat",
                          msg = "Object must be a Seurat object")
  
  assertthat::assert_that(rlang::is_integerish(n_dims),
                          msg = "n_dims must be a whole number")
  
  assertthat::assert_that(is.numeric(clustering_resolution),
                          msg = "Clustering resolution must be numeric")
  
  if(!is.null(grouping_variable)){
    assertthat::assert_that(assertthat::is.string(grouping_variable),
                            grouping_variable %in% colnames(object[[]]),
                            msg = "Grouping variable must be a string, and be present in the Seurat object metadata") 
  }
  
  assertthat::assert_that(is.logical(dim_reduction),
                          msg = "dim_reduction must be TRUE or FALSE (No quotation marks)")

  
  
  #########
  ### Setup
  #########
  
  qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.most_abundant", "complexity"#,
                  # "log10complexity"
                  )
  cell_cycle_metrics <- c("S.Score", "G2M.Score")
  all_qc_metrics <- c(qc_metrics, cell_cycle_metrics)
  
  s.genes <- Seurat::cc.genes.updated.2019$s.genes
  g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
  
  
  if(!is.null(grouping_variable)){
    group_var_sym <- rlang::sym(grouping_variable)
  }
  
  n_dims <- as.integer(n_dims)
  
  basic_plot_dir <- file.path(plot_out_dir, "qc_plots-basic")
  clustered_plot_dir <- file.path(plot_out_dir, "qc_plot-clustered")
  dir.create(path = plot_out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(path = basic_plot_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(path = clustered_plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  
  
  ########################  
  ### Calculate QC metrics
  ########################
  
  # Percent mitochondrial
  message("Calculating percentage of mitochondrial reads")
  object[["percent.mt"]] <- PercentageFeatureSet(object = object, pattern = "^MT")
  
  # Percent ribosomal
  message("Calculating percentage of ribosomal reads")
  ribo_genes <- rownames(object) %>% grep(pattern = "^RP[SL]", value = TRUE)
  object[["percent.ribo"]] <- PercentageFeatureSet(object, features = ribo_genes)
  
  # Percent haemoglobin
  message("Calculating percentage of heomoglobin reads")
  hb_genes <- rownames(object) %>% grep(pattern = "^HB", value = TRUE)
  hb_genes <- hb_genes[!hb_genes %in% c("HBEGF", "HBS1L", "HBP1")]
  object[["percent.hb"]] <- PercentageFeatureSet(object, features = hb_genes)
  
  # Percent most abundant gene
  message("Determining most abundant gene per cell")
  max_expression <-
    LayerData(object = object, layer = "counts") %>%
    as_tibble(rownames = "gene") %>%
    pivot_longer(-gene, names_to = "barcode", values_to = "counts") %>%
    group_by(barcode) %>%
    slice_max(counts, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    rename(max_counts = counts,
           max_gene = gene)
  
  object <- add_metadata(object = object, metadata = max_expression)
  
  object[[]] <-
    object[[]] %>%
    mutate(percent.most_abundant = max_counts / nCount_RNA * 100 )
  
  
  # Transcriptome complexity
  # I won't model the expected complexity (nFeature_RNA ~ nCount_RNA) as it isn't linear and it is skewed by outliers
  # It's close to linear once logged and with outliers removed, but outlier removal (eg. < 100 genes) isn't something I want to implement automatically
  # I'll just make violin scatter plots
  message("Calculating transcriptome complexity")
  object[[]] <-
    object[[]] %>%
    mutate(complexity = nFeature_RNA / nCount_RNA,
           log10complexity = log10(complexity + 1))
  # ^ I've calculated this differently to https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html#Plotting_complexity
  
  
  
  ####################
  ### Initial QC plots
  ####################
  
  # If a grouping variable is provided, just plot the number of cells in each group
  if(!is.null(grouping_variable)){
    p <-
      object[[]] %>%
      group_by(!!group_var_sym) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      ggplot(aes(x = !!group_var_sym, y = n)) +
      geom_col() +
      # labs(subtitle = "Top 15 most abundant genes by frequency") +
      theme_bw() +
      guides(x = guide_axis(angle = 65)) +
      theme(axis.title.x = element_blank())
    print(p)
    ggsave(filename = file.path(basic_plot_dir, "01-Bar-cells_per_group.png"),
           plot = p)
  }
  
  # Bar plot of max gene (how many cells for each max gene)
  message("Producing plot of top most abundant genes")
  p <-
    object[[]] %>%
    group_by(max_gene) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    slice_max(n, n = 20) %>%
    ggplot(aes(x = reorder(max_gene, desc(n)), y = n)) +
    geom_col() +
    labs(subtitle = "Top 15 most abundant genes by frequency") +
    theme_bw() +
    guides(x = guide_axis(angle = 65)) +
    theme(axis.title.x = element_blank())
  print(p)
  ggsave(filename = file.path(basic_plot_dir, "01-Bar-most_abundant_gene.png"),
         plot = p)
  
  
  # QC plots: ungrouped 
  message("Producing ungrouped QC plots")
  p <- 
    object[[]] %>%
    as_tibble(rownames = "barcode") %>%
    select(barcode, all_of(qc_metrics)) %>%
    pivot_longer(-barcode, names_to = "metric", values_to = "value") %>%
    ggplot(aes(x = value, y = metric)) +
    geom_density_ridges2(jittered_points = TRUE, quantile_lines = TRUE, quantiles = 2, vline_colour = "dodgerblue", vline_size = 1,
                         position = position_raincloud(height = 0.4, adjust_vlines = FALSE),
                         scale = 0.4, point_size = 0.4, alpha = 1, fill = "grey90") +
    coord_flip() +
    facet_wrap(~ metric, scales = "free") +
    theme_bw()
  # Used suppressMessages to hide "Picking joint bandwidth of 0.012" for every faceted plot
  suppressMessages(print(p))
  suppressMessages(
    ggsave(filename = file.path(basic_plot_dir, "02-Violin-ungrouped_QC_metrics.png"),
           plot = p,
           width = 11, height = 9)
  )
  
  
  
  # QC plots: grouped by user-defined variable
  if (!is.null(grouping_variable)) {
    message("Producing QC plots grouped by ", grouping_variable)
    
    p <- 
      object[[]] %>%
      as_tibble(rownames = "barcode") %>%
      select(barcode, !!group_var_sym, all_of(qc_metrics)) %>%
      pivot_longer(-c(barcode, !!group_var_sym), names_to = "metric", values_to = "values") %>%
      ggplot(aes(x = values, y = !!group_var_sym)) +
      geom_density_ridges2(jittered_points = TRUE, quantile_lines = TRUE, quantiles = 2, vline_colour = "dodgerblue", vline_size = 1,
                           position = position_raincloud(height = 0.4, adjust_vlines = FALSE),
                           scale = 0.4, point_size = 0.3, alpha = 1, fill = "grey90") +
      coord_flip() +
      theme_bw() + facet_wrap(~ metric, scales = "free_y") +
      guides(x = guide_axis(angle = 65))
    
    suppressMessages(print(p))
    ggsave(filename = file.path(basic_plot_dir, paste0("03-Violin-QC_metric_by_", grouping_variable, ".png")),
           plot = p,
           width = 13, height = 9)
  } else {
    message("No grouping variable provided. Skipping grouped plot.")
  }
  
  
  
  ####################################################
  ### Normalisation, dimensional reduction, clustering
  ####################################################
  
  if(dim_reduction){
    # Run a non-optimised data normalisation, scaling, dimension reduction pipeline
    message("Running normalisation and dimensional reduction algorithms")
    message("- Clustering resolution: ", clustering_resolution)
    message("- Using ", n_dims, " dimensions for neighbourhood and UMAP analyses")
    object <-
      object %>%
      NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE) %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 3000, verbose = FALSE) %>%
      ScaleData(features = rownames(object), verbose = FALSE) %>%
      RunPCA(verbose = FALSE) %>% # Defaults to features = VariableFeatures(object)), but better not to specify as these are not yet stored in that object
      FindNeighbors(dims = 1:n_dims, verbose = FALSE) %>%
      FindClusters(resolution = clustering_resolution, verbose = FALSE) %>%
      RunUMAP(dims = 1:n_dims, verbose = FALSE)
    
    n_clusters <- object[[]] %>% pull(seurat_clusters) %>% unique() %>% length()
    message("- Found ", n_clusters, " clusters")
    
    
    
    ####################
    ### Score cell cycle
    ####################
    message("Calculating cell cycle scores")
    object <-
      CellCycleScoring(object = object,
                       s.features = s.genes,
                       g2m.features = g2m.genes)
    
    
    
    #######################
    ### Dim reduction plots
    #######################

    # PCA plot by user-specified group
    if(!is.null(grouping_variable)){
      message("Plotting ", grouping_variable, " onto PCA")
      p <-
        object %>%
        DimPlot(reduction = "pca", group.by = grouping_variable)
      print(p)
      ggsave(filename = file.path(clustered_plot_dir, paste0("01-PCA-", grouping_variable, ".png")),
             plot = p)
    }

    # PCA plot showing QC metrics
    message("Plotting QC metrics onto PCA plot")
    p <-
      object %>%
      FeaturePlot(reduction = "pca",
                  features = all_qc_metrics)
    print(p)
    ggsave(filename = file.path(clustered_plot_dir, "02-PCA-QC_metrics.png"),
           plot = p, scale = 2)

    # PCA plot showing clusters
    message("Plotting clusters onto PCA")
    p <-
      object %>%
      DimPlot(reduction = "pca", label = TRUE, label.box = TRUE, repel = TRUE)
    print(p)
    ggsave(filename = file.path(clustered_plot_dir, "03-PCA-clusters.png"),
           plot = p)


    # UMAP by user-specified group
    if(!is.null(grouping_variable)){
      message("Plotting ", grouping_variable, " onto UMAP")

      plots <-
        lapply(object[[]] %>% pull(group_var_sym) %>% unique(), function(x){
          p <-
            object@reductions$umap@cell.embeddings %>%
            as_tibble(rownames = "barcode") %>%
            left_join(y = object[[]] %>% as_tibble(rownames = "barcode"),
                      by = "barcode") %>%
            mutate(assignment = ifelse(!!group_var_sym == x, !!group_var_sym, "other")) %>%
            ggplot(aes(x = umap_1,
                       y = umap_2,
                       colour = assignment,
                       alpha = assignment)) +
            geom_point(size = 0.2) +
            scale_alpha_manual(values = c(0.4, 0.9)) +
            scale_colour_manual(values = c("grey", "darkorange2")) +
            labs(subtitle = x) +
            theme_bw() +
            theme(plot.subtitle = element_text(colour = "darkorange2", face = "bold")) +
            NoLegend()
        }) %>%
        patchwork::wrap_plots()

      print(plots)
      ggsave(filename = file.path(clustered_plot_dir, paste0("04-UMAP-", grouping_variable, ".png")),
             plot = plots,
             width = 15, height = 9)



      ###########
      #
      #   p <-
      #     object %>%
      #     DimPlot(reduction = "umap", group.by = grouping_variable)
      #   print(p)
      #   ggsave(filename = file.path(clustered_plot_dir, paste0("04-UMAP-", grouping_variable, ".png")),
      #          plot = p)
    }

    ##############

    # UMAP showing clusters
    message("Plotting clusters onto UMAP")
    p <-
      object %>%
      DimPlot(reduction = "umap", label = TRUE, label.box = TRUE, repel = TRUE)
    print(p)
    ggsave(filename = file.path(clustered_plot_dir, "05-UMAP-clusters.png"),
           plot = p)

    # UMAP
    message("Plotting QC metrics onto UMAP")
    p <-
      object %>%
      FeaturePlot(reduction = "umap",
                  features = all_qc_metrics)
    print(p)
    ggsave(filename = file.path(clustered_plot_dir, "06-UMAP-QC_metrics.png"),
           plot = p,
           # scale = 2,
           width = 15, height = 11)


    # QC metrics by cluster
    message("Producing QC plots grouped by cluster")
    p <-
      object[[]] %>%
      as_tibble(rownames = "barcode") %>%
      select(barcode, seurat_clusters, all_of(all_qc_metrics)) %>%
      pivot_longer(-c(barcode, seurat_clusters), names_to = "metric", values_to = "values") %>%
      ggplot(aes(x = values, y = seurat_clusters)) +
      geom_density_ridges2(jittered_points = TRUE, quantile_lines = TRUE, quantiles = 2, vline_colour = "dodgerblue", vline_size = 1,
                           position = position_raincloud(height = 0.4, adjust_vlines = FALSE),
                           scale = 0.4, point_size = 0.5, alpha = 1, fill = "grey90") +
      coord_flip() +
      theme_bw() + facet_wrap(~ metric, scales = "free_y") +
      guides(x = guide_axis(angle = 65))
    suppressMessages(print(p))
    suppressMessages(
      ggsave(filename = file.path(clustered_plot_dir, "07-Violin-QC_metrics_by_cluster.png"),
             plot = p,
             width = 15, height = 9)
    )
  }
  
  return(object)
}




# Function to create QC metrics for FILTERED data, then plot them -------------------------------------------------

filtered_qc <- function(object, clustering_resolution, grouping_variable = NULL, plot_out_dir){
  
  ##########
  ### Checks
  ##########
  
  assertthat::assert_that(class(object) == "Seurat",
                          msg = "Object must be a Seurat object")
  
  assertthat::assert_that(is.numeric(clustering_resolution),
                          msg = "Clustering resolution must be numeric")
  
  if(!is.null(grouping_variable)){
    assertthat::assert_that(assertthat::is.string(grouping_variable),
                            grouping_variable %in% colnames(object[[]]),
                            msg = "Grouping variable must be a string, and be present in the Seurat object metadata") 
  }
  
  assertthat::assert_that(!is.null(object@reductions),
                          msg = "Must run RunUMAP() prior to running this function")
  
  assertthat::assert_that(max(grepl(pattern = "cluster", x = object[[]] %>% colnames())) == 1,
                          msg = "Seurat object must contain cluster information prior to running this function")
  
  
  
  #########
  ### Setup
  #########
  
  qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.most_abundant", "complexity"#,
                  # "log10complexity"
  )
  cell_cycle_metrics <- c("S.Score", "G2M.Score")
  all_qc_metrics <- c(qc_metrics, cell_cycle_metrics)
  
  s.genes <- Seurat::cc.genes.updated.2019$s.genes
  g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
  
  
  if(!is.null(grouping_variable)){
    group_var_sym <- rlang::sym(grouping_variable)
  }
  
  
  clustered_plot_dir <- file.path(plot_out_dir, "qc_plot-clustered")
  dir.create(path = plot_out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(path = clustered_plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  
  
  ####################
  ### Initial QC plots
  ####################
  
  # If a grouping variable is provided, just plot the number of cells in each group
  if(!is.null(grouping_variable)){
    p <-
      object[[]] %>%
      group_by(!!group_var_sym) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      ggplot(aes(x = !!group_var_sym, y = n)) +
      geom_col() +
      theme_bw() +
      guides(x = guide_axis(angle = 65)) +
      theme(axis.title.x = element_blank())
    print(p)
    ggsave(filename = file.path(clustered_plot_dir, "01-Bar-cells_per_group.png"),
           plot = p)
  }
  
  
  # QC plots: grouped by user-defined variable
  if (!is.null(grouping_variable)) {
    message("Producing QC plots grouped by ", grouping_variable)
    
    p <- 
      object[[]] %>%
      as_tibble(rownames = "barcode") %>%
      select(barcode, !!group_var_sym, all_of(all_qc_metrics)) %>%
      pivot_longer(-c(barcode, !!group_var_sym), names_to = "metric", values_to = "values") %>%
      ggplot(aes(x = values, y = !!group_var_sym)) +
      geom_density_ridges2(jittered_points = TRUE, quantile_lines = TRUE, quantiles = 2, vline_colour = "dodgerblue", vline_size = 1,
                           position = position_raincloud(height = 0.4, adjust_vlines = FALSE),
                           scale = 0.4, point_size = 0.3, alpha = 1, fill = "grey90") +
      coord_flip() +
      theme_bw() + facet_wrap(~ metric, scales = "free_y") +
      guides(x = guide_axis(angle = 65))
    
    suppressMessages(print(p))
    ggsave(filename = file.path(clustered_plot_dir, paste0("02-Violin-QC_metric_by_", grouping_variable, ".png")),
           plot = p,
           width = 13, height = 9)
  } else {
    message("No grouping variable provided. Skipping grouped plot.")
  }
  
  
  
  #######################
  ### Dim reduction plots
  #######################
  
  # PCA plot by user-specified group
  if(!is.null(grouping_variable)){
    message("Plotting ", grouping_variable, " onto PCA")
    p <-
      object %>%
      DimPlot(reduction = "pca", group.by = grouping_variable)
    print(p)
    ggsave(filename = file.path(clustered_plot_dir, paste0("03-PCA-", grouping_variable, ".png")),
           plot = p)
  }
  
  # PCA plot showing QC metrics
  message("Plotting QC metrics onto PCA plot")
  p <-
    object %>%
    FeaturePlot(reduction = "pca",
                features = all_qc_metrics)
  print(p)
  ggsave(filename = file.path(clustered_plot_dir, "04-PCA-QC_metrics.png"),
         plot = p, scale = 2)
  
  # PCA plot showing clusters
  message("Plotting clusters onto PCA")
  p <-
    object %>%
    DimPlot(reduction = "pca", label = TRUE, label.box = TRUE, repel = TRUE)
  print(p)
  ggsave(filename = file.path(clustered_plot_dir, "05-PCA-clusters.png"),
         plot = p)
  
  
  # UMAP by user-specified group
  if(!is.null(grouping_variable)){
    message("Plotting ", grouping_variable, " onto UMAP")
    
    plots <-
      lapply(object[[]] %>% pull(group_var_sym) %>% unique(), function(x){
        p <-
          object@reductions$umap@cell.embeddings %>%
          as_tibble(rownames = "barcode") %>%
          left_join(y = object[[]] %>% as_tibble(rownames = "barcode"),
                    by = "barcode") %>%
          mutate(assignment = ifelse(!!group_var_sym == x, !!group_var_sym, "other")) %>%
          ggplot(aes(x = umap_1,
                     y = umap_2,
                     colour = assignment,
                     alpha = assignment)) +
          geom_point(size = 0.2) +
          scale_alpha_manual(values = c(0.4, 0.9)) +
          scale_colour_manual(values = c("grey", "darkorange2")) +
          labs(subtitle = x) +
          theme_bw() +
          theme(plot.subtitle = element_text(colour = "darkorange2", face = "bold")) +
          NoLegend()
      }) %>%
      patchwork::wrap_plots()
    
    print(plots)
    ggsave(filename = file.path(clustered_plot_dir, paste0("06-UMAP-", grouping_variable, ".png")),
           plot = plots,
           width = 15, height = 9)
    
    # UMAP showing clusters
    message("Plotting clusters onto UMAP")
    p <-
      object %>%
      DimPlot(reduction = "umap", label = TRUE, label.box = TRUE, repel = TRUE)
    print(p)
    ggsave(filename = file.path(clustered_plot_dir, "07-UMAP-clusters.png"),
           plot = p)
    
    # UMAP
    message("Plotting QC metrics onto UMAP")
    p <-
      object %>%
      FeaturePlot(reduction = "umap",
                  features = all_qc_metrics)
    print(p)
    ggsave(filename = file.path(clustered_plot_dir, "08-UMAP-QC_metrics.png"),
           plot = p,
           # scale = 2,
           width = 15, height = 11)
    
    
    # QC metrics by cluster
    message("Producing QC plots grouped by cluster")
    p <-
      object[[]] %>%
      as_tibble(rownames = "barcode") %>%
      select(barcode, seurat_clusters, all_of(all_qc_metrics)) %>%
      pivot_longer(-c(barcode, seurat_clusters), names_to = "metric", values_to = "values") %>%
      ggplot(aes(x = values, y = seurat_clusters)) +
      geom_density_ridges2(jittered_points = TRUE, quantile_lines = TRUE, quantiles = 2, vline_colour = "dodgerblue", vline_size = 1,
                           position = position_raincloud(height = 0.4, adjust_vlines = FALSE),
                           scale = 0.4, point_size = 0.5, alpha = 1, fill = "grey90") +
      coord_flip() +
      theme_bw() + facet_wrap(~ metric, scales = "free_y") +
      guides(x = guide_axis(angle = 65))
    suppressMessages(print(p))
    suppressMessages(
      ggsave(filename = file.path(clustered_plot_dir, "09-Violin-QC_metrics_by_cluster.png"),
             plot = p,
             width = 15, height = 9)
    )
  }
  
  # return(object)
  invisible(NULL)
}
