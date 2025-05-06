# ----------------------------------------------------------------
# CNVkit Clustering - CNV Profile Analysis with pheatmap
# ----------------------------------------------------------------

rm(list=ls())

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(dendextend)

# ---- Centralized Configuration ----
config <- list(
  dirs = list(
    res = "D:/CNVkit/model/with_ref/model_res/",
    out = "D:/CNVkit/model/with_ref/model_imgs/"
  )
)

# ---- Function to import call.cns files ----
import_segments <- function() {
  samples <- list.files(config$dirs$res, pattern = ".*-t.*$") %>% 
    setdiff(c("multi_sample", "nfr", "nfrmulti_sample"))
  
  message("Found samples: ", paste(samples, collapse = ", "))
  
  map_dfr(samples, function(sample) {
    file_path <- file.path(config$dirs$res, sample, paste0(sample, "_call.cns"))
    
    if(!file.exists(file_path)) {
      message("Missing file: ", file_path)
      return(NULL)
    }
    
    message("Importing: ", file_path)
    
    tryCatch({
      read_tsv(file_path, show_col_types = FALSE) %>%
        filter(!is.na(cn)) %>%               # Exclude NA values
        mutate(sample_id = sample,
               cn = as.numeric(cn))           # Convert cn to numeric
    }, error = function(e) {
      message("Error importing ", file_path, ": ", e$message)
      return(NULL)
    })
  })
}

# ---- Prepare data for pheatmap ----
prepare_cnv_matrix <- function(segments_data) {
  # Immediately filter to keep only numerical chromosomes from 1 to 22
  segments_data <- segments_data %>%
    dplyr::filter(grepl("^chr[0-9]+$", chromosome)) %>%
    dplyr::mutate(
      # Extract chromosome number by removing "chr"
      chrom_num = as.numeric(gsub("chr", "", chromosome))
    ) %>%
    # Filter only chromosomes from 1 to 22
    dplyr::filter(chrom_num >= 1 & chrom_num <= 22) %>%
    # Sort by chromosome number and position
    dplyr::arrange(chrom_num, start)
  
  # Data aggregation by chromosome for a more compact visualization
  segments_by_chrom <- segments_data %>%
    dplyr::group_by(chromosome, sample_id) %>%
    dplyr::summarize(
      mean_cn = mean(cn, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Pivot data to create chromosome x sample matrix
  cnv_matrix <- segments_by_chrom %>%
    tidyr::pivot_wider(
      names_from = sample_id,
      values_from = mean_cn,
      values_fill = 2  # Fill NA with neutral value (2 copies)
    ) %>%
    tibble::column_to_rownames("chromosome")
  
  # Sort rows (chromosomes) in numerical order
  chrom_order <- paste0("chr", 1:22)
  cnv_matrix <- cnv_matrix[chrom_order, ]
  
  # Ensure it's a numeric matrix
  cnv_matrix <- as.matrix(cnv_matrix)
  
  return(cnv_matrix)
}

remove_outlier_by_small_cluster <- function(cnv_matrix, distance = "manhattan", k = 3) {
  dist_matrix <- dist(t(cnv_matrix), method = distance)
  hclust_obj <- hclust(dist_matrix, method = "ward.D2")
  
  clusters <- cutree(hclust_obj, k = k)
  cluster_sizes <- table(clusters)
  
  # Take smallest cluster → assumed outlier
  outlier_cluster <- as.integer(names(cluster_sizes)[which.min(cluster_sizes)])
  outlier_samples <- names(clusters[clusters == outlier_cluster])
  
  message("Samples identified as outliers in the smallest cluster (k = ", k, "): ", 
          paste(outlier_samples, collapse = ", "))
  
  cnv_matrix[, !colnames(cnv_matrix) %in% outlier_samples, drop = FALSE]
}

remove_singletons_from_dendrogram <- function(mat, n = 2, method = "complete", distance = "manhattan") {
  # Create dendrogram
  dist_mat <- dist(t(mat), method = distance)
  hc <- hclust(dist_mat, method = method)
  dend <- as.dendrogram(hc)
  
  # Function to check if a node is a singleton (single leaf)
  is_singleton <- function(node) {
    return(is.leaf(node))
  }
  
  # Function to get cluster size (number of leaves)
  get_cluster_size <- function(node) {
    if (is.leaf(node)) return(1)
    sum(sapply(node, get_cluster_size))
  }
  
  removed_samples <- c()
  
  # Recursive function to find and remove singletons in smallest clusters
  process_node <- function(node, depth = 0) {
    if (length(removed_samples) >= n || is.leaf(node)) {
      return(NULL)
    }
    
    # If not a leaf, it has exactly 2 children in the binary dendrogram
    left_child <- node[[1]]
    right_child <- node[[2]]
    
    # Calculate cluster sizes
    left_size <- get_cluster_size(left_child)
    right_size <- get_cluster_size(right_child)
    
    # Determine which side is smaller
    if (left_size <= right_size) {
      smaller_child <- left_child
      larger_child <- right_child
    } else {
      smaller_child <- right_child
      larger_child <- left_child
    }
    
    # If the smaller cluster is a singleton, remove it
    if (is_singleton(smaller_child) && !(labels(smaller_child) %in% removed_samples)) {
      removed_samples <<- c(removed_samples, labels(smaller_child))
      message("Level ", depth, ": Removing singleton ", labels(smaller_child), 
              " (sibling cluster contains ", get_cluster_size(larger_child), " samples)")
      
      # Continue with larger sibling node if needed
      if (length(removed_samples) < n && !is.leaf(larger_child)) {
        process_node(larger_child, depth + 1)
      }
    }
    # If smaller cluster is not a singleton, continue search inside
    else if (!is_singleton(smaller_child)) {
      process_node(smaller_child, depth + 1)
      
      # If we haven't reached the limit yet, also check the larger cluster
      if (length(removed_samples) < n) {
        process_node(larger_child, depth + 1)
      }
    }
    # If smaller cluster is not a singleton and we haven't reached the limit
    else if (length(removed_samples) < n) {
      process_node(larger_child, depth + 1)
    }
  }
  
  # Start process from main node
  process_node(dend, 0)
  
  message("Singleton samples removed: ", paste(removed_samples, collapse = ", "))
  return(mat[, !colnames(mat) %in% removed_samples, drop = FALSE])
}


# ---- Custom pheatmap function ----
plot_cnv_pheatmap <- function(cnv_matrix) {
  if (is.null(cnv_matrix)) return(NULL)
  
  # Define color palette (RdBu - red/blue)
  color_scale <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100)
  
  # Create breaks centered on 2
  breaks <- c(seq(0, 2, length.out = 51), seq(2.04, 6, length.out = 50))
  # Chromosomes on ROWS (already ordered) and samples on COLUMNS (clustered)
  p <- pheatmap(
    mat = cnv_matrix,
    color = color_scale,
    breaks = breaks,
    clustering_distance_cols = "manhattan",  # Distance for samples (columns)
    clustering_method = "ward.D2",
    cluster_cols = TRUE,   # Do clustering on samples (columns)
    cluster_rows = FALSE,  # Don't cluster chromosomes (rows)
    show_rownames = TRUE,  # Show chromosome names
    show_colnames = TRUE,  # Show sample names
    #cutree_cols = 2,       # Separate heatmap 
    fontsize_row = 10,     # Font size for chromosomes
    fontsize_col = 8,      # Font size for samples
    main = "CNV profiles in tumor Samples"
  )
  
  # Save heatmap to file separately
  output_file <- file.path(config$dirs$out, "cnv_clustering_pheatmap_man.png")
  png(output_file, width = 1200, height = 800, res = 120)
  print(p)
  dev.off()
  
  message("Heatmap saved in: ", output_file)
  
  return(p)
}

# ---- Main function ----
main <- function() {
  message("Starting CNV clustering analysis with pheatmap...")
  
  # Import segment data
  message("Importing segment data...")
  segments_data <- import_segments()
  
  if (is.null(segments_data) || nrow(segments_data) == 0) {
    stop("No segment data found. Check directory paths.")
  }
  
  message("Imported ", nrow(segments_data), " segments from ", 
          length(unique(segments_data$sample_id)), " samples")
  
  message("Preparing CNV matrix for visualization...")
  cnv_matrix <- prepare_cnv_matrix(segments_data)
  
  # ---- FILTER outliers ----
  message("Removing singleton outliers (one at a time per branch)...")
  cnv_matrix_filtered <- remove_singletons_from_dendrogram(
    cnv_matrix,
    n = 0,
    distance = "manhattan"
  )
  
  message("Creating heatmap with pheatmap...")
  plot_cnv_pheatmap(cnv_matrix_filtered)
  
  message("CNV clustering analysis with pheatmap completed")
}


# Execute the analysis
main()


