

num_downsamples = 5
clusterings = c('recursive', 'seurat_equivalent')

results_df <- data.frame(
  Reference = character(),
  Method = character(),
  Resolution = numeric(),
  Metric = character(),
  Value = numeric(),
  Number_of_Clusters = numeric(),
  Downsample_Num1 = numeric(),
  Downsample_Num2 = numeric(),
  stringsAsFactors = FALSE
)

# current_dir = ""
# if (.Platform$OS.type == "windows") {
#   current_dir = paste0(current_dir, 'Y:')
# }
# current_dir = paste0(current_dir, '/qiu-lab/Michael Recursive Clustering/')
current_dir <- getwd()

SetDirRead <- function() {
  setwd_path <- current_dir
  setwd_path <- paste0(setwd_path, 'Downsample ', downsample_index, '/', reference, ' Downsample ', downsample_index)
  setwd(setwd_path)
  print(setwd_path)
}

SetDirWrite <- function() {
  setwd_path = current_dir
  setwd_path = paste0(setwd_path, 'Default')
  setwd(setwd_path)
  print(setwd_path)
}

for (c_resolution in seq(0.015, 0.15, 0.005)) {
  for (reference in c("PBMC", "adipose", "tonsil", "fetus")) {

    cell_ids <- vector("list", num_downsamples) # contains the list of cell ids present, for downsample dataset key
    cluster_list_cell_ids <- vector("list", num_downsamples) # contains the list of cell ids present in a cluster, for a downsample and clustering type key
    cluster_assignments <- vector("list", num_downsamples) # contains list of cluster assignments indexed by downsampled cell num, for a downsample and clustering type key

    ref_name <- paste0('human_', reference, '_integrated_')
    for (downsample_index in 1:num_downsamples) {
      SetDirRead()

      cell_ids[[downsample_index]] <- readRDS(paste0('cell_ids_', ref_name, ' Downsample ', downsample_index, '.rds'))
      sorted_indices <- order(cell_ids[[downsample_index]])
      cell_ids[[downsample_index]] = cell_ids[[downsample_index]][sorted_indices]

      cluster_list_cell_ids[[downsample_index]] = vector("list", length(clusterings))
      cluster_assignments[[downsample_index]] = vector("list", length(clusterings))
      names(cluster_list_cell_ids[[downsample_index]]) <- clusterings
      names(cluster_assignments[[downsample_index]]) <- clusterings

      for (clustering in clusterings) {
        cluster_assignments[[downsample_index]][[clustering]] <- as.character(readRDS(paste0(clustering, '_cluster_assignments_', ref_name, c_resolution, '_', 0, ' Downsample ', downsample_index, '.rds')))
        cluster_assignments[[downsample_index]][[clustering]] = cluster_assignments[[downsample_index]][[clustering]][sorted_indices]
        unique_clusters <- unique(cluster_assignments[[downsample_index]][[clustering]])
        cluster_list_cell_ids[[downsample_index]][[clustering]] = vector("list", length(unique_clusters))
        names(cluster_list_cell_ids[[downsample_index]][[clustering]]) <- unique_clusters

        for (cluster in unique_clusters) {
          cluster_list_cell_ids[[downsample_index]][[clustering]][[as.character(cluster)]] <- cell_ids[[downsample_index]][cluster_assignments[[downsample_index]][[clustering]] == cluster]
        }
      }
    }

    num_cells_full <- length(unique(unlist(cell_ids)))

    cell_id_location <- vector("list", num_cells_full) # the list of downsampled dataset numbers which contain a cell_id, for a cell_id key

    # Populate cell_id_location with downsampled dataset numbers for each cell_id
    for (downsample_index in 1:num_downsamples) {
      for (cell_id in cell_ids[[downsample_index]]) {
        cell_id_location[[cell_id]] <- c(cell_id_location[[cell_id]], downsample_index)
      }
    }

    cluster_assignments_full <- vector("list", num_downsamples) # contains list of cluster assignments indexed by full dataset cell_id, for a downsample and clustering type key

    for (downsample_index in 1:num_downsamples) {
      cluster_assignments_full[[downsample_index]] = vector("list", length(clusterings))
      names(cluster_assignments_full[[downsample_index]]) <- clusterings

      for (clustering in clusterings) {
        cluster_assignments_full[[downsample_index]][[clustering]] <- vector("list", num_cells_full)

        cell_indices <- cell_ids[[downsample_index]]
        cluster_indices <- cluster_assignments[[downsample_index]][[clustering]]

        # Assign cluster values for cells that are present in the downsampled dataset
        cluster_assignments_full[[downsample_index]][[clustering]][cell_indices] <- cluster_indices
      }
    }

    cluster_comparisons <- list() # The results from comparing a cluster from one downsample dataset and a cluster from another downsample dataset, for both comparison directions
    cluster_shared_cells = list() # The cells from a cluster from downsample key 1 that are present anywhere in downsample key 2

    # Comparing a cell X from downsample key 1 to the same cell X from downsample key 2, what percentage
    # of cell X's cluster members in downsample key 1 are present in cell X's cluster in downsample key 2,
    # considering only cells which are present in both datasets.
    cell_consistency_values <- list()

    for (downsample_num1 in 1:5) {
      for (downsample_num2 in 1:5) {
        if (downsample_num1 != downsample_num2) {
          downsample_num_key1 <- as.character(downsample_num1)
          downsample_num_key2 <- as.character(downsample_num2)
          for (clustering in clusterings) {
            cell_consistency_values[[downsample_num_key1]][[downsample_num_key2]][[clustering]] <- vector("list", num_cells_full)
          }
        }
      }
    }

    # Function to compare cluster consistency
    compare_cluster_consistency <- function(cluster1_cells, cluster2_cells) {
      shared_cells <- length(intersect(cluster1_cells, cluster2_cells))
      consistency_value1 <- (shared_cells - 1) / length(cluster1_cells)
      consistency_value2 <- (shared_cells - 1) / length(cluster2_cells)
      return(c(consistency_value1, consistency_value2))
    }

    for (cell_id in 1:num_cells_full) {
      if (cell_id %% 10 == 0) {
        progress <- (cell_id / num_cells_full) * 100
        cat(sprintf("\rProgress: %.2f%%", progress))
      }
      datasets_with_cell <- cell_id_location[[cell_id]]

      # Iterate over each pair of downsampled dataset numbers containing the cell
      for (first_index in 1:(length(datasets_with_cell) - 1)) {
        for (second_index in (first_index + 1):length(datasets_with_cell)) {
          downsample_num1 <- datasets_with_cell[first_index]
          downsample_num2 <- datasets_with_cell[second_index]
          downsample_num_key1 <- as.character(downsample_num1)
          downsample_num_key2 <- as.character(downsample_num2)

          for (clustering in clusterings) {
            cluster1 <- cluster_assignments_full[[downsample_num1]][[clustering]][[cell_id]]
            cluster2 <- cluster_assignments_full[[downsample_num2]][[clustering]][[cell_id]]

            if (is.null(cluster_comparisons[[clustering]][[downsample_num_key1]][[downsample_num_key2]][[cluster1]][[cluster2]])) {
              if (is.null(cluster_shared_cells[[clustering]][[downsample_num_key1]][[downsample_num_key2]][[cluster1]])) {
                cluster1_cells <- cluster_list_cell_ids[[downsample_num1]][[clustering]][[cluster1]]
                cluster1_cells <- intersect(cluster1_cells, cell_ids[[downsample_num2]])
                cluster_shared_cells[[clustering]][[downsample_num_key1]][[downsample_num_key2]][[cluster1]] = cluster1_cells
              } else {
                cluster1_cells = cluster_shared_cells[[clustering]][[downsample_num_key1]][[downsample_num_key2]][[cluster1]]
              }
              if (is.null(cluster_shared_cells[[clustering]][[downsample_num_key2]][[downsample_num_key1]][[cluster2]])) {
                cluster2_cells <- cluster_list_cell_ids[[downsample_num2]][[clustering]][[cluster2]]
                cluster2_cells <- intersect(cluster2_cells, cell_ids[[downsample_num1]])
                cluster_shared_cells[[clustering]][[downsample_num_key2]][[downsample_num_key1]][[cluster2]] = cluster2_cells
              } else {
                cluster2_cells = cluster_shared_cells[[clustering]][[downsample_num_key2]][[downsample_num_key1]][[cluster2]]
              }

              if (!(cell_id %in% cluster1_cells) || !(cell_id %in% cluster2_cells)) {
                stop()
              }

              consistency_values <- compare_cluster_consistency(cluster1_cells, cluster2_cells)
              cluster_comparisons[[clustering]][[downsample_num_key1]][[downsample_num_key2]][[cluster1]][[cluster2]] <- consistency_values
            } else {
              consistency_values <- cluster_comparisons[[clustering]][[downsample_num_key1]][[downsample_num_key2]][[cluster1]][[cluster2]]
            }
            if (is.null(consistency_values)) {
              stop()
            }
            if (min(consistency_values) < 0) {
              stop()
            }
            cell_consistency_values[[downsample_num_key1]][[downsample_num_key2]][[clustering]][[cell_id]] <- consistency_values[1]
            cell_consistency_values[[downsample_num_key2]][[downsample_num_key1]][[clustering]][[cell_id]] <- consistency_values[2]

          }
        }
      }
    }
    cat("\n")  # Move to the next line after completing the loop

    for (downsample_num1 in 1:5) {
      for (downsample_num2 in 1:5) {
        if (downsample_num1 != downsample_num2) {
          downsample_num_key1 <- as.character(downsample_num1)
          downsample_num_key2 <- as.character(downsample_num2)
          for (clustering in clusterings) {
            mean_value = mean(unlist(cell_consistency_values[[downsample_num_key1]][[downsample_num_key2]][[clustering]]))
            # Determine the method name based on clustering
            if (clustering == 'recursive') {
              method <- 'Recursive'
            } else {
              method <- 'Single-Pass'
            }

            # Compute the average number of clusters between the two downsampling results
            num_clusters1 <- length(unique(cluster_assignments_full[[downsample_num1]][[clustering]]))
            num_clusters2 <- length(unique(cluster_assignments_full[[downsample_num2]][[clustering]]))
            num_clusters_avg <- (num_clusters1 + num_clusters2) / 2

            # Append the metric to results_df
            results_df <- rbind(results_df, data.frame(
              Reference = paste0(toupper(substr(reference, 1, 1)), substr(reference, 2, nchar(reference))),
              Method = method,
              Resolution = c_resolution,
              Metric = 'Shared Cluster Members',
              Value = mean_value,
              Number_of_Clusters = num_clusters_avg,
              Downsample_Num1 = downsample_num1,
              Downsample_Num2 = downsample_num2
            ))
          }
        }
      }
    }
  }
}
colnames(results_df) <- gsub("_", " ", colnames(results_df))


# Write the results_df to a CSV file
write.csv(results_df, 'shared_cluster_members.csv', row.names = FALSE)