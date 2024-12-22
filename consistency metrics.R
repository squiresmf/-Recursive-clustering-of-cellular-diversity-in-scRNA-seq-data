
library(mclust)
library(aricode)
library(dplyr)
library(ggplot2)

# current_dir = ""
# if (.Platform$OS.type == "windows") {
#   current_dir = paste0(current_dir, 'Y:')
# }
# current_dir = paste0(current_dir, '/qiu-lab/Michael Recursive Clustering/')
current_dir <- paste0(getwd(), '/')

SetDirRead <- function() {
  setwd_path = current_dir
  setwd_path <- paste0(setwd_path, 'Downsample ', downsample_index, '/', reference, ' Downsample ', downsample_index)
  setwd(setwd_path)
}

SetDirWrite <- function() {
  setwd_path = current_dir
  setwd_path <- paste0(setwd_path, 'Default')
  if (!dir.exists(setwd_path)) {
    dir.create(setwd_path, recursive = TRUE)
  }
  setwd(setwd_path)
}

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

num_downsamples = 5
clusterings = c('recursive', 'seurat_equivalent')
mean_cluster_consistency_values = list()

for (c_resolution in seq(0.015, 0.15, 0.005)) {
  for (reference in c("PBMC", "adipose", "tonsil", "fetus")) {
    print(paste(reference, c_resolution))

    cell_ids <- vector("list", num_downsamples) # contains the list of cell ids present, for downsample dataset key
    cluster_assignments <- vector("list", num_downsamples) # contains list of cluster assignments indexed by downsampled cell num, for a downsample and clustering type key

    ref_name <- paste0('human_', reference, '_integrated_')
    for (downsample_index in 1:num_downsamples) {
      SetDirRead()

      cell_ids[[downsample_index]] <- readRDS(paste0('cell_ids_', ref_name, ' Downsample ', downsample_index, '.rds'))
      sorted_indices <- order(cell_ids[[downsample_index]])
      cell_ids[[downsample_index]] = cell_ids[[downsample_index]][sorted_indices]

      cluster_assignments[[downsample_index]] = vector("list", length(clusterings))
      names(cluster_assignments[[downsample_index]]) <- clusterings

      for (clustering in clusterings) {
        cluster_assignments[[downsample_index]][[clustering]] <- as.character(readRDS(paste0(clustering, '_cluster_assignments_', ref_name, c_resolution, '_', 0, ' Downsample ', downsample_index, '.rds')))
        cluster_assignments[[downsample_index]][[clustering]] = cluster_assignments[[downsample_index]][[clustering]][sorted_indices]
      }
    }

    num_cells_full <- length(unique(unlist(cell_ids)))
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

    for (downsample_num1 in 1:4) {
      for (downsample_num2 in (downsample_num1 + 1):5) {
        downsample_num_key1 <- as.character(downsample_num1)
        downsample_num_key2 <- as.character(downsample_num2)
        cell_ids_shared = intersect(cell_ids[[downsample_num1]], cell_ids[[downsample_num2]])
        for (clustering in clusterings) {
          cluster_assignments_shared_cell_id1 = unlist(cluster_assignments_full[[downsample_num1]][[clustering]][cell_ids_shared])
          cluster_assignments_shared_cell_id2 = unlist(cluster_assignments_full[[downsample_num2]][[clustering]][cell_ids_shared])

          num_clusters_avg = (length(unique(cluster_assignments_shared_cell_id1)) + length(unique(cluster_assignments_shared_cell_id2)))/2
          ari <- adjustedRandIndex(cluster_assignments_shared_cell_id1, cluster_assignments_shared_cell_id2)
          nmi <- NMI(cluster_assignments_shared_cell_id1, cluster_assignments_shared_cell_id2)

          metrics <- list(
            'ARI' = ari,
            'NMI' = nmi
          )

          if (clustering == 'recursive') {
            method = 'Recursive'
          } else {
            method = 'Single-Pass'
          }

          for (metric_type in names(metrics)) {
            results_df <- rbind(results_df, data.frame(
              Reference = paste0(toupper(substr(reference, 1, 1)), substr(reference, 2, nchar(reference))),
              Method = method,
              Resolution = c_resolution,
              Metric = metric_type,
              Value = metrics[[metric_type]],
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

SetDirWrite()
write.csv(results_df, 'consistency_metrics_df.csv', row.names = FALSE)


