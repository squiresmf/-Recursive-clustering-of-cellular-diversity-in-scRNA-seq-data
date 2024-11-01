# library(sabre)
library(mclust)
library(aricode)
library(dplyr)
library(ggplot2)

algorithm = 'Louvain'
# algorithm = 'Leiden'
# # algorithm = 'SLM'
# # HVGs = 1000
HVGs = 2000
# # HVGs = 3000

# current_dir = ""
# if (.Platform$OS.type == "windows") {
#   current_dir = paste0(current_dir, 'Y:')
# }
# current_dir = paste0(current_dir, '/qiu-lab/Michael Recursive Clustering/')
current_dir <- getwd()

SetDirRead <- function() {
  setwd_path = current_dir
  if (name_extra == "") {
    setwd_path = paste0(setwd_path, 'Default/', reference)
  } else {
    setwd_path = paste0(setwd_path, trimws(name_extra), '/', reference, name_extra)
  }
  setwd(setwd_path)
}

SetDirWrite <- function() {
  setwd_path = current_dir
  if (name_extra == "") {
    setwd_path = paste0(setwd_path, 'Default')
  } else {
    setwd_path = paste0(setwd_path, trimws(name_extra))
  }
  setwd(setwd_path)
}

#
# for (algorithm in c('Louvain', 'Leiden', 'SLM')) {
#   # for (algorithm in c('Louvain')) {
#   if (algorithm == 'Louvain') {
#     HVGs_range = c(2000, 3000, 1000)
#     # HVGs_range = c(3000, 1000)
#     # HVGs_range = c(2000)
#   } else {
#     HVGs_range = c(2000)
#   }
#   for (HVGs in HVGs_range) {
fvf.nfeatures = HVGs

name_extra = ""
if (algorithm != 'Louvain') {
  name_extra = paste(name_extra, algorithm)
}
if (HVGs != 2000) {
  name_extra = paste(name_extra, HVGs, 'HVG')
  name_extra = paste(name_extra, fvf.nfeatures, 'fvf.nfeatures')
}

results_df <- data.frame(
  Reference = character(),
  Clustering_Method = character(),
  Annotation_Level = character(),
  Resolution = numeric(),
  Metric = character(),
  Value = numeric(),
  Number_of_Clusters = numeric(),
  stringsAsFactors = FALSE
)

mean_cell_type_cluster_purity <- function(clusters, cell_celltype) {
  cluster_celltype_count_table = as.data.frame.matrix(table(cell_cluster_assignments, cell_celltype))
  clusters = rownames(cluster_celltype_count_table)
  cluster_purity = apply(cluster_celltype_count_table, 1, max) / rowSums(cluster_celltype_count_table)
  cluster_majority_celltype = setNames(object = colnames(cluster_celltype_count_table)[max.col(cluster_celltype_count_table)], nm = clusters)
  celltypes = sort(unique(cell_celltype))
  average_purity <- tapply(cluster_purity, cluster_majority_celltype, mean)

  cell_type_cluster_purity = list()
  for (celltype in celltypes) {
    if (celltype %in% names(average_purity)) {
      cell_type_cluster_purity[[celltype]] = average_purity[[celltype]]
    } else {
      cell_type_cluster_purity[[celltype]] = 0
    }
  }
  return(mean(unlist(cell_type_cluster_purity)))
}

get_f1_results <- function(clusters, cell_celltype) {
  cell_clusters_majority_celltype = vector("character", length(cell_celltype))
  cluster_celltype_count_table = as.data.frame.matrix(table(cell_cluster_assignments, cell_celltype))
  clusters = rownames(cluster_celltype_count_table)
  cluster_majority_celltype = setNames(object = colnames(cluster_celltype_count_table)[max.col(cluster_celltype_count_table)], nm = clusters)
  for (cluster in clusters) {
    cell_clusters_majority_celltype[which(cell_cluster_assignments == cluster)] = cluster_majority_celltype[[cluster]]
  }

  celltypes = sort(unique(cell_celltype))
  predicted_vs_label_celltype_table <- as.data.frame.matrix(table(
    factor(cell_clusters_majority_celltype, levels = celltypes),
    factor(cell_celltype, levels = celltypes)
  ))

  true_positives <- diag(as.matrix(predicted_vs_label_celltype_table))
  false_positives <- rowSums(predicted_vs_label_celltype_table) - true_positives
  false_negatives <- colSums(predicted_vs_label_celltype_table) - true_positives

  precision <- ifelse(true_positives != 0, yes = true_positives / (true_positives + false_positives), no = 0)
  recall = ifelse(true_positives != 0, yes = true_positives / (true_positives + false_negatives), no = 0)
  f1 = ifelse(true_positives != 0, yes = (2 * precision * recall) / (precision + recall), no = 0)

  return(list(f1 = mean(f1), precision = mean(precision), recall = mean(recall)))
}

clusterings = c('recursive', 'seurat_equivalent')
for (reference in c("PBMC", "adipose", "tonsil", "fetus")) {
  if (name_extra == "") {
    if (reference %in% c("PBMC", "adipose", "tonsil")) {
      res_range = c(seq(0.005, 0.05, 0.002), seq(0.015, 0.15, 0.005))
      res_range = sort(unique(round(res_range, 4)))
    } else {
      res_range = c(seq(0.0025, 0.01, 0.0002), seq(0.015, 0.15, 0.005))
      res_range = sort(unique(round(res_range, 4)))
    }
  } else {
    res_range = seq(0.015, 0.15, 0.005)
    res_range = round(res_range, 4)
  }

  ref_name <- paste0('human_', reference, '_integrated_')
  SetDirRead()
  getwd()

  annotation_levels = list()
  annotation_levels[[2]] = if (reference == "tonsil") 3 else 2
  annotation_levels[[1]] = annotation_levels[[2]] - 1
  paste0('reference_human_', reference, '_integrated.rds')
  for (annotation_level in 1:2) {
    if (!file.exists(paste0(reference, ' celltype.l', annotation_level, '.rds'))) {
      integrated_data = readRDS(paste0('reference_human_', reference, '_integrated.rds'))
      cell_celltype = integrated_data@meta.data[[paste0('celltype.l', annotation_levels[[annotation_level]])]]
      saveRDS(cell_celltype, paste0(reference, ' celltype.l', annotation_level, '.rds'))
      cell_celltype_df <- data.frame(celltype = unlist(cell_celltype), stringsAsFactors = FALSE)
      write.csv(cell_celltype_df, paste0(reference, ' celltype.l', annotation_level, '.csv'), row.names = FALSE)
    } else {
      cell_celltype = readRDS(paste0(reference, ' celltype.l', annotation_level, '.rds'))
    }

    SetDirWrite()
    if (!file.exists(paste0(reference, ' celltype.l', annotation_level, ' unique.rds'))) {
      celltypes = sort(unique(cell_celltype))
      saveRDS(celltypes, paste0(reference, ' celltype.l', annotation_level, ' unique.rds'))
    }

    SetDirRead()

    for (c_resolution in res_range) {
      print(paste(reference, c_resolution, annotation_level))
      for (clustering in clusterings) {
        cell_cluster_assignments <- as.character(readRDS(paste0(clustering, '_cluster_assignments_', ref_name, c_resolution, '_', 0, name_extra, '.rds')))
        if (!file.exists(paste0(clustering, '_cluster_assignments_', ref_name, c_resolution, '_', 0, '.csv'))) {
          cell_cluster_assignments_df <- data.frame(cluster_assignment = unlist(cell_cluster_assignments), stringsAsFactors = FALSE)
          write.csv(cell_cluster_assignments_df, paste0(clustering, '_cluster_assignments_', ref_name, c_resolution, '_', 0, '.csv'), row.names = FALSE)
        }
        num_clusters = length(unique(cell_cluster_assignments))
        MCTCP = mean_cell_type_cluster_purity(cell_cluster_assignments, cell_celltype)
        f1_results = get_f1_results(cell_cluster_assignments, cell_celltype)

        metrics <- list(
          'Mean Cell Type Cluster Purity' = MCTCP,
          'F1-score' = f1_results[['f1']],
          'Precision' = f1_results[['precision']],
          'Recall' = f1_results[['recall']]
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
            Annotation_Level = paste('Annotation Level', annotation_level),
            Resolution = c_resolution,
            Metric = metric_type,
            Value = metrics[[metric_type]],
            Number_of_Clusters = num_clusters
          ))
        }

      }
    }
  }
}

colnames(results_df) <- gsub("_", " ", colnames(results_df))

SetDirWrite()

write.csv(results_df, paste0('metrics_df', name_extra, '.csv'), row.names = FALSE)
# }
# }