

library(Seurat)
library(ggplot2)
library(cluster)
library(dplyr)
library(stringr)
library(gridExtra)
library(magick)

cluster_list = c('cRc0c2', 'cRc3c0', 'cRc3c1', 'cRc4c0', 'cRc4c1', 'cRc4c2', 'cRc5', 'cRc8', 'cRc10', 'cRc11c0')

# Modify Seurat object metadata to chnage the cluster naming format from cRc1c2c3 to 1-2-3
modify_meta_data <- function(seurat_obj, target_columns, remove_n_chars = 3,
                             pattern_to_replace = "c", replacement = "-",
                             case_sensitive = TRUE, handle_short_strings = "empty") {

  for (col in target_columns) {
    if (!(col %in% colnames(seurat_obj@meta.data))) {
      next
    }

    is_factor <- is.factor(seurat_obj@meta.data[[col]])

    modify_string <- function(x) {
      x_char <- as.character(x)

      if (nchar(x_char) > remove_n_chars) {
        x_mod <- substr(x_char, remove_n_chars + 1, nchar(x_char))
      } else {
        if (handle_short_strings == "empty") {
          x_mod <- ""
        } else if (handle_short_strings == "original") {
          x_mod <- x_char
        } else {
          x_mod <- ""
        }
      }

      if (case_sensitive) {
        x_mod <- gsub(pattern_to_replace, replacement, x_mod)
      } else {
        x_mod <- gsub(pattern_to_replace, replacement, x_mod, ignore.case = TRUE)
      }
      return(x_mod)
    }

    if (is_factor) {
      original_levels <- levels(seurat_obj@meta.data[[col]])
      modified_levels <- sapply(original_levels, modify_string, USE.NAMES = FALSE)

      seurat_obj@meta.data[[col]] <- factor(
        seurat_obj@meta.data[[col]],
        levels = original_levels,
        labels = modified_levels
      )

    } else {
      seurat_obj@meta.data[[col]] <- sapply(
        seurat_obj@meta.data[[col]],
        modify_string,
        USE.NAMES = FALSE
      )
    }
  }
  return(seurat_obj)
}
#
# cell_counts_entropy_cutoff = list()
#
#
# calculate_silhouette_score <- function(seurat_object, cluster_field) {
#   embedding_data <- Embeddings(seurat_object, reduction = "pca")
#   cluster_labels <- seurat_object@meta.data[[cluster_field]]
#   distance_matrix <- dist(embedding_data)
#   sil <- silhouette(as.numeric(as.factor(cluster_labels)), distance_matrix)
#   avg_sil_width <- mean(sil[, "sil_width"])
#   return(avg_sil_width)
# }
#
# # Initialize the silhouette scores dataframe
# sil_scores_df <- data.frame(
#   parent_cluster = character(),
#   sil_parent = numeric(),
#   sil_root = numeric(),
#   sil_diff = numeric(),
#   sil_parent_celltype = numeric(),   # New column for celltype silhouette
#   sil_root_celltype = numeric(),     # New column for celltype silhouette
#   sil_diff_celltype = numeric(),
#   stringsAsFactors = FALSE
# )

umap_list = list()

# Iterate over each leaf cluster
# for (parent_cluster_name in temp_m$Cluster) {
for (parent_cluster_name in cluster_list) {
  split_parts <- unlist(strsplit(parent_cluster_name, split = "c"))
  parent_cluster_level <- length(split_parts) - 3
  print(parent_cluster_name)

  # child_clusters_names <- dataset_child_clusters$analysis[[parent_cluster_name]]
  child_clusters_names <- unique(recursive_integrated[[parent_cluster_name]]@meta.data[[paste0("level_", parent_cluster_level + 1)]])

  child_clusters_names <- sapply(vector, function(x) gsub("c", "-", substring(child_clusters_names, 4)))


  # Subset the root Seurat object to include only child cluster cells
  parent_cluster_root_subset = recursive_integrated[["cR"]]
  parent_cluster_root_subset = modify_meta_data(parent_cluster_root_subset, paste0("level_", parent_cluster_level + 1))
  current_idents <- Idents(parent_cluster_root_subset)
  Idents(parent_cluster_root_subset) <- paste0("level_", parent_cluster_level + 1)
  parent_cluster_root_subset <- subset(parent_cluster_root_subset, idents = child_clusters_names)
  Idents(parent_cluster_root_subset) <- current_idents
  Idents(parent_cluster_root_subset) <- "seurat_clusters"

  parent_cluster_root_subset <- RunUMAP(parent_cluster_root_subset, reduction = "pca", dims = 1:30, assay = 'integrated', verbose = FALSE)

  parent_cluster = recursive_integrated[[parent_cluster_name]]
  parent_cluster = modify_meta_data(parent_cluster, paste0("level_", parent_cluster_level + 1))
  current_idents <- Idents(parent_cluster)
  Idents(parent_cluster) <- paste0("level_", parent_cluster_level + 1)
  parent_cluster <- subset(parent_cluster, idents = child_clusters_names)
  Idents(parent_cluster) <- current_idents
  Idents(parent_cluster) <- "seurat_clusters"

  parent_cluster <- RunUMAP(parent_cluster, reduction = "pca", dims = 1:30, assay = 'integrated', verbose = FALSE)

  # cell_counts_entropy_cutoff[[parent_cluster_name]] = parent_cluster@assays$integrated@data@Dim[[2]]

  umap_list[[parent_cluster_name]][['recursive']] <- DimPlot(parent_cluster, group.by = paste0("level_", parent_cluster_level + 1), label = TRUE) +
    ggtitle(paste0("Features selected from cluster ", gsub("c", "-", substring(parent_cluster_name, 4)))) +
    theme(plot.title = element_text(
      hjust = 0.5))

  umap_list[[parent_cluster_name]][['root']] <- DimPlot(parent_cluster_root_subset, group.by = paste0("level_", parent_cluster_level + 1), label = TRUE) +
    ggtitle(paste0("Features selected from full dataset")) +
    theme(plot.title = element_text(
      hjust = 0.5))

}

graphics.off()

combined_plots <- list()

label_counter <- 1

for (parent in names(umap_list)) {
  recursive_plot <- umap_list[[parent]][['recursive']]
  root_plot <- umap_list[[parent]][['root']]

  shared_title <- paste0("Recursive cluster ", gsub("c", "-", substring(parent, 4)), " subclusters")

  paired_grob <- arrangeGrob(
    recursive_plot, root_plot,
    ncol = 2,
    top = textGrob(shared_title, gp = gpar(fontsize = 18, fontface = "bold"), hjust = 0.5)
  )

  label <- paste0(' (', letters[label_counter], ')')
  label_grob <- textGrob(
    label,
    x = unit(0, "npc") + unit(5, "points"),
    y = unit(1, "npc") - unit(5, "points"),
    just = c("left", "top"),
    gp = gpar(fontsize = 20, fontface = "bold")
  )

  paired_grob <- gTree(children = gList(
    paired_grob,
    label_grob
  ))

  paired_grob <- gTree(children = gList(
    rectGrob(gp = gpar(col = "#3b3a3a", fill = NA, lwd = 2)),
    paired_grob
  ))

  combined_plots <- c(combined_plots, list(paired_grob))

  label_counter <- label_counter + 1
}

num_plots <- length(combined_plots)
num_rows <- ceiling(num_plots / 2)

grid_layout <- arrangeGrob(
  grobs = combined_plots,
  ncol = 2
)

output_filename <- paste(reference, c_resolution, "UMAPs")

ggsave(
  filename = paste0(output_filename, '.png'),
  plot = grid_layout,
  width = 8.5 * 2,
  height = 11 * 2 / 8 * num_rows * 1.5,
  dpi = 300,
  limitsize = FALSE
)

pp1 <- image_trim(image_read(paste0(output_filename, '.png')))
image_write(pp1, paste0(output_filename, '.pdf'), format = 'pdf')
pp1 <- image_scale(pp1, "50%")
image_write(pp1, paste0(output_filename, '.eps'), format = 'eps')
