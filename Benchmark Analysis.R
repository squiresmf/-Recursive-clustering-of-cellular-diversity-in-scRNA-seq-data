# Code for generating benchmark mean cell type cluster purity metric plots for the four scRNA-seq benchmark datasets
# This code will reproduce figure 3 if main has been run for all four datasets for all 28 resolution settings
# Load the libraries from main first before running


assay = 'integrated'

plot_list = list()

annotation_levels = list()
annotation_levels[['human_PBMC']] = 2
annotation_levels[['human_adipose']] = 2
annotation_levels[['human_tonsil']] = 3
annotation_levels[['human_fetus']] = 2

for (reference_name in c("human_PBMC", "human_adipose", "human_tonsil", "human_fetus")) {
  ref_name = paste0(reference_name, '_', assay, '_')

  results_filename = paste0(ref_name, 'full_results', '.rds')
  if (!file.exists(results_filename)) {
    c_resolution_list = vector(mode = 'numeric')
    DEG_cutoff_list = vector(mode = 'numeric')
    annotation_level_list = vector(mode = 'numeric')
    level_list = vector(mode = 'numeric')
    num_clusters_list = vector(mode = 'numeric')
    result_type_list = vector(mode = 'character')
    values_iterative_list = vector(mode = 'numeric')
    values_seurat_list = vector(mode = 'numeric')
    values_list = vector(mode = 'numeric')
    for (c_resolution in seq(0.015, 0.15, 0.0005)) {
      for (DEG_cutoff in c(0)) {
        equivalent_clustering_filename = paste0("seurat_equivalent_clustering_reference", '_', ref_name)
        equivalent_clustering_filename = paste0(equivalent_clustering_filename, c_resolution, '_', DEG_cutoff, '.rds')
        cluster_purity_filename = paste0('results_', ref_name, c_resolution, '_', DEG_cutoff, '.rds')
        if (file.exists(equivalent_clustering_filename) & file.exists(cluster_purity_filename)) {
          equivalent_clustering = readRDS(equivalent_clustering_filename)
          cluster_purity_results = readRDS(cluster_purity_filename)
          if (c_resolution == 0.015 & reference_name == 'human_PBMC'){
            annotation_1_iterative = mean(cluster_purity_results$`annotation level 1`$`level 2`$Iterative$`Cell Type Mean Cluster Purity`)
            annotation_1_seurat = mean(cluster_purity_results$`annotation level 1`$`level 2`$`Seurat Equivalent`$`Cell Type Mean Cluster Purity`)
            annotation_2_iterative = mean(cluster_purity_results$`annotation level 2`$`level 2`$Iterative$`Cell Type Mean Cluster Purity`)
            annotation_2_seurat = mean(cluster_purity_results$`annotation level 2`$`level 2`$`Seurat Equivalent`$`Cell Type Mean Cluster Purity`)
            print(annotation_1_seurat/annotation_1_iterative - 1)
            print(annotation_2_iterative/annotation_2_seurat - 1)
          }
          for (annotation_level in seq(annotation_levels[[reference_name]])) {
            for (level in seq(3)) {
              num_clusters = length(unique(equivalent_clustering[[paste0(level)]]))
              for (result_type in c('Cell Type Mean Cluster Purity', 'Cluster Purity')) {
                value_iterative = mean(cluster_purity_results[[paste0('annotation level ', annotation_level)]][[paste0('level ', level)]][['Iterative']][[result_type]])
                value_seurat = mean(cluster_purity_results[[paste0('annotation level ', annotation_level)]][[paste0('level ', level)]][['Seurat Equivalent']][[result_type]])
                value = value_iterative - value_seurat
                c_resolution_list = c(c_resolution_list, c_resolution)
                DEG_cutoff_list = c(DEG_cutoff_list, DEG_cutoff)
                annotation_level_list = c(annotation_level_list, annotation_level)
                level_list = c(level_list, level)
                num_clusters_list = c(num_clusters_list, num_clusters)
                result_type_list = c(result_type_list, result_type)
                values_iterative_list = c(values_iterative_list, value_iterative)
                values_seurat_list = c(values_seurat_list, value_seurat)
                values_list = c(values_list, value)
                print(paste('', annotation_level, level, c_resolution, DEG_cutoff, num_clusters, value, sep = " "))
              }
            }
          }
        }
      }
    }
    for (level in seq(3)) {
      print(length(unique(equivalent_clustering[[paste0(level)]])))
    }
    results_df = cbind.data.frame(c_resolution_list, DEG_cutoff_list, annotation_level_list, level_list, num_clusters_list, result_type_list, values_iterative_list, values_seurat_list, values_list)
    colnames(results_df) = cbind("c_resolution", "DEG_cutoff", "annotation_level", "level", "num_clusters", "result_type", "values_iterative", "values_seurat", "values")
    saveRDS(results_df, results_filename)

  }
}

results_file_list = list()
assay = 'integrated'
reference_names = c("human_PBMC", "human_adipose", "human_tonsil", "human_fetus")
reference_names_pretty = c("PBMC", "Adipose", "Tonsil", "Fetus")
reference_names_cell_types = list()

for (analysis_type in c('no DEG cutoff')) {
  for (annotation_level in c(2)){
    for (dataset_number in seq(length(reference_names))) {
      reference_name = reference_names[[dataset_number]]
      reference_name_pretty = reference_names_pretty[[dataset_number]]
      results_filename = paste0(reference_name, '_', assay, '_full_results.rds')
      results_file_list[[reference_name]] = readRDS(results_filename)
      results_file_list[[reference_name]]$reference = reference_name_pretty
      annotation_level_final = annotation_levels[[reference_name]] - (2 - annotation_level)

      results_file_list[[reference_name]] = results_file_list[[reference_name]][(results_file_list[[reference_name]]$annotation_level <= annotation_level_final &
        results_file_list[[reference_name]]$annotation_level >= (annotation_level_final - 1) &
        results_file_list[[reference_name]]$result_type == 'Cell Type Mean Cluster Purity' &
        results_file_list[[reference_name]]$DEG_cutoff == 0 &
        results_file_list[[reference_name]]$level == 2),]
        results_file_list[[reference_name]]$c_resolution = round(results_file_list[[reference_name]]$c_resolution, 10)
        results_file_list[[reference_name]] = results_file_list[[reference_name]][results_file_list[[reference_name]]$c_resolution %in% seq(0.015, 0.15, 0.005),]

      results_file_list[[reference_name]]$annotation_level = results_file_list[[reference_name]]$annotation_level + (2 - annotation_level_final)

      names(results_file_list[[reference_name]]) = c("Resolution", "DEG Cutoff", "Annotation Level", "Cluster Tree Depth", "Number of Clusters", "Metric", "Recursive", "Single-Pass", "Difference", "Reference")
    }

    base_colors = c("#669bcc", "#ccb166", "#CC6677")

    results_subset_df = bind_rows(results_file_list)

    results_subset_df = melt(results_subset_df, id.vars = c("Resolution", "DEG Cutoff", "Annotation Level", "Cluster Tree Depth", "Number of Clusters", "Metric", "Reference"), measure.vars = c("Recursive", "Single-Pass", "Difference"))
    names(results_subset_df) = c("Resolution", "DEG Cutoff", "Annotation Level", "Cluster Tree Depth", "Number of Clusters", "Metric", "Reference", "Method", "Value")

    results_subset_df["Metric"][results_subset_df["Metric"] == "Cell Type Mean Cluster Purity"] <- "Cell Type Mean\nCluster Purity"
    results_subset_df$Reference = factor(results_subset_df$Reference, levels = reference_names_pretty)

    results_subset_df = results_subset_df[results_subset_df$Method != 'Difference', ]
    results_subset_df$`Annotation Level` = paste0('Annotation Level ', results_subset_df$`Annotation Level`, '\n')

    aes_string2 <- function(...) {
      args <- lapply(list(...), function(x) sprintf("`%s`", x))
      do.call(aes_string, args)
    }

    scale = 0.925
    results_plot = ggplot(results_subset_df, aes_string2(x = "Number of Clusters", y = 'Value')) +
      geom_point(aes_string(color = 'Method'), alpha = 0.3, size = 2 * scale) +
      scale_size_area() +
      scale_y_continuous(limits = c(0, 1), breaks = c(0, .30, .60, .90), expand = c(0.001, 0)) +
      scale_color_manual(values = c("#669bcc", "#ccb166", "#CC6677")) +
      geom_hline(yintercept = 0, color = "dark grey", linetype = 1, size = 0.5) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "dark grey")) +
      facet_grid(reformulate('Reference', 'Annotation Level'), switch = "y", scales = "free") +
      theme(strip.background = element_rect(colour = "white", fill = "white"), strip.text.y.left = element_text(angle = 90), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12), strip.placement = "outside",
            panel.spacing = unit(10, "pt"), plot.title = element_text(hjust = 0.5), axis.title.y = element_text(vjust = -15, angle = 90), axis.text.x = element_text(size = 9),
            axis.text.y = element_text(size = 10)) +
      labs(y = "Mean Cell Type Cluster Purity")

    print(results_plot)
    ggsave(paste0(analysis_type, annotation_level, '.pdf'), results_plot, width = 8, height = 6, device = 'pdf')
    ggsave(paste0(analysis_type, annotation_level, '.png'), results_plot, width = 8*scale, height = 6*scale, device = 'png', dpi = 400)
  }
}

pp1 <- image_trim(image_read(paste0(analysis_type, annotation_level, '.png')))
image_write(pp1, paste0(analysis_type, annotation_level, '.eps'), format = 'eps')
