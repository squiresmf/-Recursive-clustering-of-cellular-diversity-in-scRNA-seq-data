# This code generates heatmap plots for recursive and seurat equivalent cluster purity results for two annotation levels
# This code will reproduce figure 2, if main has been run with c_resolution = 0.02 and reference_name = "human_PBMC/"
# Load the libraries from main first before running

library(showtext)
library(ggpubr)

font_add_google("PT Sans Narrow", family = "narrow")
showtext_auto()
showtext_opts(dpi = 400)


plot_list = list()

# Code for generating the cluster purity analysis for benchmark analysis
if (!is_helmsley) {
  values = vector(mode = 'list')
  annotation_level_name_vector = vector(mode = 'character')
  level_name_vector = vector(mode = 'character')
  type_vector = vector(mode = 'character')
  values_list = vector(mode = 'list')
  values_df_list = vector(mode = 'list')
  values_vector_list = vector(mode = 'list')
  mean_values_list = vector(mode = 'list')
  for (value_type in c('mean', 't_value', 'w_value')) {
    values_vector_list[[value_type]] = vector(mode = 'numeric')
  }
  for (annotation_level in seq(2)) {
    for (level in 2) {
      if (paste0('level_', level, '_seurat_equivalent') %in% names(recursive_integrated$cR@meta.data)) {
        for (clustering in c('Iterative', 'Seurat Equivalent')) {
          if (clustering == 'Iterative') {
            cluster_level_name = paste0('level_', level)
          } else {
            cluster_level_name = paste0('level_', level, '_seurat_equivalent')
          }
          annotation_title = 'celltype.l'


          dataset = recursive_integrated$cR
          cluster_level = cluster_level_name
          annotation_level2 = paste0(annotation_title, annotation_level)
          title = paste('Level', level, clustering, 'Clustering, Annotation Level', annotation_level2, ref_name, c_resolution, DEG_cutoff)
          subset_cluster_df = data.frame(cluster = as.vector(dataset@meta.data[[cluster_level]]), l1 = dataset@meta.data[[annotation_level2]])

          subset_grouped <- suppressMessages(suppressWarnings(list(l1 = dcast(subset_cluster_df, cluster ~ l1, fun.aggregate = length))))

          subset_correspondence_row_normalized <- list(l1 = NULL)

          for (level2 in names(subset_correspondence_row_normalized)) {
            subset_grouped_matrix = data.matrix(subset_grouped[[level2]][tail(names(subset_grouped[[level2]]), length(names(subset_grouped[[level2]])) - 1)])
            rownames(subset_grouped_matrix) <- as.vector(subset_grouped[[level2]][['cluster']])

            subset_grouped_matrix_row_normalized = subset_grouped_matrix / rowSums(subset_grouped_matrix)

            row_order_names <- vector()
            row = 0
            for (column_idx in seq(ncol(subset_grouped_matrix_row_normalized))) {
              temp_reordered <- subset_grouped_matrix_row_normalized[order(subset_grouped_matrix_row_normalized[, column_idx], decreasing = TRUE),]
              for (row_idx in seq(nrow(temp_reordered))) {
                argmax = length(integer(which.max(temp_reordered[row_idx,])))
                if (argmax == column_idx) {
                  row = row + 1
                  row_order_names[[row]] = rownames(temp_reordered)[[row_idx]]
                }
              }
            }

            subset_grouped_matrix_row_normalized = subset_grouped_matrix_row_normalized[rev(row_order_names),]
            subset_grouped_df_final <- as.data.frame(as.table(subset_grouped_matrix_row_normalized))
            subset_correspondence_row_normalized[[level2]] <- subset_grouped_df_final
            colnames(subset_correspondence_row_normalized[[level2]]) = c("Cluster", "Type", "Proportion")

            cell_names = aggregate(Proportion ~ Type, data = subset_correspondence_row_normalized[[level2]], FUN = max)
            cell_names$Proportion = 0

            out_temp = merge(aggregate(Proportion ~ Cluster, data = subset_correspondence_row_normalized[[level2]], FUN = max), subset_correspondence_row_normalized[[level2]], by = c("Proportion", "Cluster"))
            out = list()
            out_temp2 = aggregate(Proportion ~ Type, data = rbind(cell_names, aggregate(Proportion ~ Type, data = out_temp, FUN = mean)), FUN = max)

            out[['Cell Type Mean Cluster Purity']] = out_temp2$Proportion
            out[['Cluster Purity']] = out_temp$Proportion

            color_values = seq(100) / 100
            palette = colorRampPalette(c('white', 'orange', 'red'))(length(color_values))

            if (clustering == 'Iterative') {
              for (i in seq_along(levels(subset_correspondence_row_normalized$l1$Cluster))) {
                levels(subset_correspondence_row_normalized$l1$Cluster)[i] = sub('cRc', '', levels(subset_correspondence_row_normalized$l1$Cluster)[i])
                levels(subset_correspondence_row_normalized$l1$Cluster)[i] = sub('c', '-', levels(subset_correspondence_row_normalized$l1$Cluster)[i])
                levels(subset_correspondence_row_normalized$l1$Cluster)[i] = sub('c', '-', levels(subset_correspondence_row_normalized$l1$Cluster)[i])
              }
            }

            if (clustering == 'Iterative') {
              title_plot = paste0('Recursive Clustering, Annotation Level ', annotation_level)
            } else {
              title_plot = paste0('Single-Pass Clustering, Annotation Level ', annotation_level)
            }

            plot <- ggplot(subset_correspondence_row_normalized[[level2]], aes(x = Type, y = Cluster, fill = Proportion)) +
              geom_tile() +
              geom_text(aes(label = round(Proportion, 1), color = ifelse(Proportion < 0.05, "light grey", 'black'), family = "narrow", fontface = "bold"), size = 3.5) +
              scale_fill_gradient2(midpoint = 0.5, low = "white", mid = palette[20], high = palette[70], space = "Lab") +
              scale_color_identity() +
              theme(axis.text.x = element_text(family = "narrow", angle = 45, hjust = 1, face = "bold"),
                    axis.text.y = element_text(family = "narrow", face = "bold"), #
                    axis.title.x = element_text(family = "narrow"), #, face = "bold"),
                    axis.title.y = element_text(family = "narrow"), #, face = "bold"),
                    plot.title = element_text(hjust = 0.5),
                    legend.title = element_text(family = "narrow", face = "bold"),
                    legend.text = element_text(family = "narrow", size = 10, face = "bold"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank()
              )
            #ggtitle(title_plot)

            if (annotation_level == 2) {
              plot = plot + theme(axis.title.y = element_blank())
            }


            #if (clustering == 'Iterative' & annotation_level == 1) {
            #  legend = cowplot::get_legend(plot)
            #  legend_name = paste0('cluster purity legend')
            #  ggsave(paste0(legend_name, '.png'), legend, device = "png", dpi = 400)
            #  image_write(image_trim(image_read(paste0(legend_name, '.png'), strip = TRUE)), paste0(legend_name, '.png'), format = 'png')
            #}
            #
            #plot = plot + theme(legend.position = "none")

            if (clustering == 'Seurat Equivalent') {
              plot = plot + theme(plot.margin = unit(c(0.85, 0, 0.2, 0.85), 'lines'))
            } else {
              plot = plot + theme(plot.margin = unit(c(0.85, 0, 0.2, 0.85), 'lines'))
            }

            if (annotation_level == 1) {
              plot = plot + labs(x = "\n\nCell Type")
              if (clustering == 'Seurat Equivalent') {
                plot = plot + labs(y = "Cluster\n")
                suppressMessages(ggsave(filename = paste0(title, '.png'), width = 2.75, height = 6, plot = plot, device = 'png', dpi = 400))
              } else {
                suppressMessages(ggsave(filename = paste0(title, '.png'), width = 2.75, height = 6, plot = plot, device = 'png', dpi = 400))
              }

            } else {
              plot = plot + labs(x = "Cell Type")
              if (clustering == 'Seurat Equivalent') {
                plot = plot + theme(plot.margin = unit(c(0.85, 0, 0.2, 1.5), 'lines'))
                suppressMessages(ggsave(filename = paste0(title, '.png'), width = 8.25, height = 6, plot = plot, device = 'png', dpi = 400))
              } else {
                suppressMessages(ggsave(filename = paste0(title, '.png'), width = 8.25, height = 6, plot = plot, device = 'png', dpi = 400))
              }
            }
          }


          plot_list[[clustering]][[annotation_level]] = plot


          values[[paste0('annotation level ', annotation_level)]][[paste0('level ', level)]][[clustering]] = out
        }
        for (result_type in c('Cell Type Mean Cluster Purity', 'Cluster Purity')) {
          if (result_type == 'Cell Type Mean Cluster Purity') {
            paired = TRUE
          } else {
            paired = FALSE
          }
          print(paste0(paste0('annotation level ', annotation_level, ' level ', level, ' ', result_type)))
          values_list[['mean']] = mean(values[[paste0('annotation level ', annotation_level)]][[paste0('level ', level)]][['Iterative']][[result_type]]) - mean(values[[paste0('annotation level ', annotation_level)]][[paste0('level ', level)]][['Seurat Equivalent']][[result_type]])
          values_list[['t_value']] = t.test(values[[paste0('annotation level ', annotation_level)]][[paste0('level ', level)]][['Iterative']][[result_type]], values[[paste0('annotation level ', annotation_level)]][[paste0('level ', level)]][['Seurat Equivalent']][[result_type]], paired = paired)$p.value
          values_list[['w_value']] = wilcox.test(values[[paste0('annotation level ', annotation_level)]][[paste0('level ', level)]][['Iterative']][[result_type]], values[[paste0('annotation level ', annotation_level)]][[paste0('level ', level)]][['Seurat Equivalent']][[result_type]], paired = paired)$p.value
          values_list[['Iterative mean']] = mean(values[[paste0('annotation level ', annotation_level)]][[paste0('level ', level)]][['Iterative']][[result_type]])
          values_list[['Seurat Equivalent mean']] = mean(values[[paste0('annotation level ', annotation_level)]][[paste0('level ', level)]][['Seurat Equivalent']][[result_type]])
          annotation_level_name_vector = c(annotation_level_name_vector, paste0('annotation level ', annotation_level))
          level_name_vector = c(level_name_vector, paste0('level ', level))
          type_vector = c(type_vector, result_type)
          for (value_type in c('mean', 't_value', 'w_value')) {
            values_vector_list[[value_type]] = c(values_vector_list[[value_type]], values_list[[value_type]])
          }
          for (value_type in c('Iterative mean', 'Seurat Equivalent mean')) {
            mean_values_list[[value_type]][[paste0('annotation level ', annotation_level)]][[result_type]] = values_list[[value_type]]
          }
        }
      }
    }
  }
}

temp = ggarrange(plot_list[['Iterative']][[1]],
                 plot_list[['Iterative']][[2]],
                 plot_list[['Seurat Equivalent']][[1]],
                 plot_list[['Seurat Equivalent']][[2]],
                 widths = c(2.75, 8.25),
                 legend = 'right',
                 common.legend = TRUE
)

#                         PBMC Cluster Purity Analysis\n
temp2 = annotate_figure(temp,
                        top = text_grob("Annotation Level 1                                                                                Annotation Level 2"
                          , size = 14, family = "narrow", face = "bold", hjust = 0.72),

                        left = text_grob("Single-Pass Clustering                                                                                                  Recursive Clustering",
                                         rot = 90, size = 14, family = "narrow", face = "bold", hjust = 0.475)
)

suppressMessages(ggsave(filename = paste0('PBMC example heatmaps.png'), width = 10, height = 11, plot = temp2, device = 'png', dpi = 400))
pp1 <- image_read(paste0('PBMC example heatmaps.png'))
pp1 <- image_scale(pp1, "50%")
image_write(pp1, paste0('PBMC example heatmaps.eps'), format = 'eps')




