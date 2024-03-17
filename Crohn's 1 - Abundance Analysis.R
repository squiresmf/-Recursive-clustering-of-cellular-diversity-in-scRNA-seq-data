# Code for calculating patient abundance values and related plots for Crohn's scRNA-seq data analysis
# The Crohn's dataset is not currently publicly available

patient_names = names(paths)
num_cells_patient = list()
for (patient_name in patient_names) {
  num_cells_patient[[patient_name]] = sum(integrated_data@meta.data$patient == patient_name)
}

pairs = list(list("control", "CD"), list("control", "treatment naïve-CD"), list("CD", "treatment naïve-CD"))
#pairs = list(list("control", "CD"), list("control", "treatment naïve-CD"), list("CD", "treatment naïve-CD"),
#             list("control", c("CD", "treatment naïve-CD")), list("CD", c("control", "treatment naïve-CD")), list("treatment naïve-CD", c("control", "CD")))
var_names = c('CTRL vs CD', 'CTRL vs TN-CD', 'TN-CD vs CD')
#var_names = c('CTRL vs CD', 'CTRL vs TN-CD', 'TN-CD vs CD', 'CTRL vs CD+TN-CD', 'CD vs CTRL+TN-CD', 'TN-CD vs CTRL+CD')

tests = c('t', 'wilcox')

patient_ids = names(paths)
disease_phenotypes = patient_metadata$'Disease phenotype'

for (dataset_name in names(recursive_integrated)) {
  print(dataset_name)
  cluster_nums = as.numeric(levels(recursive_integrated[[dataset_name]]@meta.data$seurat_clusters))

  recursive_integrated[[dataset_name]]@misc[['patient_fracs']] = list()
  recursive_integrated[[dataset_name]]@misc[['total_patient_fracs']] = list()

  for (cluster in cluster_nums) {
    cluster_idxs = seq(recursive_integrated[[dataset_name]]$RNA@
                         counts@
                         Dim[[2]])[recursive_integrated[[dataset_name]]@meta.data$seurat_clusters == cluster]

    cluster_patient_fracs = list()
    cluster_total_patient_fracs = list()
    cluster_patients = recursive_integrated[[dataset_name]]@meta.data$patient[cluster_idxs]

    for (patient_name in patient_names) {
      cluster_patient_fracs[[patient_name]] = length(cluster_patients[cluster_patients == patient_name]) / length(cluster_patients)
      cluster_total_patient_fracs[[patient_name]] = length(cluster_patients[cluster_patients == patient_name]) / num_cells_patient[[patient_name]]
    }

    cluster_name = paste0('c', cluster)
    recursive_integrated[[dataset_name]]@misc[['patient_fracs']][[cluster_name]] = unlist(cluster_patient_fracs)
    recursive_integrated[[dataset_name]]@misc[['total_patient_fracs']][[cluster_name]] = unlist(cluster_total_patient_fracs)

  }
}


fracs_by_level = list()
fracs_by_level_df = list()
abundance_box_plots = list()
p_value_tables = list()
p_value_table_plots = list()
is_significant = list()

#for (cluster_type in c('iterative', 'seurat_equivalent')) {
for (cluster_type in c('iterative')) {
  for (metric in c('patient_total')) {
    for (level in seq(0, m_level)) {
      level_name = paste0('level ', level)
      if (cluster_type == 'iterative') {
        sub_cluster_names = dataset_level_clusters[['analysis']][['cR']][[paste0(level)]]
      } else if (cluster_type == 'seurat_equivalent') {
        sub_cluster_names = names(seurat_equivalent_cluster_total_patient_fracs_by_level[[paste0(level)]])
      }
      for (sub_cluster_name in sub_cluster_names) {
        sub_cluster_name_split = unlist(strsplit(sub_cluster_name, split = 'c'))
        sub_cluster_parent_name = paste(sub_cluster_name_split[1:(length(sub_cluster_name_split) - 1)], collapse = 'c')
        sub_cluster_num = paste0('c', as.numeric(sub_cluster_name_split[[length(sub_cluster_name_split)]]))
        if (cluster_type == 'iterative') {
          if (metric == 'patient_total') {
            fracs_by_level[[cluster_type]][[metric]][[level_name]][[sub_cluster_name]] = recursive_integrated[[sub_cluster_parent_name]]@misc$total_patient_fracs[[sub_cluster_num]]
          } else if (metric == 'patient_cluster') {
            fracs_by_level[[cluster_type]][[metric]][[level_name]][[sub_cluster_name]] = recursive_integrated[[sub_cluster_parent_name]]@misc$patient_fracs[[sub_cluster_num]]
          }
          fracs_by_level_df[[cluster_type]][[metric]][[level_name]] = as.data.frame(do.call(cbind, fracs_by_level[[cluster_type]][[metric]][[level_name]]))
        } else if (cluster_type == 'seurat_equivalent') {
          fracs_by_level[[cluster_type]][[metric]][[paste0(level)]][[sub_cluster_name]] = seurat_equivalent_cluster_total_patient_fracs_by_level[[paste0(level)]][[sub_cluster_name]]
          fracs_by_level_df[[cluster_type]][[metric]][[level_name]] = as.data.frame(do.call(cbind, fracs_by_level[[cluster_type]][[metric]][[paste0(level)]]))
        }
      }

      fracs_by_level_df[[cluster_type]][[metric]][[level_name]][['Phenotype']] = disease_phenotypes
      temp_fracs_melted = melt(fracs_by_level_df[[cluster_type]][[metric]][[level_name]], id = "Phenotype")
      temp_fracs_melted = rename(temp_fracs_melted, Cluster = variable, Abundance = value)

      boxplot_temp_fracs_melted = temp_fracs_melted
      for (normalization in c('unnormalized', 'normalized')) {
        if (normalization == 'normalized') {
          for (boxplot_sub_cluster in sub_cluster_names) {
            boxplot_temp_fracs_melted$Abundance[boxplot_temp_fracs_melted$Cluster == boxplot_sub_cluster] = boxplot_temp_fracs_melted$Abundance[boxplot_temp_fracs_melted$Cluster == boxplot_sub_cluster] / max(boxplot_temp_fracs_melted$Abundance[boxplot_temp_fracs_melted$Cluster == boxplot_sub_cluster])
          }
        }
        if (use_annotations) {
          boxplot_df = GetAnnotatedClusterNamesDF(df = boxplot_temp_fracs_melted, recursive_cluster_list = recursive_cluster_list, cluster_annotation_list = cluster_annotation_list, cluster_type = cluster_type)
        } else {
          boxplot_df = boxplot_temp_fracs_melted
        }
        boxplot_df$Phenotype <- factor(boxplot_df$Phenotype, levels = c("CD", "treatment naïve-CD", "control"))
        plot = ggplot(boxplot_df, aes(x = Cluster, y = Abundance, color = Phenotype)) +
          geom_boxplot(outlier.shape = NA) +
          geom_point(position = position_jitterdodge()) +
          theme_bw() +
          scale_color_manual(values = darken(c("#CC6677", "#DDCC77", "#88CCEE"), amount = 0.3), guide = guide_legend(reverse = TRUE)) +
          ggtitle(str_to_title(str_replace(paste0(metric, ', ', level_name), pattern = '_', replacement = ' '))) +
          theme(axis.text.x = element_text(angle = 0, vjust = 0.35, hjust = 0, size = 14), axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), title = element_text(size = 15)) +
          # theme(axis.text.x = element_text(angle = 0, vjust = 0.35, hjust = 0, size = 14), axis.text.y = element_text(size = 12, family = "mono"), axis.title=element_text(size=14,face="bold"), title = element_text(size=15)) +
          scale_x_discrete(limits = rev) +
          coord_flip() +
          theme(panel.grid = element_blank()) +
          geom_vline(xintercept = seq(0.5, length(unique(boxplot_df$Cluster)), by = 1), color = "black", size = 0.75, alpha = .5)
        if (normalization == 'normalized') {
          plot = plot + ylab('Normalized Abundance')
        }
        abundance_box_plots[[cluster_type]][[metric]][[level_name]][[normalization]] = plot
      }
      temp_combined_m = NULL
      for (test in tests) {
        print(paste0(cluster_type, ' clustering ', metric, ' level ', level, ' ', test, ' test'))
        temp_p_value_table = matrix(, nrow = length(sub_cluster_names), ncol = length(var_names), dimnames = list(sub_cluster_names, var_names))
        for (i in seq(length(sub_cluster_names))) {
          for (j in seq(length(var_names))) {
            temp_fracs_melted_for_test = temp_fracs_melted
            x = temp_fracs_melted_for_test$Abundance[temp_fracs_melted_for_test$Phenotype %in% pairs[[j]][[1]] & temp_fracs_melted_for_test$Cluster == sub_cluster_names[[i]]]
            y = temp_fracs_melted_for_test$Abundance[temp_fracs_melted_for_test$Phenotype %in% pairs[[j]][[2]] & temp_fracs_melted_for_test$Cluster == sub_cluster_names[[i]]]
            if (sum(x) + sum(y) == 0) {
              p = 1
            } else {
              if (min(length(x), length(y)) == 0) {
                p = 1
              } else {
                if (test == 't') {
                  p = t.test(x, y)$p.value
                } else if (test == 'wilcox') {
                  p = suppressWarnings(wilcox.test(x, y)$p.value)
                }
              }
            }
            if (is.nan(p)) {
              stop()
            }
            temp_p_value_table[sub_cluster_names[[i]], var_names[[j]]] = p
          }
        }
        p_value_tables[[cluster_type]][[metric]][[level_name]][[test]] = temp_p_value_table


        for (alpha in c(0.05, 0.01)) {

          m = temp_p_value_table

          m = melt(t(as.data.frame(m)))
          m = rename(m, Comparison = Var1, Cluster = Var2)
          m['P.value'] = m['value']
          m['significant'] = as.vector(ifelse(m['P.value'] < alpha, 1, 0))
          m['significant'][is.na(m['significant'])] = 0
          if (adjusting == 'unadjusted') {
            if (is.null(temp_combined_m)) {
              temp_combined_m = m[c('Comparison', 'Cluster')]
            }
            temp_combined_m[[test]] = m[['P.value']]
          }
          is_significant[[cluster_type]][[metric]][[level_name]][[test]][[paste0(alpha)]][[adjusting]] = m[c('Comparison', 'Cluster', 'significant')]
          if (use_annotations) {
            m = GetAnnotatedClusterNamesDF(df = m, recursive_cluster_list = recursive_cluster_list, cluster_annotation_list = cluster_annotation_list, cluster_type = cluster_type)
          }
          p_value_table_plots[[cluster_type]][[metric]][[level_name]][[test]][[paste0(alpha)]][[adjusting]] =
            GetTestHeatmap(mat = m, significance = alpha, x_var_name = 'Cluster', y_var_name = 'Comparison', fill_name = 'P.value', title_list = list(metric, level_name, paste0(test, ' Test')), text = level <= 1, rounding = 4) +
              coord_flip() +
              scale_x_discrete(limits = rev) +
              scale_y_discrete() +
              theme(axis.text.x = element_text(angle = 0, vjust = 0.35, size = 14), axis.text.y = element_text(size = 12, family = "mono"), axis.title = element_text(size = 14, face = "bold"), title = element_text(size = 15))
        }
      }
    }
    temp_combined_m[['P.value']] = rowMeans(subset(temp_combined_m, select = tests))
    temp_p_value_column_combined = NULL
    for (test in tests) {
      temp_p_value_column = as.list(temp_combined_m[[test]])
      for (i in seq_along(temp_p_value_column)) {
        temp_p_value_column[[i]] = PadRight(value = temp_p_value_column[[i]], width = 6, side = "right", pad = "0")
        if (!is.null(temp_p_value_column_combined)) {
          temp_p_value_column_combined[[i]] = paste0(temp_p_value_column_combined[[i]], ',\n', temp_p_value_column[[i]], ' ')
        }
      }
      if (is.null(temp_p_value_column_combined)) {
        temp_p_value_column_combined = temp_p_value_column
      }
    }
    temp_combined_m[['value']] = unlist(temp_p_value_column_combined)
    temp_combined_m[['True P']] = temp_combined_m$P.value
    p_value_table_plots[[cluster_type]][[metric]][[level_name]][['combined']] = list()
    for (alpha in c(0.05, 0.01)) {
      temp_p_value_column_color = as.list(temp_combined_m[['True P']])
      for (test in tests) {
        temp_p_value_column = as.list(temp_combined_m[[test]])
        for (i in seq_along(temp_p_value_column)) {
          if (is.nan(temp_p_value_column[[i]])) {
            temp_p_value_column_color[[i]] = 1
          } else if (temp_p_value_column[[i]] > alpha) {
            temp_p_value_column_color[[i]] = temp_p_value_column_color[[i]] + alpha
          }
        }
      }
      temp_combined_m$P.value = unlist(temp_p_value_column_color)
      temp_combined_m['significant'] = as.vector(ifelse(temp_combined_m['P.value'] < alpha, 1, 0))
      is_significant[[cluster_type]][[metric]][[level_name]][['combined']][[paste0(alpha)]] = temp_combined_m[c('Comparison', 'Cluster', 'significant')]
      temp_combined_m_for_plot = temp_combined_m[temp_combined_m[['Comparison']] %in% c('CTRL vs CD', 'CTRL vs TN-CD', 'TN-CD vs CD'),]
      if (use_annotations) {
        temp_combined_m_for_plot = GetAnnotatedClusterNamesDF(df = temp_combined_m_for_plot, recursive_cluster_list = recursive_cluster_list, cluster_annotation_list = cluster_annotation_list, cluster_type = cluster_type)
      }
      temp_combined_m_for_plot$Comparison = factor(temp_combined_m_for_plot$Comparison, ordered = TRUE, levels = c('CTRL vs TN-CD', 'CTRL vs CD', 'TN-CD vs CD'))
      p_value_table_plots[[cluster_type]][[metric]][[level_name]][['combined']][[paste0(alpha)]] =
        GetTestHeatmap(mat = temp_combined_m_for_plot, significance = alpha, x_var_name = 'Cluster', y_var_name = 'Comparison', fill_name = 'P.value', title_list = list(metric, level_name, paste0(paste(tests, collapse = ' Test, '), ' Test')), text = level <= 1) +
          coord_flip() +
          scale_x_discrete(limits = rev) +
          scale_y_discrete() +
          theme(axis.text.x = element_text(angle = 0, vjust = 0.35, size = 14), axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), title = element_text(size = 15))
    }
  }
}


print('AbundanceAnalysis')