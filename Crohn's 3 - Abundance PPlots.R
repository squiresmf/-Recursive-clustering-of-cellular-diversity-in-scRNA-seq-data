# Code for creating heirarchically structured heat map plots for patient abundance Crohn's scRNA-seq data analysis
# The Crohn's dataset is not currently publicly available

# Abundance
annotation_colors = list()
annotation_colors[['CTRL vs TN-CD']][['cRc7c0']] = "#ff00a0" # higher vs parents

annotation_colors[['CTRL vs CD']][['cRc3c0']] = "#55ff00" # higher vs children
# annotation_colors[['CTRL vs CD']][['cRc0c1c1']] = "#ff00a0" # higher vs parents
# annotation_colors[['CTRL vs CD']][['cRc4c3c1']] = "#ff00a0" # higher vs parents
annotation_colors[['CTRL vs CD']][['cRc3']] = "#ff8c00" # consistent
annotation_colors[['CTRL vs CD']][['cRc6']] = "#ff8c00" # consistent

annotation_colors[['TN-CD vs CD']][['cRc3c0']] = "#55ff00" # higher vs children
# annotation_colors[['TN-CD vs CD']][['cRc0c1']] = "#ff00a0" # higher vs parents
annotation_colors[['TN-CD vs CD']][['cRc3']] = "#ff8c00" # consistent
annotation_colors[['TN-CD vs CD']][['cRc8']] = "#ff8c00" # consistent


p_value_table_plots_final = list()
abundance_tree_plot_list = list()

#for (level in seq(m_level)){
for (level in 2){
  level_name = paste0('level ', level)
  # for (test in c('t', 'wilcox', 'combined')){
  for (test in 'wilcox'){
    temp_p_value_table_melt = setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("Var1", "Var2", "value", "comparison"))
    for (i in seq_along(var_names)){
      temp_p_value_table = melt(t(p_value_tables_final[[level_name]][[test]][[var_names[i]]]))
      temp_p_value_table['comparison'] = var_names[i]
      temp_p_value_table_melt = rbind(temp_p_value_table_melt, temp_p_value_table)
    }
    #for (alpha in c(0.05, 0.01)){
    for (alpha in 0.05){
      m = temp_p_value_table_melt
      m = rename(m, Level = Var1, Cluster = Var2)
      m['P.value'] = m['value']

      temp_m = m[m[['comparison']] %in% c('CTRL vs CD', 'CTRL vs TN-CD', 'TN-CD vs CD'), ]
      temp_m_heatmap = temp_m
      if (use_annotations){
        temp_m_heatmap = GetAnnotatedClusterNamesDF(df = temp_m, recursive_cluster_list = recursive_cluster_list, cluster_annotation_list = cluster_annotation_list, cluster_type = 'iterative')
      }
      temp_m['Level'] = substring(unlist(temp_m['Level']), 7, 7)
      # p_value_table_plots_final[[level_name]][[test]][[paste0(alpha)]] =
      #   GetTestHeatmap(mat = temp_m_heatmap, significance = alpha, x_var_name = 'Cluster', y_var_name = 'Level', fill_name = 'P.value', title_list = list(metric, 'Difference in Abundance', paste0(test, ' Test')), text = level <= 1, wrap_var_name = 'comparison', rounding = 4) +
      #     scale_x_discrete(limits=rev) +
      #     facet_wrap(~comparison, scales = 'fixed', ncol = 6) +
      #     coord_flip() +
      #     scale_y_discrete() +
      #     theme(axis.text.x = element_text(angle = 0, vjust = 0.35, size = 14), axis.text.y = element_text(size = 12, family = "mono"), axis.title=element_text(size=14,face="bold"), title = element_text(size=15))
      #print(p_value_table_plots_final[[level_name]][[test]][[paste0(alpha)]])

      GetRadialTreeFromDataframe(df = temp_m, tree_comparisons = c('CTRL vs TN-CD', 'CTRL vs CD', 'TN-CD vs CD'), val_name = 'P.value', tree_level = level, gravities = c('NorthEast', 'NorthEast', 'NorthEast'), annotation_colors = annotation_colors)

    }
  }
}
print('AbundancePPlots')

ggsave(paste0('Abundance Level 0', '.png'), abundance_box_plots$iterative$patient_total$`level 0`$unnormalized +
  theme(legend.text = element_text(size = 16), axis.title.y = element_text(angle = 0, hjust = -1)) +
  ggtitle('Patient Level Abundance, Level 0'), device = "png", width = 14, height = 9, dpi = 400)
ggsave(paste0('Abundance Level 0 Heatmap', '.png'), p_value_table_plots$iterative$patient_total$`level 0`$wilcox$`0.05`$unadjusted +
  ggtitle('Wilcox Test of Phenotypic Difference\nIn Patient Level Abundance, Level 0'), device = "png", width = 14, height = 9, dpi = 400)

