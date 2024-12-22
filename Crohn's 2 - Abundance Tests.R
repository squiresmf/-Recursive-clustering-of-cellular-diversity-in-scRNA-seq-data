# Code for performing statistical tests for patient abundance values and related plots for Crohn's scRNA-seq data analysis

for (cluster_type in c('iterative')) {
  for (metric in c('patient_total')) {
    for (level in seq(0, m_level)){
      p_value_tables[[cluster_type]]$patient_total[[paste0('level ', level)]][['combined']] = (p_value_tables[[cluster_type]]$patient_total[[paste0('level ', level)]][['t']] + p_value_tables[[cluster_type]]$patient_total[[paste0('level ', level)]][['wilcox']])/2
    }
  }
}


tests = c('t', 'wilcox', 'combined')
var_names_list = list()
var_names = c('CTRL vs CD', 'CTRL vs TN-CD', 'TN-CD vs CD')
p_value_tables_final = list()
p_value_tables_temp = list()
p_value_table_plots_temp = list()

for (level in seq(m_level)){
  level_name = paste0('level ', level)
  for (test in tests){
    for (level2 in seq(level, 0, -1)){
      level2_name = paste0('level ', level2)
      temp_p_value_table = as.data.frame(p_value_tables$iterative$patient_total[[paste0('level ', level)]][[test]][, var_names])
      for (cluster in rownames(temp_p_value_table)){
        cluster_name_split = unlist(strsplit(cluster, split = 'c'))
        cluster_upper_level_name = paste(cluster_name_split[1:(length(cluster_name_split) - max((length(cluster_name_split) - 2) - (level2+1), 0))], collapse = 'c')
        temp_p_value_table[cluster, var_names] = p_value_tables$iterative$patient_total[[paste0('level ', level2)]][[test]][cluster_upper_level_name, var_names]
      }
      p_value_tables_temp[[level_name]][[test]][[level2_name]] = temp_p_value_table
    }
    for (var_name in var_names){
      temp_p_value_table = temp_p_value_table[, FALSE]
      for (level2_name in rev(names(p_value_tables_temp[[level_name]][[test]]))){
        temp_p_value_table[level2_name] = p_value_tables_temp[[level_name]][[test]][[level2_name]][[var_name]]
      }
      p_value_tables_final[[level_name]][[test]][[var_name]] = temp_p_value_table
      for (alpha in c(0.05)){
        m = temp_p_value_table
        m = melt(t(as.data.frame(m)))
        m = rename(m, Level = Var1, Cluster = Var2)
        m['P.value'] = m['value']
        # p_value_table_plots_temp[[level_name]][[test]][[var_name]][[paste0(alpha)]] =
        #   GetTestHeatmap(mat = m, significance = alpha, x_var_name = 'Cluster', y_var_name = 'Level', fill_name = 'P.value', title_list = list(metric, var_name, paste0(test, ' Test')), text = level <= 1, rounding = 4, angle = 90)
      }
    }
  }
}

print('AbundanceTests')


