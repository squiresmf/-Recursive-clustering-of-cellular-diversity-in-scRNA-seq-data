# Code for creating heirarchically structured heat map plots for DEG Crohn's scRNA-seq data analysis
# The Crohn's dataset is not currently publicly available


# DEG Analysis
annotation_colors = list()
# annotation_colors[["control vs treatment naive-CD"]][['cRc0c1c1']] = "#ff00a0" # higher vs parents
annotation_colors[["control vs treatment naive-CD"]][['cRc2c1c0']] = "#ff00a0" # higher vs parents
annotation_colors[["control vs treatment naive-CD"]][['cRc7c0']] = "#ff00a0" # higher vs parents
annotation_colors[["control vs treatment naive-CD"]][['cRc1']] = "#ff8c00" # consistent
annotation_colors[["control vs treatment naive-CD"]][['cRc2']] = "#ff8c00" # consistent

annotation_colors[["control vs CD"]][['cRc6']] = "#55ff00" # higher vs children
annotation_colors[["control vs CD"]][['cRc0c1c1']] = "#ff00a0" # higher vs parents
# annotation_colors[["control vs CD"]][['cRc0c1c1']] = "#ff00a0" # higher vs parents

# annotation_colors[["treatment naive-CD vs CD"]][['cRc4c3']] = "#55ff00" # higher vs children
annotation_colors[["treatment naive-CD vs CD"]][['cRc1']] = "#55ff00" # higher vs children
# annotation_colors[["treatment naive-CD vs CD"]][['cRc2c2']] = "#ff00a0" # higher vs parents


phenotype_test_DEGs_list = list()
phenotype_test_DEGs_unfiltered = list()
n_trials = 300

cluster_comparison_gene_count = readRDS(seed_cluster_comparison_gene_count_filename)


for (seed_num in seq(0, n_trials)) {
  seed_cluster_comparison_gene_count_filename = paste0(cluster_comparison_gene_count_filename, seed_num, '.rds')
  if (seed_num == 0) {
    phenotype_test_DEGs_unfiltered = readRDS(seed_cluster_comparison_gene_count_filename)
  } else {
    phenotype_test_DEGs_list[[seed_num]] <- readRDS(seed_cluster_comparison_gene_count_filename)
  }
}

phenotype_test_DEGs_distribution = list()
for (cluster_name in names(phenotype_test_DEGs_list[[seed_num]])) {
  print(cluster_name)
  for (comparison_name in names(phenotype_test_DEGs_list[[seed_num]][[cluster_name]])) {
    for (seed_num in seq(n_trials)) {
      phenotype_test_DEGs_distribution[[cluster_name]][[comparison_name]] = c(phenotype_test_DEGs_distribution[[cluster_name]][[comparison_name]],
                                                                              phenotype_test_DEGs_list[[seed_num]][[cluster_name]][[comparison_name]])
    }
  }
}

phenotype_test_DEGs_df_lists = list()
phenotype_test_DEGs_df_lists[['Cluster']] = list()
phenotype_test_DEGs_df_lists[['DEGs']] = list()
phenotype_test_DEGs_df_lists[['DEGs Filtered']] = list()
phenotype_test_DEGs_df_lists[['comparison']] = list()
phenotype_test_DEGs_df_lists[['P.value']] = list()

for (cluster_name in names(phenotype_test_DEGs_unfiltered)) {
  print(cluster_name)
  for (comparison_name in names(phenotype_test_DEGs_unfiltered[[cluster_name]])) {
    true_count = phenotype_test_DEGs_unfiltered[[cluster_name]][[comparison_name]]
    cdf_values = unlist(phenotype_test_DEGs_distribution[[cluster_name]][[comparison_name]])
    empirical_p = 1 - ecdf(cdf_values)(true_count-100)
    if (empirical_p <= 0.05) {
      phenotype_test_DEGs_df_lists[['DEGs Filtered']] = c(phenotype_test_DEGs_df_lists[['DEGs Filtered']], (true_count/mean(cdf_values)))
      # phenotype_test_DEGs_df_lists[['DEGs Filtered']] = c(phenotype_test_DEGs_df_lists[['DEGs Filtered']], true_count)
    } else {
      phenotype_test_DEGs_df_lists[['DEGs Filtered']] = c(phenotype_test_DEGs_df_lists[['DEGs Filtered']], 0)
    }
    phenotype_test_DEGs_df_lists[['P.value']] = c(phenotype_test_DEGs_df_lists[['P.value']], empirical_p)
    phenotype_test_DEGs_df_lists[['comparison']] = c(phenotype_test_DEGs_df_lists[['comparison']], comparison_name)
    phenotype_test_DEGs_df_lists[['Cluster']] = c(phenotype_test_DEGs_df_lists[['Cluster']], cluster_name)
    phenotype_test_DEGs_df_lists[['DEGs']] = c(phenotype_test_DEGs_df_lists[['DEGs']], true_count)
  }
}


phenotype_test_DEGs_df = as.data.frame(do.call(cbind, phenotype_test_DEGs_df_lists))
GetRadialTreeFromDataframe(df = phenotype_test_DEGs_df, tree_comparisons = c("control vs treatment naive-CD", "control vs CD", "treatment naive-CD vs CD"), val_name = 'P.value', tree_level = 2, gravities = c('NorthEast', 'NorthEast', 'NorthEast'), title_extra = '_empirical')
GetRadialTreeFromDataframe(df = phenotype_test_DEGs_df, tree_comparisons = c("control vs treatment naive-CD", "control vs CD", "treatment naive-CD vs CD"), val_name = 'DEGs', tree_level = 2, gravities = c('NorthEast', 'NorthEast', 'NorthEast'), title_extra = '_empirical_nofilter')
phenotype_test_DEGs_df = phenotype_test_DEGs_df[c('Cluster', 'DEGs Filtered', 'comparison')]
colnames(phenotype_test_DEGs_df) = c('Cluster', 'DEGs', 'comparison')
GetRadialTreeFromDataframe(df = phenotype_test_DEGs_df, tree_comparisons = c("control vs treatment naive-CD", "control vs CD", "treatment naive-CD vs CD"), val_name = 'DEGs', tree_level = 2, gravities = c('NorthEast', 'NorthEast', 'NorthEast'), title_extra = '_empirical', annotation_colors = annotation_colors)

print('DEG Analysis Tree Map')

