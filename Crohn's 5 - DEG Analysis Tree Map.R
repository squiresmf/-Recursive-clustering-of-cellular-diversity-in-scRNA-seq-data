# Code for creating heirarchically structured heat map plots for DEG Crohn's scRNA-seq data analysis
# Also produces tables for supplementary section
library(dplyr)
library(gridExtra)
library(ggplot2)

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


phenotype_test_DEGs_distribution <- list()
# Aggregate DEGs distribution across all trials
for (cluster_name in names(phenotype_test_DEGs_list[[1]])) {
  print(cluster_name)
  comparisons <- names(phenotype_test_DEGs_list[[1]][[cluster_name]])
  for (comparison_name in comparisons) {
    # Iterate over each trial to collect DEGs
    for (trial_num in seq_len(n_trials)) {
      # Append the DEGs from the current trial to the distribution list
      current_deg <- phenotype_test_DEGs_list[[trial_num]][[cluster_name]][[comparison_name]]
      phenotype_test_DEGs_distribution[[cluster_name]][[comparison_name]] <- c(
        phenotype_test_DEGs_distribution[[cluster_name]][[comparison_name]],
        current_deg
      )
    }
  }
}

# Initialize lists to store data frame columns
phenotype_test_DEGs_df_lists <- list()

degs_cutoff <- quantile(unlist(phenotype_test_DEGs_unfiltered, recursive = TRUE), probs = 0.8, na.rm = TRUE)[["80%"]]
# Loop over each cluster in the unfiltered DEGs
for (cluster_name in names(phenotype_test_DEGs_unfiltered)) {
  print(cluster_name)
  comparisons <- names(phenotype_test_DEGs_unfiltered[[cluster_name]])
  for (comparison_name in comparisons) {
    # Retrieve the true count of DEGs for the current cluster and comparison
    true_count <- phenotype_test_DEGs_unfiltered[[cluster_name]][[comparison_name]]
    # Retrieve the distribution of DEGs from all trials
    cdf_values <- phenotype_test_DEGs_distribution[[cluster_name]][[comparison_name]]
    # Calculate the empirical p-value
    empirical_p <- 1 - ecdf(cdf_values)(true_count)
    null_distribution_mean_DEGs = mean(cdf_values)
    null_distribution_max_DEGs = max(cdf_values)
    # Apply filtering criteria
    if (empirical_p < 0.05 && true_count >= degs_cutoff) {
      degs_filtered <- true_count / null_distribution_mean_DEGs
    } else {
      degs_filtered <- 0
    }
    phenotype_test_DEGs_df_lists[['DEGs Filtered']] <- c(phenotype_test_DEGs_df_lists[['DEGs Filtered']], degs_filtered)
    phenotype_test_DEGs_df_lists[['Null Distribution Mean DEG Count']] <- c(phenotype_test_DEGs_df_lists[['Null Distribution Mean DEG Count']], null_distribution_mean_DEGs)
    phenotype_test_DEGs_df_lists[['P value']] <- c(phenotype_test_DEGs_df_lists[['P value']], empirical_p)
    phenotype_test_DEGs_df_lists[['comparison']] <- c(phenotype_test_DEGs_df_lists[['comparison']], comparison_name)
    phenotype_test_DEGs_df_lists[['Cluster']] <- c(phenotype_test_DEGs_df_lists[['Cluster']], cluster_name)
    phenotype_test_DEGs_df_lists[['DEG Count']] <- c(phenotype_test_DEGs_df_lists[['DEG Count']], true_count)
  }
}

phenotype_test_DEGs_df <- data.frame(
  Cluster = phenotype_test_DEGs_df_lists[['Cluster']],
  `DEG Count` = phenotype_test_DEGs_df_lists[['DEG Count']],
  `Null Distribution Mean DEG Count` = phenotype_test_DEGs_df_lists[['Null Distribution Mean DEG Count']],
  `P value` = phenotype_test_DEGs_df_lists[['P value']],
  `DEGs Filtered` = phenotype_test_DEGs_df_lists[['DEGs Filtered']],
  comparison = phenotype_test_DEGs_df_lists[['comparison']],
  stringsAsFactors = FALSE,
  check.names = FALSE
)
phenotype_test_DEGs_df_table = phenotype_test_DEGs_df[c('Cluster', 'DEG Count', 'Null Distribution Mean DEG Count', 'P value', 'comparison')]

phenotype_test_DEGs_df_table[['Fold Change']] = phenotype_test_DEGs_df_table[['DEG Count']] / phenotype_test_DEGs_df_table[['Null Distribution Mean DEG Count']]
phenotype_test_DEGs_df_table[['Fold Change']] = sprintf("%.1f", phenotype_test_DEGs_df_table[['Fold Change']])
phenotype_test_DEGs_df_table[['Null Distribution Mean DEG Count']] = sprintf("%.2f", phenotype_test_DEGs_df_table[['Null Distribution Mean DEG Count']])
phenotype_test_DEGs_df_table[['P value']] = sprintf("%.3f", phenotype_test_DEGs_df_table[['P value']])
phenotype_test_DEGs_df_table = phenotype_test_DEGs_df_table[, c('Cluster', 'DEG Count', 'Null Distribution Mean DEG Count', 'Fold Change', 'P value', 'comparison')]
phenotype_test_DEGs_df_table[['Cluster']] = gsub("c", "-", sub("^cRc", "", phenotype_test_DEGs_df_table[['Cluster']]))

colnames(phenotype_test_DEGs_df_table) = c('Cluster\nName', 'Observed\nDEG Count', 'Expected\nFalse-Positive\nDEG Count', 'DEG Count\nFold Change\n(Observed/FP)', 'Empirical\nP value', 'comparison')

write.csv(phenotype_test_DEGs_df_table, 'phenotype_test_DEGs_df.csv', row.names = FALSE)
# GetRadialTreeFromDataframe(df = phenotype_test_DEGs_df, tree_comparisons = c("control vs treatment naive-CD", "control vs CD", "treatment naive-CD vs CD"), val_name = 'P.value', tree_level = 2, gravities = c('NorthEast', 'NorthEast', 'NorthEast'), title_extra = '_empirical')
# GetRadialTreeFromDataframe(df = phenotype_test_DEGs_df, tree_comparisons = c("control vs treatment naive-CD", "control vs CD", "treatment naive-CD vs CD"), val_name = 'DEG Count', tree_level = 2, gravities = c('NorthEast', 'NorthEast', 'NorthEast'), title_extra = '_empirical_nofilter')
phenotype_test_DEGs_df = phenotype_test_DEGs_df[c('Cluster', 'DEGs Filtered', 'comparison')]
colnames(phenotype_test_DEGs_df) = c('Cluster', 'DEG Count', 'comparison')
# GetRadialTreeFromDataframe(df = phenotype_test_DEGs_df, tree_comparisons = c("control vs treatment naive-CD", "control vs CD", "treatment naive-CD vs CD"), val_name = 'DEG Count', tree_level = 2, gravities = c('NorthEast', 'NorthEast', 'NorthEast'), title_extra = '_empirical', annotation_colors = annotation_colors)

# Split the dataframe by 'comparison'
comparison_levels <- unique(phenotype_test_DEGs_df_table$comparison)
for (comp in comparison_levels) {

  # Filter the dataframe for the current comparison
  df_comp <- phenotype_test_DEGs_df_table %>%
    filter(comparison == comp) %>%
    select(-comparison)  # Remove the 'comparison' column

  # Determine the split point
  n_rows <- nrow(df_comp)
  split_point <- ceiling(n_rows / 2)

  # Split the dataframe into two halves
  df_half1 <- df_comp[1:split_point,]
  df_half2 <- df_comp[(split_point + 1):n_rows,]

  custom_theme <- ttheme_minimal(
    base_size = 10,
    core = list(
      padding = unit(c(2, 2), "mm"),
      fg_params = list(cex = 0.8)
    ),
    colhead = list(
      padding = unit(c(4, 4), "mm"),
      fg_params = list(cex = 0.9)
    )
  )

  table1 <- tableGrob(df_half1, rows = NULL, theme = custom_theme)
  table2 <- tableGrob(df_half2, rows = NULL, theme = custom_theme)

  if (comp == "control vs treatment naive-CD"){
    comp = "Control vs Treatment Naive-CD"
  }
  if (comp == "control vs CD"){
    comp = "Control vs CD"
  }
  if (comp == "treatment naive-CD vs CD"){
    comp = "Treatment Naive-CD vs CD"
  }

  title <- textGrob(paste("Comparison:", comp),
                    gp = gpar(fontsize = 16, fontface = "bold"))

  # Arrange the title and the two tables side by side
  combined <- arrangeGrob(
    title,
    arrangeGrob(table1, table2, ncol = 2, widths = c(1, 1)),
    nrow = 2,
    heights = c(0.02, 1)
  )

  extra = 1.05
  ggsave(paste0("comparison_tables ", comp, ".pdf"), combined, width = 8.5 * extra, height = 11 * extra)
}

print('DEG Analysis Tree Map')

