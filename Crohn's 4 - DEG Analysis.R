# Code for for Crohn's scRNA-seq data DEG analysis

patient_names = samples

# patient_metadata <- read.csv(file = 'patient_metadata helmsley.csv', fileEncoding = "UTF-8-BOM", check.names = FALSE)

# # Using split to create the named list
# phenotype_patients = split(patient_metadata$`Sample ID`, patient_metadata$`Disease phenotype`)

# level_ordered_cluster_names = readRDS('level_ordered_cluster_names.rds')

phenotype_patients = split(samples, patient_disease)

cluster_patient_average_gene_expression_filename = paste0("cluster_patient_average_gene_expression_", assay, "_", c_resolution, '_', DEG_cutoff, '_', entropy_cutoff, '_')
cluster_comparison_gene_count_filename = paste0("cluster_comparison_gene_count_", assay, "_", c_resolution, '_', DEG_cutoff, '_', entropy_cutoff, '_')
cluster_comparison_patient_count_filename = paste0("cluster_comparison_patient_count_", assay, "_", c_resolution, '_', DEG_cutoff, '_', entropy_cutoff, '_')

for (seed_num in seq(0, 300)) {
  #for (seed_num in 0) {
  seed_cluster_patient_average_gene_expression_filename = paste0(cluster_patient_average_gene_expression_filename, seed_num, '.rds')
  if (!file.exists(seed_cluster_patient_average_gene_expression_filename)) {
    cluster_patient_average_gene_expression = list()
    for (cluster in level_ordered_cluster_names[2:length(level_ordered_cluster_names)]) {
      cluster_split = strsplit(cluster, 'c')
      level = length(cluster_split[[1]]) - 3
      level_name = paste0('level_', level)
      if (level < 3) {
        print(cluster)
        patient_ids = recursive_integrated$cR@meta.data$patient

        #shuffle patient ids for cluster
        if (seed_num > 0) {
          patient_ids = sample(patient_ids)
        }

        for (patient in patient_names) {
          cluster_bool = recursive_integrated$cR@meta.data[[level_name]] == cluster
          patient_bool = patient_ids == patient
          cluster_patient_gene_counts_by_cell = recursive_integrated$cR@assays$RNA@data[recursive_integrated$cR@assays$integrated@var.features, cluster_bool & patient_bool]
          if (sum(cluster_bool & patient_bool) >= 3) {
            cluster_patient_mean_gene_counts = unname(rowMeans(cluster_patient_gene_counts_by_cell))
            cluster_patient_average_gene_expression[[cluster]][[patient]] = cluster_patient_mean_gene_counts
          } else {
            print(paste0(cluster, ' ', patient))
            cluster_patient_average_gene_expression[[cluster]][[patient]] = rep(NA, 2000)
          }
        }
        cluster_patient_average_gene_expression[[cluster]] = as.data.frame(do.call(rbind, cluster_patient_average_gene_expression[[cluster]]))
        colnames(cluster_patient_average_gene_expression[[cluster]]) = recursive_integrated$cR@assays$integrated@var.features
      }
    }
    saveRDS(cluster_patient_average_gene_expression, file = seed_cluster_patient_average_gene_expression_filename)
  } else {
    cluster_patient_average_gene_expression = readRDS(seed_cluster_patient_average_gene_expression_filename)
  }

  cluster_comparison_patient_count = list()

  seed_cluster_comparison_gene_count_filename = paste0(cluster_comparison_gene_count_filename, seed_num, '.rds')
  seed_cluster_comparison_patient_count_filename = paste0(cluster_comparison_patient_count_filename, seed_num, '.rds')

  if (!file.exists(seed_cluster_comparison_gene_count_filename)) {
    cluster_comparison_gene_count = list()
    for (cluster in level_ordered_cluster_names[2:length(level_ordered_cluster_names)]) {
      cluster_split = strsplit(cluster, 'c')
      level = length(cluster_split[[1]]) - 3
      level_name = paste0('level_', level)
      if (level < 3) {
        print(cluster)
        for (comparison in list("control vs treatment naive-CD", "control vs CD", "treatment naive-CD vs CD")) {
          comparison_split = strsplit(comparison, ' vs ')[[1]]
          cluster_comparison_gene_count[[cluster]][[comparison]] = 0
          for (gene in recursive_integrated$cR@assays$integrated@var.features) {
            phenotype_patient_expression_list = list()
            for (phenotype_idx in seq(2)) {
              phenotype_patient_expression_list[[phenotype_idx]] = cluster_patient_average_gene_expression[[cluster]][phenotype_patients[[comparison_split[[phenotype_idx]]]], gene]
              phenotype_patient_expression_list[[phenotype_idx]] = phenotype_patient_expression_list[[phenotype_idx]][!is.na(phenotype_patient_expression_list[[phenotype_idx]])]
            }
            if (length(phenotype_patient_expression_list[[1]]) >= 1 & length(phenotype_patient_expression_list[[2]] >= 1)) {
              p = suppressWarnings(wilcox.test(phenotype_patient_expression_list[[1]], phenotype_patient_expression_list[[2]])$p.value)
            } else {
              p = 1
            }
            if (!is.nan(p)) {
              if (p < 0.05) {
                cluster_comparison_gene_count[[cluster]][[comparison]] = cluster_comparison_gene_count[[cluster]][[comparison]] + 1
              }
            }
          }
          if (seed_num == 0) {
            if (length(phenotype_patient_expression_list[[1]]) >= 1 & length(phenotype_patient_expression_list[[2]] >= 1)) {
              cluster_comparison_patient_count[[cluster]][[comparison]] = length(phenotype_patient_expression_list[[1]]) + length(phenotype_patient_expression_list[[2]])
            } else {
              cluster_comparison_patient_count[[cluster]][[comparison]] = 0
            }
          }
        }
      }
    }
    saveRDS(cluster_comparison_gene_count, file = seed_cluster_comparison_gene_count_filename)
    saveRDS(cluster_comparison_patient_count, file = seed_cluster_comparison_patient_count_filename)
  } else {
    cluster_comparison_gene_count = readRDS(seed_cluster_comparison_gene_count_filename)
    if (seed_num == 0) {
      cluster_comparison_patient_count = readRDS(seed_cluster_comparison_patient_count_filename)
    }
  }

  # Function to recursively extract names, sub-names, and values
  flatten_list <- function(x, prefix = "") {
    rows = do.call(rbind, lapply(names(x), function(name) {
      item = x[[name]]
      if (is.list(item)) {
        flatten_list(item, paste(prefix, name, sep = "."))
      } else {
        data.frame(main_name = prefix, sub_name = name, value = item, stringsAsFactors = FALSE)
      }
    }))
    # Remove the dot from main_name
    rows$main_name = sub("\\.", "", rows$main_name)
    rows
  }

  df_for_multilevel = flatten_list(cluster_comparison_gene_count)
  colnames(df_for_multilevel) = c('Cluster', 'comparison', 'DEGs')
  if (seed_num <= 1) {
    GetRadialTreeFromDataframe(df = df_for_multilevel, tree_comparisons = c("control vs treatment naive-CD", "control vs CD", "treatment naive-CD vs CD"), val_name = 'DEGs', tree_level = 2, gravities = c('NorthEast', 'NorthEast', 'NorthEast'), title_extra = seed_num)
  }
  if (seed_num == 0) {
    df_for_multilevel_patient = flatten_list(cluster_comparison_patient_count)
    colnames(df_for_multilevel_patient) = c('Cluster', 'comparison', 'Samples')
    GetRadialTreeFromDataframe(df = df_for_multilevel_patient, tree_comparisons = c("control vs treatment naive-CD", "control vs CD", "treatment naive-CD vs CD"), val_name = 'Samples', tree_level = 2, gravities = c('NorthEast', 'NorthEast', 'NorthEast'))
  }
}