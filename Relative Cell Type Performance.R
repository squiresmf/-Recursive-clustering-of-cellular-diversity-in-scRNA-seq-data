library(dplyr)
library(stringr)
library(colorspace)
library(ggplot2)
library(magick)

# Code for generating scatter plots of cell type cluster purity for recursive verus single-pass method results for each scRNA-seq dataset
# This code will reproduce figure 4 as well as the frequently better captured analysis if main has been run for all four datasets for all 28 resolution settings
# Load the libraries from main first before running

# current_dir = ""
# if (.Platform$OS.type == "windows") {
#   current_dir = paste0(current_dir, 'Y:')
# }
# current_dir = paste0(current_dir, '/qiu-lab/Michael Recursive Clustering/')
current_dir <- paste0(getwd(), '/')

SetDirRead <- function() {
  setwd_path = current_dir
  setwd_path = paste0(setwd_path, 'Default/', reference)
  setwd(setwd_path)
  print(setwd_path)
}

SetDirWrite <- function() {
  setwd_path = current_dir
  setwd_path = paste0(setwd_path, 'Default')
  setwd(setwd_path)
  print(setwd_path)
}

assay = 'integrated'

annotation_levels = list()
annotation_levels[['human_PBMC']] = 2
annotation_levels[['human_adipose']] = 2
annotation_levels[['human_tonsil']] = 3
annotation_levels[['human_fetus']] = 2

cell_names_list = list()
cell_nums_list = list()
cell_type_list = list()

if (!file.exists('cell_names_list.rds')) {
  for (reference in c("PBMC", "adipose", "tonsil", "fetus")) {
    reference_name = paste0("human_", reference)

    SetDirRead()
    getwd()
    reference_integrated = readRDS(paste0('reference_', reference_name, '_integrated.rds'))
    for (annotation_level in seq(2)) {
      annotation_level_final = annotation_levels[[reference_name]] - (2 - annotation_level)
      cell_names_list[[reference_name]][[annotation_level]] = sort(unique(reference_integrated@meta.data[[paste0('celltype.l', annotation_level_final)]]))
    }
    for (cell_type_l2 in cell_names_list[[reference_name]][[2]]) {
      cell_type_l1 = reference_integrated@meta.data[[paste0('celltype.l', annotation_level_final - 1)]][[match(cell_type_l2, reference_integrated@meta.data[[paste0('celltype.l', annotation_level_final)]])]]
      cell_type_list[[reference_name]][[cell_type_l2]] = cell_type_l1
      for (cell_type_idx in seq_along(c(cell_type_l1, cell_type_l2))) {
        cell_type = c(cell_type_l1, cell_type_l2)[[cell_type_idx]]
        num_cells_cell_type = sum(reference_integrated@meta.data[[paste0('celltype.l', cell_type_idx)]] == cell_type)
        num_cells_total = length(reference_integrated@meta.data[[paste0('celltype.l1')]])
        cell_nums_list[[reference_name]][[paste0('celltype.l', cell_type_idx)]][[cell_type]] = num_cells_cell_type / num_cells_total
      }
    }
  }
  SetDirWrite()
  saveRDS(cell_nums_list, 'cell_nums_list.rds')
  saveRDS(cell_type_list, 'cell_type_list.rds')
  saveRDS(cell_names_list, 'cell_names_list.rds')
} else {
  SetDirWrite()
  cell_names_list = readRDS('cell_names_list.rds')
  cell_type_list = readRDS('cell_type_list.rds')
  cell_nums_list = readRDS('cell_nums_list.rds')
}

cell_nums_list_sorted = list()
cell_nums_list_ranked = list()
cell_nums_list_ranked_rev = list()
for (reference_name in names(cell_nums_list)) {
  for (annotation_level_name in names(cell_nums_list[[reference_name]])) {
    cell_nums_list_sorted[[reference_name]][[annotation_level_name]] = cell_nums_list[[reference_name]][[annotation_level_name]][order(unlist(cell_nums_list[[reference_name]][[annotation_level_name]]), decreasing = TRUE)]
    for (cell_type in names(cell_nums_list_sorted[[reference_name]][[annotation_level_name]])) {
      cell_nums_list_ranked[[reference_name]][[annotation_level_name]][[cell_type]] = match(cell_type, names(cell_nums_list_sorted[[reference_name]][[annotation_level_name]]))
      cell_nums_list_ranked_rev[[reference_name]][[annotation_level_name]][[cell_type]] = match(cell_type, rev(names(cell_nums_list_sorted[[reference_name]][[annotation_level_name]])))
    }
  }
}

results_df_list = list()
cell_type_purity_results_list = list()

reference_names = c("human_PBMC", "human_adipose", "human_tonsil", "human_fetus")
reference_names_pretty = c("PBMC", "Adipose", "Tonsil", "Fetus")

for (reference in c("PBMC", "adipose", "tonsil", "fetus")) {
  reference_name = paste0("human_", reference)

  # for (reference_name in c("human_PBMC")) {
  for (c_resolution in seq(0.015, 0.15, 0.005)) {
    # for (c_resolution in 0.02) {
    for (DEG_cutoff in c(0)) {
      cluster_purity_filename = paste0('results_', reference_name, '_integrated_', c_resolution, '_', DEG_cutoff, '.rds')
      if (DEG_cutoff == 0) {
        SetDirRead()
        cluster_purity_results = readRDS(cluster_purity_filename)
        for (annotation_level in seq(2)) {
          annotation_level_final = annotation_levels[[reference_name]] - (2 - annotation_level)
          annotation_level_name = paste0('Annotation Level ', annotation_level)
          for (level in 2) {
            for (result_type in 'Cell Type Mean Cluster Purity') {
              value_iterative = cluster_purity_results[[paste0('annotation level ', annotation_level_final)]][[paste0('level ', level)]][['Iterative']][[result_type]]
              value_seurat = cluster_purity_results[[paste0('annotation level ', annotation_level_final)]][[paste0('level ', level)]][['Seurat Equivalent']][[result_type]]
              results_df = cbind.data.frame(list(value_iterative = value_iterative, value_seurat = value_seurat))
              rownames(results_df) = cell_names_list[[reference_name]][[annotation_level]]

              results_df_list[[reference_name]][[annotation_level_name]][[paste(c_resolution)]] = results_df
              results_df_list[[reference_name]][[annotation_level_name]][[paste(c_resolution)]]['Reference'] = reference_names_pretty[[match(reference_name, reference_names)]]
              results_df_list[[reference_name]][[annotation_level_name]][[paste(c_resolution)]]['Annotation'] = paste0('Annotation Level ', annotation_level, '\n\n\n')
              results_df_list[[reference_name]][[annotation_level_name]][[paste(c_resolution)]]['Resolution'] = paste(c_resolution)

              if (c_resolution == 0.015) {
                cell_type_purity_results_list[[reference_name]][[annotation_level_name]][['iterative_better']] = results_df
                cell_type_purity_results_list[[reference_name]][[annotation_level_name]][['seurat_better']] = results_df
              }
              cell_type_purity_results_list[[reference_name]][[annotation_level_name]][['iterative_better']][[paste(c_resolution)]] = (results_df$value_iterative > 0.75) & (results_df$value_iterative > 1.5 * results_df$value_seurat)
              cell_type_purity_results_list[[reference_name]][[annotation_level_name]][['seurat_better']][[paste(c_resolution)]] = (results_df$value_seurat > 0.75) & (results_df$value_seurat > 1.5 * results_df$value_iterative)
              if (c_resolution == 0.015) {
                cell_type_purity_results_list[[reference_name]][[annotation_level_name]][['iterative_better']] = cell_type_purity_results_list[[reference_name]][[annotation_level]][['iterative_better']]['0.015']
                cell_type_purity_results_list[[reference_name]][[annotation_level_name]][['seurat_better']] = cell_type_purity_results_list[[reference_name]][[annotation_level]][['seurat_better']]['0.015']
              }
            }
          }
        }
      }
    }
  }
}
results_df_all = bind_rows(unlist(unlist(results_df_list, recursive = FALSE), recursive = FALSE))


print(results_df)
better_results_list = list()
for (annotation_level in seq(2)) {
  annotation_level_name = paste0('Annotation Level ', annotation_level)
  for (reference_name in names(cell_type_purity_results_list)) {
    for (better_type in c('iterative_better', 'seurat_better')) {
      #print(paste('annotation_level', annotation_level, 'reference_name', reference_name, 'better_type', better_type))
      better_results = rowMeans(cell_type_purity_results_list[[reference_name]][[annotation_level_name]][[better_type]])
      better_results = better_results[rev(order(better_results))]
      #better_results = better_results[better_results > 0]
      better_results_list[[annotation_level_name]][[reference_name]][[better_type]] = better_results
    }
  }
}

cell_type_l1_vec = vector()
cell_type_l2_vec = vector()
cell_type_nums_1_vec = vector()
cell_type_nums_2_vec = vector()
cell_type_rank_1_vec = vector()
cell_type_rank_2_vec = vector()
cell_type_rank_rev_1_vec = vector()
cell_type_rank_rev_2_vec = vector()
better_value_1_vec = vector()
better_value_2_vec = vector()
worse_value_1_vec = vector()
worse_value_2_vec = vector()
reference_name_vec = vector()
better_type_vec = vector()
top_5_vec = vector()

for (reference_name in names(cell_type_list)) {
  for (better_type in c('iterative_better', 'seurat_better')) {
    for (cell_type_l2 in names(better_results_list[['Annotation Level 2']][[reference_name]][[better_type]])) {
      if (length(better_results_list[[paste0('Annotation Level ', annotation_level)]][[reference_name]][[better_type]]) > 0) {
        cell_type_l1 = cell_type_list[[reference_name]][[cell_type_l2]]

        cell_type_nums = list()
        cell_type_rank = list()
        cell_type_rank_rev = list()
        better_value = list()
        worse_value = list()

        for (annotation_level in seq(2)) {
          cell_type_nums[[annotation_level]] = cell_nums_list[[reference_name]][[paste0('celltype.l', annotation_level)]][[c(cell_type_l1, cell_type_l2)[[annotation_level]]]]
          cell_type_rank[[annotation_level]] = paste0(cell_nums_list_ranked[[reference_name]][[paste0('celltype.l', annotation_level)]][[c(cell_type_l1, cell_type_l2)[[annotation_level]]]],
                                                      '/',
                                                      max(unlist(cell_nums_list_ranked[[reference_name]][[paste0('celltype.l', annotation_level)]])))
          cell_type_rank_rev[[annotation_level]] = paste0(cell_nums_list_ranked_rev[[reference_name]][[paste0('celltype.l', annotation_level)]][[c(cell_type_l1, cell_type_l2)[[annotation_level]]]],
                                                          '/',
                                                          max(unlist(cell_nums_list_ranked_rev[[reference_name]][[paste0('celltype.l', annotation_level)]])))
          better_value[[annotation_level]] = better_results_list[[paste0('Annotation Level ', annotation_level)]][[reference_name]][[better_type]][[c(cell_type_l1, cell_type_l2)[[annotation_level]]]]
          worse_value[[annotation_level]] = better_results_list[[paste0('Annotation Level ', annotation_level)]][[reference_name]][[c('iterative_better', 'seurat_better')[[3 - match(better_type, c('iterative_better', 'seurat_better'))]]]][[c(cell_type_l1, cell_type_l2)[[annotation_level]]]]
        }
        cell_type_l1_vec = c(cell_type_l1_vec, cell_type_l1)
        cell_type_l2_vec = c(cell_type_l2_vec, cell_type_l2)
        cell_type_nums_1_vec = c(cell_type_nums_1_vec, cell_type_nums[[1]])
        cell_type_nums_2_vec = c(cell_type_nums_2_vec, cell_type_nums[[2]])
        cell_type_rank_1_vec = c(cell_type_rank_1_vec, cell_type_rank[[1]])
        cell_type_rank_2_vec = c(cell_type_rank_2_vec, cell_type_rank[[2]])
        cell_type_rank_rev_1_vec = c(cell_type_rank_rev_1_vec, cell_type_rank_rev[[1]])
        cell_type_rank_rev_2_vec = c(cell_type_rank_rev_2_vec, cell_type_rank_rev[[2]])
        better_value_1_vec = c(better_value_1_vec, better_value[[1]])
        better_value_2_vec = c(better_value_2_vec, better_value[[2]])
        worse_value_1_vec = c(worse_value_1_vec, worse_value[[1]])
        worse_value_2_vec = c(worse_value_2_vec, worse_value[[2]])
        reference_name_vec = c(reference_name_vec, reference_name)
        better_type_vec = c(better_type_vec, better_type)
        if (cell_type_l2 %in% names(better_results_list[['Annotation Level 2']][[reference_name]][[better_type]])[1:min(5, length(names(better_results_list[['Annotation Level 2']][[reference_name]][[better_type]])))]) {
          top_5_vec = c(top_5_vec, TRUE)
        } else {
          top_5_vec = c(top_5_vec, FALSE)
        }
      }
    }
  }
}

better_results_df_2 = cbind.data.frame(reference_name_vec, better_type_vec, cell_type_l2_vec, better_value_2_vec, worse_value_2_vec, cell_type_nums_2_vec, cell_type_rank_2_vec, cell_type_rank_rev_2_vec, cell_type_l1_vec, better_value_1_vec, worse_value_1_vec, cell_type_nums_1_vec, cell_type_rank_1_vec, cell_type_rank_rev_1_vec, top_5_vec)
colnames(better_results_df_2) = c('reference_name', 'better_type', 'cell_type_l2', 'better_value_2', 'worse_value_2', 'cell_type_nums_2', 'cell_type_rank_2', 'cell_type_rank_rev_2', 'cell_type_l1', 'better_value_1', 'worse_value_1', 'cell_type_nums_1', 'cell_type_rank_1', 'cell_type_rank_rev_1', 'top_5')
better_results_df_2_final = better_results_df_2[(better_results_df_2['better_value_2'] > 0.25) & (better_results_df_2['better_value_2'] >= better_results_df_2['worse_value_2'] * 2),]

frequently_better_captured_count = list()

for (reference_name in names(cell_type_list)) {
  for (better_type in c('iterative_better', 'seurat_better')) {
    frequently_better_captured_count[['Annotation Level 2']][[better_type]][[reference_name]] = paste0(sum((better_results_df_2_final[['reference_name']] == reference_name) & (better_results_df_2_final[['better_type']] == better_type)), ' of the ', length(cell_names_list[[reference_name]][[2]]))
  }
}

better_results_df_2_final = better_results_df_2_final[better_results_df_2_final[['top_5']],]

cell_type_l1_vec = vector()
cell_type_nums_1_vec = vector()
cell_type_rank_1_vec = vector()
cell_type_rank_rev_1_vec = vector()
better_value_1_vec = vector()
worse_value_1_vec = vector()
reference_name_vec = vector()
better_type_vec = vector()
top_5_vec = vector()

for (reference_name in names(cell_type_list)) {
  for (better_type in c('iterative_better', 'seurat_better')) {
    for (cell_type_l1 in names(better_results_list[['Annotation Level 1']][[reference_name]][[better_type]])) {

      if (length(better_results_list[[paste0('Annotation Level ', annotation_level)]][[reference_name]][[better_type]]) > 0) {

        cell_type_nums = list()
        cell_type_rank = list()
        cell_type_rank_rev = list()
        better_value = list()
        worse_value = list()

        for (annotation_level in seq(1)) {
          cell_type_nums[[annotation_level]] = cell_nums_list[[reference_name]][[paste0('celltype.l', annotation_level)]][[c(cell_type_l1, cell_type_l2)[[annotation_level]]]]
          cell_type_rank[[annotation_level]] = paste0(cell_nums_list_ranked[[reference_name]][[paste0('celltype.l', annotation_level)]][[c(cell_type_l1, cell_type_l2)[[annotation_level]]]],
                                                      '/',
                                                      max(unlist(cell_nums_list_ranked[[reference_name]][[paste0('celltype.l', annotation_level)]])))
          cell_type_rank_rev[[annotation_level]] = paste0(cell_nums_list_ranked_rev[[reference_name]][[paste0('celltype.l', annotation_level)]][[c(cell_type_l1, cell_type_l2)[[annotation_level]]]],
                                                          '/',
                                                          max(unlist(cell_nums_list_ranked_rev[[reference_name]][[paste0('celltype.l', annotation_level)]])))
          better_value[[annotation_level]] = better_results_list[[paste0('Annotation Level ', annotation_level)]][[reference_name]][[better_type]][[c(cell_type_l1, cell_type_l2)[[annotation_level]]]]
          worse_value[[annotation_level]] = better_results_list[[paste0('Annotation Level ', annotation_level)]][[reference_name]][[c('iterative_better', 'seurat_better')[[3 - match(better_type, c('iterative_better', 'seurat_better'))]]]][[c(cell_type_l1, cell_type_l2)[[annotation_level]]]]
        }
        cell_type_l1_vec = c(cell_type_l1_vec, cell_type_l1)
        cell_type_nums_1_vec = c(cell_type_nums_1_vec, cell_type_nums[[1]])
        cell_type_rank_1_vec = c(cell_type_rank_1_vec, cell_type_rank[[1]])
        cell_type_rank_rev_1_vec = c(cell_type_rank_rev_1_vec, cell_type_rank_rev[[1]])
        better_value_1_vec = c(better_value_1_vec, better_value[[1]])
        worse_value_1_vec = c(worse_value_1_vec, worse_value[[1]])
        reference_name_vec = c(reference_name_vec, reference_name)
        better_type_vec = c(better_type_vec, better_type)
        if (cell_type_l1 %in% names(better_results_list[['Annotation Level 1']][[reference_name]][[better_type]])[1:min(5, length(names(better_results_list[['Annotation Level 2']][[reference_name]][[better_type]])))]) {
          top_5_vec = c(top_5_vec, TRUE)
        } else {
          top_5_vec = c(top_5_vec, FALSE)
        }
      }
    }
  }
}

better_results_df_1 = cbind.data.frame(reference_name_vec, better_type_vec, cell_type_l1_vec, better_value_1_vec, worse_value_1_vec, cell_type_nums_1_vec, cell_type_rank_1_vec, cell_type_rank_rev_1_vec, top_5_vec)
colnames(better_results_df_1) = c('reference_name', 'better_type', 'cell_type_l1', 'better_value_1', 'worse_value_1', 'cell_type_nums_1', 'cell_type_rank_1', 'cell_type_rank_rev_1', 'top_5')
better_results_df_1_final = better_results_df_1[(better_results_df_1['better_value_1'] > 0.25) & (better_results_df_1['better_value_1'] >= better_results_df_1['worse_value_1'] * 2),]

for (reference_name in names(cell_type_list)) {
  for (better_type in c('iterative_better', 'seurat_better')) {
    frequently_better_captured_count[['Annotation Level 1']][[better_type]][[reference_name]] = paste0(sum((better_results_df_1_final[['reference_name']] == reference_name) & (better_results_df_1_final[['better_type']] == better_type)), ' of the ', length(cell_names_list[[reference_name]][[1]]))
  }
}

better_results_df_1_final = better_results_df_1_final[better_results_df_1_final[['top_5']],]

better_results_df_final = list()
better_results_df_final[['1']] = better_results_df_1_final
better_results_df_final[['2']] = better_results_df_2_final

for (reference_name in names(cell_type_list)) {
  for (annotation_level in names(better_results_df_final)) {
    for (better_type in c('iterative_better', 'seurat_better')) {
      cell_types = better_results_df_final[[annotation_level]][(better_results_df_final[[annotation_level]][['reference_name']] == reference_name) & (better_results_df_final[[annotation_level]][['better_type']] == better_type),]
      cell_types = cell_types[, c(paste0('cell_type_l', annotation_level), paste0('better_value_', annotation_level))]
      cell_types[[paste0('better_value_', annotation_level)]] = round(cell_types[[paste0('better_value_', annotation_level)]], 2) * 100
      cols <- c(paste0('cell_type_l', annotation_level), paste0('better_value_', annotation_level))
      cell_types$x <- apply(cell_types[, cols], 1, paste, collapse = " (")
      if (nrow(cell_types) > 0) {
        cell_types$x <- paste0(cell_types$x, "%)")
      }
      cell_types <- cell_types[, !(names(cell_types) %in% cols)]
      cell_types = paste(unlist(cell_types), collapse = ', ')
      reference_name_split = 1
      if (reference_name == 'human_fetus') {
        reference_name_split = 'Fetal'
      } else if (reference_name == 'human_PBMC') {
        reference_name_split = (strsplit(reference_name, split = '_')[[1]][[2]])
      } else {
        reference_name_split = str_to_title((strsplit(reference_name, split = '_')[[1]][[2]]))
      }
      if (better_type == 'iterative_better') {
        method_name = 'recursive'
      } else {
        method_name = 'single-pass'
      }
      count_better = strtoi(str_trim(substr(frequently_better_captured_count[[paste0('Annotation Level ', annotation_level)]][[better_type]][[reference_name]], start = 1, stop = 2)))
      if (count_better > 5) {
        extra_message = ', where the five most frequently better captured are'
      } else {
        extra_message = ''
      }
      if (annotation_level == 1) {
        annotation_level_word = 'first'
      } else {
        annotation_level_word = 'second'
      }
      if (count_better > 0) {
        print(paste0('For the ', reference_name_split, ' reference, at the ', annotation_level_word, ' annotation level, the ', method_name, ' method frequently better captures ', frequently_better_captured_count[[paste0('Annotation Level ', annotation_level)]][[better_type]][[reference_name]], ' cell types', extra_message, ': ', cell_types, sep = ', '))
      } else {
        print(paste0('For the ', reference_name_split, ' reference, at the ', annotation_level_word, ' annotation level, the ', method_name, ' method does not frequently better capture any cell types.'))
      }
    }
  }
}

# Custom function to format labels
format_labels <- function(x) {
  sapply(x, function(num) {
    if (is.na(num)) {
      return(NA)
    } else if (num == floor(num)) {
      return(as.character(as.integer(num)))
    } else {
      return(as.character(num))
    }
  })
}


segments_data <- data.frame(x = c(1, 2, 3), xend = c(4, 5, 6),
                            y = c(1, 2, 3), yend = c(4, 5, 6),
                            Key = c("Recursive Better Captured ", "Equivalent", "Single-Pass Better Captured"))

results_df_all$Reference = factor(results_df_all$Reference, levels = reference_names_pretty)
bbox_color1 = lighten("#669bcc", amount = 0.5)
bbox_color2 = "#FFD700"
middle_line = lighten("#CC6677", 0.3)
size = 0.6
linetype = 'longdash'


plot = ggplot(results_df_all, aes(x = value_seurat, y = value_iterative)) +
  geom_point(alpha = 0.15, color = darken("#A9A9A9", amount = 0.9)) +
  theme_minimal(base_size = 16) +  # Increased base text size
  geom_segment(data = segments_data, aes(x = x, y = y, xend = xend, yend = yend,
                                         color = Key, linetype = Key), size = 0.6) +
  scale_color_manual(values = c("Recursive Better Captured             " = bbox_color1, "Equivalent             " = middle_line, "Single-Pass Better Captured             " = bbox_color2)) +
  scale_linetype_manual(values = c("Recursive Better Captured             " = "longdash", "Equivalent             " = "longdash", "Single-Pass Better Captured             " = "longdash")) +
  scale_x_continuous(limits = c(0, 1), labels = format_labels) +  # Apply custom function to x-axis
  scale_y_continuous(limits = c(0, 1), labels = format_labels) +    # Apply custom function to y-axis
  geom_segment(aes(x = 0.5, y = 0.75, xend = 2 / 3, yend = 1), linetype = linetype, color = bbox_color1, size = size) +
  geom_segment(aes(x = 0, y = 0.75, xend = 0.5, yend = 0.75), linetype = linetype, color = bbox_color1, size = size) +
  geom_segment(aes(y = 0.5, x = 0.75, yend = 2 / 3, xend = 1), linetype = linetype, color = bbox_color2, size = size) +
  geom_segment(aes(y = 0, x = 0.75, yend = 0.5, xend = 0.75), linetype = linetype, color = bbox_color2, size = size) +
  geom_segment(aes(y = 0, x = 0, yend = 1, xend = 1), linetype = linetype, color = middle_line, size = size) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'), plot.background = element_rect(fill = 'white', colour = 'white'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(color = "grey", fill = NA, size = 1),
        strip.text.y = element_text(angle = 0, vjust = 0, size = 13, margin = margin(r = -11)),  # Increased facet label size
        strip.placement = "outside",
        axis.title.y = element_text(size = 13, vjust = -15, angle = 90),
        legend.key.width = unit(1.5, "cm"),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        legend.margin = margin(0, 20, 0, 0),
        legend.title = element_blank(),
        legend.position = 'bottom',
        legend.spacing.x = unit(0.05, "cm"),
        axis.title.x = element_text(size = 13, margin = margin(t = 5, unit = "pt"), vjust = 0),  # Increased axis text size
        axis.text = element_text(size = 13),  # Increased axis text size
        legend.text = element_text(size = 13),  # Increased legend text size
        panel.spacing = unit(0.35, "lines")
  ) +
  labs(x = "Single-Pass Cell Type Purity", y = "Recursive Cell Type Purity") +
  facet_grid(Annotation ~ Reference, switch = "y")

SetDirWrite()
ggsave('cell type mean cluster purity plot.png', plot, dpi = 300, width = 9, height = 5, device = 'png')
#suppressMessages(ggsave(filename = paste0('PBMC example heatmaps.png'), width = 10, height = 11, plot = temp2, device = 'png', dpi = 400))
pp1 <- image_trim(image_read(paste0('cell type mean cluster purity plot.png')))
#pp1 <- image_scale(pp1, "50%")
image_write(pp1, paste0('cell type mean cluster purity plot.eps'), format = 'eps')
image_write(pp1, paste0('cell type mean cluster purity plot.pdf'), format = 'pdf')

