library(magick)
library(ggplot2)


algorithm = 'Louvain'
# algorithm = 'Leiden'
# algorithm = 'SLM'
HVGs = 2000
# HVGs = 1000
# HVGs = 3000

evaluation_types = 'cell type'
# evaluation_types = 'consistency'

is_res_range_smalls = FALSE
# is_res_range_smalls = c(TRUE, TRUE)

# current_dir = ""
# if (.Platform$OS.type == "windows") {
#   current_dir = paste0(current_dir, 'Y:')
# }
# current_dir = paste0(current_dir, '/qiu-lab/Michael Recursive Clustering/')
current_dir <- getwd()

celltype_counts <- data.frame(
  Reference = character(),
  `Annotation Level` = character(),
  Count = integer(),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cell_types = list()
for (reference in c("PBMC", "Adipose", "Tonsil", "Fetus")) {
  for (level in c(1, 2)) {
    file_name <- paste0(ifelse(reference == "PBMC", reference, tolower(reference)), " celltype.l", level, " unique.rds")
    cell_types[[reference]][[level]] <- readRDS(file_name)
    count <- length(cell_types[[reference]][[level]])

    annotation_label <- paste("Annotation Level", level)

    celltype_counts <- rbind(
      celltype_counts,
      data.frame(
        Reference = reference,
        `Annotation Level` = annotation_label,
        Count = count,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    )
  }
}

celltype_counts$Count <- sprintf("%-14s", paste(celltype_counts$Count, 'cell types'))

#
# for (algorithm in c('Louvain', 'Leiden', 'SLM')) {
#   # for (algorithm in c('Louvain')) {
#   if (algorithm == 'Louvain') {
#     HVGs_range = c(2000, 3000, 1000)
#     # HVGs_range = c(2000)
#   } else {
#     HVGs_range = c(2000)
#   }
#   for (HVGs in HVGs_range) {

fvf.nfeatures = HVGs

name_extra = ""
if (algorithm != 'Louvain') {
  name_extra = paste(name_extra, algorithm)
}
if (HVGs != 2000) {
  name_extra = paste(name_extra, HVGs, 'HVG')
  name_extra = paste(name_extra, fvf.nfeatures, 'fvf.nfeatures')
}

setwd_path = current_dir
if (name_extra == "") {
  setwd_path = paste0(setwd_path, 'Default')
} else {
  setwd_path = paste0(setwd_path, trimws(name_extra))
}
setwd(setwd_path)
print(setwd_path)

#
# if (name_extra == "") {
#   evaluation_types = c('cell type', 'consistency')
# } else {
#   evaluation_types = c('cell type')
# }
# for (evaluation_type in evaluation_types) {
#   if (evaluation_type == 'cell type' && name_extra == "") {
#     is_res_range_smalls = c(FALSE, TRUE, TRUE)
#   } else {
#     is_res_range_smalls = FALSE
#   }

for (evaluation_type in evaluation_types) {
  if (evaluation_type == 'cell type') {
    if (name_extra == "") {
      results_df1 <- read.csv('metrics_df.csv')
      results_df2 <- read.csv('balanced_metrics_df.csv')
      results_df <- rbind(results_df1, results_df2)
      colnames(results_df) <- gsub("\\.", " ", colnames(results_df))
    } else {
      results_df <- read.csv(paste0('metrics_df', name_extra, '.csv'))
      colnames(results_df) <- gsub("\\.", " ", colnames(results_df))
    }

  } else {
    results_df1 <- read.csv('consistency_metrics_df.csv')
    results_df2 <- read.csv('shared_cluster_members.csv')
    results_df <- rbind(results_df1, results_df2)
    colnames(results_df) <- gsub("\\.", " ", colnames(results_df))
  }

  results_df$Reference = factor(results_df$Reference, levels = c("PBMC", "Adipose", "Tonsil", "Fetus"))
  celltype_counts$Reference = factor(celltype_counts$Reference, levels = c("PBMC", "Adipose", "Tonsil", "Fetus"))

  if (evaluation_type == 'cell type') {
    metric_types = c('Mean Cell Type Cluster Purity', 'Precision', 'Recall', 'F1-score')
    if (name_extra == "") {
      metric_types = c(metric_types, 'Balanced Completeness', 'Balanced Homogeneity', 'Balanced V-Measure')
    }
  }
  if (evaluation_type == 'consistency') {
    metric_types = c('ARI', 'NMI', 'Shared Cluster Members')
  }

  res_range_smalls_num = 1
  for (is_res_range_small in is_res_range_smalls) {

    name_extra = ""
    if (is_res_range_small) {
      name_extra = paste(name_extra, "small")
    }

    for (metric_type in metric_types) {
      print(metric_type)
      results_subset_df = results_df[results_df$Metric == metric_type,]

      if (metric_type == 'Mean Cell Type Cluster Purity') {
        metric_type = 'Mean Cell Type Purity'
      }

      if (!is_res_range_small) {
        results_subset_df <- results_subset_df[results_subset_df$Resolution %in% round(seq(0.015, 0.15, 0.005), 4),]
      }  else {
        # small_rows = (results_subset_df$Reference %in% c("PBMC", "Adipose", "Tonsil") & results_subset_df$Resolution %in% round(seq(0.005, 0.05, 0.001), 4)) |
        #   (!results_subset_df$Reference %in% c("PBMC", "Adipose", "Tonsil") & results_subset_df$Resolution %in% round(seq(0.0025, 0.01, 0.0002), 4))
        small_rows = (results_subset_df$`Number of Clusters` < 50)
        results_subset_df[!small_rows & results_subset_df$Method == 'Recursive', 'Method'] = '\nRecursive\n50+ Clusters\n'
        results_subset_df[!small_rows & results_subset_df$Method == 'Single-Pass', 'Method'] = '\nSingle-Pass\n50+ Clusters\n'
        results_subset_df[small_rows & results_subset_df$Method == 'Recursive', 'Method'] = '\nRecursive\n<50 Clusters\n'
        results_subset_df[small_rows & results_subset_df$Method == 'Single-Pass', 'Method'] = '\nSingle-Pass\n<50 Clusters\n'
      }

      aes_string2 <- function(...) {
        args <- lapply(list(...), function(x) sprintf("`%s`", x))
        do.call(aes_string, args)
      }

      if (metric_type == 'NMI') {
        alpha = 0.2
        size = 1.75
      } else {
        alpha = 0.3
        size = 2
      }
      scale_factor = 0.925

      if (!is_res_range_small) {
        plot_scale_color_manual = scale_color_manual(values = c("#669bcc", "#ccb166"))
      } else {
        if (res_range_smalls_num == 1) {
          plot_scale_color_manual = scale_color_manual(values = c("#aaaeb3", "#669bcc", "#b3b1aa", "#ccb166"))
        }
        if (res_range_smalls_num == 2) {
          plot_scale_color_manual = scale_color_manual(values = c("#669bcc", "#aaaeb3", "#ccb166", "#b3b1aa"))
        }
      }

      results_plot = ggplot(results_subset_df, aes_string2(x = "Number of Clusters", y = 'Value')) +
        geom_point(aes_string(color = 'Method'), alpha = alpha, size = size * scale_factor) +
        scale_size_area() +
        scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.30, 0.60, 0.90), expand = c(0.001, 0)) +
        plot_scale_color_manual +
        geom_hline(yintercept = 0, color = "dark grey", linetype = 1, size = 0.5) +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "dark grey")
        ) +
        theme(
          strip.background = element_rect(colour = "white", fill = "white"),
          strip.text.y.left = element_text(angle = 90),
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = 12),
          strip.placement = "outside",
          panel.spacing = unit(10, "pt"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 10)
        ) +
        labs(y = metric_type)

      if (evaluation_type == 'cell type') {
        custom_labeller <- labeller(
          `Annotation Level` = function(x) paste0(x, "\n")
        )
        results_plot <- results_plot +
          facet_grid(
            `Annotation Level` ~ Reference,
            switch = "y",
            scales = "free",
            labeller = custom_labeller
          )
      }
      if (evaluation_type == 'consistency') {
        results_plot = results_plot +
          facet_grid(~Reference, switch = "y", scales = "free")
      }
      if (is_res_range_small) {
        results_plot = results_plot + theme(
          legend.background = element_blank(),  # Removes the background of the entire legend
          legend.key = element_blank()          # Removes the background behind the legend symbols
        )       # Removes the background behind the legend symbols
        # # scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60)) +
        # # expand_limits(x = c(7, 53))
      }

      if (evaluation_type == 'cell type') {
        results_plot = results_plot +
          geom_text(
            data = celltype_counts,
            aes(x = Inf, y = Inf, label = Count),
            hjust = 1.45,
            vjust = 29.8,
            inherit.aes = FALSE,
            size = 2.5,
            color = "#1a1a1a",
            fontface = "italic"
          ) +
          theme(
            axis.title.y = element_text(vjust = -15, angle = 90),
          )
      } else {
        results_plot = results_plot +
          theme(axis.title.y = element_text(vjust = 2, angle = 90),
          )
      }
      graphics.off()
      # print(results_plot)

      if (!is_res_range_small) {
        name = paste0(metric_type, name_extra)
      } else {
        if (res_range_smalls_num == 1) {
          name = paste0(metric_type, name_extra, ' big colors')
        }
        if (res_range_smalls_num == 2) {
          name = paste0(metric_type, name_extra, ' small colors')
        }
      }
      ggsave(paste0(name, '.png'), results_plot, width = 8 * scale_factor, height = 6 * scale_factor, device = 'png', dpi = 400)
      pp1 <- image_trim(image_read(paste0(name, '.png')))
      image_write(pp1, paste0(name, '.eps'), format = 'eps')
      image_write(pp1, paste0(name, '.pdf'), format = 'pdf')
    }
    if (is_res_range_small && res_range_smalls_num == 1) {
      res_range_smalls_num = 2
    }
  }
}
#   }
# }