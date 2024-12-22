library(Seurat)
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)
library(stringr)
library(magick)
library(scales)

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
}

SetDirWrite <- function() {
  setwd_path = current_dir
  setwd_path = paste0(setwd_path, 'Default/')
  setwd(setwd_path)
}

# A helper function to reduce repetitive code for generating Seurat clustering plots
GetClusterPlot <- function(dataset, title, group.by = "seurat_clusters", label = TRUE, shuffle_colors = FALSE, label_alpha = 0.7, label.size = 2.75) {
  set.seed(0)
  if (length(group.by) > 1) {
    plot_title = paste0(title, ' ', group.by[[length(group.by)]])
  } else {
    plot_title = paste0(title, ' ', group.by)
  }

  if (shuffle_colors) {
    num_clusters = suppressWarnings(length(as.numeric(levels(as.factor(dataset@meta.data[[group.by]])))))
    temp_plot = DimPlot(dataset, group.by = group.by, label = label, label.size = label.size) + NoLegend()
    plot_colors = ggplot_build(temp_plot)$data[[1]]$colour
    plot_cluster_nums = ggplot_build(temp_plot)$data[[1]]$group
    new_colors = sample(plot_colors[match(seq(num_clusters), plot_cluster_nums)])
    seurat_clustering_umap_plot = DimPlot(dataset, group.by = group.by, label = label, cols = new_colors, label.size = label.size, label.color = alpha("black", label_alpha)) +
      NoLegend() +
      ggtitle(plot_title)
  } else {
    seurat_clustering_umap_plot = DimPlot(dataset, group.by = group.by, label = label, label.size = label.size, label.color = alpha("black", label_alpha)) +
      NoLegend() +
      ggtitle(plot_title)
  }

  if (label) {
    seurat_clustering_umap_plot = seurat_clustering_umap_plot + NoLegend()
  }

  # Store the plot in the misc slot with a unique name
  plot_key = paste0(paste0(group.by, collapse = '_'), '_UMAP_', title)
  dataset@misc[[plot_key]] = seurat_clustering_umap_plot
  return(dataset)
}

# Modify Seurat object metadata to change the cluster naming format from cRc1c2c3 to 1-2-3
modify_meta_data <- function(seurat_obj, target_columns, remove_n_chars = 3,
                             pattern_to_replace = "c", replacement = "-",
                             case_sensitive = TRUE, handle_short_strings = "empty") {

  for (col in target_columns) {
    if (!(col %in% colnames(seurat_obj@meta.data))) {
      next
    }

    is_factor <- is.factor(seurat_obj@meta.data[[col]])

    modify_string <- function(x) {
      x_char <- as.character(x)

      if (nchar(x_char) > remove_n_chars) {
        x_mod <- substr(x_char, remove_n_chars + 1, nchar(x_char))
      } else {
        if (handle_short_strings == "empty") {
          x_mod <- ""
        } else if (handle_short_strings == "original") {
          x_mod <- x_char
        } else {
          x_mod <- ""
        }
      }

      if (case_sensitive) {
        x_mod <- gsub(pattern_to_replace, replacement, x_mod)
      } else {
        x_mod <- gsub(pattern_to_replace, replacement, x_mod, ignore.case = TRUE)
      }
      return(x_mod)
    }

    if (is_factor) {
      original_levels <- levels(seurat_obj@meta.data[[col]])
      modified_levels <- sapply(original_levels, modify_string, USE.NAMES = FALSE)

      seurat_obj@meta.data[[col]] <- factor(
        seurat_obj@meta.data[[col]],
        levels = original_levels,
        labels = modified_levels
      )

    } else {
      seurat_obj@meta.data[[col]] <- sapply(
        seurat_obj@meta.data[[col]],
        modify_string,
        USE.NAMES = FALSE
      )
    }
  }
  return(seurat_obj)
}

# Define resolutions for different references
resolution_for_100_plus_clusters = list()
resolution_for_100_plus_clusters['PBMC'] = 0.1
resolution_for_100_plus_clusters['adipose'] = 0.095
resolution_for_100_plus_clusters['tonsil'] = 0.075
resolution_for_100_plus_clusters['fetus'] = 0.03


reference_dataset_list = list()
for (reference in c("PBMC", "adipose", "tonsil", "fetus")) {
  ref_name <- paste0('human_', reference, '_integrated_')
  SetDirRead()
  reference_dataset_list[reference] = readRDS(paste0('cR with metadata reference_human_', reference, '_integrated_', resolution_for_100_plus_clusters[reference], '_', 0, '.rds'))
  seurat_equivalent_clusters = read.csv(paste0('seurat_equivalent_cluster_assignments_', ref_name, resolution_for_100_plus_clusters[reference], '_0.csv'), header = TRUE, sep = ",")
  seurat_equivalent_clusters = seurat_equivalent_clusters$cluster_assignment
  reference_dataset_list[[reference]] = AddMetaData(object = reference_dataset_list[[reference]], metadata = seurat_equivalent_clusters, col.name = 'seurat_equivalent_clusters')
  reference_dataset_list[[reference]] = modify_meta_data(reference_dataset_list[[reference]], "level_2")
  print(length(unique(reference_dataset_list[[reference]]@meta.data$level_2)))
  print(length(unique(reference_dataset_list[[reference]]@meta.data$seurat_equivalent_clusters)))
}

cell_shorthands_l2_tonsil <- list(
  "Mono/Macro" = "Mon/Mac",
  "DC" = "DC",
  "Granulocytes" = "Granul",
  "Mast" = "Mast",
  "cycling myeloid" = "cycl myeloid",
  "epithelial" = "Epi",
  "PDC" = "PDC",
  "FDC" = "FDC",
  "cycling FDC" = "cycling FDC",
  "NK" = "NK",
  "ILC" = "ILC",
  "Naive CD8 T" = "Naive CD8 T",
  "CD8 T" = "CD8 T",
  "DN" = "DN",
  "CD4 T" = "CD4 T",
  "Naive CD4 T" = "Naive CD4 T",
  "cycling T" = "cycling T",
  "PC" = "PC",
  "MBC" = "MBC",
  "NBC" = "NBC",
  "Activated NBC" = "Activated NBC",
  "GCBC" = "GCBC",
  "preB/T" = "preB/T"
)


cell_shorthands_l3_tonsil <- list(
  "ITGAX ZEB2 macrophages" = "ITGAXMac",
  "IL7R DC" = "IL7R_DC",
  "Neutrophil Granulocytes" = "NeuGran",
  "DC2" = "DC2",
  "Monocytes" = "Mono",
  "DC1 mature" = "DC1Mat",
  "SELENOP FUCA1 PTGDS macrophages" = "SELMac",
  "DC5" = "DC5",
  "IL7R MMP12 macrophages" = "IL7RMac",
  "Mast cells" = "MastC",
  "aDC1" = "aDC1",
  "C1Q HLA macrophages" = "C1QMac",
  "Cycling" = "Cycl",
  "aDC2" = "aDC2",
  "DC3" = "DC3",
  "M1 Macrophages" = "M1_Mac",
  "aDC3" = "aDC3",
  "DC1 precursor" = "DC1Pre",
  "DC4" = "DC4",
  "Crypt" = "Crypt",
  "Surface epithelium" = "SurfEpi",
  "FDCSP epithelium" = "FDCSPEpi",
  "Basal cells" = "Basal",
  "PDC" = "PDC",
  "IFN1+ PDC" = "IFN1PDC",
  "LZ FDC" = "LZFDC",
  "CD14+CD55+ FDC" = "CD14FDC",
  "DZ FDC" = "DZFDC",
  "MRC" = "MRC",
  "unknown" = "Unk",
  "cycling FDC" = "CyclFDC",
  "CD16-CD56- NK" = "16-56-NK",
  "CD16+CD56- NK" = "16+56-NK",
  "NKp44+ ILC3" = "NKp44+",
  "ILC1" = "ILC1",
  "NKp44- ILC3" = "NKp44-",
  "CD16-CD56+ NK" = "16-56+NK",
  "Naive CD8 T" = "NaiCD8T",
  "TCRVδ+ gd T" = "TCRgdT",
  "SCM CD8 T" = "SCMCD8T",
  "IFN CD8 T" = "IFNCD8T",
  "CXCR6+ RM CD8 T" = "CXCCD8T",
  "CM CD8 T" = "CMCD8T",
  "CD8 Tf" = "CD8Tf",
  "RM CD8 T" = "RMCD8T",
  "DC recruiters CD8 T" = "DCrCD8T",
  "MAIT" = "MAIT",
  "DN" = "DN",
  "CD56+ gd T" = "CD56gdT",
  "Nksig CD8 T" = "NksCD8T",
  "Tfh-LZ-GC" = "TfhLZGC",
  "Naive" = "Naive",
  "GC-Tfh-SAP" = "GCTfhS",
  "CM Pre-non-Tfh" = "CMPnTfh",
  "GC-Tfh-0X40" = "GCTX40",
  "CM PreTfh" = "CMP_Tfh",
  "T-Trans-Mem" = "TTransM",
  "T-Eff-Mem" = "TEffM",
  "GC-Tf-regs" = "GCTfRg",
  "Tfh-Mem" = "TfhMem",
  "T-helper" = "ThHelp",
  "Eff-Tregs" = "EffTrg",
  "Tfh T:B border" = "TfhTB",
  "non-GC-Tf-regs" = "nonGCT",
  "cycling T" = "CyclT",
  "IgG+ PC precursor" = "IgGPCP",
  "MBC derived early PC precursor" = "MBCePC",
  "DZ migratory PC precursor" = "DZmPCP",
  "PB" = "PB",
  "IgM+ PC precursor" = "IgMPCP",
  "MBC derived IgA+ PC" = "MBCIgA",
  "PB committed early PC precursor" = "PBcPCP",
  "preMature IgG+ PC" = "preIgG",
  "Mature IgM+ PC" = "MatIgM",
  "MBC derived IgG+ PC" = "MBCIgG",
  "Transitional PB" = "TransPB",
  "Mature IgG+ PC" = "MatIgG",
  "preMature IgM+ PC" = "preIgM",
  "Short lived IgM+ PC" = "ShortIgM",
  "Mature IgA+ PC" = "MatIgA",
  "IgM+ early PC precursor" = "IgMePC",
  "Early PC precursor" = "EarlyPC",
  "IgD PC precursor" = "IgDPCP",
  "MBC derived PC precursor" = "MBCPCP",
  "PC committed Light Zone GCBC" = "PCLZGC",
  "class switch MBC" = "cSwMBC",
  "Dark Zone GCBC" = "DZGCBC",
  "ncsMBC FCRL4+" = "ncsF4M",
  "NBC IFN-activated" = "NBCIFN",
  "GC-commited MIR155HG+" = "GCM155",
  "csMBC" = "csMBC",
  "NBC main" = "NBCMain",
  "ncsMBC" = "ncsMBC",
  "csMBC FCRL4+" = "csF4MBC",
  "GC-commited MIR155HG-" = "GCM155m",
  "preGC_1" = "preGC1",
  "NBC first step activation" = "NBCFSA",
  "NBC S100A+" = "NBCS100",
  "GC LZ Noproli" = "GCLZNp",
  "FCRL5+ MBC" = "FCRL5M",
  "Proliferative NBC" = "ProNBC",
  "preGC_2" = "preGC2",
  "Early GC-commited MYC+" = "EGCMYC",
  "GC DZ Noproli" = "GCDZNp",
  "NBC CD229+" = "NBCD229",
  "Early MBC" = "EarlyMBC",
  "GC-commited metabolic activation" = "GCMetA",
  "NBC/OA first step activation" = "NBCOA1",
  "DZ_G2M_HistoneHigh" = "DZG2MH",
  "DZ_Sphase" = "DZSph",
  "DZ/LZ" = "DZLZ",
  "DZ-nonproliferative_FOXP1hi" = "DZnpF1",
  "LZ" = "LZ",
  "LZ-BCL2A1 neg" = "LZBCLN",
  "PC-precursors" = "PCPre",
  "LZ-DZ-re-entry early commitment" = "LZDZEC",
  "DZ_Sphase_HistoneHigh" = "DZSphH",
  "DZ-cell cycle exit" = "DZCCEx",
  "DZ_G2M_CCNBHigh" = "DZG2MC",
  "LZ-proliferative_BCL2A1neg" = "LZpBCL",
  "LZ-proliferative_BCL2A1pos" = "LZpBCLP",
  "DZ-nonproliferative" = "DZNonP",
  "MBC-like_nonproli" = "MBCnp",
  "MBC-like_proli3" = "MBCp3",
  "MBC-like_proli1" = "MBCp1",
  "MBC-like_proli2" = "MBCp2",
  "MBC-like_FCRL4+" = "MBCF4",
  "Outer surface" = "OutSurf",
  "VEGFA+" = "VEGFA",
  "FRC" = "FRC",
  "preB" = "preB",
  "preT" = "preT",
  "myeloid_multiome" = "MyeloM",
  "PDC_multiome" = "PDcM",
  "FDC_multiome" = "FDcM",
  "Light Zone GCBC" = "LZGCBC",
  "epithelial_multiome" = "EpiM"
)


cell_shorthands_l1_fetus <- list(
  "Acinar cells" = "Acin",
  "Adrenocortical cells" = "AdrC",
  "AFP_ALB positive cells" = "AFP",
  "Amacrine cells" = "Amac",
  "Antigen presenting cells" = "APC",
  "Astrocytes" = "Astro",
  "Bipolar cells" = "Bipol",
  "Bronchiolar and alveolar epithelial cells" = "BrAl",
  "Cardiomyocytes" = "CMyo",
  "CCL19_CCL21 positive cells" = "CCL19",
  "Chromaffin cells" = "Chrom",
  "Ciliated epithelial cells" = "Cili",
  "CLC_IL5RA positive cells" = "CLC",
  "Corneal and conjunctival epithelial cells" = "CorCon",
  "CSH1_CSH2 positive cells" = "CSH1",
  "Ductal cells" = "Duct",
  "ELF3_AGBL2 positive cells" = "ELF3",
  "Endocardial cells" = "Endoc",
  "ENS glia" = "ENS_G",
  "ENS neurons" = "ENS_N",
  "Epicardial fat cells" = "EpiF",
  "Erythroblasts" = "Eryth",
  "Excitatory neurons" = "ExNe",
  "Extravillous trophoblasts" = "ExTro",
  "Ganglion cells" = "Gangl",
  "Goblet cells" = "Gob",
  "Granule neurons" = "Gran",
  "Hematopoietic stem cells" = "HSC",
  "Hepatoblasts" = "Hepat",
  "Horizontal cells" = "Horiz",
  "IGFBP1_DKK1 positive cells" = "IGFBP",
  "Inhibitory interneurons" = "InhInt",
  "Inhibitory neurons" = "InhNeu",
  "Intestinal epithelial cells" = "Intest",
  "Islet endocrine cells" = "Islet",
  "Lens fibre cells" = "Lens",
  "Limbic system neurons" = "LimbN",
  "Lymphatic endothelial cells" = "LymEnd",
  "Lymphoid cells" = "Lymph",
  "Megakaryocytes" = "Mega",
  "Mesangial cells" = "Mesang",
  "Mesothelial cells" = "Mesot",
  "Metanephric cells" = "Meta",
  "Microglia" = "Micro",
  "MUC13_DMBT1 positive cells" = "MUC13",
  "Myeloid cells" = "Myelo",
  "Neuroendocrine cells" = "NeuroE",
  "Oligodendrocytes" = "Oligo",
  "PAEP_MECOM positive cells" = "PAEP",
  "Parietal and chief cells" = "ParCh",
  "PDE11A_FAM19A2 positive cells" = "PDE11",
  "PDE1C_ACSM3 positive cells" = "PDE1C",
  "Photoreceptor cells" = "Phot",
  "Purkinje neurons" = "Purk",
  "Retinal pigment cells" = "RetPig",
  "Retinal progenitors and Muller glia" = "RetPro",
  "SATB2_LRRC7 positive cells" = "SATB2",
  "Satellite cells" = "Sat",
  "Schwann cells" = "Schw",
  "Skeletal muscle cells" = "SkelM",
  "SKOR2_NPSR1 positive cells" = "SKOR2",
  "SLC24A4_PEX5L positive cells" = "SLC24",
  "SLC26A4_PAEP positive cells" = "SLC26",
  "Smooth muscle cells" = "SmMus",
  "Squamous epithelial cells" = "Squam",
  "STC2_TLX1 positive cells" = "STC2",
  "Stellate cells" = "Stell",
  "Stromal cells" = "Strm",
  "Sympathoblasts" = "Symp",
  "Syncytiotrophoblasts and villous cytotrophoblasts" = "Sync",
  "Thymic epithelial cells" = "Thym",
  "Thymocytes" = "Thymo",
  "Trophoblast giant cells" = "TroG",
  "Unipolar brush cells" = "UnipB",
  "Ureteric bud cells" = "UretB",
  "Vascular endothelial cells" = "VasEnd",
  "Visceral neurons" = "VisN"
)

cell_shorthands_l2_fetus = list()
for (ctype in sort(unique(reference_dataset_list$fetus@meta.data$celltype.l2))) {
  ctype_vec = strsplit(ctype, "-")[[1]]
  cell_shorthands_l2_fetus[[ctype]] = paste0(substr(ctype_vec[[1]], 1, 2), '-', cell_shorthands_l1_fetus[[ctype_vec[[2]]]])
}

reference_dataset_list$PBMC@meta.data$celltype.short.l1 <- reference_dataset_list$PBMC@meta.data$celltype.l1
reference_dataset_list$PBMC@meta.data$celltype.short.l2 <- reference_dataset_list$PBMC@meta.data$celltype.l2
reference_dataset_list$adipose@meta.data$celltype.short.l1 <- reference_dataset_list$adipose@meta.data$celltype.l1
reference_dataset_list$adipose@meta.data$celltype.short.l2 <- reference_dataset_list$adipose@meta.data$celltype.l2
reference_dataset_list$tonsil@meta.data$celltype.short.l2 <- sapply(reference_dataset_list$tonsil@meta.data$celltype.l2, function(name) cell_shorthands_l2_tonsil[[name]])
reference_dataset_list$tonsil@meta.data$celltype.short.l3 <- sapply(reference_dataset_list$tonsil@meta.data$celltype.l3, function(name) cell_shorthands_l3_tonsil[[name]])
reference_dataset_list$fetus@meta.data$celltype.short.l1 <- sapply(reference_dataset_list$fetus@meta.data$celltype.l1, function(name) cell_shorthands_l1_fetus[[name]])
reference_dataset_list$fetus@meta.data$celltype.short.l2 <- sapply(reference_dataset_list$fetus@meta.data$celltype.l2, function(name) cell_shorthands_l2_fetus[[name]])

SetDirWrite()
for (reference in c("PBMC", "adipose", "tonsil", "fetus")) {
  annotation_levels = list()
  annotation_levels[[2]] = if (reference == "tonsil") 3 else 2
  annotation_levels[[1]] = annotation_levels[[2]] - 1
  do_label = TRUE
  reference_dataset_list[[reference]] <- GetClusterPlot(reference_dataset_list[[reference]], '', group.by = "seurat_equivalent_clusters", shuffle_colors = TRUE, label = do_label)
  reference_dataset_list[[reference]] <- GetClusterPlot(reference_dataset_list[[reference]], '', group.by = "level_2", shuffle_colors = TRUE, label = do_label)
  reference_dataset_list[[reference]] <- GetClusterPlot(reference_dataset_list[[reference]], '', group.by = paste0('celltype.short.l', annotation_levels[[1]]), shuffle_colors = TRUE, label = do_label)
  reference_dataset_list[[reference]] <- GetClusterPlot(reference_dataset_list[[reference]], '', group.by = paste0('celltype.short.l', annotation_levels[[2]]), shuffle_colors = TRUE, label = do_label)


  # Retrieve the four plots from the misc slot
  plot_list = list()
  plot_list[['plot_level_2']] <- reference_dataset_list[[reference]]@misc$level_2_UMAP_ + ggtitle('(a) Recursive Clusters')
  plot_list[['plot_seurat_eq']] <- reference_dataset_list[[reference]]@misc$seurat_equivalent_clusters_UMAP_ + ggtitle('(b) Single-Pass Clusters')
  plot_list[['plot_celltype_l1']] <- reference_dataset_list[[reference]]@misc[[paste0('celltype.short.l', annotation_levels[[1]], '_UMAP_')]] + ggtitle('(c) Cell Type Annotational Level 1')
  plot_list[['plot_celltype_l2']] <- reference_dataset_list[[reference]]@misc[[paste0('celltype.short.l', annotation_levels[[2]], '_UMAP_')]] + ggtitle('(d) Cell Type Annotational Level 2')

  for (plot_name in names(plot_list)) {
    plot_list[[plot_name]] = plot_list[[plot_name]] +
      theme(plot.title = element_text(size = 12)) +
      theme(
        axis.title = element_blank(),       # Remove axis titles
        axis.text = element_blank(),        # Remove axis tick numbers
        axis.ticks = element_blank())       # Remove axis ticks
  }


  reference_title = if (reference == "PBMC") reference  else (paste0(toupper(substr(reference, 1, 1)), substr(reference, 2, nchar(reference))))
  # Arrange the four plots in a 2x2 grid
  combined_plot = grid.arrange(
    plot_list[['plot_level_2']], plot_list[['plot_seurat_eq']],
    plot_list[['plot_celltype_l1']], plot_list[['plot_celltype_l2']],
    ncol = 2, nrow = 2,
    top = textGrob(paste("UMAP Plots for", reference_title), gp = gpar(fontsize = 14, fontface = "bold"))
  )

  output_filename = paste0("UMAP_plots_", reference)
  ggsave(
    filename = paste0(output_filename, ".png"),
    plot = combined_plot,
    width = 8.5,
    height = 11,
    dpi = 300,
    limitsize = FALSE
  )

  pp1 <- image_trim(image_read(paste0(output_filename, '.png')))
  image_write(pp1, paste0(output_filename, '.pdf'), format = 'pdf')
  # pp1 <- image_scale(pp1, "50%")
  # image_write(pp1, paste0(output_filename, '.eps'), format = 'eps')
}
