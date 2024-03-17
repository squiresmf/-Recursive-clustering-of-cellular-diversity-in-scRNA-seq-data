# The main code for applying recursive clustering and single-pass equivalent clustering
# To run, choose a recursive clustering resolution parameter c_resolution and a dataset reference_name

#####################################################################################
c_resolution = 0.02 # Recursive Clustering resolution parameter

# Specify which reference data to use
reference_name = "human_PBMC/"
# reference_name = "human_adipose/"
# reference_name = "human_tonsil/"
# reference_name = "human_fetus/"

is_helmsley = FALSE #If FALSE, use a reference dataset, if TRUE, use the Chron's disease dataset
#####################################################################################
do_integrate = FALSE
assay = 'integrated'
#####################################################################################
# These are the all libraries needed to run all code across all files
library(Matrix)
library(Seurat)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(assertthat)
library(stringr)
library(colorspace)
library(dplyr)
library(data.tree)
library(networkD3)
library(webshot)
library(grid)
library(magick)
library("png")
library(stringi)
library(SeuratDisk)

# Other parameters
DEG_cutoff = 0
if (is_helmsley) {
  entropy_cutoff = 0.5
} else {
  entropy_cutoff = 0
}

do_lower_clusters_higher_UMAP = FALSE
m_level = 3
use_annotations=FALSE

options(future.globals.maxSize = 20000 * 1024^2)

# A function to make strings the same length by padding with spaces
PadRight <- function(value, width, side = "right", pad = " ") {
  if (is.numeric(value)) {
    if (is.nan(value)) {
      return("NaN")
    }
    if (value - floor(value) > 0) {
      value = round(value, digits = width - 2)
    }
    if (value - floor(value) == 0) {
      pad = " "
    }
    string = paste0(value)
    string = substring(string, first = 1, last = width)
    if (grepl(pattern = 'e', x = string, fixed = TRUE)) {
      pad = " "
    }
  } else {
    string = paste0(value)
  }
  string = str_pad(string, width = width, side = side, pad = pad)
  return(string)
}

paths = list()
data = list()

if (is_helmsley) {
  paths[['GCA1']] = 'raw_count_data/first_batch_13_libraries_CellRanger_filtered/GCA1_prot_01788_5GEX_C7/'
  paths[['GCA2']] = 'raw_count_data/first_batch_13_libraries_CellRanger_filtered/GCA2_prot_21394_5GEX_C5/'
  paths[['GCA3']] = 'raw_count_data/first_batch_13_libraries_CellRanger_filtered/GCA3_prot_21460_5GEX_D5/'
  paths[['GCA5']] = 'raw_count_data/first_batch_13_libraries_CellRanger_filtered/GCA5_prot_21478_5GEX_H5/'
  paths[['GCA7']] = 'raw_count_data/second_batch_8_libraries_CellRanger_filtered/GCA7/'
  paths[['GCA8']] = 'raw_count_data/second_batch_8_libraries_CellRanger_filtered/GCA8/'
  paths[['GCA9']] = 'raw_count_data/second_batch_8_libraries_CellRanger_filtered/GCA9/'
  paths[['GCA10']] = 'raw_count_data/second_batch_8_libraries_CellRanger_filtered/GCA10/'
  paths[['GCA11']] = 'raw_count_data/second_batch_8_libraries_CellRanger_filtered/GCA11/'
  paths[['GCA12']] = 'raw_count_data/second_batch_8_libraries_CellRanger_filtered/GCA12/'
  paths[['GCA13']] = 'raw_count_data/second_batch_8_libraries_CellRanger_filtered/GCA13/'
  paths[['GCA14']] = 'raw_count_data/second_batch_8_libraries_CellRanger_filtered/GCA14/'
  paths[['GCA15']] = 'raw_count_data/third_batch_8_libraries_CellRanger_filtered/GCA15/filtered_feature_bc_matrix/'
  paths[['GCA16']] = 'raw_count_data/third_batch_8_libraries_CellRanger_filtered/GCA16/filtered_feature_bc_matrix/'
  paths[['GCA17']] = 'raw_count_data/third_batch_8_libraries_CellRanger_filtered/GCA17/filtered_feature_bc_matrix/'
  paths[['GCA18']] = 'raw_count_data/third_batch_8_libraries_CellRanger_filtered/GCA18/filtered_feature_bc_matrix/'
  paths[['GCA19']] = 'raw_count_data/third_batch_8_libraries_CellRanger_filtered/GCA19/filtered_feature_bc_matrix/'
  paths[['GCA20']] = 'raw_count_data/third_batch_8_libraries_CellRanger_filtered/GCA20/filtered_feature_bc_matrix/'
  paths[['GCA21']] = 'raw_count_data/third_batch_8_libraries_CellRanger_filtered/GCA21/filtered_feature_bc_matrix/'
  paths[['GCA22']] = 'raw_count_data/third_batch_8_libraries_CellRanger_filtered/GCA22/filtered_feature_bc_matrix/'
  paths[['GCA23']] = 'raw_count_data/batch_4/Helmsley Batch 4/p21063-s001_GCA23_5GEX_C12/filtered_feature_bc_matrix/'
  paths[['GCA24']] = 'raw_count_data/batch_4/Helmsley Batch 4/p21063-s002_GCA24_5GEX_D12/filtered_feature_bc_matrix/'
  paths[['GCA25']] = 'raw_count_data/batch_4/Helmsley Batch 4/p21063-s003_GCA25_5GEX_E12/filtered_feature_bc_matrix/'
  paths[['GCA26']] = 'raw_count_data/batch_4/Helmsley Batch 4/p21063-s004_GCA26_5GEX_F12/filtered_feature_bc_matrix/'
  paths[['GCA27']] = 'raw_count_data/batch_4/Helmsley Batch 4/p21063-s005_GCA27_5GEX_G12/filtered_feature_bc_matrix/'
  paths[['GCA28']] = 'raw_count_data/batch_4/Helmsley Batch 4/p21063-s006_GCA28_5GEX_H12/filtered_feature_bc_matrix/'
  paths[['GCA29']] = 'raw_count_data/batch_4/Helmsley Batch 4/p21063-s007_GCA29_5GEX_E2/filtered_feature_bc_matrix/'
}
disease_phenotype_names = c("control", "treatment naïve-CD", "CD")

if (is_helmsley) {
  patient_metadata <- read.csv(file = 'patient_metadata helmsley.csv', fileEncoding = "UTF-8-BOM", check.names = FALSE)
  row.names(patient_metadata) <- patient_metadata$'Sample ID'
  data_filename = 'data.rds'
  ref_name = paste0('helmsley', '_', assay, '_')
} else {
  ref_name = paste0(strsplit(reference_name, '/')[[1]], '_', assay, '_')
  data_filename = paste0(strsplit(reference_name, '/')[[1]], '_', 'RNA', '_', 'ref_data', '.rds')
}

if (!is_helmsley) {
  if (reference_name == "human_PBMC/") {
    annotation_levels = 3
  }
  if (reference_name == "human_tonsil/") {
    annotation_levels = 3
  }
  if (reference_name == "human_fetus/") {
    annotation_levels = 2
  }
  if (reference_name == "human_adipose/") {
    annotation_levels = 2
  }
}


# Create the initial reference dataset seurat objects from output of Azimuth snakemake files
if (!file.exists(data_filename)) {
  if (is_helmsley) {
    patient_ids = names(paths)
    for (patient_id in patient_ids) {
      print(patient_id)
      data[[patient_id]] <- CreateSeuratObject(counts = Read10X(data.dir = paths[[patient_id]], gene.column = 2))
    }
  } else {
    if (reference_name == "human_PBMC/") {
      data <- LoadH5Seurat("azimuth-references/human_pbmc/data/pbmc_multimodal.h5seurat")
      temp = CreateSeuratObject(data@assays$SCT@counts, assay = "RNA")
      data@assays$RNA = temp@assays$RNA
      data = AddMetaData(object = data, metadata = data@meta.data$orig.ident, col.name = 'patient')
      for (metadata_name in names(data@meta.data)) {
        if (!(metadata_name %in% c('patient', 'nCount_RNA', 'nFeature_RNA', 'celltype.l1', 'celltype.l2', 'celltype.l3'))) {
          data@meta.data[[metadata_name]] = NULL
        }
      }
    }
    if (reference_name == "human_tonsil/") {
      data2 <- readRDS('azimuth-references/human_tonsil/data/obj.rds')
      data = AddMetaData(object = data, metadata = data@meta.data$donor_id, col.name = 'patient')
      data = AddMetaData(object = data, metadata = data@meta.data$annotation_level_1, col.name = 'celltype.l1')
      data = AddMetaData(object = data, metadata = data@meta.data$annotation_figure_1, col.name = 'celltype.l2')
      data = AddMetaData(object = data, metadata = data@meta.data$annotation_20220215, col.name = 'celltype.l3')
      for (metadata_name in names(data@meta.data)) {
        if (!(metadata_name %in% c('patient', 'nCount_RNA', 'nFeature_RNA', 'celltype.l1', 'celltype.l2', 'celltype.l3'))) {
          data@meta.data[[metadata_name]] = NULL
        }
      }
    }
    if (reference_name == "human_fetus/") {
      data3 <- readRDS('azimuth-references/human_fetus/seurat_objects/fullref.Rds')
      data = AddMetaData(object = data, metadata = data@meta.data$Fetus_id, col.name = 'patient')
      data = AddMetaData(object = data, metadata = data@meta.data$annotation.l1, col.name = 'celltype.l1')
      data = AddMetaData(object = data, metadata = data@meta.data$annotation.l2, col.name = 'celltype.l2')
      for (metadata_name in names(data@meta.data)) {
        if (!(metadata_name %in% c('patient', 'nCount_RNA', 'nFeature_RNA', 'celltype.l1', 'celltype.l2', 'celltype.l3'))) {
          data@meta.data[[metadata_name]] = NULL
        }
      }
    }
    if (reference_name == "human_adipose/") {
      data4 <- readRDS('azimuth-references/human_adipose/full_reference.Rds')
      data = AddMetaData(object = data, metadata = data@meta.data$Number, col.name = 'patient')
      for (metadata_name in names(data@meta.data)) {
        if (!(metadata_name %in% c('patient', 'nCount_RNA', 'nFeature_RNA', 'celltype.l1', 'celltype.l2'))) {
          data@meta.data[[metadata_name]] = NULL
        }
      }
    }
    for (a in names(data@assays)) {
      if (a != 'RNA') {
        print(paste0('removing ', a))
        data@assays[[a]] = NULL
      }
    }
    DefaultAssay(data) = 'RNA'
    data@graphs = list()
    data@reductions = list()
    data@commands = list()
    data <- SplitObject(object = data, split.by = "patient")
    patient_ids = names(data)
  }

  for (patient_id in patient_ids) {
    cell_id = data[[patient_id]][['RNA']]@counts@Dimnames[[2]]
    cell_num = seq(ncol(data[[patient_id]]))
    data[[patient_id]] = AddMetaData(object = data[[patient_id]], metadata = cell_id, col.name = 'cell_id')
    data[[patient_id]] = AddMetaData(object = data[[patient_id]], metadata = cell_num, col.name = 'cell_num')

    data[[patient_id]] <- NormalizeData(data[[patient_id]])
    data[[patient_id]] <- FindVariableFeatures(data[[patient_id]])

    if (is_helmsley) {
      # add patient metadata
      patient = rep(patient_metadata[patient_id,]$'Sample ID', ncol(data[[patient_id]]))
      disease_phenotype = rep(patient_metadata[patient_id,]$'Disease phenotype', ncol(data[[patient_id]]))
      gender = rep(patient_metadata[patient_id,]$Sex, ncol(data[[patient_id]]))
      batch = rep(patient_metadata[patient_id,]$'Batch #', ncol(data[[patient_id]]))
      age = rep(patient_metadata[patient_id,]$Age, ncol(data[[patient_id]]))
      race = rep(patient_metadata[patient_id,]$Race, ncol(data[[patient_id]]))

      data[[patient_id]] = AddMetaData(object = data[[patient_id]], metadata = disease_phenotype, col.name = 'disease_phenotype')
      data[[patient_id]] = AddMetaData(object = data[[patient_id]], metadata = gender, col.name = 'gender')
      data[[patient_id]] = AddMetaData(object = data[[patient_id]], metadata = batch, col.name = 'batch')
      data[[patient_id]] = AddMetaData(object = data[[patient_id]], metadata = age, col.name = 'age')
      data[[patient_id]] = AddMetaData(object = data[[patient_id]], metadata = race, col.name = 'race')
      data[[patient_id]] = AddMetaData(object = data[[patient_id]], metadata = patient, col.name = 'patient')
    }
  }
  saveRDS(data, file = data_filename)
} else {
  data = readRDS(data_filename)
}


patient_ids = names(data)
for (patient_id in patient_ids) {
  DefaultAssay(data[[patient_id]]) = 'RNA'
}

# The seurat clustering pipeline:
# - Center and scale the raw single cell RNA sequencing counts matrix
# - Perform principal components analysis (PCA) to reduce data dimensions
# - Generate a shared nearest neighbor (SNN) graph based on the Euclidean distances of PCA embeddings of each cell
# - Apply Louvain clustering on the SNN graph.

FindSeuratClusters <- function(dataset, assay, resolution, rerun, random_seed = 0, savetime = FALSE, cluster_name = NULL, method_name = NULL) {
  if (savetime) {
    start.time <- Sys.time()
  }
  if (length(dataset[[assay]]@scale.data) <= 1 | rerun) {
    print('scaling data')
    dataset <- ScaleData(dataset, assay = assay, verbose = FALSE)
  }
  if (length(dataset@reductions$pca) == 0 | rerun) {
    print('running pca')
    dataset = RunPCA(dataset, assay = assay, verbose = FALSE)
  }
  if (length(dataset@reductions$umap) == 0 | rerun) {
    print('running umap')
    dataset <- RunUMAP(dataset, reduction = "pca", dims = 1:30, assay = assay, verbose = FALSE)
  }
  if (length(dataset@graphs[[paste0(assay, '_nn')]]) == 0 | rerun) {
    print('finding neighbors')
    dataset <- FindNeighbors(dataset, reduction = "pca", dims = 1:30, assay = assay, verbose = FALSE)
  }
  print('finding clusters')
  dataset <- FindClusters(dataset, resolution = resolution, verbose = FALSE, random.seed = random_seed)

  if (savetime) {
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    num_clusters = suppressWarnings(length(as.numeric(levels(as.factor(dataset@meta.data[["seurat_clusters"]])))))
    time_list = list(time = time.taken, ref_name = ref_name, assay = assay, cluster_name = cluster_name, resolution = resolution, method_name = method_name, num_clusters = num_clusters, random_seed = random_seed)
    saveRDS(time_list, paste0(paste(ref_name, assay, cluster_name, resolution, method_name, random_seed, num_clusters), '.rds'))
  }
  return(dataset)
}

# A helper function to reduce repetitive code for generating seurat clustering plots
GetClusterPlot <- function(dataset, title, group.by = "seurat_clusters", no_legend = TRUE, shuffle_colors = FALSE) {
  #print('getting cluster stats')
  if (length(group.by) > 1) {
    plot_title = paste0(title, ' ', group.by[[length(group.by)]])
  } else {
    plot_title = paste0(title, ' ', group.by)
  }
  if (shuffle_colors) {
    num_clusters = suppressWarnings(length(as.numeric(levels(as.factor(dataset@meta.data[[group.by]])))))
    temp_plot = DimPlot(dataset, group.by = group.by, label = no_legend)
    plot_colors = ggplot_build(temp_plot)$data[[1]]$colour
    plot_cluster_nums = ggplot_build(temp_plot)$data[[1]]$group
    new_colors = sample(plot_colors[match(seq(num_clusters), plot_cluster_nums)])
    seurat_clustering_umap_plot = DimPlot(dataset, group.by = group.by, label = no_legend, cols = new_colors) + ggtitle(plot_title)
  } else {
    seurat_clustering_umap_plot = DimPlot(dataset, group.by = group.by, label = no_legend) + ggtitle(plot_title)
  }
  if (no_legend) {
    seurat_clustering_umap_plot = seurat_clustering_umap_plot + NoLegend()
  }
  #print(paste0(paste0(group.by, collapse = '_'), '_UMAP'))
  dataset@misc[[paste0(paste0(group.by, collapse = '_'), '_UMAP_', title)]] = seurat_clustering_umap_plot
  return(dataset)
}


# Calculates shannon entropy of a distribution
GetEntropy <- function(frac_list) {
  shannon_entropy = 0
  for (i in seq(length(frac_list))) {
    if (frac_list[[i]] != 0) {
      shannon_entropy = shannon_entropy - frac_list[[i]] * log(frac_list[[i]], 2)
    }
  }
  return(shannon_entropy / log(length(frac_list), 2))
}

# This is the recursive clustering function.  It inlcudes code for recursively applying the default seurat pipeline,
# as well as recursively applying the default seurat pipeline alongside recursively re-doing seurat integration.
# There is also code for using number of differentially expressed genes between subclusters as a cluster merging criterion.
# And there is an option to include the entropy cutoff criterion so that computation time is reduced for clustering
# where we only care about high-entropy clusters and subclusters.
RecursiveClustering <- function(parent_data, cluster_resolution, max_level, do_integrate, assay, DEG_cutoff, entropy_cutoff = 0, cluster_num = 'cR', cluster_sequence = 'cR', parent_integrated = NULL, level = 0, min_cells = 50) {
  print(cluster_sequence)
  if (level > 0) {
    Idents(parent_integrated[[level]]) <- "seurat_clusters"
    parent_integrated_cluster_subset = subset(parent_integrated[[level]], idents = cluster_num)
    Idents(parent_integrated_cluster_subset) <- "patient"

    # Getting a patient-wise list of cells present in a cluster for the purpose of performing integration recursively
    if (do_integrate) {
      child_data = vector(mode = 'list')
      npcs = 50
      for (patient_id in names(parent_data)) {
        #print(patient_id)
        if (patient_id %in% parent_integrated_cluster_subset@meta.data$patient) {
          parent_integrated_cluster_subset_patient_subset = subset(parent_integrated_cluster_subset, idents = patient_id)
          temp_child_data = parent_data[[patient_id]]
          Idents(temp_child_data) <- "cell_id"
          temp_child_data = subset(temp_child_data, idents = parent_integrated_cluster_subset_patient_subset@meta.data$cell_id)
          if (temp_child_data[[assay]]@counts@Dim[[2]] > 30) {
            child_data[[patient_id]] = temp_child_data
            npcs = min(npcs, child_data[[patient_id]][[assay]]@counts@Dim[[2]] - 1)
          }
        }
      }
    }
  } else if (do_integrate) {
    child_data = parent_data
    npcs = 50
    parent_integrated = vector(mode = 'list')
  }
  if (is_helmsley) {
    integrated_name = paste0("helmsley_", assay, '_', cluster_resolution, '_', DEG_cutoff, '_', entropy_cutoff, '_', cluster_sequence, ".rds")
  } else {
    integrated_name = paste0("reference_", ref_name, cluster_resolution, '_', DEG_cutoff, '_', cluster_sequence, ".rds")
  }
  if (!file.exists(integrated_name)) {
    # Code for re-performing seurat integration pipeline on subclusters
    if (do_integrate) {
      if (length(child_data) < 2) {
        #out = vector(mode = 'list')
        #out[[cluster_sequence]] = parent_integrated_cluster_subset
        #return(out)
        return()
      }

      child_integration_features <- SelectIntegrationFeatures(object.list = child_data)
      if (assay == 'SCT') {
        child_data <- PrepSCTIntegration(object.list = child_data, anchor.features = child_integration_features)
        normalization.method = "SCT"
      } else {
        normalization.method = "LogNormalize"
      }
      for (patient_id in names(child_data)) {
        child_data[[patient_id]] <- ScaleData(child_data[[patient_id]], features = child_integration_features, assay = assay, verbose = FALSE)
        child_data[[patient_id]] <- RunPCA(child_data[[patient_id]], features = child_integration_features, assay = assay, npcs = npcs, verbose = FALSE)
      }
      child_anchors <- FindIntegrationAnchors(object.list = child_data, anchor.features = child_integration_features, normalization.method = normalization.method, reduction = 'rpca', verbose = FALSE)
      for (k.weight in seq(100, 10, -5)) {
        child_integrated = tryCatch(
          expr = {
            IntegrateData(anchorset = child_anchors, features = child_integration_features, verbose = FALSE, normalization.method = normalization.method, k.weight = k.weight)
          },
          error = function(e) {
            print(paste0("k.weight value of ", k.weight, " may be too large, let's try a smaller value"))
            return(paste0(e))
          }
        )
        if (is.string(child_integrated)) {
          next
        }
        break
      }
      if (is.string(child_integrated)) {
        print(child_integrated)
        print('Integration failed.  Returning as single cluster.')
        #out = vector(mode = 'list')
        #out[[cluster_sequence]] = parent_integrated_cluster_subset
        #return(out)
        return()
      }
    } else {
      # The seurat clustering pipeline:
      # - Center and scale the raw single cell RNA sequencing counts matrix
      # - Perform principal components analysis (PCA) to reduce data dimensions
      # - Generate a shared nearest neighbor (SNN) graph based on the Euclidean distances of PCA embeddings of each cell
      # - Apply Louvain clustering on the SNN graph.
      if (level == 0) {
        child_integrated = FindVariableFeatures(parent_data)
        child_integrated <- FindSeuratClusters(dataset = child_integrated, assay = assay, resolution = cluster_resolution, rerun = TRUE, savetime = TRUE, cluster_name = cluster_sequence, method_name = 'recursive')
        # Number of cells is used as a stopping criterion for recursive clustering
      } else if (parent_integrated_cluster_subset[[assay]]@
        data@
        Dim[[2]] > min_cells) {
        child_integrated = FindVariableFeatures(parent_integrated_cluster_subset)
        child_integrated <- FindSeuratClusters(dataset = child_integrated, assay = assay, resolution = cluster_resolution, rerun = TRUE, savetime = TRUE, cluster_name = cluster_sequence, method_name = 'recursive')
      } else {
        #out = vector(mode = 'list')
        #out[[cluster_sequence]] = parent_integrated_cluster_subset
        #return(out)
        return()
      }
    }

    # code for using number of DEG between child clusters and min number of cells as a cluster merging criterion.  We
    # Assess the DEGs between all pairs of clusters and merge the pair with the fewest DEGs if it is below the cutoff.
    # We repeat this process until all pairwise DEGs are above the cutoff.  Then we repeat this process for clusters
    # with fewer than min_cells, merging small clusters with the lowest min DEGs across all other child clusters first.
    if ((level > 0) & (DEG_cutoff > 0)) {
      while (TRUE) {
        merge = FALSE
        Idents(child_integrated) <- "seurat_clusters"
        cluster_nums = as.numeric(levels(child_integrated@meta.data$seurat_clusters))
        c1s = vector()
        c2s = vector()
        nDEGs = vector()
        lowest_DE = list(nDEG = DEG_cutoff, cluster1 = NULL, cluster2 = NULL)
        num_clusters = length(cluster_nums)
        for (cluster1 in cluster_nums) {
          for (cluster2 in cluster_nums) {
            if (cluster2 > cluster1) {
              markers = FindMarkers(child_integrated, ident.1 = cluster1, ident.2 = cluster2, verbose = FALSE, min.cells.feature = 1, min.cells.group = 1, assay = 'integrated')
              markers = markers[markers$p_val_adj < 0.05,]
              nDEG = nrow(markers)
              c1s[[length(c1s) + 1]] = cluster1
              c2s[[length(c2s) + 1]] = cluster2
              nDEGs[[length(nDEGs) + 1]] = nDEG
              if (nDEG < lowest_DE$nDEG) {
                lowest_DE$nDEG = nDEG
                lowest_DE$cluster1 = cluster1
                lowest_DE$cluster2 = cluster2
              }
              print(paste0('c', cluster1, ' c', cluster2, ' ', nDEG))
            }
          }
        }
        if (lowest_DE$nDEG < DEG_cutoff) {
          merge = TRUE
          merge_c1 = lowest_DE$cluster1
          merge_c2 = lowest_DE$cluster2
        } else {
          for (i in order(nDEGs)) {
            if (length(child_integrated@meta.data$seurat_clusters[child_integrated@meta.data$seurat_clusters == c2s[[i]]]) < min_cells) {
              merge = TRUE
              merge_c1 = c1s[[i]]
              merge_c2 = c2s[[i]]
              break
            }
          }
        }
        if (merge) {
          child_integrated@meta.data$seurat_clusters[child_integrated@meta.data$seurat_clusters == merge_c2] = merge_c1
          if (merge_c2 < max(cluster_nums)) {
            for (cluster in seq(merge_c2 + 1, max(cluster_nums))) {
              child_integrated@meta.data$seurat_clusters[child_integrated@meta.data$seurat_clusters == cluster] = cluster - 1
            }
          }
          child_integrated@meta.data$seurat_clusters = factor(as.numeric(child_integrated@meta.data$seurat_clusters) - 1)
        }
        if (length(levels(child_integrated@meta.data$seurat_clusters)) == num_clusters) {
          break
        }
      }
    }
    subset_cell_idx = seq(child_integrated[[assay]]@data@Dim[[2]])
    child_integrated = AddMetaData(object = child_integrated, metadata = subset_cell_idx, col.name = cluster_sequence)
    saveRDS(child_integrated, file = integrated_name)
  } else {
    child_integrated <- readRDS(integrated_name)
  }
  # Only for recursively integrating, we map the cells in the integrated dataset back to the patient-wise list of
  # data so we can reperform integration on the child clusters in the next round of clustering
  if (do_integrate) {
    subset_cell_idx = seq(child_integrated[[assay]]@data@Dim[[2]])
    start_idx = 1
    for (patient_id in names(child_data)) {
      num_cells = child_data[[patient_id]][[assay]]@data@Dim[[2]]
      child_data[[patient_id]] <- AddMetaData(object = child_data[[patient_id]], metadata = subset_cell_idx[start_idx:(start_idx + num_cells - 1)], col.name = cluster_sequence)
      start_idx = start_idx + num_cells
    }
  } else {
    child_data = NULL
  }
  #out = vector(mode = 'list')
  #out[[cluster_sequence]] = child_integrated

  # We can include entropy cutoff criterion so that we don't find subclusters of low entropy clusters in case we aren't
  # interested in low entropy clusters
  if ((level > 0) & (entropy_cutoff > 0)) {
    cluster_nums_entropy = as.numeric(levels(child_integrated@meta.data$seurat_clusters))
    cluster_nums = list()

    for (cluster1 in cluster_nums_entropy) {
      cluster_patient_fracs = list()
      cluster_idxs = seq(child_integrated$RNA@data@Dim[[2]])[child_integrated@meta.data$seurat_clusters == cluster1]
      cluster_patients = child_integrated@meta.data$patient[cluster_idxs]

      for (patient_name in names(paths)) {
        cluster_patient_fracs[[patient_name]] = length(cluster_patients[cluster_patients == patient_name]) / length(cluster_patients)
      }

      patient_entropy = GetEntropy(cluster_patient_fracs)
      print(paste(cluster1, patient_entropy))
      if (patient_entropy > entropy_cutoff) {
        cluster_nums = c(cluster_nums, cluster1)
      }
    }
  } else {
    cluster_nums = as.numeric(levels(child_integrated@meta.data$seurat_clusters))
  }

  if (length(as.numeric(levels(child_integrated@meta.data$seurat_clusters))) == 1 |
    length(cluster_nums) == 0 |
    level == max_level) {
    #return(out)
    return()
  } else {
    parent_integrated[[level + 1]] = child_integrated
    print(cluster_nums)
    for (cluster_num in cluster_nums) {
      #out[[paste0(cluster_sequence, 'c', cluster_num)]] = RecursiveClustering(parent_data = child_data, cluster_resolution = cluster_resolution, cluster_num = cluster_num, cluster_sequence = paste0(cluster_sequence, 'c', cluster_num), parent_integrated = parent_integrated, level = level + 1, max_level = max_level, do_integrate = do_integrate, min_cells = min_cells, DEG_cutoff = DEG_cutoff, entropy_cutoff = entropy_cutoff, assay = assay)
      RecursiveClustering(parent_data = child_data, cluster_resolution = cluster_resolution, cluster_num = cluster_num, cluster_sequence = paste0(cluster_sequence, 'c', cluster_num), parent_integrated = parent_integrated, level = level + 1, max_level = max_level, do_integrate = do_integrate, min_cells = min_cells, DEG_cutoff = DEG_cutoff, entropy_cutoff = entropy_cutoff, assay = assay)
    }
  }
  #return(out)
  return()
}

if (is_helmsley) {
  prefix = paste0("helmsley_", assay, '_', c_resolution, '_', DEG_cutoff, '_', entropy_cutoff)
} else {
  prefix = paste0("reference_", ref_name, c_resolution, '_', DEG_cutoff)
}
done_title = paste0('done integrating ', prefix, '.rds')

# We need to have the integrated dataset
if (assay == 'integrated') {
  if (!is_helmsley) {
    integrated_data_filename = paste0('reference_', strsplit(reference_name, '/')[[1]], '_integrated.rds')
  } else {
    integrated_data_filename = paste0("helmsley_integrated.rds")
  }
}

if (!file.exists(integrated_data_filename)) {
  ref_name = paste0(strsplit(reference_name, '/')[[1]], '_', 'RNA', '_')
  RecursiveClustering(parent_data = data, cluster_resolution = c_resolution, max_level = 0, do_integrate = TRUE, DEG_cutoff = DEG_cutoff, entropy_cutoff = entropy_cutoff, assay = 'RNA')
  if (is_helmsley) {
  integrated_name = paste0("helmsley_", 'RNA', '_', c_resolution, '_', DEG_cutoff, '_', entropy_cutoff, '_', 'cR', ".rds")
  } else {
    integrated_name = paste0("reference_", paste0(strsplit(reference_name, '/')[[1]], '_', 'RNA', '_'), c_resolution, '_', DEG_cutoff, '_', 'cR', ".rds")
  }
  file.rename(integrated_name, integrated_data_filename)
  ref_name = paste0(strsplit(reference_name, '/')[[1]], '_', assay, '_')
}

integrated_data = readRDS(integrated_data_filename)

if (!is_helmsley) {
  print(paste('resolution', c_resolution, 'DEG cutoff', DEG_cutoff, sep = ' '))
}

if (!file.exists(done_title)) {
  if (is_helmsley) {
    RecursiveClustering(parent_data = integrated_data, cluster_resolution = c_resolution, max_level = m_level, do_integrate = do_integrate, DEG_cutoff = DEG_cutoff, entropy_cutoff = entropy_cutoff, assay = assay)
  } else {
    RecursiveClustering(parent_data = integrated_data, cluster_resolution = c_resolution, max_level = m_level, do_integrate = do_integrate, DEG_cutoff = DEG_cutoff, assay = assay)
  }
}

done = 'done'
saveRDS(done, done_title)


# We read in all saved seurat object for each cluster in the heirarchy which has subclusters into a flat list
# This way we don't need to re-run the recursive clustering pipeline once the pipeline has finished.
files <- (Sys.glob(paste0(prefix, "*.rds")))
recursive_integrated = vector(mode = 'list')
recursive_integrated_extra = vector(mode = 'list')
for (file in files) {
  dataset_name = substr(file, nchar(prefix) + 2, nchar(file) - 4)
  print(dataset_name)
  temp_integrated = readRDS(file)
  cluster_nums = as.numeric(levels(temp_integrated@meta.data$seurat_clusters))
  num_clusters = length(cluster_nums)
  if (num_clusters > 1) {
    recursive_integrated[[dataset_name]] = temp_integrated
  } else {
    recursive_integrated_extra[[dataset_name]] = temp_integrated
  }
}

# Some basic plotting for the full dataset and each cluster in the heirarchy which has subclusters
for (dataset_name in names(recursive_integrated)) {
  print(dataset_name)
  if (is_helmsley | do_lower_clusters_higher_UMAP) {
    recursive_integrated[[dataset_name]] <- GetClusterPlot(recursive_integrated[[dataset_name]], title = paste0(dataset_name, ' ', c_resolution), no_legend = FALSE)
  }
  if (is_helmsley) {
    recursive_integrated[[dataset_name]] <- GetClusterPlot(recursive_integrated[[dataset_name]], dataset_name, group.by = "batch", no_legend = FALSE)
    recursive_integrated[[dataset_name]] <- GetClusterPlot(recursive_integrated[[dataset_name]], dataset_name, group.by = "disease_phenotype", no_legend = FALSE)
    recursive_integrated[[dataset_name]] <- GetClusterPlot(recursive_integrated[[dataset_name]], dataset_name, group.by = "patient", no_legend = FALSE)
  } else {
    recursive_integrated[[dataset_name]] <- GetClusterPlot(recursive_integrated[[dataset_name]], dataset_name, group.by = 'celltype.l1', no_legend = FALSE)
    recursive_integrated[[dataset_name]] <- GetClusterPlot(recursive_integrated[[dataset_name]], dataset_name, group.by = 'celltype.l2', no_legend = FALSE)
  }
}


# Create a dictionary where the key is the name of each cluster in the heirarchy which has subclusters, and the
# value is a list of subcluster names for that cluster
# In case we have an entropy cutoff, we keep two seperate lists, as some labeling tasks will be expecting all clusters in the heirarchy,
# while the analysis tasks we only want clusters who meet the entropy cutoff
dataset_child_clusters = list()
for (type in c('analysis', 'labeling')) {
  for (dataset_name in names(recursive_integrated)) {
    cluster_nums = as.numeric(levels(recursive_integrated[[dataset_name]]@meta.data$seurat_clusters))
    clusters = list()
    for (i in seq(length(cluster_nums))) {
      cluster = cluster_nums[[i]]
      if (entropy_cutoff > 0 & type == 'analysis') {
        cluster_patient_fracs = list()
        cluster_idxs = seq(recursive_integrated[[dataset_name]]$RNA@data@Dim[[2]])[recursive_integrated[[dataset_name]]@meta.data$seurat_clusters == cluster]
        cluster_patients = recursive_integrated[[dataset_name]]@meta.data$patient[cluster_idxs]
        for (patient_name in names(paths)) {
          cluster_patient_fracs[[patient_name]] = length(cluster_patients[cluster_patients == patient_name]) / length(cluster_patients)
        }

        patient_entropy = GetEntropy(cluster_patient_fracs)

        if (patient_entropy > entropy_cutoff) {
          clusters[[i]] = paste0(dataset_name, 'c', cluster)
        }
      } else {
        clusters[[i]] = paste0(dataset_name, 'c', cluster)
      }
    }
    clusters = unlist(clusters)
    dataset_child_clusters[[type]][[dataset_name]] = clusters
  }
}

for (cluster in names(dataset_child_clusters$analysis)) {
  if (length(dataset_child_clusters$analysis[[cluster]]) == 1) {
    dataset_child_clusters$analysis[[cluster]] = NULL
  }
}


# Create a dictionary where the key is the name of each cluster in the heirarchy which has subclusters, and the
# value is a dictionary of the list of subclusters at each level in the heirarchy contained within the parent cluster
# In case we have an entropy cutoff, we keep two seperate lists, as some labeling tasks will be expecting all clusters in the heirarchy,
# while the analysis tasks we only want clusters who meet the entropy cutoff
# We also label the cells in each seurat object with the names of  subclusters
# at each level in the heirarchy contained within the parent cluster
dataset_level_clusters = list()
dataset_sub_cluster_idxs = list()
for (type in c('analysis', 'labeling')) {
  for (dataset_name in names(dataset_child_clusters[[type]])) {
    print("")
    print(dataset_name)
    dataset_level = str_count(dataset_name, "c") - 1
    dataset_level_clusters[[type]][[dataset_name]] = list()
    dataset_sub_cluster_idxs[[dataset_name]] = list()
    dataset_n_cells = recursive_integrated[[dataset_name]][[assay]]@
      data@
      Dim[[2]]
    for (level in seq(dataset_level, m_level)) {
      dataset_sub_cluster_level_ids = rep("NULL", dataset_n_cells)
      print(paste0('level ', level))
      if (level == dataset_level) {
        dataset_level_clusters[[type]][[dataset_name]][[paste0(level)]] = dataset_child_clusters[[type]][[dataset_name]]
      } else {
        dataset_level_clusters[[type]][[dataset_name]][[paste0(level)]] = as.list(dataset_level_clusters[[type]][[dataset_name]][[paste0(level - 1)]])
        for (i in seq(length(dataset_level_clusters[[type]][[dataset_name]][[paste0(level)]]))) {
          child_dataset_name = dataset_level_clusters[[type]][[dataset_name]][[paste0(level)]][[i]]
          if (!is.null(dataset_child_clusters[[type]][[child_dataset_name]])) {
            dataset_level_clusters[[type]][[dataset_name]][[paste0(level)]][[i]] = dataset_child_clusters[[type]][[child_dataset_name]]
          }
        }
        dataset_level_clusters[[type]][[dataset_name]][[paste0(level)]] = unlist(dataset_level_clusters[[type]][[dataset_name]][[paste0(level)]])
      }
      print(dataset_level_clusters[[type]][[dataset_name]][[paste0(level)]])
      for (sub_cluster_name in dataset_level_clusters[[type]][[dataset_name]][[paste0(level)]]) {
        sub_cluster_name_split = unlist(strsplit(sub_cluster_name, split = 'c'))
        sub_cluster_parent_name = paste(sub_cluster_name_split[1:(length(sub_cluster_name_split) - 1)], collapse = 'c')
        sub_cluster_num = as.numeric(sub_cluster_name_split[[length(sub_cluster_name_split)]])
        print(sub_cluster_name)
        print(sub_cluster_parent_name)
        print(sub_cluster_num)
        if (!(sub_cluster_name %in% names(dataset_sub_cluster_idxs[[dataset_name]]))) {
          dataset_sub_cluster_idxs[[dataset_name]][[sub_cluster_name]] = recursive_integrated[[sub_cluster_parent_name]]@meta.data[[dataset_name]][recursive_integrated[[sub_cluster_parent_name]]@meta.data$seurat_clusters == sub_cluster_num]
        }
        idxs = dataset_sub_cluster_idxs[[dataset_name]][[sub_cluster_name]]
        dataset_sub_cluster_level_ids[idxs] = sub_cluster_name
      }
      if (type == 'labeling') {
        recursive_integrated[[dataset_name]] = AddMetaData(object = recursive_integrated[[dataset_name]], metadata = dataset_sub_cluster_level_ids, col.name = paste0('level_', level))
      }
    }
  }
}

# Determine how many clusters are present as leaf nodes if we consider each level as a stopping criterion
num_clusters_per_level = list()
for (level in seq(0, m_level)) {
  num_clusters_per_level[[paste0(level)]] = length(dataset_level_clusters[['analysis']][['cR']][[paste0(level)]])
}
num_clusters_title = paste0('num_clusters_', prefix, '.rds')
saveRDS(num_clusters_per_level, num_clusters_title)
num_clusters_per_level = readRDS(num_clusters_title)
print(num_clusters_per_level)


# This code is for generating plots for each cluster in the heirarchy which has subclusters where we can see what
# the subclusters would have looked like if they were visualized in UMAP space using the PCA embeddings of upper-level
# (parent, grandparent, etc.) clusters in the heirarchy
if (do_lower_clusters_higher_UMAP) {
  for (dataset_name in names(recursive_integrated)) {
    print(dataset_name)
    for (level in seq(0, m_level)) {
      if (paste0('level_', level) %in% names(recursive_integrated[[dataset_name]]@meta.data)) {
        recursive_integrated[[dataset_name]] <- GetClusterPlot(recursive_integrated[[dataset_name]], group.by = paste0('level_', level), title = dataset_name, no_legend = FALSE, shuffle_colors = TRUE)
      }
    }
  }

  for (dataset_name in names(recursive_integrated)) {
    higher_level_plots = list()
    sub_cluster_name_split = unlist(strsplit(dataset_name, split = 'c'))
    level = length(sub_cluster_name_split) - 2
    if (level > 0) {
      for (back_level in seq(level)) {
        upper_level = level - 1
        upper_level_name = paste0('level_', upper_level)
        sub_cluster_parent_name = paste(sub_cluster_name_split[1:(length(sub_cluster_name_split) - back_level)], collapse = 'c')
        # print(paste(dataset_name, sub_cluster_parent_name, sep = ' '))
        parent_integrated = recursive_integrated[[sub_cluster_parent_name]]
        Idents(parent_integrated) = upper_level_name
        parent_integrated_cluster_subset = subset(parent_integrated, idents = dataset_name)
        if (suppressWarnings(!all(recursive_integrated[[dataset_name]]@meta.data$cR == parent_integrated_cluster_subset@meta.data$cR))) {
          # print('some cells are different')
          print(paste(dataset_name, sub_cluster_parent_name, sep = ' '))
          if (do_lower_clusters_higher_UMAP) {
            Idents(parent_integrated_cluster_subset) = 'cR'
            parent_integrated_cluster_subset = subset(parent_integrated_cluster_subset, idents = recursive_integrated[[dataset_name]]@meta.data$cR)
          }
        }
        if (do_lower_clusters_higher_UMAP) {
          parent_integrated_cluster_subset@meta.data$seurat_clusters = recursive_integrated[[dataset_name]]@meta.data$seurat_clusters
          higher_title = paste0(sub_cluster_parent_name, ' Subset Equal to ', dataset_name)
          parent_integrated_cluster_subset = GetClusterPlot(parent_integrated_cluster_subset, title = higher_title, no_legend = FALSE)
          higher_level_plots[[sub_cluster_parent_name]] = parent_integrated_cluster_subset@misc[[paste0('seurat_clusters', '_UMAP_', higher_title)]]
        }
      }
      if (do_lower_clusters_higher_UMAP) {
        recursive_integrated[[dataset_name]]@misc$higher_level_plots = higher_level_plots
      }
    }
  }
}

assay_equivalent = 'integrated'

# For benchmark datasets, here we generate a high-resolution seurat clustering where the number of clusters generated
# via a single round of clustering matches the number of (leaf node) clusters generated in a corresponding recursive clustering
if (!is_helmsley) {
  m_level_equivalent = m_level

  if (is_helmsley) {
    equivalent_clustering_filename = paste0("seurat_equivalent_clustering_helmsley_", assay, "_")
  } else {
    equivalent_clustering_filename = paste0("seurat_equivalent_clustering_reference", '_', ref_name)
  }
  equivalent_clustering_filename = paste0(equivalent_clustering_filename, c_resolution, '_', DEG_cutoff, '.rds')

  # We first start with the resolution parameter used recursive clustering, then we keep doubling the resolution parameter
  # Until we produce more clusters than in the recursive approach
  if (!file.exists(equivalent_clustering_filename)) {
    num_clusters = 0
    num_clusters_per_resolution = list()
    num_clusters_per_resolution[[paste0(c_resolution)]] = num_clusters_per_level[[paste0(0)]]
    seurat_equivalent_clustering = list()
    seurat_equivalent_clustering_resolution = list()
    current_resolution = c_resolution
    while (num_clusters < num_clusters_per_level[[paste0(m_level_equivalent)]]) {
      temp_integrated = FindSeuratClusters(dataset = recursive_integrated$cR, assay = assay_equivalent, resolution = current_resolution, rerun = FALSE, savetime = TRUE, cluster_name = 'cR', method_name = 'seurat_equivalent')
      num_clusters = length(as.numeric(levels(temp_integrated@meta.data$seurat_clusters)))
      print(paste0(current_resolution, ' ', num_clusters))
      num_clusters_per_resolution[[paste0(current_resolution)]] = num_clusters
      for (level in seq(0, m_level_equivalent)) {
        if (num_clusters == num_clusters_per_level[[paste0(level)]]) {
          seurat_equivalent_clustering[[paste0(level)]] = temp_integrated@meta.data$seurat_clusters
          seurat_equivalent_clustering_resolution[[paste0(level)]] = current_resolution
        }
      }
      current_resolution = 2 * current_resolution
    }

    print(num_clusters_per_level)

    # Next, we use a pseudo binary search type approach to narrow in on a resolution parameter value which can generate
    # the same number of clusters as the recursive clustering approach.  For the louvain clustering which is used, the
    # number of clusters generated is not quite a monotonic function of resolution parameter value, and some stochasticity
    # is involved in the number of clusters generated by the louvain approach, so we use a different random seed at each
    # iteration so that an equivalent number of clusters can be generated in spite of the lack of monotonicity of the
    # number of clusters as a function of resolution
    for (level in seq(0, m_level_equivalent)) {
      random_seed = -1
      while (!(paste0(level) %in% names(seurat_equivalent_clustering))) {
        random_seed = random_seed + 1
        lower_bound = min(as.numeric(names(num_clusters_per_resolution)))
        upper_bound = max(as.numeric(names(num_clusters_per_resolution)))
        for (res in names(num_clusters_per_resolution)) {
          res = as.numeric(res)
          if (res > lower_bound && num_clusters_per_resolution[[paste0(res)]] < num_clusters_per_level[[paste0(level)]]) {
            lower_bound = res
          }
          if (res < upper_bound && num_clusters_per_resolution[[paste0(res)]] > num_clusters_per_level[[paste0(level)]]) {
            upper_bound = res
          }
        }
        print(paste0('level ', level, ' num clusters ', num_clusters_per_level[[paste0(level)]], ' lower bound ', lower_bound, ' upper bound ', upper_bound, ' lower bound clusters ', num_clusters_per_resolution[[paste0(lower_bound)]], ' upper bound clusters ', num_clusters_per_resolution[[paste0(upper_bound)]]))
        current_resolution = upper_bound - (upper_bound - lower_bound) / 2
        temp_integrated <- FindSeuratClusters(dataset = recursive_integrated$cR, assay = assay_equivalent, resolution = current_resolution, rerun = FALSE, random_seed = random_seed, savetime = TRUE, cluster_name = 'cR', method_name = 'seurat_equivalent')

        num_clusters = length(as.numeric(levels(temp_integrated@meta.data$seurat_clusters)))
        print(paste0(current_resolution, ' ', num_clusters))
        num_clusters_per_resolution[[paste0(current_resolution)]] = num_clusters
        for (level2 in seq(m_level_equivalent)) {
          if (num_clusters == num_clusters_per_level[[paste0(level2)]]) {
            seurat_equivalent_clustering[[paste0(level2)]] = temp_integrated@meta.data$seurat_clusters
            seurat_equivalent_clustering_resolution[[paste0(level)]] = current_resolution
          }
        }
      }
    }
    temp = seurat_equivalent_clustering
    seurat_equivalent_clustering = list()
    for (level in seq(0, m_level_equivalent)) {
      seurat_equivalent_clustering[[paste0(level)]] = temp[[paste0(level)]]
    }
    temp = NULL
    print(length(as.numeric(levels(temp_integrated@meta.data$seurat_clusters))))
    saveRDS(seurat_equivalent_clustering, file = equivalent_clustering_filename)
  } else {
    seurat_equivalent_clustering = readRDS(equivalent_clustering_filename)
  }


  # Create cluster plots for the generated high-resolution equivalent clustering
  for (dataset_name in names(recursive_integrated)) {
    print(dataset_name)
    for (level in seq(0, m_level)) {
      if (paste0('level_', level, '_seurat_equivalent') %in% names(recursive_integrated[[dataset_name]]@meta.data)) {
        recursive_integrated[[dataset_name]] <- GetClusterPlot(recursive_integrated[[dataset_name]], group.by = paste0('level_', level, '_seurat_equivalent'), title = dataset_name, no_legend = FALSE, shuffle_colors = TRUE)
      }
    }
  }
  # Label the base dataset with the generated high-resolution equivalent clustering labels for each cell
  for (level in seq(0, m_level_equivalent)) {
    recursive_integrated[['cR']] = AddMetaData(object = recursive_integrated[['cR']], metadata = seurat_equivalent_clustering[[paste0(level)]], col.name = paste0('level_', level, '_seurat_equivalent'))
  }
}

# A recursive algorithm for finding a depth-first search-like ordered list of cluster names in the cluster heirarchy
GetRecursiveClusterNames <- function(recursive_data, cluster_sequence = 'cR', cluster_list = NULL) {
  cluster_list = c(cluster_list, cluster_sequence)
  if (is.null(dataset_level_clusters$analysis[[cluster_sequence]])) {
    return(cluster_list)
  }
  cluster_names = dataset_level_clusters$analysis[[cluster_sequence]][[names(dataset_level_clusters$analysis[[cluster_sequence]])[[1]]]]
  for (cluster_name in cluster_names) {
    cluster_list = GetRecursiveClusterNames(recursive_data, cluster_sequence = cluster_name, cluster_list = cluster_list)
  }
  return(cluster_list)
}

# Generate a depth-first search-like ordered list of cluster names in the cluster heirarchy
recursive_cluster_list = GetRecursiveClusterNames(recursive_data = recursive_integrated)

# Generate a breadth-first search-like ordered list of cluster names in the cluster heirarchy
level_ordered_cluster_names = c('cR', as.vector(unlist(dataset_level_clusters[['analysis']][['cR']])))
level_ordered_cluster_names = level_ordered_cluster_names[!duplicated(level_ordered_cluster_names)]


# Get a list of seurat objects for all clusters in the heirarchy for which we performed an additional iteration
# of the seurat clustering pipeline but found no additional clusters
recursive_unintegrated = list()
for (cluster_name in level_ordered_cluster_names) {
  sub_cluster_name_split = unlist(strsplit(cluster_name, split = 'c'))
  upper_level_cluster_name = paste(sub_cluster_name_split[1:(length(sub_cluster_name_split) - 1)], collapse = 'c')
  cluster_num = tail(sub_cluster_name_split, 1)
  if (!(cluster_name %in% names(recursive_integrated))) {
    if (!(cluster_name %in% names(recursive_integrated_extra))) {
      print(cluster_name)
      parent_integrated = recursive_integrated[[upper_level_cluster_name]]
      Idents(parent_integrated) <- "seurat_clusters"
      child_integrated = subset(parent_integrated, idents = cluster_num)
      temp = parent_integrated@assays[[assay_equivalent]]@var.features
      DefaultAssay(child_integrated) = assay
      child_integrated = FindVariableFeatures(child_integrated, verbose = FALSE)
      print(1 - length(intersect(temp, child_integrated@assays[[assay]]@var.features)) / length(temp))
      recursive_unintegrated[[cluster_name]] = child_integrated
    }
  }
}


# Code for generating the heatmaps for cluster statistical significance in the Abundance and DEG based patient-wise
# comparison of phenotypes by cluster analysis
GetTestHeatmap <- function(mat, x_var_name, y_var_name, fill_name, text = TRUE, significance = NULL, wrap_var_name = NULL, title_list = NULL, rounding = NULL, angle = NULL, font_size = 3.88) {
  if (!is.null(rounding)) {
    g_text = geom_text(aes(label = round(value, rounding)), size = font_size)
  } else {
    g_text = geom_text(aes(label = value), size = font_size)
  }
  plot <- ggplot(mat, aes_string(x = x_var_name, y = y_var_name, fill = fill_name)) +
    geom_tile() +
    ylim(rev(levels(as.factor(mat[[y_var_name]]))))
  base_color = '#669bcc'
  upper_color = toupper('#FFFFFF')
  middle_color = toupper(lighten(base_color, amount = 0.9))
  lower_color = toupper(darken(base_color, amount = 0.4))
  if (fill_name == 'P.value') {
    my_breaks = c(significance / 100, significance / 10, significance)
    plot = plot + scale_fill_gradient2(midpoint = 0.005, low = lower_color, mid = middle_color, high = upper_color, space = "Lab", trans = "log10", limits = c(min(mat$P.value), significance), labels = my_breaks, breaks = my_breaks, na.value = "white")
  } else {
    plot = plot + scale_fill_gradient2(low = "white", high = "darkgreen", space = "Lab", na.value = "white")
  }
  if (!is.null(angle)) {
    plot = plot + theme(axis.text.x = element_text(angle = angle, vjust = 0.35, hjust = 0))
  }

  if (!is.null(wrap_var_name)) {
    if (wrap_var_name == 'comparison') {
      plot = plot + facet_wrap(~comparison, scales = 'free_y', ncol = length(unique(mat['comparison'])))
    }
    if (wrap_var_name == 'type') {
      plot = plot + facet_wrap(~type, scales = 'free_y', ncol = length(unique(mat['type'])))
    }
  }
  if (text) {
    plot = plot + g_text
  }
  if (!is.null(title)) {
    plot = plot + ggtitle(str_to_title(str_replace(paste(title_list, collapse = ', '), pattern = '_', replacement = ' ')))
  }
  return(plot)
}


# Code for generating the cluster purity plots for benchmark analysis
GetCorrespondencePlot <- function(dataset, cluster_level, annotation_level, title) {
  subset_cluster_df = data.frame(cluster = as.vector(dataset@meta.data[[cluster_level]]), l1 = dataset@meta.data[[annotation_level]])

  subset_grouped <- suppressMessages(suppressWarnings(list(l1 = dcast(subset_cluster_df, cluster ~ l1, fun.aggregate = length))))

  subset_correspondence_row_normalized <- list(l1 = NULL)

  for (level in names(subset_correspondence_row_normalized)) {
    subset_grouped_matrix = data.matrix(subset_grouped[[level]][tail(names(subset_grouped[[level]]), length(names(subset_grouped[[level]])) - 1)])
    rownames(subset_grouped_matrix) <- as.vector(subset_grouped[[level]][['cluster']])

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
    subset_correspondence_row_normalized[[level]] <- subset_grouped_df_final
    colnames(subset_correspondence_row_normalized[[level]]) = c("Cluster", "Type", "Proportion")
    cell_names = aggregate(Proportion ~ Type, data = subset_correspondence_row_normalized[[level]], FUN = max)
    cell_names$Proportion = 0

    out_temp = merge(aggregate(Proportion ~ Cluster, data = subset_correspondence_row_normalized[[level]], FUN = max), subset_correspondence_row_normalized[[level]], by = c("Proportion", "Cluster"))
    out = list()
    out_temp2 = aggregate(Proportion ~ Type, data = rbind(cell_names, aggregate(Proportion ~ Type, data = out_temp, FUN = mean)), FUN = max)

    out[['Cell Type Mean Cluster Purity']] = out_temp2$Proportion
    out[['Cluster Purity']] = out_temp$Proportion

    color_values = seq(100) / 100
    palette = colorRampPalette(c('white', 'orange', 'red'))(length(color_values))

    plot <- ggplot(subset_correspondence_row_normalized[[level]], aes(x = Type, y = Cluster, fill = Proportion)) +
      geom_tile() +
      geom_text(aes(label = round(Proportion, 2), color = ifelse(Proportion < 0.005, "light grey", 'black')), size = 3.5) +
      scale_fill_gradient2(midpoint = 0.5, low = "white", mid = palette[20], high = palette[70], space = "Lab") +
      scale_color_identity() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(title)
    #suppressMessages(ggsave(filename = paste0(title, '.png'), width = 11, height = 9.92, plot = plot, device = 'png', dpi = 400))
    print(plot)
    suppressMessages(ggsave(filename = paste0(title, '.png'), width = 11.7, height = 10.75, plot = plot, device = 'png', dpi = 400))
  }
  return(out)
}

# Code for generating the cluster purity analysis for benchmark analysis
if (!is_helmsley) {
  values = vector(mode = 'list')
  annotation_level_name_vector = vector(mode = 'character')
  level_name_vector = vector(mode = 'character')
  type_vector = vector(mode = 'character')
  values_list = vector(mode = 'list')
  values_df_list = vector(mode = 'list')
  values_vector_list = vector(mode = 'list')
  for (value_type in c('mean', 't_value', 'w_value')) {
    values_vector_list[[value_type]] = vector(mode = 'numeric')
  }
  for (annotation_level in seq(annotation_levels)) {
    for (level in seq(m_level)) {
      if (paste0('level_', level, '_seurat_equivalent') %in% names(recursive_integrated$cR@meta.data)) {
        for (clustering in c('Iterative', 'Seurat Equivalent')) {
          if (clustering == 'Iterative') {
            cluster_level_name = paste0('level_', level)
          } else {
            cluster_level_name = paste0('level_', level, '_seurat_equivalent')
          }
          annotation_title = 'celltype.l'
          values[[paste0('annotation level ', annotation_level)]][[paste0('level ', level)]][[clustering]] = GetCorrespondencePlot(recursive_integrated$cR, cluster_level = cluster_level_name, annotation_level = paste0(annotation_title, annotation_level), title =
            paste('Level', level, clustering, 'Clustering, Annotation Level', annotation_level, ref_name, c_resolution, DEG_cutoff))
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
          annotation_level_name_vector = c(annotation_level_name_vector, paste0('annotation level ', annotation_level))
          level_name_vector = c(level_name_vector, paste0('level ', level))
          type_vector = c(type_vector, result_type)
          for (value_type in c('mean', 't_value', 'w_value')) {
            values_vector_list[[value_type]] = c(values_vector_list[[value_type]], values_list[[value_type]])
          }
        }
      }
    }
  }

  for (value_type in c('mean', 't_value', 'w_value')) {
    values_df_list[[value_type]] = data.frame(Annotation = annotation_level_name_vector, Level = level_name_vector, type = type_vector, value = values_vector_list[[value_type]])
    results_heatmap = GetTestHeatmap(mat = values_df_list[[value_type]], x_var_name = 'Annotation', y_var_name = 'Level', fill_name = 'value', text = TRUE, significance = NULL, wrap_var_name = 'type', title_list = paste0(str_to_title(str_replace(value_type, pattern = '_', replacement = ' ')), ' Difference Between Iterative vs High-Resolution Clustering ', ref_name, c_resolution, ' ', DEG_cutoff), rounding = 3, angle = NULL, font_size = 3.88) +
      scale_fill_gradient2(low = 'darkblue', mid = 'white', high = "darkgreen") +
      theme(panel.background = element_blank(), plot.background = element_blank())
    if (value_type == 'mean') {
      print(results_heatmap)
    }
    ggsave(paste0('results_', value_type, '_', ref_name, c_resolution, '_', DEG_cutoff, '.png'), plot = results_heatmap, device = 'png')
  }
  saveRDS(values, paste0('results_', ref_name, c_resolution, '_', DEG_cutoff, '.rds'))

}
num_clusters_title = paste0('num_clusters_', prefix, '.rds')
num_clusters_per_level = readRDS(num_clusters_title)
print(num_clusters_per_level)

print('main')
