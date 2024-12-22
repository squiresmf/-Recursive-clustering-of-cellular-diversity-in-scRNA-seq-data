# Code for creating heirarchically structured heat map plots for PCA Similarity analysis
# This code produces figure 5.

shared_variable_individual_plots = list()
shared_variable_final_plots = list()

for (cluster_name in level_ordered_cluster_names[2:length(level_ordered_cluster_names)]) {
  level = length(unlist(strsplit(cluster_name, split = 'c'))) - 3
  if (level == 2) {
    print(cluster_name)
    if (cluster_name %in% names(recursive_integrated)) {
      recursive_integrated[[cluster_name]] = RunPCA(recursive_integrated[[cluster_name]], assay = 'integrated', npcs = 50, verbose = FALSE)
    } else if (cluster_name %in% names(recursive_integrated_extra)) {
      recursive_integrated_extra[[cluster_name]] = RunPCA(recursive_integrated_extra[[cluster_name]], assay = 'integrated', npcs = 50, verbose = FALSE)
    } else {
      new_pca_results = tryCatch(
        expr = {
          RunPCA(recursive_unintegrated[[cluster_name]], assay = 'integrated', npcs = 50, verbose = FALSE)
        },
        error = function(e) {
          print(paste0("PCA not possible due to cluster size.  Retaining parent PCA results."))
          return('error')
        }
      )
      if (!is.string(new_pca_results)) {
        recursive_unintegrated[[cluster_name]] = new_pca_results
      }
    }
  }
}

root_eigenvectors = recursive_integrated[['cR']]@reductions$pca@feature.loadings
root_eigenvectors_df = as.data.frame(root_eigenvectors)
root_eigenvectors_df$id <- 1:nrow(root_eigenvectors_df)

for (max_level in 2) {
  i = 0
  m_rows = list()
  for (cluster_name in level_ordered_cluster_names[2:length(level_ordered_cluster_names)]) {
    print(cluster_name)
    sub_cluster_name_split = unlist(strsplit(cluster_name, split = 'c'))
    upper_level_cluster_name = paste(sub_cluster_name_split[1:(length(sub_cluster_name_split) - 1)], collapse = 'c')
    cluster_num = tail(sub_cluster_name_split, 1)
    if (cluster_name %in% names(recursive_integrated)) {
      cluster_eigenvectors = recursive_integrated[[cluster_name]]@reductions$pca@feature.loadings
    } else if (cluster_name %in% names(recursive_integrated_extra)) {
      cluster_eigenvectors = recursive_integrated_extra[[cluster_name]]@reductions$pca@feature.loadings
    } else {
      cluster_eigenvectors = recursive_unintegrated[[cluster_name]]@reductions$pca@feature.loadings
    }

    print(nrow(cluster_eigenvectors))

    # Some feature loading matrices are missing genes, we add those genes back with 0 values
    cluster_eigenvectors_df = as.data.frame(cluster_eigenvectors)
    cluster_eigenvectors_df = merge(x = root_eigenvectors_df, y = cluster_eigenvectors_df, by = 0, all.x = TRUE)
    cluster_eigenvectors_df = cluster_eigenvectors_df[order(cluster_eigenvectors_df$id),]
    cluster_eigenvectors_df[is.na(cluster_eigenvectors_df)] = 0
    cluster_eigenvectors_df = cluster_eigenvectors_df[, grepl(".y", names(cluster_eigenvectors_df))]
    colnames(cluster_eigenvectors_df) = gsub('.y', '', colnames(cluster_eigenvectors_df))
    rownames(cluster_eigenvectors_df) = rownames(root_eigenvectors_df)
    cluster_eigenvectors = as.matrix(cluster_eigenvectors_df)

    pca_similarity_factor = sum(diag(t(root_eigenvectors) %*%
                                       cluster_eigenvectors %*%
                                       t(cluster_eigenvectors) %*%
                                       root_eigenvectors)) / ncol(cluster_eigenvectors_df)

    i = i + 1
    m_rows[[i]] = list(cluster_name, pca_similarity_factor)
  }
  m = as.data.frame(do.call(rbind, m_rows))
  colnames(m) = c("Cluster", "PCA Similarity")
  m['comparison'] = 'cluster vs root'
  if (entropy_cutoff == 0.5) {
    GetRadialTreeFromDataframe(df = m, tree_comparisons = 'cluster vs root', val_name = 'PCA Similarity', tree_level = 2, gravities = "NorthEast")
  } else {
    GetRadialTreeFromDataframe(df = m, tree_comparisons = 'cluster vs root', val_name = 'PCA Similarity', tree_level = 2, gravities = "NorthEast", dimension = 150)
  }
}

temp = as.matrix(recursive_integrated$cR@assays$integrated@data)
identities = sample(seq(ncol(temp)))
halfway = floor(length(identities) / 2)
half1 = temp[, identities[1:halfway]]
half2 = temp[, identities[(halfway + 1):length(identities)]]
npcs = 50
pca.results1 <- irlba::irlba(A = t(half1), nv = npcs)
pca.results2 <- irlba::irlba(A = t(half2), nv = npcs)

pca_similarity_factor = sum(diag(t(pca.results1$v) %*%
                                   pca.results2$v %*%
                                   t(pca.results2$v) %*%
                                   pca.results1$v)) / npcs


print(nrow(temp))
temp2 = sample(c(as.matrix(recursive_integrated$cR@assays$integrated@data)))
temp2 <- matrix(temp2, nrow = nrow(temp), byrow = TRUE)
identities = sample(seq(ncol(temp2)))
halfway = floor(length(identities) / 2)
half1 = temp2[, identities[1:halfway]]
half2 = temp2[, identities[(halfway + 1):length(identities)]]
npcs = 50
pca.results1 <- irlba::irlba(A = t(half1), nv = npcs)
pca.results2 <- irlba::irlba(A = t(half2), nv = npcs)

pca_similarity_factor2 = sum(diag(t(pca.results1$v) %*%
                                    pca.results2$v %*%
                                    t(pca.results2$v) %*%
                                    pca.results1$v)) / npcs


print('SharedVariables')




