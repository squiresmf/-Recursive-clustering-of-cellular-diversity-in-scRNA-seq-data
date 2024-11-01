# Code for generating hierarchical clustering tree plots.
# This code will reproduce figure 1 if main has been run with c_resolution = 0.02 and reference_name = "human_PBMC/"
# This code also includes functions used to generate Crohn's scRNA-seq data analysis figures 5 and 6
# graphics.off()

radial_tree_list = list()
radial_tree_list_no_color = list()

# par(mar = c(1, 1, 1, 1))
# par("mar")

# A function to replace the node name in cluster with just the parent cluster subcluster id rather than the full id
# of the cluster
ReplaceNodeName <- function(tree) {
  tree$name = tail(unlist(strsplit(tree$name, 'c')), 1)
  if (!is.null(tree$children)) {
    for (child in seq_along(tree$children)) {
      tree$children[[child]] = ReplaceNodeName(tree$children[[child]])
    }
  }
  return(tree)
}

# Create the tree data structure which will be used in future analysis to fill in with color values
for (max_level in seq(0, m_level)) {
  tree_list = list()
  tree_list[['cR']] = Node$new("cR")
  for (cluster_name in recursive_cluster_list[2:length(recursive_cluster_list)]) {
    cluster_name_split = unlist(strsplit(cluster_name, split = 'c'))
    if (length(cluster_name_split) - 3 <= max_level) {
      cluster_upper_level_name = paste(cluster_name_split[1:(length(cluster_name_split) - 1)], collapse = 'c')
      tree_list[[cluster_name]] = tree_list[[cluster_upper_level_name]]$AddChild(cluster_name) }
  }

  nested_tree_list <- ToListExplicit(tree_list[['cR']], unname = TRUE)
  radial_tree_list_no_color[[paste0(max_level)]] = ReplaceNodeName(nested_tree_list)
  radial_tree = radialNetwork(radial_tree_list_no_color[[paste0(max_level)]], fontSize = 12 - max_level, fontFamily = "verdana", linkColour = "#bbb", nodeColour = '#FFFFFF', textColour =, nodeStroke = '#808080')
  tree_name = paste0('radial_tree', max_level)
  saveNetwork(radial_tree, paste0(tree_name, '.html'), selfcontained = FALSE)
  webshot(paste0(tree_name, '.html'), paste0(tree_name, '.png'), vwidth = (max_level + 2) * 100, vheight = (max_level + 2) * 100, zoom = 3)
  pp <- readPNG(paste0(tree_name, '.png'))
  plot.new()
  radial_tree_list[[paste0(max_level)]] = rasterImage(pp, 0, 0, 1, 1)
  # print(radial_tree_list[[paste0(max_level)]])
  pp1 <- image_read(paste0(tree_name, '.png'))
  image_write(pp1, paste0(tree_name, '.eps'), format = 'eps')
  image_write(pp1, paste0(tree_name, '.pdf'), format = 'pdf')
}


# This code generates heirarchical tree plot visualizations for abunance, DEG, and PCA similarity analysis.
# It also creates a seperate plot with the same heatmap color properties and retreives the heatmap color scale from that
# plot and writes it overtop of a image capture of the heirarcical tree

GetRadialTreeFromDataframe <- function(df, tree_comparisons, val_name, tree_level, gravities, title_extra = '', dimension = 100, annotation_colors = NULL) {
  radial_tree_m = df
  if ((val_name == 'P.value') & !(grepl('empirical', title_extra, fixed = TRUE))) {
    cluster_names_true = character(length = nrow(df))
    clusters = unlist(as.character(df$Cluster))
    cluster_levels = unlist(as.character(df$Level))
    for (i in seq(nrow(df))) {
      cluster_level = cluster_levels[[i]]
      cluster = clusters[[i]]
      cluster_split = strsplit(cluster, 'c')
      deepest_level = length(cluster_split[[1]])
      cluster_true = paste(cluster_split[[1]][1:min(as.integer(cluster_level) + 3, deepest_level)], collapse = 'c')
      cluster_names_true[[i]] = cluster_true
    }
    radial_tree_m['Cluster'] = cluster_names_true
  }
  radial_tree_m = radial_tree_m[c('Cluster', 'comparison', val_name)]
  radial_tree_m = unique(radial_tree_m)
  x_labels = intersect(level_ordered_cluster_names, unique(unlist(radial_tree_m['Cluster'])))
  radial_tree_m = radial_tree_m[order(factor(unlist(radial_tree_m['Cluster']), levels = x_labels)),]

  base_colors = c("#669bcc", "#ccb166", "#CC6677")
  for (tree_comparison_idx in seq_along(tree_comparisons)) {
    tree_comparison = tree_comparisons[[tree_comparison_idx]]
    base_color = '#669bcc'
    if (val_name == 'PCA Similarity') {
      upper_color = toupper(lighten(base_color, amount = 0.9))
    } else {
      upper_color = toupper('#FFFFFF')
    }
    middle_color = toupper(lighten(base_color, amount = 0.6))
    lower_color = toupper(darken(base_color, amount = 0.8))
    radial_tree_m_comparison = radial_tree_m[radial_tree_m[['comparison']] == tree_comparison,]
    color_values = radial_tree_m_comparison[val_name]
    color_values = as.numeric(unlist(color_values))
    if (val_name == 'P.value') {
      # All non-significant p-values are shaded white
      color_values[color_values > 0.05] = 1
      # The scale for p-value is limited to a minimum of 0.0015, so all p-values lower will have the same color as 0.0015
      color_values[color_values < 0.0015] = 0.0015
      color_values = -log10(color_values)
    }
    if (val_name == 'DEGs') {
      max_percentile = 0.99
      # The scale for DEG is limited to a maximum of the max_percentile percentile of the DEG values, so all DEG values higher
      # than max_percentile percentile will have the same color as the max_percentile percentile value
      max_color = quantile(color_values, max_percentile)[[1]]
      color_values[color_values > max_color] = max_color
    }
    if (val_name == 'PCA Similarity') {
      palette = colorRampPalette(c(lower_color, middle_color, upper_color))(length(color_values))
    } else {
      palette = colorRampPalette(c(upper_color, middle_color, lower_color))(length(color_values))
    }
    colors_node = palette[cut(color_values, length(color_values))]
    colors_stroke = colors_node
    if (val_name == 'DEGs') {
      colors_stroke[color_values < 0.20 * max(color_values)] = 'D3D3D3'
    }

    colors_node = c('#FFFFFF', colors_node)
    colors_stroke = c('#FFFFFF', colors_stroke)

    if (val_name == 'PCA Similarity') {
      colors_stroke[[1]] = '#FFFFFF'
    } else {
      colors_stroke[[1]] = 'D3D3D3'
      colors_stroke[colors_node == upper_color] = 'D3D3D3'
    }

    if (!is.null(annotation_colors)) {
      for (annotation_cluster_name in names(annotation_colors[[tree_comparison]])) {
        annotation_index = which(radial_tree_m_comparison$Cluster == annotation_cluster_name)
        colors_stroke[[annotation_index + 1]] = annotation_colors[[tree_comparison]][[annotation_cluster_name]]
      }
    }

    colors_node_JS = JS(paste0('function(d, i) { return ', paste0('["', paste(colors_node, collapse = '", "'), '"]'), '[i]; }'))  # Code from https://rpubs.com/securl/NetworkDiagramR
    colors_stroke_JS = JS(paste0('function(d, i) { return ', paste0('["', paste(colors_stroke, collapse = '", "'), '"]'), '[i]; }'))  # Code from https://rpubs.com/securl/NetworkDiagramR
    radial_tree = radialNetwork(radial_tree_list_no_color[[paste0(tree_level)]], fontSize = 12 - tree_level, fontFamily = "verdana", linkColour = "#bbb", nodeColour = colors_node_JS, textColour = '232b2b', nodeStroke = colors_stroke_JS)
    tree_name = paste0(stri_trans_general(str = paste('radial_tree', tree_level, tree_comparison, val_name, sep = ' '), id = "Latin-ASCII"), title_extra)
    saveNetwork(radial_tree, paste0(tree_name, '.html'), selfcontained = FALSE)
    webshot(paste0(tree_name, '.html'), paste0(tree_name, '.png'), vwidth = (tree_level + 2) * dimension, vheight = (level + 2) * dimension, zoom = 3)
    image_write(image_trim(image_read(paste0(tree_name, '.png'), strip = TRUE)), paste0(tree_name, '.png'), format = 'png')

    if (val_name == 'P.value') {
      color_values = 10^(-color_values)
    }

    plot_df = data.frame(values = color_values)
    ordered_legend_colors = colors_node[2:length(colors_node)][order(color_values)]

    legend_start_color = ordered_legend_colors[[1]]
    legend_middle_color = colors_node[2:length(colors_node)][[which.min(abs(color_values - (max(color_values) - (max(color_values) - min(color_values)) / 2)))]]
    if (val_name == 'P.value') {
      legend_end_color = tail(ordered_legend_colors[sort(color_values) < 0.05], 1)
    } else {
      legend_end_color = tail(ordered_legend_colors, 1)
    }
    plot = ggplot(plot_df, aes(values, values, color = values)) +
      geom_point() +
      theme(legend.text = element_text(size = 18)) +
      theme(legend.title = element_text(size = 18))

    if (val_name == 'P.value') {
      plot = plot + scale_color_gradientn(colours = c(legend_start_color, legend_middle_color, legend_end_color), limits = c(min(color_values), 0.05), name = 'P Value', trans = "log10",
                                          breaks = c(0.05, 0.015, 0.005, 0.0015), labels = c(0.05, 0.015, 0.005, 0.0015))
    } else if (val_name == 'PCA Similarity') {
      plot = plot + scale_color_gradientn(colours = c(legend_start_color, legend_middle_color, legend_end_color), limits = c(min(color_values), ceiling(max(color_values) * 10)) / 10, name = val_name)
    } else if (val_name == 'DEGs') {
      plot = plot + scale_color_gradientn(colours = c(legend_start_color, legend_middle_color, legend_end_color), limits = c(min(color_values), ceiling(max(color_values) * 10)) / 10, name = 'DEG Fold\nChange')
    } else {
      plot = plot + scale_color_gradientn(colours = c(legend_start_color, legend_middle_color, legend_end_color), limits = c(min(color_values), max(color_values)), name = val_name)
    }
    legend = cowplot::get_legend(plot)

    legend_name = paste0(tree_name, ' ', 'legend')
    ggsave(paste0(legend_name, '.png'), legend, device = "png", dpi = 400)
    image_write(image_trim(image_read(paste0(legend_name, '.png'), strip = TRUE)), paste0(legend_name, '.png'), format = 'png')

    pp1 <- image_read(paste0(tree_name, '.png'))
    pp2 <- image_read(paste0(legend_name, '.png'))
    if (val_name == 'PCA Similarity') {
      pp2 <- image_scale(pp2, "25%x")
    } else {
      pp2 <- image_scale(pp2, "28%x")
    }
    img_with_inset <- pp1 %>% image_composite(
      pp2,
      operator = "Atop",
      gravity = gravities[[tree_comparison_idx]],
    )
    img_with_inset = image_transparent(image = img_with_inset, color = '#FFFFFF', fuzz = 0)
    plot.new()
    rasterImage(img_with_inset, 0, 0, 1, 1)
    print(paste0(tree_name, '.png'))
    image_write(img_with_inset, paste0(tree_name, '.png'), format = 'png')
    image_write(img_with_inset, paste0(tree_name, '.eps'), format = 'eps')
    image_write(img_with_inset, paste0(tree_name, '.pdf'), format = 'pdf')
  }
}

print('TreePlot')

