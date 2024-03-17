# This code was able to install Suerat v4.2 as well as other libraries needed to reproduce the results of the 'Recursive clustering of cellular diversity in scRNA-seq' paper
# This code worked with R version 4.0.5 on Windows 10 on 3/7/24

install.packages('igraph', type = 'binary')
install.packages(c('remotes', 'assertthat', 'magick', 'jsonlite', 'renv', 'showtext'))
remotes::install_version('Matrix', version = '1.5.0')
remotes::install_version('SeuratObject', version = '4.1.2')

library(renv)
library(jsonlite)

# The following uses renv to install only Seurat and SeuratDisk versions corresponding to the renv .lock file
packages_to_install <- c('Seurat', 'SeuratDisk')
lock_content <- fromJSON("renv.lock", flatten = TRUE)
lock_content$Packages <- lock_content$Packages[names(lock_content$Packages) %in% packages_to_install]
modified_lock_json <- toJSON(lock_content, auto_unbox = TRUE, pretty = TRUE)
temp_lock_path <- tempfile(fileext = ".lock")
writeLines(modified_lock_json, temp_lock_path)
renv::restore(lockfile = temp_lock_path)
unlink(temp_lock_path)

install.packages(c('ggpubr', 'data.tree', 'networkD3', 'webshot'))

remotes::install_version('Matrix', version = '1.5.1')
