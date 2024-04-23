# Script for processing the Bill et al data for comparison with TimeLapse-seq. 
# Output delivers only the macrophage subpopulations 

library(tidyverse)
library(Matrix)

# Set directories
in_dir <- 'CHANGE/PATH/TO/GSE234933'


# Read annotation files and raw data
mac_notes <- read_delim(file.path(in_dir, 'MGH_HNSCC_cell_annotation.txt'))
samp_notes <- read_delim(file.path(in_dir, 'MGH_HNSCC_sample_annotation.txt'))

test <- read_rds(file.path(in_dir, 'GSE234933_MGH_HNSCC_gex_raw_counts', 'HN1_gex_raw_counts.rds'))


# Filter for macrophages --------------------------------------------------

mac_notes <- mac_notes %>% 
      filter(major_state == 'Macrophages')

# Iterate through the samples and attempt to build a big ass master matrix
files <- list.files(file.path(in_dir, 'GSE234933_MGH_HNSCC_gex_raw_counts'))

data_all <- matrix()

for (file in files) {
      if (nrow(data_all) == 1) {
            data_all <- read_rds(file.path(in_dir, 'GSE234933_MGH_HNSCC_gex_raw_counts', file))
            
      } else {
            next_matrix <- read_rds(file.path(in_dir, 'GSE234933_MGH_HNSCC_gex_raw_counts', file))
            data_all <- cbind(data_all, next_matrix)
      }
 }

data_all <- data_all[,mac_notes$sample_barcode]

saveRDS(object = data_all, file = file.path(in_dir, 'analysis_out', 'macrophage_counts.rds'))


