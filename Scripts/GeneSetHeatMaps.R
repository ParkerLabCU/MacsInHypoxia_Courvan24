# Script for making heatmaps associated with specific gene sets. 

library(tidyverse)

# Load data for hypoxic conmparisons. This data can be found preprocessed in Supplementary Table 1

M0M0H_combined <- read_csv('PATH/TO/COMBINED/OUTPUT')
M1M1H_combined <- read_csv('PATH/TO/COMBINED/OUTPUT')

# Load the grp and filter the datasets

hallmarkHypoxia <- read_delim('PATH/TO/GRP/HALLMARK_HYPOXIA.v2023.2.Hs.grp', delim = ' ', skip = 1)

M0M0H_filtered <- M0M0H_combined %>% 
  filter(gene_ID %in% hallmarkHypoxia$`#`) %>% 
  select(gene_ID, log2FoldChange, L2FC_kdeg, ksyn) %>% 
  arrange(-log2FoldChange) %>% 
  column_to_rownames(var = 'gene_ID') %>% 
  as.matrix()

M1M1H_filtered <- M1M1H_combined %>% 
  filter(gene_ID %in% hallmarkHypoxia$`#`) %>% 
  select(gene_ID, log2FoldChange, L2FC_kdeg, ksyn) %>% 
  arrange(-log2FoldChange) %>% 
  column_to_rownames(var = 'gene_ID') %>% 
  as.matrix()

# Plot for hypoxia

library(pheatmap)
library(RColorBrewer)

breaks <- seq(-2, 2, by = .1)

pheatmap(M0M0H_filtered, color = colorRampPalette(c('#087E8B', 'white','#FF5A5F'))(length(breaks)),
         breaks = breaks,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = NA,
         # cellwidth = 20,
         # cellheight = 5,
         show_colnames = FALSE,
         show_rownames = FALSE,
         legend = FALSE,
         # cutree_rows = 6,
         treeheight_row = 0,
         filename = 'PATH/FOR/OUTPUT/M0M0H_hypoxia.png',
         fontsize = 5,
         width = 2,
         height = 4
         
)


pheatmap(M1M1H_filtered, color = colorRampPalette(c('#087E8B', 'white','#FF5A5F'))(length(breaks)),
         breaks = breaks,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = NA,
         # cellwidth = 20,
         # cellheight = 5,
         show_colnames = FALSE,
         show_rownames = FALSE,
         legend = FALSE,
         # cutree_rows = 6,
         treeheight_row = 0,
         filename = 'PATH/FOR/OUTPUT/M1M1_hypoxia.png',
         fontsize = 5,
         width = 2,
         height = 4
         
)

# Plots for ECM

ReactomeECM <- read_delim('PATH/TO/GRP/REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION.v2023.2.Hs.grp', delim = ' ', skip = 1)

M0M0H_filtered <- M0M0H_combined %>% 
  filter(gene_ID %in% ReactomeECM$`#`) %>% 
  select(gene_ID, log2FoldChange, L2FC_kdeg, ksyn) %>% 
  arrange(-log2FoldChange) %>% 
  column_to_rownames(var = 'gene_ID') %>% 
  as.matrix()

M1M1H_filtered <- M1M1H_combined %>% 
  filter(gene_ID %in% ReactomeECM$`#`) %>% 
  select(gene_ID, log2FoldChange, L2FC_kdeg, ksyn) %>% 
  arrange(-log2FoldChange) %>% 
  column_to_rownames(var = 'gene_ID') %>% 
  as.matrix()


pheatmap(M0M0H_filtered, color = colorRampPalette(c('#087E8B', 'white','#FF5A5F'))(length(breaks)),
         breaks = breaks,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = NA,
         # cellwidth = 20,
         # cellheight = 5,
         show_colnames = FALSE,
         show_rownames = FALSE,
         legend = FALSE,
         # cutree_rows = 6,
         treeheight_row = 0,
         filename = 'PATH/FOR/OUTPUT/M0M0H_ECM.png',
         fontsize = 5,
         width = 2,
         height = 4
         
)


pheatmap(M1M1H_filtered, color = colorRampPalette(c('#087E8B', 'white','#FF5A5F'))(length(breaks)),
         breaks = breaks,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = NA,
         # cellwidth = 20,
         # cellheight = 5,
         show_colnames = FALSE,
         show_rownames = FALSE,
         legend = FALSE,
         # cutree_rows = 6,
         treeheight_row = 0,
         filename = 'PATH/FOR/OUTPUT/M1M1_ECM.png',
         fontsize = 5,
         width = 2,
         height = 4
         
)

# Finally make plots for ribosome 

Ribosome <- read_delim('PATH/TO/GRP/KEGG_RIBOSOME.v2023.2.Hs.grp', delim = ' ', skip = 1)

M0M0H_filtered <- M0M0H_combined %>% 
  filter(gene_ID %in% Ribosome$`#`) %>% 
  select(gene_ID, log2FoldChange, L2FC_kdeg, ksyn) %>% 
  arrange(-log2FoldChange) %>% 
  column_to_rownames(var = 'gene_ID') %>% 
  as.matrix()

M1M1H_filtered <- M1M1H_combined %>% 
  filter(gene_ID %in% Ribosome$`#`) %>% 
  select(gene_ID, log2FoldChange, L2FC_kdeg, ksyn) %>% 
  arrange(-log2FoldChange) %>% 
  column_to_rownames(var = 'gene_ID') %>% 
  as.matrix()


pheatmap(M0M0H_filtered, color = colorRampPalette(c('#087E8B', 'white','#FF5A5F'))(length(breaks)),
         breaks = breaks,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = NA,
         # cellwidth = 20,
         # cellheight = 5,
         show_colnames = FALSE,
         show_rownames = FALSE,
         legend = FALSE,
         # cutree_rows = 6,
         treeheight_row = 0,
         filename = 'PATH/FOR/OUTPUT/M0M0H_ribosome.png',
         fontsize = 5,
         width = 2,
         height = 4
         
)


pheatmap(M1M1H_filtered, color = colorRampPalette(c('#087E8B', 'white','#FF5A5F'))(length(breaks)),
         breaks = breaks,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = NA,
         # cellwidth = 20,
         # cellheight = 5,
         show_colnames = FALSE,
         show_rownames = FALSE,
         legend = FALSE,
         # cutree_rows = 6,
         treeheight_row = 0,
         filename = 'PATH/FOR/OUTPUT/M1M1_ribosome.png',
         fontsize = 5,
         width = 2,
         height = 4
         
)