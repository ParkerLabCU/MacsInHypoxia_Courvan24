# This script generates GSEA plots for PCA of TAMs compared to timelapse samples

library(tidyverse)

# Load undirected GSEA results for dimension 1 and 2 

Dim1_gsea_neg <- read_delim('PATH/TO/GSEA/OUTPUT') %>% 
      top_n(10, wt = -NES)
Dim1_gsea_pos <- read_delim('PATH/TO/GSEA/OUTPUT') %>% 
      top_n(10, wt = NES)


Dim2_gsea_neg <- read_delim('PATH/TO/GSEA/OUTPUT') %>% 
      top_n(10, wt = -NES)
Dim2_gsea_pos <- read_delim('PATH/TO/GSEA/OUTPUT') %>% 
      top_n(10, wt = NES)


# Combine tables and modify as necessary 

Dim1_gsea <- bind_rows(Dim1_gsea_pos, Dim1_gsea_neg) %>% 
      mutate(
            NAME = gsub('_', ' ', NAME),
            NAME = gsub('REACTOME', '', NAME),
            NAME = gsub('KEGG', '', NAME),
            
            `FDR q-val` = `FDR q-val` + 10^-9) %>% 
      filter(`FDR q-val` < 0.25)

Dim2_gsea <- bind_rows(Dim2_gsea_pos, Dim2_gsea_neg) %>% 
      mutate(
            NAME = gsub('_', ' ', NAME),
            NAME = gsub('REACTOME', '', NAME),
            NAME = gsub('KEGG', '', NAME),
            NAME = gsub('RESPIRATORY ELECTRON TRANSPORT ATP SYNTHESIS BY CHEMIOSMOTIC COUPLING AND HEAT PRODUCTION BY UNCOUPLING PROTEINS', 'ELECTRON TRANSPORT COUPLING AND UNCOUPLING', NAME),
            
            `FDR q-val` = `FDR q-val` + 10^-9) %>% 
      filter(`FDR q-val` < 0.25)

# Generate plots for undirected GSEA analysis

point_color <- '#75DDDD'

Dim1_plot <- ggplot(data = Dim1_gsea) +
      geom_point(aes(x = NES, y = reorder(NAME, NES), size = -log10(`FDR q-val`)), color = point_color, alpha = 0.75) +
      theme_light() +
      expand_limits(x = c(min(Dim1_gsea$NES*1.1), max(Dim1_gsea$NES*1.1))) +
      scale_size(range = c(1,4)) +
      
      theme(
            # axis.text.x = element_blank(),
            axis.text = element_text(family = 'helvetica', size = 7.5, color = 'black'),
            axis.title = element_blank(),
            # axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(lineend = 'round'),
            axis.line = element_line(lineend = 'round'),
            legend.title = element_blank(),
            # legend.text = element_blank(),
            legend.key.size = unit(3.5,'mm')
            
      ) 

print(Dim1_plot)

Dim2_plot <- ggplot(data = Dim2_gsea) +
      geom_point(aes(x = NES, y = reorder(NAME, NES), size = -log10(`FDR q-val`)), color = point_color, alpha = 0.75) +
      theme_light() +
      expand_limits(x = c(min(Dim2_gsea$NES*1.1), max(Dim2_gsea$NES*1.1))) +
      scale_size(range = c(1,4)) +
      
      theme(
            # axis.text.x = element_blank(),
            axis.text = element_text(family = 'helvetica', size = 7.5, color = 'black'),
            axis.title = element_blank(),
            # axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(lineend = 'round'),
            axis.line = element_line(lineend = 'round'),
            legend.title = element_blank(),
            legend.text = element_blank(),
            legend.key.size = unit(3.5,'mm')
            
      ) 

print(Dim2_plot)

# save output

ggsave('PATH/FOR/OUTPUT/Dim1_gsea.png',
       plot = Dim1_plot, dpi = 1000, height = 30, width = 70, units = 'mm', scale = 2)

ggsave('PATH/FOR/OUTPUT/Dim2_gsea.png',
       plot = Dim2_plot, dpi = 1000, height = 30, width = 70, units = 'mm', scale = 2)

