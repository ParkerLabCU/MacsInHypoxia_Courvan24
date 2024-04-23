# Generate Enrichment charts for hypoxic RNA stability in both activation states

library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)

# Load in the GSEA tables from differential stability in hypoxia 

M0M0H_UP <- read_delim('PATH/TO/GSEA/OUTPUT') %>% 
      top_n(10, wt = NES)

M0M0H_DOWN <- read_delim('PATH/TO/GSEA/OUTPUT') %>% 
      top_n(10, wt = -NES)

M1M1H_UP <- read_delim('PATH/TO/GSEA/OUTPUT') %>% 
      top_n(10, wt = NES)

M1M1H_DOWN <- read_delim('PATH/TO/GSEA/OUTPUT') %>% 
      top_n(10, wt = -NES)

# Combine the top 10 enriched and depleted terms and trim the pathway names. Filter FDR < 0.25
M0M0H <- bind_rows(M0M0H_DOWN, M0M0H_UP) %>% 
      mutate(
            NAME = gsub('_', ' ', NAME),
            NAME = gsub('REACTOME', '', NAME),
            NAME = gsub('KEGG', '', NAME)) %>% 
      filter(`FDR q-val` < 0.25)


M1M1H <- bind_rows(M1M1H_UP, M1M1H_DOWN) %>% 
      mutate(
            NAME = gsub('_', ' ', NAME),
            NAME = gsub('REACTOME', '', NAME),
            NAME = gsub('KEGG', '', NAME)) %>% 
      filter(`FDR q-val` < 0.25)

# Generate enrichment plots for each term with the normalized enrichment score and FDR plotted. 

point_color <- '#306BAC'

M0M0H_kdeg_GSEA <- ggplot(data = M0M0H) +
      geom_point(aes(x = NES, y = reorder(NAME, NES), size = 1/`FDR q-val`), color = point_color, alpha = 0.75) +
      theme_light() +
      expand_limits(x = c(0, max(M0M0H$NES*1.1))) +
      scale_size("FDR", breaks = c(5, 20), range = c(1,2), labels = c(0.20, 0.05)) +

      theme(
            # axis.text.x = element_blank(),
            axis.text = element_text(family = 'helvetica', size = 7.5, color = 'black'),
            axis.title = element_blank(),
            # axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(lineend = 'round'),
            axis.line = element_line(lineend = 'round'),
            legend.title = element_blank(),
            legend.text = element_blank(),
      
      ) 

print(M0M0H_kdeg_GSEA)

point_color <- '#CA7DF9'

M1M1H_kdeg_GSEA <- ggplot(data = M1M1H) +
      geom_point(aes(x = NES, y = reorder(NAME, NES), size = 1/`FDR q-val`), color = point_color, alpha = 0.75) +
      theme_light() +
      expand_limits(x = c(min(M1M1H$NES*1.1), max(M1M1H$NES*1.1))) +
      scale_size("FDR", breaks = c(5, 20, 40, 60), range = c(1,5), labels = c(0.20, 0.05, 0.025, 0.017)) +
      theme(
            # axis.text.x = element_blank(),
            axis.text = element_text(family = 'helvetica', size = 7.5, color = 'black'),
            axis.title = element_blank(),
            # axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(lineend = 'round'),
            axis.line = element_line(lineend = 'round'),
            legend.title = element_blank(),
            legend.text = element_blank(),
      ) 

print(M1M1H_kdeg_GSEA)

# Save plots

ggsave('PATH/FOR/OUTPUT/M0M0H_kdegGSEA.png',
       plot = M0M0H_kdeg_GSEA, dpi = 1000, scale = 2, width = 65, height = 30, units = 'mm')

ggsave('PATH/FOR/OUTPUT/M1M1H_kdegGSEA.png',
       plot = M1M1H_kdeg_GSEA, dpi = 1000, scale = 2, width = 65, height = 30, units = 'mm')
