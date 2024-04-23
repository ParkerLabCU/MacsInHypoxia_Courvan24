# this script plots contour lines of the ksyn vs kdeg plots. 


library(tidyverse)
library(ggplot2)
library(ggrepel)

# Read in combined files (this data is found in Supplementary Table 1)
M0M1_combined <- read_csv('PATH/TO/COMBINED/FILES')

M0M0H_combined <- read_csv('PATH/TO/COMBINED/FILES')

M1M1H_combined <- read_csv('PATH/TO/COMBINED/FILES')


# Combine data and compute parameters

M0M1_combined <- M0M1_combined %>% 
      mutate(diffexpress = case_when(log2FoldChange > 0 ~ 'up',
                                     log2FoldChange < 0 ~ 'down'),
             treatment = 'M0M1')

M0M0H_combined <- M0M0H_combined %>% 
      mutate(diffexpress = case_when(log2FoldChange > 0 ~ 'up',
                                     log2FoldChange < 0 ~ 'down'),
             treatment = 'M0M0H')
      

M1M1H_combined <- M1M1H_combined %>% 
      mutate(diffexpress = case_when(log2FoldChange > 0 ~ 'up',
                                     log2FoldChange < 0 ~ 'down'),
             treatment = 'M1M1H')
      
# Make a master dataframe

Data_all <- bind_rows(M0M1_combined, M0M0H_combined, M1M1H_combined) %>% 
      mutate(treatment = factor(treatment, levels = c('M0M1', 'M0M0H', 'M1M1H')))

# plot countours for all conditions in the same plot. 

color <- c('#F87060', '#306BAC', '#CA7DF9')

contourPlot <- ggplot(Data_all, aes(x = ksyn, y = L2FC_kdeg)) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      geom_density_2d(aes(color = treatment), linewidth = 0.2) +
      scale_color_manual(values = color) +
      stat_density2d(aes(alpha = after_stat(level), fill = treatment), geom = 'polygon') +
      scale_fill_manual(values = color) +
      scale_alpha_continuous(range = c(0.1, 0.2)) +
      theme_light() +
      theme(legend.position = 'none',
            axis.title = element_blank(),
            # axis.text = element_blank(),
            )

print(contourPlot)

ggsave('PATH/FOR/OUTPUT/countours.png', 
       dpi = 1000, units = 'mm', height = 40, width = 40, scale = 2)
