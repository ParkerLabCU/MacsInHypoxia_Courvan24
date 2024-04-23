# Script for plotting individual enrichment scores for specific terms in both datasets

library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)

# Load in the data 

M0M0H_UP <- read_delim('PATH/TO/GSEA/OUTPUT')

M0M0H_DOWN <- read_delim('PATH/TO/GSEA/OUTPUT')

M1M1H_UP <- read_delim('PATH/TO/GSEA/OUTPUT')

M1M1H_DOWN <- read_delim('PATH/TO/GSEA/OUTPUT')

M0M0H <- bind_rows(M0M0H_UP, M0M0H_DOWN) %>% 
      mutate(treatment = 'Resting')

M1M1H <- bind_rows(M1M1H_UP, M1M1H_DOWN) %>% 
      mutate(treatment = 'Inflammed')


# combine the dataframes, trim the labels and add nominal p value 
combined <- bind_rows(M0M0H, M1M1H) %>% 
      mutate(
            NAME = gsub('_', ' ', NAME),
            NAME = gsub('\\..*', '', NAME),
            `NOM p-val` = `NOM p-val` + 10^-9,
            sig = case_when(`NOM p-val` > 0.05 ~ 'not',
                            `NOM p-val` < 0.05 ~ 'sig')
      ) %>% 
  dplyr::arrange(desc(NAME)) %>% 
  mutate(
    plot_order = row_number(NAME)
  )


# make the plot 

GSEA_plot <- ggplot(data = combined) +
      scale_shape_manual(values = c(8,16)) +
      scale_fill_manual(values = c('#CA7DF9','#306BAC')) +
      scale_color_manual(values = c('#CA7DF9','#306BAC')) +
      scale_alpha(range = c(0.2, 0.8)) +
      geom_bar(aes(x = NES, y = reorder(NAME, -plot_order), fill = treatment, color = treatment), stat = 'identity', position = position_dodge(), width = 0.5, alpha = 0.7) +
      geom_point(aes(x = NES, y = reorder(NAME, -plot_order), fill = treatment, shape = sig), size = 2, position = position_dodge(width = .5)) +
      geom_vline(xintercept = 0) +
      theme_light() +
      theme(
            axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = 'none',
            
      )

print(GSEA_plot)

ggsave('PATH/FOR/OUTPUT/GSEA_combined.png', 
       plot = GSEA_plot, dpi = 1000, units = 'mm', height = 60, width = 40, scale = 2)

