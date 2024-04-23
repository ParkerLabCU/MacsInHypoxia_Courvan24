# This script takes the output from enrichr on the TRRUST TF enrichment and plots the pvalue and odds ratio for each condition
library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)


#set directories and load in the data

out_dr <- 'PATH/FOR/OUTPUT'

#set plot dimensions
plot_height <- 35
plot_width <- 60

# Set color
point_color <- '#75DDDD'

#load M0M1 data
M0M1ksynUP <- read.table('PATH/TO/ENRICHR/OUTUPUT', sep = '\t', header = TRUE) %>% 
      separate_wider_delim(cols = Term, delim = ' ', names = c('TF', 'species')) %>% 
      filter(species == 'human')

# Plot TRRUST TFs p value and odds ratio for M0M1 comparison
M0M1_TFs_UP <- ggplot(data = M0M1ksynUP[1:10,]) +
      geom_point(aes(x = -log10(P.value), y = reorder(TF, -log10(P.value)), size = Odds.Ratio), color = point_color) +
      theme_light() +
      expand_limits(x = 0) +
      scale_size("Odds Ratio", range = c(1,5)) +
      theme(
            # axis.text.x = element_blank(),
            axis.text = element_text(family = 'helvetica', size = 10, color = 'black'),
            axis.title = element_blank(),
            # axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(lineend = 'round'),
            axis.line = element_line(lineend = 'round'),
            legend.title = element_text(family = 'helvetica', size = 10),
            legend.box.margin = margin(t=0,r=0,b=0,l=-10),
      ) 


print(M0M1_TFs_UP)
ggsave(file.path(out_dr,'M0M1_TFs_UP.png'), height = plot_height, width = plot_width, units = 'mm', dpi = 1000, scale = 2)

# Load data for M0M0H comparison
M0M0H_ksynUP <- read.table('PATH/TO/ENRICHR/OUTUPUT', sep = '\t', header = TRUE) %>% 
      separate_wider_delim(cols = Term, delim = ' ', names = c('TF', 'species')) %>% 
      filter(species == 'human')


# Plot TRRUST TFs p value and odds ratio for M0M0H comparison
M0M0H_TFs_UP <- ggplot(data = M0M0H_ksynUP[1:10,]) +
      geom_point(aes(x = -log10(P.value), y = reorder(TF, -log10(P.value)), size = Odds.Ratio), color = point_color) +
      theme_light() +
      expand_limits(x = 0) +
      scale_size("Odds Ratio", range = c(1,5)) +
      theme(
            # axis.text.x = element_blank(),
            axis.text = element_text(family = 'helvetica', size = 10, color = 'black'),
            axis.title = element_blank(),
            # axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(lineend = 'round'),
            axis.line = element_line(lineend = 'round'),
            legend.title = element_text(family = 'helvetica', size = 10),
            legend.box.margin = margin(t=0,r=0,b=0,l=-10),
      ) 


print(M0M0H_TFs_UP)
ggsave(file.path(out_dr,'M0M0H_TFs_UP.png'), height = plot_height, width = plot_width, units = 'mm', dpi = 1000, scale = 2)

# Load data for M1M1H comparison
M1M1H_ksynUP <- read.table('PATH/TO/ENRICHR/OUTUPUT', sep = '\t', header = TRUE) %>% 
      separate_wider_delim(cols = Term, delim = ' ', names = c('TF', 'species')) %>% 
      filter(species == 'human')


# Plot TRRUST TFs p value and odds ratio for M1M1H comparison
M1M1H_TFs_UP <- ggplot(data = M1M1H_ksynUP[1:12,]) +
      geom_point(aes(x = -log10(P.value), y = reorder(TF, -log10(P.value)), size = Odds.Ratio), color = point_color) +
      theme_light() +
      expand_limits(x = 0) +
      scale_size("Odds Ratio", range = c(1,5)) +
      theme(
            # axis.text.x = element_blank(),
            axis.text = element_text(family = 'helvetica', size = 10, color = 'black'),
            axis.title = element_blank(),
            # axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(lineend = 'round'),
            axis.line = element_line(lineend = 'round'),
            legend.title = element_text(family = 'helvetica', size = 10),
            legend.box.margin = margin(t=0,r=0,b=0,l=-10),
      ) 


print(M1M1H_TFs_UP)
ggsave(file.path(out_dr,'M1M1H_TFs_UP.png'), height = plot_height, width = plot_width, units = 'mm', dpi = 1000, scale = 2)
