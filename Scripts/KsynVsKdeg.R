# This script generates quadrant plots of l2fc kdeg vs l2fc ksyn

library(tidyverse)
library(ggplot2)
library(ggrepel)

# Read in combined files (this data is found in Supplementary Table 1)
M0M1_combined <- read_csv('PATH/TO/COMBINED/FILES')

M0M0H_combined <- read_csv('PATH/TO/COMBINED/FILES')

M1M1H_combined <- read_csv('PATH/TO/COMBINED/FILES')


# read some configs for which genes to label 
M0M1_labels <- read_delim('PATH/TO/LABELS/M0M1_plottedGenes.txt', ' ')
M0M0H_labels <- read_delim('PATH/TO/LABELS/M0M0H_plottedGenes.txt', ' ')
M1M1H_labels <- read_delim('PATH/TO/LABELS/M1M1H_plottedGenes.txt', ' ')

# Combine data and compute parameters

M0M1_combined <- M0M1_combined %>% 
      mutate(diffexpress = case_when(log2FoldChange > 0 ~ 'up',
                                     log2FoldChange < 0 ~ 'down'),
             label = case_when(gene_ID %in% M0M1_labels$to_label ~ TRUE,
                               TRUE ~ FALSE))
                 
M0M0H_combined <- M0M0H_combined %>% 
      mutate(diffexpress = case_when(log2FoldChange > 0 ~ 'up',
                                     log2FoldChange < 0 ~ 'down'),
             label = case_when(gene_ID %in% M0M0H_labels$to_label ~ TRUE,
                               TRUE ~ FALSE)
      )

M1M1H_combined <- M1M1H_combined %>% 
      mutate(diffexpress = case_when(log2FoldChange > 0 ~ 'up',
                                     log2FoldChange < 0 ~ 'down'),
             label = case_when(gene_ID %in% M1M1H_labels$to_label ~ TRUE,
                               TRUE ~ FALSE)
      )

# set parameters and plot for M0M1

colors = c('#F9C80E', '#F87060')

M0M1_quad <- ggplot(data = M0M1_combined) +
      scale_color_manual(values = colors) +
      scale_alpha_continuous(range = c(0.01, 1)) +
      scale_size(range = c(0.1,4)) + 
      geom_hline(yintercept = 0, color = 'gray') +
      geom_vline(xintercept = 0, color = 'gray') +
      geom_abline(intercept = c(0,0), slope = 1, color = '#3C3C3C', linetype = 'dashed') +
      geom_point(aes(x=ksyn, y = L2FC_kdeg, color = diffexpress, alpha = abs(log2FoldChange), size = abs(log2FoldChange))) +
      geom_text_repel(data = M0M1_combined[M0M1_combined$label,], aes(x=ksyn, y = L2FC_kdeg, label=gene_ID), max.overlaps = 20) +
      geom_point(data = M0M1_combined[M0M1_combined$label,], aes(x=ksyn, y = L2FC_kdeg,), size = 0.1)+
      coord_equal(xlim = c(min(M0M1_combined$L2FC_kdeg, M0M1_combined$ksyn), max(M0M1_combined$L2FC_kdeg, M0M1_combined$ksyn)), 
                  ylim = c(min(M0M1_combined$L2FC_kdeg, M0M1_combined$ksyn), max(M0M1_combined$L2FC_kdeg, M0M1_combined$ksyn))) +
      theme_classic() +
      theme(legend.position = 'none',
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.line = element_line(lineend = 'round'),
            axis.ticks = element_line(lineend= 'round', color = 'black')
            
      )

print(M0M1_quad)

# set parameters and plot for M0M0H

colors = c('#F9C80E', '#306BAC')

M0M0H_quad <- ggplot(data = M0M0H_combined) +
      scale_color_manual(values = colors) +
      scale_alpha_continuous(range = c(0.01, 1)) +
      scale_size(range = c(0.1,3)) + 
      geom_hline(yintercept = 0, color = 'gray') +
      geom_vline(xintercept = 0, color = 'gray') +
      geom_abline(intercept = c(0,0), slope = 1, color = '#3C3C3C', linetype = 'dashed') +
      geom_point(aes(x=ksyn, y = L2FC_kdeg, color = diffexpress, alpha = abs(log2FoldChange), size = abs(log2FoldChange))) +
      geom_text_repel(data = M0M0H_combined[M0M0H_combined$label,], aes(x=ksyn, y = L2FC_kdeg, label=gene_ID), max.overlaps = 20) +
      geom_point(data = M0M0H_combined[M0M0H_combined$label,], aes(x=ksyn, y = L2FC_kdeg,), size = 0.1)+
      coord_equal(xlim = c(min(M0M0H_combined$L2FC_kdeg, M0M0H_combined$ksyn), max(M0M0H_combined$L2FC_kdeg, M0M0H_combined$ksyn)), 
                  ylim = c(min(M0M0H_combined$L2FC_kdeg, M0M0H_combined$ksyn), max(M0M0H_combined$L2FC_kdeg, M0M0H_combined$ksyn))) +
      theme_classic() +
      theme(legend.position = 'none',
            axis.title = element_blank(),
            # axis.text = element_blank(),
            axis.line = element_line(lineend = 'round'),
            axis.ticks = element_line(lineend= 'round', color = 'black')
      )

print(M0M0H_quad)

# set parameters and plot for M1M1H

colors = c('#F87060', '#CA7DF9')

M1M1H_quad <- ggplot(data = M1M1H_combined) +
      scale_color_manual(values = colors) +
      scale_alpha_continuous(range = c(0.01, 1)) +
      scale_size(range = c(.1,3)) + 
      geom_hline(yintercept = 0, color = 'gray') +
      geom_vline(xintercept = 0, color = 'gray') +
      geom_abline(intercept = c(0,0), slope = 1, color = '#3C3C3C', linetype = 'dashed') +
      geom_point(aes(x=ksyn, y = L2FC_kdeg, color = diffexpress, alpha = abs(log2FoldChange), size = abs(log2FoldChange))) +
      geom_text_repel(data = M1M1H_combined[M1M1H_combined$label,], aes(x=ksyn, y = L2FC_kdeg, label=gene_ID), max.overlaps = 20) +
      geom_point(data = M1M1H_combined[M1M1H_combined$label,], aes(x=ksyn, y = L2FC_kdeg,), size = 0.1)+
      coord_equal(xlim = c(min(M1M1H_combined$L2FC_kdeg, M1M1H_combined$ksyn), max(M1M1H_combined$L2FC_kdeg, M1M1H_combined$ksyn)), 
                  ylim = c(min(M1M1H_combined$L2FC_kdeg, M1M1H_combined$ksyn), max(M1M1H_combined$L2FC_kdeg, M1M1H_combined$ksyn))) +
      theme_classic() +
      theme(legend.position = 'none',
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.line = element_line(lineend = 'round'),
            axis.ticks = element_line(lineend = 'round', color = 'black')
      )

print(M1M1H_quad)

#Save plots

out_dir <- 'PATH/FOR/OUTPUT'

ggsave(file.path(out_dir, 'M0M1_quad.png'), plot = M0M1_quad, units = 'mm', width = 45, height = 45, dpi = 1000, scale = 2)
ggsave(file.path(out_dir, 'M0M0H_quad.png'), plot = M0M0H_quad, units = 'mm', width = 45, height = 45, dpi = 1000, scale = 2)
ggsave(file.path(out_dir, 'M1M1H_quad.png'), plot = M1M1H_quad, units = 'mm', width = 45, height = 45, dpi = 1000, scale = 2)
