# Script for comparing changes in degradation rate for different RNA types


library(tidyverse)
library(rtracklayer)
library(plyranges)
library(rstatix)

# Load GTF

GTF <- rtracklayer::import('PATH/TO/GTF/GRCh38.gtf') 

# filter for genes
GTF <- GTF %>% 
      filter(type == 'gene')

# extract geneIDs and biotype

types <- GTF %>% 
      select(gene_id, gene_biotype) %>% 
      as.data.frame()

# Load the kdeg data and merge it with labels for each comparison 
kdeg_M0M1 <- read.csv('PATH/TO/BAKR/OUTPUT/BReffects.csv') %>% 
      filter(Exp_ID == 2)

kdeg_M0M0H <- read.csv('PATH/TO/BAKR/OUTPUT/BReffects.csv')

kdeg_M1M1H <- read.csv('PATH/TO/BAKR/OUTPUT/BReffects.csv')

kdeg_M0M0H <- inner_join(kdeg_M0M0H, types, by = join_by(XF == gene_id)) %>% 
      mutate(
            type = case_when(
                  gene_biotype == 'lncRNA' ~ 'lncRNA',
                  gene_biotype == 'protein_coding' ~ 'mRNA',
                  gene_biotype == 'transcribed_pseudogene' ~ 'pseudogene',
                  TRUE ~ 'snRNA'
            ),
            type = factor(type, levels = c('mRNA', 'lncRNA', 'snRNA', 'pseudogene')),
            treatment = 'M0M0H'
            
            
      )

kdeg_M1M1H <- inner_join(kdeg_M1M1H, types, by = join_by(XF == gene_id)) %>% 
      mutate(
            type = case_when(
                  gene_biotype == 'lncRNA' ~ 'lncRNA',
                  gene_biotype == 'protein_coding' ~ 'mRNA',
                  gene_biotype == 'transcribed_pseudogene' ~ 'pseudogene',
                  TRUE ~ 'snRNA'
            ),
            type = factor(type, levels = c('mRNA', 'lncRNA', 'snRNA', 'pseudogene')),
            treatment = 'M1M1H'
            
            
      )

kdeg_M0M1 <- inner_join(kdeg_M0M1, types, by = join_by(XF == gene_id)) %>% 
      mutate(
            type = case_when(
                  gene_biotype == 'lncRNA' ~ 'lncRNA',
                  gene_biotype == 'protein_coding' ~ 'mRNA',
                  gene_biotype == 'transcribed_pseudogene' ~ 'pseudogene',
                  TRUE ~ 'snRNA'
            ),
            type = factor(type, levels = c('mRNA', 'lncRNA', 'snRNA', 'pseudogene')),
            treatment = 'M0M1'
            
            
      )

# plot out the differences between samples 

typePlot_M0M0H <- ggplot(data = kdeg_M0M0H) +
      scale_alpha_discrete(range = c(1.0, 0.5)) +

      geom_jitter(aes(x= type, y = L2FC_kdeg), size = 0.05, alpha = 0.1) +
      geom_boxplot(aes(x = type, y = L2FC_kdeg, alpha = type),
                   fill = '#306BAC', outliers = FALSE) +
      theme_classic() +
      scale_fill_manual(values = colors) +
      theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            legend.position = 'none',
            axis.line = element_line(lineend = 'round'),
            axis.ticks = element_line(lineend = 'round', color = 'black')
      )

print(typePlot_M0M0H)

typePlot_M1M1H <- ggplot(data = kdeg_M1M1H) +
      scale_alpha_discrete(range = c(1.0, 0.5)) +

      geom_jitter(aes(x= type, y = L2FC_kdeg), size = 0.05, alpha = 0.1) +
      geom_boxplot(aes(x = type, y = L2FC_kdeg, alpha = type),
                   fill = '#CA7DF9', outliers = FALSE) +
      theme_classic() +
      scale_fill_manual(values = colors) +
      theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            legend.position = 'none',
            axis.line = element_line(lineend = 'round'),
            axis.ticks = element_line(lineend = 'round', color = 'black')
      )

print(typePlot_M1M1H)

typePlot_M0M1 <- ggplot(data = kdeg_M0M1) +
      scale_alpha_discrete(range = c(1.0, 0.5)) +
      
      geom_jitter(aes(x= type, y = L2FC_kdeg), size = 0.05, alpha = 0.1) +
      geom_boxplot(aes(x = type, y = L2FC_kdeg, alpha = type),
                   fill = '#F87060', outliers = FALSE) +
      theme_classic() +
      scale_fill_manual(values = colors) +
      theme(
            axis.title = element_blank(),
            # axis.text = element_blank(),
            legend.position = 'none',
            axis.line = element_line(lineend = 'round'),
            axis.ticks = element_line(lineend = 'round', color = 'black')
      )

print(typePlot_M0M1)

# save output

ggsave('PATH/FOR/OUTPUT/RNA_type_M0M0H.png', 
       plot = typePlot_M0M0H, dpi = 1000, width = 20, height = 30, units = 'mm', scale = 2)

ggsave('PATH/FOR/OUTPUT/RNA_type_M1M1H.png', 
       plot = typePlot_M1M1H, dpi = 1000, width = 20, height = 30, units = 'mm', scale = 2)

ggsave('PATH/FOR/OUTPUT/RNA_type_M0M1.png', 
       plot = typePlot_M0M1, dpi = 1000, width = 20, height = 30, units = 'mm', scale = 2)
