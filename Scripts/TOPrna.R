# Script for describing correlation between TOP motif RNAs and degradation in hypoxia

library(tidyverse)
library(readxl)
library(rtracklayer)
library(plyranges)

# Load GTF and filter to mRNAs
GTF <- rtracklayer::import('PATH/TO/GTF') %>% 
      filter(type == 'gene') %>% 
      filter(gene_biotype == 'protein_coding') %>% 
      select(gene_id, gene_biotype) %>% 
      as.data.frame()


# Load top score and L2FCkdeg data. TOP scores from https://doi.org/10.1073/pnas.1912864117

M0M0H_kdeg <- read_csv('PATH/TO/BAKR/OUTPUT/BReffects.csv') %>% 
      mutate(treatment = 'M0M0H') %>% 
      filter(XF %in% GTF$gene_id)

M1M1H_kdeg <- read_csv('PATH/TO/BAKR/OUTPUT/BReffects.csv') %>% 
      mutate(treatment = 'M1M1H')  %>% 
      filter(XF %in% GTF$gene_id)


kdeg <- bind_rows(M0M0H_kdeg, M1M1H_kdeg)

SuperTops <- read_excel('PATH/TO/TOPS/pnas.1912864117.sd07.xlsx')

TopScore <- read_excel('PATH/TO/TOPS/pnas.1912864117.sd05.xlsx')

# Merge the data 

combined <- kdeg %>% 
      mutate(isSuper = case_when(XF %in% SuperTops$gene ~ TRUE,
                                 TRUE ~ FALSE))


#compute mann-whitney

M0M0H_mw <- wilcox_test(combined[combined$treatment == 'M0M0H',], L2FC_kdeg ~ isSuper)

M1M1H_mw <- wilcox_test(combined[combined$treatment == 'M1M1H',], L2FC_kdeg ~ isSuper)

# plot results 

colors <- c('grey', '#75DDDD')

TOPplot <- ggplot(data = combined) +
      scale_fill_manual(values = colors) +
      geom_point(aes(x = treatment, y = L2FC_kdeg, fill = isSuper), position = position_jitterdodge(jitter.width = 0.2), size = 0.05, alpha = 0.1,) + 
      geom_boxplot(aes(x = treatment, y = L2FC_kdeg, fill = isSuper), outliers = FALSE, alpha = 0.5) +
      theme_classic() +
      theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            legend.position = 'none',
            axis.line = element_line(lineend = 'round'),
            axis.ticks = element_line(lineend = 'round', color = 'black')
      )
      

print(TOPplot)

ggsave('PATH/FOR/OUTPUT/TOPplot.png',
       dpi = 3000, width = 20, height = 30, units = 'mm', scale = 2)
