# This script generates plots of the kdegs from different conditions showing the range of possible options in macrophages. 

library(tidyverse)
library(dplyr)

out_dir <- 'PATH/FOR/OUTPUT'


Kdegs_rest <- read_csv('PATH/FOR/BAKR/OUTPUT/BReffects.csv') %>% 
      dplyr::rename(restHypox = L2FC_kdeg)

Kdegs_active <- read.csv('PATH/FOR/BAKR/OUTPUT/BReffects.csv') %>%
      filter(Exp_ID == 2) %>% 
      dplyr::rename(M1 = L2FC_kdeg)
      
Kdegs_activeHypox <- read.csv('PATH/FOR/BAKR/OUTPUT/BReffects.csv')  %>% 
      dplyr::rename(M1Hypox = L2FC_kdeg)



L2FC_frame = inner_join(Kdegs_rest, Kdegs_active, by = join_by(XF == XF)) %>% 
      inner_join(Kdegs_activeHypox, by = join_by(XF == XF)) %>% 
      select(XF, restHypox, M1, M1Hypox) %>% 
      pivot_longer(cols = c(restHypox, M1, M1Hypox), names_to= 'treatment', values_to = 'L2FC') %>% 
      mutate(treatment = factor(treatment, levels = c('M1Hypox', 'restHypox', 'M1')))

my_colors = c('#CA7DF9', '#306BAC', '#F87060')
plot <- ggplot(L2FC_frame) +
      scale_fill_manual(values=my_colors) +
      geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 0.3) +
      geom_violin(aes(y=treatment, x = L2FC, fill = treatment), draw_quantiles = 0.5, alpha = 0.8, trim = TRUE, linewidth = 0.3) +
      theme_classic() +
      theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.line = element_line(lineend = 'round'),
            axis.ticks = element_line(lineend = 'round', color = 'black'),
            legend.position = 'none'
      )
      

print(plot)

ggsave(file.path(out_dir,'L2FC_plot.png'), width = 30, height = 20, units = 'mm', dpi = 1000, scale = 2)


library(rstatix)      

kruskal <- kruskal.test(data = L2FC_frame, L2FC ~ treatment)

dunn <- dunn_test(L2FC_frame, L2FC ~ treatment)

M1_onesample <- L2FC_frame[L2FC_frame$treatment == 'M1',] %>%  wilcox_test(L2FC ~ 1, mu = 0)

M0H_onesample <- L2FC_frame[L2FC_frame$treatment == 'restHypox',] %>%  wilcox_test(L2FC ~ 1, mu = 0)

M1H_onesample <- L2FC_frame[L2FC_frame$treatment == 'M1Hypox',] %>%  wilcox_test(L2FC ~ 1, mu = 0)

