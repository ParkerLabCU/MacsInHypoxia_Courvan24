# Script plots CS ratios to PCs and calculates correlations for those parameters

library(tidyverse)

#load cs and PCA data 

PCA <- read_csv('PATH/TO/PCA/OUTPUT') %>% 
      column_to_rownames(var = 'subPop')

# compute correlation matrix 
pearsons = cor(PCA, method = 'pearson')

# make plots for dim 1 and dim 2

dim1Plot <- ggplot(data = PCA, aes(x = Dim.1, y = cs)) +
      geom_point() +
      geom_smooth(method = 'lm', color = '#75DDDD') +
      theme_classic()+
      theme(
            axis.title = element_blank(),
            axis.line = element_line(lineend = 'round'),
            axis.ticks = element_line(lineend = 'round'),
            axis.text = element_text(family = 'helvetica', color = 'black')
      )

print(dim1Plot)



dim2Plot <- ggplot(data = PCA, aes(x = Dim.2, y = cs)) +
      geom_point() +
      geom_smooth(method = 'lm', color = '#75DDDD') +
      theme_classic() +
      theme(
            axis.title = element_blank(),
            axis.line = element_line(lineend = 'round'),
            axis.ticks = element_line(lineend = 'round'),
            axis.text = element_text(family = 'helvetica', color = 'black')
      )

print(dim2Plot)


# save output

ggsave('PATH/FOR/OUTPUT/Dim1_cs.png',
       plot = dim1Plot, dpi = 1000, units = 'mm', height = 20, width = 20, scale = 2)


ggsave('PATH/FOR/OUTPUT/Dim2_cs.png',
       plot = dim2Plot, dpi = 1000, units = 'mm', height = 20, width = 20, scale = 2)


