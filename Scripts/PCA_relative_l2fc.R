# Script analyzes data from clusters in Bill et al for the CS ratio compared to mine and the gene sets in the paper
library(tidyverse)
library(Matrix)

# read in the data from Bill et all ---------------------------------------

in_dir <- 'PATH/TO/MACROPHAGE_COUNTS.RDS'

out_dir <- 'PATH/FOR/OUTPUT'

mac_notes <- read_delim(file.path(in_dir, 'MGH_HNSCC_cell_annotation.txt')) %>% 
      filter(major_state == 'Macrophages')

sample_notes <- read_delim(file.path(in_dir, 'MGH_HNSCC_sample_annotation.txt'))

mac_data <- read_rds(file.path(in_dir, 'analysis_out', 'macrophage_counts.rds'))



# Load TimeLapse Data -----------------------------------------------------

# Load deseq data l2fc
M0M1 <- read_csv('PATH/TO/NORMOXIA/INFLAMMED/DESEQ/OUTPUT') %>% 
      select(log2FoldChange, gene_ID) %>% 
      dplyr::rename(M0M1 = log2FoldChange)


M1M1H <- read_csv('PATH/TO/HYPOXIA/RESTING/DESEQ/OUTPUT') %>% 
      select(log2FoldChange, gene_ID) %>% 
      dplyr::rename(M1M1H = log2FoldChange)


M0M0H <- read_csv('PATH/TO/HYPOXIA/INFLAMMED/DESEQ/OUTPUT') %>% 
      select(log2FoldChange, gene_ID) %>% 
      dplyr::rename(M0M0H = log2FoldChange)

All_results <- inner_join(M0M1, M1M1H, by= join_by('gene_ID')) %>% 
      inner_join(M0M0H, by = join_by('gene_ID')) 



# Restructure the count Matrix into a dataframe of counts per clus --------

subPops <- unique(mac_notes$minor_state)
subPop_sums <- data.frame(row.names = rownames(mac_data) )


for (pop in subPops) {
      pop_cells <- filter(mac_notes, minor_state == pop)
      
      pop_dat <- mac_data[,pop_cells$sample_barcode]
      pop_sums <- rowSums(pop_dat)
      print(pop)
      subPop_sums[,pop] <- pop_sums
}

# Filter rows with less than 5000 counts total

keep <- rowSums(subPop_sums) > 5000

subPop_sums <- subPop_sums[keep,]

subPop_sums <- rownames_to_column(subPop_sums, var = 'gene_ID') %>% 
      select(gene_ID, Mac_SPP1, Mac_MT1H, Mac_MARK4, Mac_CCL18, Mac_APOE, Mac_CXCL10, Mac_F13A1)

# Compute differential expression for subclusters -------------------------
library(Thresher)
subPop_sums <- column_to_rownames(subPop_sums, var = 'gene_ID') + 1

Total_sums <- as.matrix(rowSums(subPop_sums)) 

Unitized_total <- unitize(Total_sums)

Unitized_subPop <- unitize(subPop_sums)

L2FC_subPop <- log2(as.data.frame(sweep(Unitized_subPop, 1, Unitized_total, FUN='/'))) %>% 
      as.data.frame() 

L2FC_subPop_means <- rowMeans(L2FC_subPop)

L2FC_subPop_meanCentered <- (L2FC_subPop - L2FC_subPop_means) 

L2FC_subPop_meanCentered <- L2FC_subPop_meanCentered %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'gene_ID')

# Do PCA on the new dataframe 
library(FactoMineR)
library(ggrepel)


L2FC_subPop <- L2FC_subPop_meanCentered %>% 
      filter(gene_ID %in% All_results$gene_ID) %>% 
      column_to_rownames(var = 'gene_ID')

theirData <- PCA(t(L2FC_subPop))

# Make a PCA plot of their sub populations 

dim1_2 <- theirData$ind$coord %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'subPop')

subpop_df <- subPop_sums %>% 
      t() %>% 
      as.data.frame() %>% 
      select(CXCL9, SPP1) %>% 
      mutate(cs = CXCL9/SPP1) %>% 
      select(cs) %>% 
      rownames_to_column(var = "subPop") %>% 
      inner_join(dim1_2) %>% 
      mutate(subPop = gsub('Mac_', 'TAM-', subPop))

write_csv(subpop_df, file.path(out_dir, 'PCAvalues.csv'))

subset_PCA <- ggplot(data= subpop_df) +
      scale_color_gradient(low = '#087E8B',  high = '#FF5A5F') +
      geom_hline(yintercept = 0, linetype = 'dashed', color = 'darkgray') +
      geom_vline(xintercept = 0, linetype = 'dashed', color = 'darkgray') +
      geom_point(aes(x = Dim.1, y = Dim.2, color = cs), size = 2) +
      geom_text_repel(aes(x = Dim.1, y = Dim.2, label = subPop), max.overlaps = 1, family = 'helvetica', size = 3) +
      expand_limits(x = c(min(subpop_df$Dim.1) - 10, max(subpop_df$Dim.1) + 10), y = c(min(subpop_df$Dim.2) - 10, max(subpop_df$Dim.2) + 10)) +
      theme_linedraw() +
      scale_x_continuous(breaks = seq(-200, 100, by= 50)) +
      theme(
            axis.text = element_blank(),
            axis.title = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = 'gray'),
            legend.position = 'none'
      )


print(subset_PCA)

ggsave(file.path(out_dir, 'SupPop_PCA.png'), subset_PCA,
       dpi = 1000, units = 'mm', width = 40, height = 40, scale = 2)

# plot the components histogram
comp_df <- theirData$eig  %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'comp')

my_colors <- c('#725D68', '#F7D08A', 'darkgrey', 'darkgrey', 'darkgrey', 'darkgrey')

components <- ggplot(data = comp_df) +
      scale_fill_manual(values = my_colors)+
      geom_bar(aes(x = comp, y = `percentage of variance`, fill = comp), stat = 'identity', alpha = 0.5) +
      theme_classic() + 
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.line = element_line(lineend = 'round', linewidth = 0.35),
            axis.ticks = element_line(linewidth = 0.35, lineend = 'round', color = 'black'),
            legend.position = 'none'
            )


plot(components)

ggsave(file.path(out_dir, 'component_plot.png'), components,
       dpi = 1000, units = 'mm', width = 20, height = 20, scale = 2)

# Get out a list of the genes most influential in PC2 

influential <- theirData$var$coord %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'gene_ID') %>% 
      select(gene_ID, Dim.1) %>% 
      arrange(Dim.1)
      

write_tsv(influential, file = file.path(out_dir, 'pc1.rnk.txt'), col_names = FALSE, )

influential <- theirData$var$coord %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'gene_ID') %>% 
      select(gene_ID, Dim.2) %>% 
      arrange(Dim.2)

write_csv(influential, file = file.path(out_dir, 'pc2.rnk.txt'), col_names = FALSE)


# Make plots of the loadings ----------------------------------------------

loads <- tibble(load1 = theirData$var$coord[,1] / sqrt(theirData$eig[1]),
                load2 = theirData$var$coord[,2] / sqrt(theirData$eig[2])
                
) 

load_plot <- ggplot(data = loads) +
      geom_vline(xintercept = 0, linetype = 'dashed') +
            geom_histogram(aes(load2), bins = 9, fill = '#F7D08A', color = '#3C3C3C', alpha = 0.5) +
      geom_histogram(aes(load1), bins = 9, fill = '#725D68', color = '#3C3C3C', alpha = 0.5) +
      theme_classic() +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.line = element_line(lineend = 'round', linewidth = 0.35),
            axis.ticks = element_line(linewidth = 0.35, lineend = 'round', color = 'black')
      )

print(load_plot)

ggsave(file.path(out_dir, 'loading_hist.png'), plot = load_plot,
       dpi = 1000, units = 'mm', width = 30, height = 20, scale = 2)


# project my data into PCA space  -----------------------------------------


# Load deseq data l2fc
M0M1 <- read_csv('PATH/TO/NORMOXIA/INFLAMMED/DESEQ/OUTPUT') %>% 
      select(log2FoldChange, gene_ID) %>% 
      dplyr::rename(M0M1 = log2FoldChange)


M1M1H <- read_csv('PATH/TO/HYPOXIA/RESTING/DESEQ/OUTPUT') %>% 
      select(log2FoldChange, gene_ID) %>% 
      dplyr::rename(M1M1H = log2FoldChange)


M0M0H <- read_csv('PATH/TO/HYPOXIA/INFLAMMED/DESEQ/OUTPUT') %>% 
      select(log2FoldChange, gene_ID) %>% 
      dplyr::rename(M0M0H = log2FoldChange)

All_results <- inner_join(M0M1, M1M1H, by= join_by('gene_ID')) %>% 
      inner_join(M0M0H, by = join_by('gene_ID')) %>% 
      filter(gene_ID %in% rownames(theirData$var$coord)) 

L2FC_subPop_means <- as.data.frame(L2FC_subPop_means) %>% 
      rownames_to_column(var = 'gene_ID')

All_results <- inner_join(All_results, L2FC_subPop_means, by = join_by(gene_ID)) %>% 
      rename(means = L2FC_subPop_means)

All_results <- mutate(All_results,
      M0M1 = M0M1 - means,
      M1M1H = M1M1H - means,
      M0M0H = M0M0H - means
)

# Set up a kernel for PCA space
kernel <- sweep(theirData$var$coord, 2, sqrt(theirData$eig[1:ncol(theirData$var$coord),1]), FUN = '/') %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'gene_ID') %>% 
      filter(gene_ID %in% All_results$gene_ID) %>% 
      arrange(gene_ID) %>% 
      column_to_rownames(var = 'gene_ID') %>% 
      as.matrix()

M0M1 <- All_results %>% 
      select(M0M1, gene_ID) %>%
      column_to_rownames('gene_ID') %>% 
      as.matrix()

M1M1H <- All_results %>% 
      select(M1M1H, gene_ID) %>%
      column_to_rownames('gene_ID') %>% 
      as.matrix()

M0M0H <- All_results %>% 
      select(M0M0H, gene_ID) %>%
      column_to_rownames('gene_ID') %>% 
      as.matrix()

project_M1M1H <- t(M1M1H) %*% kernel %>% 
      as.data.frame() %>% 
      # select(Dim.1, Dim.2) %>% 
      mutate(subPop = 'Hypoxia',
             hypoxia = 'yes')

project_M0M1 <- t(M0M1) %*% kernel %>% 
      as.data.frame() %>% 
      # select(Dim.1, Dim.2) %>% 
      mutate(subPop = '',
             hypoxia = 'LPS')


project_M0M0H <- t(M0M0H) %*% kernel %>% 
  as.data.frame() %>% 
  # select(Dim.1, Dim.2) %>% 
  mutate(subPop = '',
         hypoxia = 'rest_yes')

subpop_df <- subpop_df %>% 
      mutate(hypoxia = 'no')  %>% 
      select(subPop, Dim.1, Dim.2, hypoxia) %>% 
      bind_rows(project_M1M1H) %>% 
      bind_rows(project_M0M1) %>% 
      bind_rows(project_M0M0H)


subset_PCA_withMydata <- ggplot(data= subpop_df) +
      scale_color_manual(values = c('#F87060', '#725D68', '#306BAC', '#CA7DF9')) +
      geom_hline(yintercept = 0, linetype = 'dashed', color = 'darkgray') +
      geom_vline(xintercept = 0, linetype = 'dashed', color = 'darkgray') +
      geom_point(aes(x = Dim.1, y = Dim.2, color = hypoxia), size = 3) +
      # geom_text_repel(aes(x = Dim.1, y = Dim.2, label = subPop), max.overlaps = 2, family = 'helvetica', size = 3) +
      # expand_limits(x = c(min(subpop_df$Dim.1) - 10, max(subpop_df$Dim.1) + 10), y = c(min(subpop_df$Dim.2) - 10, max(subpop_df$Dim.2) + 10)) +
      theme_linedraw() +
      scale_x_continuous(breaks = seq(-200, 100, by= 50)) +
      theme(
            # axis.text = element_blank(),
            axis.title = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = 'gray'),
            legend.position = 'none'
      )

print(subset_PCA_withMydata)

ggsave(file.path(out_dir, 'PC2_plot.png'), plot = subset_PCA_withMydata, 
       dpi = 1000, units = 'mm', width = 40, height = 40, scale = 2)


