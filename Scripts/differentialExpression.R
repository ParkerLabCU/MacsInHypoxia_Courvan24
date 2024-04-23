## This script sets up DESEQ2 analysis of macrophages in each condition

# Set paths to the count matrix for each condition, adjusting for whether your output is a master file or raw HTSeq .txt

#Load libraries
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)


#Set input and output data directory 

in_dir <- "PATH/TO/COUNTS"
out_dir <- "PATH/FOR/OUTPUT"

# Load metadata table

expt_table <- read.csv(file.path(in_dir, "DESEQ2_EXPERIMENT_METADATA_FILE.csv"))

# Set up DESeq2 input data

ddHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable= expt_table,
                                      directory= in_dir,
                                      design= ~ condition)

# Filter out any row with fewer than 100 counts

keep <- rowSums(counts(ddHTSeq)) >= 100
ddHTSeq <- ddHTSeq[keep,]

# Filter out rRNA from input data if counts are not prefiltered
rRNA <- read.table("PATH/TO/RRNA/LIST", sep = "\t", header = TRUE)
keep <- !grepl(paste(rRNA$Approved.symbol, collapse = '|'),row.names(ddHTSeq))
ddHTSeq <- ddHTSeq[keep,]



# Plot PCA of output

rld <-  rlog(ddHTSeq, blind = FALSE)
plotPCA(rld, intgroup=c('condition'))

# Execute DESeq2 (as written for normoxia inflammed comparison)

ddHTSeq <- DESeq(ddHTSeq)
res_M0M1 <- results(ddHTSeq, contrast=c("condition", "M1", "M0"))



# Plot MA of comparisons
plotMA(res_M0M1)

#pull most relevant results

M0M1_out <- results(ddHTSeq, contrast=c("condition", "M1", "M0")) %>%
  as.data.frame %>%
  rownames_to_column("gene_ID")



write.csv(M0M1_out, file.path(out_dir, "M0M1_DESeq.csv"))



counts <- counts(estimate, normalize = TRUE)
pairs(counts)
