## This notebook runs bakR analysis on the 24 hr timepoint of EMC003011


# Load up bakR and read in  data

library(bakR)
library(dplyr)

# Read in the cB and filter out any samples you don't want to analyze
# Depending on the perturbation it maybe necessary to hold out the unfed control
# BakR will filter out anything it can't compare to the control and depending on what controls are available you may miss important changes

cB_24hr <- read.csv("PATH/TO/cB/FILE") %>% 
              filter(!sample %in% c("SAMPLES YOU DONT WANT TO ANALYZE"))

# Read in rRNA table

rRNA_genes <- read.table("PATH/TO/RRNA", sep = "\t", header = TRUE)

#Filter rRNA out of cB table
cB_24hr <- cB_24hr %>% 
            filter(!grepl(paste(rRNA_genes$Approved.symbol, collapse = '|'), XF))

# Load metadata csv

meta <- read.csv("PATH/TO/TIMELAPSE/METADATA", header = TRUE, row.names = "sample")

#Create the bakRData object

bakRdata <- bakRData(cB_24hr,meta)

#Run the fast model for a first pass

Fit <- bakRFit(bakRdata)

# Rerun using the hybrid Stan model

Fit <- bakRFit(Fit, HybridFit = TRUE)

# Generate a few plots from the data 
bakR::plotMA(Fit, Model = 'MLE')
bakR::plotMA(Fit, Model = 'Hybrid')


# PCA plots 
FnPCA(Fit$Fast_Fit)
FnPCA(Fit$Hybrid_Fit)

Kdegs_M0 <- Fit$Hybrid_Fit$Kdeg_df %>% 
      filter(Exp_ID == 1)

Kdegs_M1 <- Fit$Hybrid_Fit$Kdeg_df %>% 
      filter(Exp_ID == 2)

Kdegs_M1H <- Fit$Hybrid_Fit$Kdeg_df %>% 
      filter(Exp_ID == 3)

library(BSDA)

Kdegs <- inner_join(Kdegs_M1, Kdegs_M1H, by = 'XF') %>% 
      mutate(l2fc = log2(kdeg.y/kdeg.x),
             quart = as.factor(ntile(kdeg.x, 10))
      )

Kdegs_M0M1 <- inner_join(Kdegs_M0, Kdegs_M1, by = 'XF') %>% 
      mutate(l2fc = log2(kdeg.y/kdeg.x),
             quart = as.factor(ntile(kdeg.x, 10))
      )

#Save tabular output for later

out_dir <- 'PATH/FOR/OUTPUT'

write.csv(Kdegs,file.path(out_dir,'Kdegs_M1M1H.csv'))
write.csv(Kdegs_M0M1, file.path(out_dir, 'Kdegs_M0M1.csv'))
write.csv(Fit$Hybrid_Fit$Effects_df, file.path(out_dir,'BReffects_24hr.csv'))


