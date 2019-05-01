# Validation of REACTOME pathway significance

# We want to validate the selected Reactome pathways for p value significance when comparing different samples together. The initial filtering criteria used p value cutoff but that is relative to PBS. We want to calculate significance relative to STM for comparison.

# Steps involved:
# 1. Differential expression analysis for sample over STM (already done)
# 2. Run GSEA REACTOME for diffexpress files (already done)
# 3. From selected list of pathways, look at p value significance
# 4. Validate or remove pathways based on significance

#####################################################

library(here)
library(dplyr)
library(magrittr)

# Get list of pathways to validate (from blocks in Figure 3)
# SPI2 mutant 2h pathways
paths <- read.csv(here('img/stm_mutants/mut_2h_gsea_diff_heatmap_spi2.csv'))
head(paths)



######################################################

# Load in the reactome GSEA files for spi2 over stm at 2h
spi2.2h <- read.csv(here('results/DESeq2_human/GSEA/reactome_stm/GSEA_reactome_SPI2_2h_over_STM_2h.csv'))


# Filter out only the pathways from the SPI2 block in Figure 3
spi2.2h <- dplyr::filter(spi2.2h, Description %in% paths$Description)
