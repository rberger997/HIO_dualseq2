# Find pathway similarities and differences between samples

library(here)
library(dplyr)
library(ggplot2)
library(magrittr)


# Function to make heatmap:
# 1. Filter reactome enrichment from 3 samples by pathways significant in at least 1 sample.
# 2. Make cluster heatmap of NES values from filtered pathways list

pathway_filter <- function(d1, d2, d3, 
                            label1, label2, label3, p = 1){
  
  # Get list of pathways - Must be significant in at least one sample
  paths <- NULL
  
  for(i in list(d1, d2, d3)){
    a <- filter(i, pvalue <= p) %>% 
      select(Description)
    
    print(paste(nrow(a), 'significant pathways'))
    
    paths <- rbind(paths, a) %>% 
      unique()
  }
  
  print(paste(nrow(paths), 'pathways significant in at least one sample with p value cutoff at', p))
  
  
  # Get NES for each pathway from each sample
  paths1 <- paths
  
  for(i in list(d1, d2, d3)){
    a <- filter(i, Description %in% paths$Description) %>% 
      select(c(Description, NES))
    
    paths1 <- left_join(paths1, a, by = 'Description')
    
  }
  colnames(paths1)[2:4] <- c(label1,label2,label3)
  
  return(paths1)
}
  
  
path_heatmap <- function(paths1, rows = T, cutrows = 1){
  x <- as.matrix(paths1[2:4])
rownames(x) <- paths1$Description

#install.packages("matrixStats")
library(matrixStats)

x <- x[order(-rowVars(x)),]


library(pheatmap)
# Define top variance genes, extract from full data
topVarGenes <- order(rowVars(x), decreasing = TRUE)
mat <- x[topVarGenes, ]


# convert to fold over mean of all samples
#mat <- mat - rowMeans(mat)  

pheatmap(mat,
         show_rownames = rows,
         breaks = seq(-2.6, 2.7, by = .055),
         cluster_cols = F,
         cutree_rows = cutrows,
         border_color = NA)
#fontsize = 10,
#cellwidth = 20,
#cellheight = 10

}
  
  


# Load samples
gsea.dir <- here('results/DESeq2_human/GSEA/reactome/')


# 2h mutants
stm_2h <- read.csv(file.path(gsea.dir, 'GSEA_reactome_STM_2h_over_PBS_2h.csv'), stringsAsFactors = F)
spi1_2h <- read.csv(file.path(gsea.dir, 'GSEA_reactome_SPI1_2h_over_PBS_2h.csv'), stringsAsFactors = F)
spi2_2h <- read.csv(file.path(gsea.dir, 'GSEA_reactome_SPI2_2h_over_PBS_2h.csv'), stringsAsFactors = F)


# 8h mutants
stm_8h <- read.csv(file.path(gsea.dir, 'GSEA_reactome_STM_8h_over_PBS_8h.csv'), stringsAsFactors = F)
spi1_8h <- read.csv(file.path(gsea.dir, 'GSEA_reactome_SPI1_8h_over_PBS_8h.csv'), stringsAsFactors = F)
spi2_8h <- read.csv(file.path(gsea.dir, 'GSEA_reactome_SPI2_8h_over_PBS_8h.csv'), stringsAsFactors = F)



# 2h serovars
stm_2h <- read.csv(file.path(gsea.dir, 'GSEA_reactome_STM_2h_over_PBS_2h.csv'), stringsAsFactors = F)
se_2h <- read.csv(file.path(gsea.dir, 'GSEA_reactome_SE_2h_over_PBS_2h.csv'), stringsAsFactors = F)
st_2h <- read.csv(file.path(gsea.dir, 'GSEA_reactome_ST_2h_over_PBS_2h.csv'), stringsAsFactors = F)


# 8h serovars
stm_8h <- read.csv(file.path(gsea.dir, 'GSEA_reactome_STM_8h_over_PBS_8h.csv'), stringsAsFactors = F)
se_8h <- read.csv(file.path(gsea.dir, 'GSEA_reactome_SE_8h_over_PBS_8h.csv'), stringsAsFactors = F)
st_8h <- read.csv(file.path(gsea.dir, 'GSEA_reactome_ST_8h_over_PBS_8h.csv'), stringsAsFactors = F)



# Look at specific groups from the pathway heatmaps
# 2h serovars 
ser2h <- pathway_filter(stm_2h, se_2h, st_2h, p = 0.05, 
                        label1 = 'STM_2h', label2 = 'SE_2h', label3 = 'ST_2h')
head(ser2h)

# pathways up in STM and SE but down in ST
st_paths <- filter(ser2h, ST_2h < 0 & STM_2h > 0 & SE_2h > 0)
path_heatmap(st_paths)
# pathways down in STM and SE but up in ST
st_paths1 <- filter(ser2h, ST_2h > 0 & STM_2h < 0 & SE_2h < 0)
path_heatmap(st_paths1)
# Pathways up in all three
ser_paths1 <- filter(ser2h, ST_2h > 1.4 & STM_2h > 1.4 & SE_2h > 1.4)
path_heatmap(ser_paths1)




# 8h serovars
ser8h <- pathway_filter(stm_8h, se_8h, st_8h, p = 0.05,
                        label1 = 'STM_8h', 'SE_8h', 'ST_8h')

# Pathways down in STM and SE but up in ST
st_dn8h <- filter(ser8h, ST_8h < 0 & STM_8h > 0 & SE_8h > 0)
path_heatmap(st_dn8h)

# Pathways down in STM and SE but up in ST
st_up8h <- filter(ser8h, ST_8h > 0 & STM_8h < 0 & SE_8h < 0)
path_heatmap(st_up8h)

# Pathways up in all three
ser_paths8 <- filter(ser8h, ST_8h > 1.5 & STM_8h > 1.5 & SE_8h > 1.5)
path_heatmap(ser_paths8)



# STM Mutants
# 2h mutants

mut2h <- pathway_filter(stm_2h, spi1_2h, spi2_2h, p = 0.05,
                        label1 = 'STM_2h', 'SPI1_2h', 'SPI2_2h')

# Pathways up in all three
nes <- 1.4
mut_paths2 <- filter(mut2h, STM_2h > nes & SPI1_2h > nes & SPI2_2h > nes)
path_heatmap(mut_paths2)

# Pathways down in SPI-1 and up in STM and SPI2
a <- filter(mut2h, STM_2h > 0 & SPI2_2h > 0 & SPI1_2h < 0)


# Pathways up in SPI-1 and down in STM and SPI-2
b <- filter(mut2h, STM_2h < 0 & SPI1_2h > 0 & SPI2_2h < 0)
path_heatmap(b)




# 8h mutants
mut8h <- pathway_filter(stm_8h, spi1_8h, spi2_8h, p = 0.05,
                        label1 = 'STM_8h', label2 = 'SPI1_8h', label3 = 'SPI2_8h')

# Pathways up in all three
nes <- 1.4
mut_paths8 <- filter(mut8h, STM_8h > nes & SPI1_8h > nes & SPI2_8h > nes)
path_heatmap(mut_paths8)

# Pathways down in SPI-2 and up in STM
c <- filter(mut8h, STM_8h > 0 & SPI1_8h < 0)
path_heatmap(c)


# Pathways down in SPI-2 and up in STM
d <- filter(mut8h, STM_8h > 0 & SPI2_8h < 0)
path_heatmap(d)
