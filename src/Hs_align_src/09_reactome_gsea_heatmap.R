# Find pathway similarities and differences between samples

library(here)
library(dplyr)
library(ggplot2)
library(magrittr)


# Function to make heatmap:
# 1. Filter reactome enrichment from 3 samples by pathways significant in at least 1 sample.
# 2. Make cluster heatmap of NES values from filtered pathways list

pathway_heatmap <- function(d1, d2, d3, 
                            label1, label2, label3, 
                            n, rows = T, p = 1, cutrows = 1){
  
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

x <- as.matrix(paths1[2:4])
rownames(x) <- paths1$Description

#install.packages("matrixStats")
library(matrixStats)

x <- x[order(-rowVars(x)),]


library(pheatmap)
# Define top variance genes, extract from full data
topVarGenes <- head(order(rowVars(x), decreasing = TRUE), n) 
mat <- x[topVarGenes, ]


# convert to fold over mean of all samples
#mat <- mat - rowMeans(mat)  

pheatmap(mat,
         show_rownames = rows,
         breaks = seq(-2.6, 2.7, by = .055),
         cluster_cols = F,
         cutree_rows = cutrows)
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





# Make zoomed out heatmaps without labels

# Mutant heatmaps
# Original heatmap used 50 for mut2h p<0.05
mut2h <- pathway_heatmap(stm_2h, spi1_2h, spi2_2h,
                         'STM', 'SPI1','SPI2', 
                         n = 1357, rows = F, p = 1, cutrows = 2)

# original heatmap used 30 for mut8h p<0.05
mut8h <- pathway_heatmap(stm_8h, spi1_8h, spi2_8h,
                         'STM', 'SPI1','SPI2', 
                         n = 1357, rows = F, p = 1, cutrows = 2)



# serovars heatmaps
# original heatmap used 80 for ser2h p<0.05
ser2h <- pathway_heatmap(stm_2h, se_2h, st_2h,
                         'STM', 'SE','ST', n = 1357, rows = F, p = 1)

# original heatmap used 100 for ser8h p<0.05
ser8h <- pathway_heatmap(stm_8h, se_8h, st_8h,
                         'STM', 'SE','ST', n = 1357, rows = F, p = 1)


# Save heatmaps
width = 3
height = 6.5
res = 500

# Save mutants as png
png(filename = here('img/stm_mutants/mut_2h_gsea_diff_heatmap.png'),
    width = width, height = height, units = 'in', res = res)
mut2h
dev.off()


png(filename = here('img/stm_mutants/mut_8h_gsea_diff_heatmap.png'),
    width = width, height = height, units = 'in', res = res)
mut8h
dev.off()


# Save serovars as png
png(filename = here('img/serovars/ser_2h_gsea_diff_heatmap.png'),
    width = width, height = height, units = 'in', res = res)
ser2h
dev.off()


png(filename = here('img/serovars/ser_8h_gsea_diff_heatmap.png'),
    width = width, height = height, units = 'in', res = res)
ser8h
dev.off()




# Make detailed heatmaps with labels (huge images) for zoom closeup

# Mutant heatmaps
# Original heatmap used 50 for mut2h p<0.05
mut2h <- pathway_heatmap(stm_2h, spi1_2h, spi2_2h,
                         'STM', 'SPI1','SPI2', n = 1357, rows = T, p = 1)

# original heatmap used 30 for mut8h p<0.05
mut8h <- pathway_heatmap(stm_8h, spi1_8h, spi2_8h,
                         'STM', 'SPI1','SPI2', n = 1357, rows = T, p = 1)



# serovars heatmaps
# original heatmap used 80 for ser2h p<0.05
ser2h <- pathway_heatmap(stm_2h, se_2h, st_2h,
                         'STM', 'SE','ST', n = 1357, rows = T, p = 1)

# original heatmap used 100 for ser8h p<0.05
ser8h <- pathway_heatmap(stm_8h, se_8h, st_8h,
                         'STM', 'SE','ST', n = 1357, rows = T, p = 1)



# Save heatmaps

# Save mutants as png
png(filename = here('img/stm_mutants/mut_2h_gsea_diff_heatmap_detailed.png'),
    width = 11, height = 185, units = 'in', res = 200)
mut2h
dev.off()


png(filename = here('img/stm_mutants/mut_8h_gsea_diff_heatmap_detailed.png'),
    width = 11, height = 185, units = 'in', res = 200)
mut8h
dev.off()


# Save serovars as png
png(filename = here('img/serovars/ser_2h_gsea_diff_heatmap_detailed.png'),
    width = 11, height = 185, units = 'in', res = 200)
ser2h
dev.off()


png(filename = here('img/serovars/ser_8h_gsea_diff_heatmap_detailed.png'),
    width = 11, height = 185, units = 'in', res = 200)
ser8h
dev.off()





