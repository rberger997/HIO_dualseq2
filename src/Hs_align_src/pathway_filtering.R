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
         border_color = NA,
         fontsize = 10,
         cellwidth = 25,
         cellheight = 10,
         treeheight_row = 20)

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
# Remove redundant pathways (determined manually)
st_paths <- filter(ser2h, ST_2h < 0 & STM_2h > 0 & SE_2h > 0) %>% 
  filter(!Description %in% c('Signaling by MET',
                             'MAPK1/MAPK3 signaling',
                             'Visual phototransduction',
                             'Diseases of metabolism',
                             'Cardiac conduction',
                             'Formation of the cornified envelope'))
p2a <- path_heatmap(st_paths)



# pathways down in STM and SE but up in ST
st_paths1 <- filter(ser2h, ST_2h > 0 & STM_2h < 0 & SE_2h < 0)
p2b <- path_heatmap(st_paths1)


# Pathways up in all three
ser_paths1 <- filter(ser2h, ST_2h > 1.4 & STM_2h > 1.4 & SE_2h > 1.4) %>% 
  filter(!Description %in% c('Beta defensins',
                             'TNFR1-induced NFkappaB signaling pathway',
                             'Regulation of TNFR1 signaling',
                             'Other interleukin signaling',
                             'TNFR1-induced proapoptotic signaling',
                             'Signaling by Interleukins'))
ser_paths1$Description <- gsub(pattern = '\\(TNFSF\\) members ', 
                               replacement = '',
                               x = ser_paths1$Description)



p2c <- path_heatmap(ser_paths1)




# 8h serovars
ser8h <- pathway_filter(stm_8h, se_8h, st_8h, p = 0.05,
                        label1 = 'STM_8h', 'SE_8h', 'ST_8h')

# Pathways down in ST but up in STM and SE
st_dn8h <- filter(ser8h, ST_8h < 0 & STM_8h > 0 & SE_8h > 0) %>% 
  filter(!Description %in% c('Prefoldin mediated transfer of substrate  to CCT/TriC',
                             'Smooth Muscle Contraction',
                             'Integrin alphaIIb beta3 signaling',
                             'Integrin signaling',
                             'Formation of tubulin folding intermediates by CCT/TriC',
                             'Digestion'))
p8a <- path_heatmap(st_dn8h)


# Pathways down in STM and SE but up in ST
st_up8h <- filter(ser8h, ST_8h > 0 & STM_8h < 0 & SE_8h < 0)
p8b <- path_heatmap(st_up8h)



# Pathways up in all three
ser_paths8 <- filter(ser8h, ST_8h > 1.5 & STM_8h > 1.5 & SE_8h > 1.5) %>% 
  filter(!Description %in% c('TNFR1-induced NFkappaB signaling pathway',
                             'MyD88-independent TLR4 cascade',
                             'TRIF (TICAM1)-mediated TLR4 signaling',
                             'TNFR1-induced proapoptotic signaling',
                             'Signaling by Interleukins',
                             'Interleukin- 1 family signaling',
                             'Formation of the cornified envelope'))
p8c <- path_heatmap(ser_paths8)




# Save PNG files of each 2h heatmap subset
png(filename = here('img/serovars/ser_2h_gsea_diff_heatmap_subset1.png'),
      width = 12, height = 4, units = 'in', res = 300)
p2a
dev.off()
  

png(filename = here('img/serovars/ser_2h_gsea_diff_heatmap_subset2.png'),
    width = 12, height = 4, units = 'in', res = 300)
p2b
dev.off()

png(filename = here('img/serovars/ser_2h_gsea_diff_heatmap_subset3.png'),
    width = 12, height = 4, units = 'in', res = 300)
p2c
dev.off()



# Save PNG files of each 8h heatmap subset
png(filename = here('img/serovars/ser_8h_gsea_diff_heatmap_subset1.png'),
    width = 12, height = 4, units = 'in', res = 300)
p8a
dev.off()


png(filename = here('img/serovars/ser_8h_gsea_diff_heatmap_subset2.png'),
    width = 12, height = 4, units = 'in', res = 300)
p8b
dev.off()

png(filename = here('img/serovars/ser_8h_gsea_diff_heatmap_subset3.png'),
    width = 12, height = 4, units = 'in', res = 300)
p8c
dev.off()









# STM Mutants
# 2h mutants
mut2h <- pathway_filter(stm_2h, spi1_2h, spi2_2h, p = 0.05,
                        label1 = 'STM_2h', 'SPI1_2h', 'SPI2_2h')

# Pathways up in all three
nes <- 1.4
mut_paths2 <- filter(mut2h, STM_2h > nes & SPI1_2h > nes & SPI2_2h > nes) %>% 
  filter(!Description %in% c('TNFs bind their physiological receptors',
                             'Beta defensins',
                             'Diseases of Immune System',
                             'Diseases associated with the TLR signaling cascade',
                             'TNFR1-induced NFkappaB signaling pathway',
                             'Regulation of TNFR1 signaling',
                             'Other interleukin signaling',
                             'Signaling by Interleukins'))
mut_paths2$Description <- gsub(pattern = '\\(TNFSF\\) members ', 
                               replacement = '',
                               x = mut_paths2$Description)

p2a <- path_heatmap(mut_paths2)


# Pathways down in all three
down <- filter(mut2h, STM_2h < 0 & SPI2_2h < 0 & SPI1_2h < 0)
path_heatmap(down)


# Pathways down in SPI-2 and STM, up in SPI2
a <- filter(mut2h, STM_2h < 0 & SPI2_2h > 0 & SPI1_2h < 0)
path_heatmap(a)


# Pathways up in SPI-1 and STM, down in SPI2
x <- filter(mut2h, STM_2h > 0 & SPI2_2h < 0 & SPI1_2h > 0) %>% 
  filter(!Description %in% c('Olfactory Signaling Pathway',
                             'Common Pathway of Fibrin Clot Formation',
                             'Visual phototransduction',
                             'Retinoid metabolism and transport',
                             'Plasma lipoprotein assembly',
                             'Intrinsic Pathway of Fibrin Clot Formation',
                             'Plasma lipoprotein assembly, remodeling, and clearance'))

x$Description <- gsub('by Insulin-like.*', '', x$Description)
p2b <- path_heatmap(x)



# Pathways up in SPI-1 and down in STM and SPI-2
b <- filter(mut2h, STM_2h < 0 & SPI1_2h > 0 & SPI2_2h < 0) %>% 
  filter(!Description %in% c('Amplification of signal from the kinetochores',
                             'Amplification  of signal from unattached  kinetochores via a MAD2  inhibitory signal',
                             'Negative regulation of activity of TFAP2 (AP-2) family transcription factors',
                             'Activation of the TFAP2 (AP-2) family of transcription factors',
                             'Digestion',
                             'Activation of Nicotinic Acetylcholine Receptors',
                             'Postsynaptic nicotinic acetylcholine receptors',
                             'Ionotropic activity of kainate receptors',
                             'Activation of Ca permeable Kainate Receptor',
                             'Diseases of carbohydrate metabolism'))
p2c <- path_heatmap(b)




# 8h mutants
mut8h <- pathway_filter(stm_8h, spi1_8h, spi2_8h, p = 0.05,
                        label1 = 'STM_8h', label2 = 'SPI1_8h', label3 = 'SPI2_8h')

# Pathways up in all three
nes <- 1.4
mut_paths8 <- filter(mut8h, STM_8h > nes & SPI1_8h > nes & SPI2_8h > nes) %>% 
  filter(!Description  %in% c('Signaling by Interleukins'))
p8a <- path_heatmap(mut_paths8)

# Pathways down in SPI-2 and up in STM
c <- filter(mut8h, STM_8h > 0 & SPI1_8h < 0) %>% 
  filter(!Description == 'Response to metal ions')
path_heatmap(c)


# Pathways down in SPI-2 and up in STM
d <- filter(mut8h, STM_8h > 0 & SPI2_8h < 0) %>% 
  filter(!Description %in% c('Response to metal ions',
                             'Visual phototransduction'))
  
p8b <- path_heatmap(d)

# Pathways down in STM and SPI-1 and up in SPI2
d <- filter(mut8h, STM_8h < 0 & SPI2_8h > 0 & SPI1_8h < 0) %>% 
  filter(!Description %in% c('mRNA Splicing - Major Pathway',
                             'Cleavage of Growing Transcript in the Termination Region '))
p8c <- path_heatmap(d)



# Save PNG files of each 2h heatmap subset
png(filename = here('img/stm_mutants/mut_2h_gsea_diff_heatmap_subset1.png'),
    width = 12, height = 4, units = 'in', res = 300)
p2a
dev.off()


png(filename = here('img/stm_mutants/mut_2h_gsea_diff_heatmap_subset2.png'),
    width = 12, height = 4, units = 'in', res = 300)
p2b
dev.off()

png(filename = here('img/stm_mutants/mut_2h_gsea_diff_heatmap_subset3.png'),
    width = 12, height = 4, units = 'in', res = 300)
p2c
dev.off()



# Save PNG files of each 8h heatmap subset
png(filename = here('img/stm_mutants/mut_8h_gsea_diff_heatmap_subset1.png'),
    width = 12, height = 4, units = 'in', res = 300)
p8a
dev.off()


png(filename = here('img/stm_mutants/mut_8h_gsea_diff_heatmap_subset2.png'),
    width = 12, height = 4, units = 'in', res = 300)
p8b
dev.off()


png(filename = here('img/stm_mutants/mut_8h_gsea_diff_heatmap_subset3.png'),
    width = 12, height = 4, units = 'in', res = 300)
p8c
dev.off()
