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
    i$core_enrichment <- gsub('\\/', ', ', i$core_enrichment)
    
    a <- filter(i, Description %in% paths$Description) %>% 
      select(c(Description, NES, pvalue, core_enrichment))
    
    paths1 <- left_join(paths1, a, by = 'Description')
    
  }
  colnames(paths1)[2:10] <- c(label1, paste0(label1, '_pvalue'),
                             paste0(label1, '_genes'),
                             label2, paste0(label2, '_pvalue'),
                             paste0(label2, '_genes'),
                             label3, paste0(label3, '_pvalue'),
                             paste0(label3, '_genes'))
                             
                             
  paths1 <- paths1[,c(1,2,5,8,3,6,9,4,7,10)]
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
  
# Function to save heatmaps and csv data from heatmaps
save_heatmaps <- function(df, dir, name){
  
  # Save csv file in subfolder
  dir.create(file.path(dir, 'reactome_heatmap_files'))
  write.csv(df, file = paste0(dir,'reactome_heatmap_files/',name, '.csv'), row.names = F)
  
  # Save heatmap
  dev.off()
  png(filename = paste0(dir, name, '.png'),
    width = 12, height = 5, units = 'in', res = 500)
  print(path_heatmap(df))
  dev.off()
  
  print(paste(name, 'saved'))
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

# Directory for saving all serovars heatmaps
ser_dir <- here('img/serovars/')


# Subset: Pathways up in all three
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
path_heatmap(ser_paths1)

save_heatmaps(ser_paths1, 
              dir = ser_dir,
              name = 'ser_2h_gsea_diff_heatmap_all_up')



# pathways disregulated in ST compared to STM
# Subset: pathways up in STM and SE but down in ST
# Remove redundant pathways (determined manually)
st_paths <- filter(ser2h, ST_2h < 0 & STM_2h > 0 & SE_2h > 0) %>% 
  filter(STM_2h_pvalue < 0.05 | ST_2h_pvalue < 0.05) %>% 
  filter(!Description %in% c('Signaling by MET',
                             'MAPK1/MAPK3 signaling',
                             'Visual phototransduction',
                             'Diseases of metabolism',
                             'Cardiac conduction',
                             'Formation of the cornified envelope'))
path_heatmap(st_paths)



## Subset: pathways down in STM and SE but up in ST
st_paths1 <- filter(ser2h, ST_2h > 0 & STM_2h < 0 & SE_2h < 0) %>% 
  filter(STM_2h_pvalue < 0.05 | ST_2h_pvalue < 0.05)
path_heatmap(st_paths1)

# Combine pathways differentially regulated in ST
st_path_diffs_2h <- rbind(st_paths, st_paths1)


save_heatmaps(st_path_diffs_2h, 
              dir = ser_dir,
              name = 'ser_2h_gsea_diff_heatmap_st')



#######################################################
# pathways disregulated in SE compared to STM
# Subset: pathways up in STM and SE but down in SE
# Remove redundant pathways (determined manually)
se_paths <- filter(ser2h, ST_2h > 0 & STM_2h > 0 & SE_2h < 0) %>% 
  filter(STM_2h_pvalue < 0.05 | SE_2h_pvalue < 0.05)


## Subset: pathways down in STM and ST but up in SE
se_paths1 <- filter(ser2h, ST_2h < 0 & STM_2h < 0 & SE_2h > 0) %>% 
  filter(STM_2h_pvalue < 0.05 | SE_2h_pvalue < 0.05)
path_heatmap(se_paths1)

# Combine pathways differentially regulated in ST
se_path_diffs_2h <- rbind(se_paths, se_paths1)


save_heatmaps(se_path_diffs_2h, 
              dir = ser_dir,
              name = 'ser_2h_gsea_diff_heatmap_se')

#######################################################


# 8h serovars
ser8h <- pathway_filter(stm_8h, se_8h, st_8h, p = 0.05,
                        label1 = 'STM_8h', 'SE_8h', 'ST_8h')


#######################################################

# Subset: Pathways up in all three
ser_paths8 <- filter(ser8h, ST_8h > 1.5 & STM_8h > 1.5 & SE_8h > 1.5) %>% 
  filter(!Description %in% c('TNFR1-induced NFkappaB signaling pathway',
                             'MyD88-independent TLR4 cascade',
                             'TRIF (TICAM1)-mediated TLR4 signaling',
                             'TNFR1-induced proapoptotic signaling',
                             'Signaling by Interleukins',
                             'Interleukin- 1 family signaling',
                             'Formation of the cornified envelope'))
path_heatmap(ser_paths8)

save_heatmaps(ser_paths8, 
              dir = ser_dir,
              name = 'ser_8h_gsea_diff_heatmap_all_up')

#######################################################


# Pathways differentially regulated in ST
# Subset: Pathways down in ST but up in STM and SE
st_dn8h <- filter(ser8h, ST_8h < 0 & STM_8h > 0 & SE_8h > 0) %>% 
  filter(STM_8h_pvalue < 0.05 | ST_8h_pvalue < 0.05) %>% 
  filter(!Description %in% c('Prefoldin mediated transfer of substrate  to CCT/TriC',
                             'Smooth Muscle Contraction',
                             'Integrin alphaIIb beta3 signaling',
                             'Integrin signaling',
                             'Formation of tubulin folding intermediates by CCT/TriC',
                             'Digestion'))
path_heatmap(st_dn8h)


# Subset: Pathways down in STM and SE but up in ST
st_up8h <- filter(ser8h, ST_8h > 0 & STM_8h < 0 & SE_8h < 0) %>% 
  filter(STM_8h_pvalue < 0.05 | ST_8h_pvalue < 0.05)
path_heatmap(st_up8h)

# Combine ST 8h differentially regulated pathways
st_path_diffs_8h <- rbind(st_dn8h, st_up8h)

save_heatmaps(st_path_diffs_8h, 
              dir = ser_dir,
              name = 'ser_8h_gsea_diff_heatmap_st')

#######################################################

# Pathways differentially regulated in SE
# Subset: Pathways down in SE but up in STM and ST
se_dn8h <- filter(ser8h, SE_8h < 0 & STM_8h > 0 & ST_8h > 0) %>% 
  filter(STM_8h_pvalue < 0.05 | SE_8h_pvalue < 0.05) 

path_heatmap(se_dn8h)


# Subset: Pathways down in STM and ST but up in SE
se_up8h <- filter(ser8h, SE_8h > 0 & STM_8h < 0 & ST_8h < 0) %>% 
  filter(STM_8h_pvalue < 0.05 | SE_8h_pvalue < 0.05)
path_heatmap(se_up8h)

# Combine ST 8h differentially regulated pathways
se_path_diffs_8h <- rbind(se_dn8h, se_up8h)

save_heatmaps(se_path_diffs_8h, 
              dir = ser_dir,
              name = 'ser_8h_gsea_diff_heatmap_se')





#######################################################


# STM Mutants
mut_dir <- here('img/stm_mutants/')

# 2h mutants
mut2h <- pathway_filter(stm_2h, spi1_2h, spi2_2h, p = 0.05,
                        label1 = 'STM_2h', 'SPI1_2h', 'SPI2_2h')

#######################################################
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

path_heatmap(mut_paths2)

save_heatmaps(mut_paths2, 
              dir = mut_dir,
              name = 'mut_2h_gsea_diff_heatmap_all_up')


# Pathways down in all three
down <- filter(mut2h, STM_2h < 0 & SPI2_2h < 0 & SPI1_2h < 0)
path_heatmap(down)

#######################################################

# Pathways differentially regulated in SPI2
# Pathways down in SPI-1 and STM, up in SPI2
spi2_up <- filter(mut2h, STM_2h < 0 & SPI2_2h > 0 & SPI1_2h < 0) %>% 
  filter(STM_2h_pvalue < 0.05 | SPI2_2h_pvalue < 0.05)
path_heatmap(spi2_up)


# Pathways up in SPI-1 and STM, down in SPI2
spi2_dn <- filter(mut2h, STM_2h > 0 & SPI2_2h < 0 & SPI1_2h > 0) %>% 
  filter(STM_2h_pvalue < 0.05 | SPI2_2h_pvalue < 0.05) %>% 
  filter(!Description %in% c('Olfactory Signaling Pathway',
                             'Common Pathway of Fibrin Clot Formation',
                             'Visual phototransduction',
                             'Retinoid metabolism and transport',
                             'Plasma lipoprotein assembly',
                             'Intrinsic Pathway of Fibrin Clot Formation',
                             'Plasma lipoprotein assembly, remodeling, and clearance'))

spi2_dn$Description <- gsub('by Insulin-like.*', '', 
                            spi2_dn$Description)

path_heatmap(spi2_dn)


# Combine spi2 differentially regulated pathways
spi2_diff <- rbind(spi2_up, spi2_dn)

save_heatmaps(spi2_diff, 
              dir = mut_dir,
              name = 'mut_2h_gsea_diff_heatmap_spi2')


#######################################################

# Pathways differentially regulated in SPI1
# Subset: Pathways up in SPI-1 and down in STM and SPI-2
spi1_up <- filter(mut2h, STM_2h < 0 & SPI1_2h > 0 & SPI2_2h < 0) %>% 
  filter(STM_2h_pvalue < 0.05 | SPI1_2h_pvalue < 0.05) %>% 
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
path_heatmap(spi1_up)

# Subset: pathways down in SPI-1 (there aren't any)
spi1_dn <- filter(mut2h, STM_2h > 0 & SPI1_2h < 0 & SPI2_2h > 0)


# Save heatmap of SPI-1 differentially regulated pathways
save_heatmaps(spi1_up, 
              dir = mut_dir,
              name = 'mut_2h_gsea_diff_heatmap_spi1')

#######################################################

# 8h mutants
mut8h <- pathway_filter(stm_8h, spi1_8h, spi2_8h, p = 0.05,
                        label1 = 'STM_8h', label2 = 'SPI1_8h', label3 = 'SPI2_8h')

#######################################################

# Pathways up in all three
nes <- 1.4
mut_paths8 <- filter(mut8h, STM_8h > nes & SPI1_8h > nes & SPI2_8h > nes) %>% 
  filter(!Description  %in% c('Signaling by Interleukins'))
path_heatmap(mut_paths8)

save_heatmaps(mut_paths8, 
              dir = mut_dir,
              name = 'mut_8h_gsea_diff_heatmap_all_up')

#######################################################

## Pathways differentially regulated in SPI-1
# Pathways down in SPI-1
down <- filter(mut8h, STM_8h > 0 & SPI1_8h < 0 & SPI2_8h > 0) %>% 
  filter(STM_8h_pvalue < 0.05 | SPI1_8h_pvalue < 0.05) %>% 
  filter(!Description == 'Response to metal ions')
path_heatmap(down)

# Pathways up in SPI-1
up <- filter(mut8h, STM_8h < 0 & SPI1_8h > 0 & SPI2_8h < 0) %>% 
  filter(STM_8h_pvalue < 0.05 | SPI1_8h_pvalue < 0.05)
path_heatmap(up)

# Combine 
spi1_diff <- rbind(down, up)

# Save
save_heatmaps(spi1_diff, 
              dir = mut_dir,
              name = 'mut_8h_gsea_diff_heatmap_spi1')

#######################################################


## Pathways differentially regulated in SPI-2
# Pathways down in SPI-2
down <- filter(mut8h, STM_8h > 0 & SPI2_8h < 0 & SPI1_8h > 0) %>% 
  filter(STM_8h_pvalue < 0.05 | SPI2_8h_pvalue < 0.05) %>% 
  filter(!Description %in% c('Response to metal ions',
                             'Visual phototransduction'))
  
path_heatmap(down)

# Pathways up in SPI-2
up <- filter(mut8h, STM_8h < 0 & SPI2_8h > 0 & SPI1_8h < 0) %>% 
  filter(STM_8h_pvalue < 0.05 | SPI2_8h_pvalue < 0.05) %>% 
  filter(!Description %in% c('mRNA Splicing - Major Pathway',
                             'Cleavage of Growing Transcript in the Termination Region '))
path_heatmap(up)

# Combine 
spi2_diff <- rbind(down, up)

# Save
save_heatmaps(spi2_diff, 
              dir = mut_dir,
              name = 'mut_8h_gsea_diff_heatmap_spi2')

