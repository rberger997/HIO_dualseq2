# Find pathway similarities and differences between samples
library(here)
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(magrittr)



# Look at STM mutants first

# Load samples
stm2 <- read.csv(here('results/DESeq2_human/GSEA/reactome/GSEA_reactome_STM_2h_over_PBS_2h.csv'), stringsAsFactors = F) %>% 
  filter(pvalue < 0.05)
spi1_2h <- read.csv(here('results/DESeq2_human/GSEA/reactome/GSEA_reactome_SPI1_2h_over_PBS_2h.csv'), stringsAsFactors = F)%>% 
  filter(pvalue < 0.05)
spi2_2h <- read.csv(here('results/DESeq2_human/GSEA/reactome/GSEA_reactome_SPI2_2h_over_PBS_2h.csv'), stringsAsFactors = F)%>% 
  filter(pvalue < 0.05)

# Find pathways significant in all 
shared <- stm2$Description [stm2$Description%in%spi1_2h$Description [spi1_2h$Description%in%spi2_2h$Description]]


# Extract NES and pathway from each sample and combine together
x <- filter(stm2, Description %in% shared) %>% 
  select(c(Description, NES, label))

y <- filter(spi1_2h, Description %in% shared) %>% 
  select(c(Description, NES, label))

z <- filter(spi2_2h, Description %in% shared) %>% 
  select(c(Description, NES, label))

paths <- rbind(x, y) %>% 
  rbind(z)






df <- group_by(paths, Description) %>% 
  summarise(avg.NES = mean(NES),
            sd = sd(NES)) %>% 
  arrange(-avg.NES)





# Function for loading data, selecting shared significant pathways from all

avg.reactome <- function(a,b,c){
# Load samples
stm2 <- read.csv(here(paste0('results/DESeq2_human/GSEA/reactome/GSEA_reactome_',a,'.csv')), stringsAsFactors = F) %>% 
  filter(pvalue < 0.05)
spi1_2h <- read.csv(here(paste0('results/DESeq2_human/GSEA/reactome/GSEA_reactome_',b,'.csv')), stringsAsFactors = F)%>% 
  filter(pvalue < 0.05)
spi2_2h <- read.csv(here(paste0('results/DESeq2_human/GSEA/reactome/GSEA_reactome_',c,'.csv')), stringsAsFactors = F)%>% 
  filter(pvalue < 0.05)

# Find pathways significant in all 
shared <- stm2$Description [stm2$Description%in%spi1_2h$Description [spi1_2h$Description%in%spi2_2h$Description]]

# Extract NES and pathway from each sample and combine together
x <- filter(stm2, Description %in% shared) %>% 
  select(c(Description, NES, label))

y <- filter(spi1_2h, Description %in% shared) %>% 
  select(c(Description, NES, label))

z <- filter(spi2_2h, Description %in% shared) %>% 
  select(c(Description, NES, label))

paths <- rbind(x, y) %>% 
  rbind(z)

return(paths)
}



# Function for making barplots of average NES for pathways
reactome.barplot <- function(df){
# Make barplot
  ggplot(df, 
         aes(x = reorder(Description, NES), 
             y = NES, fill = as.factor(label)))+
    geom_bar(stat = 'identity', position = 'dodge')+
    coord_flip()+
    theme(panel.grid.minor = element_line(linetype = "blank"))+
    theme_bw()+
    scale_fill_discrete(name = "Injection")+
    labs(x = '')+
    ylim(-2.5, 2.6)+
    geom_hline(yintercept = 0)
  
}


# Function for removing redundant pathways from ggplot
remove.pathways <- function(df, paths){
  for(i in seq(paths)){
    df <- dplyr::filter(df, Description != paths[i])
  }
  return(df)
}


# STM mutants 2h  
mut.2h <- avg.reactome(a = 'STM_2h_over_PBS_2h',
                       b = 'SPI1_2h_over_PBS_2h',
                       c = 'SPI2_2h_over_PBS_2h')

# Make list of pathways to remove from barplot
x1 <- c('Signaling by Interleukins',
        'Interleukin-1 family signaling',
        'Other interleukin signaling',
        'Diseases of Immune System',
        'Downstream signaling events of B Cell Receptor (BCR)',
        'Diseases associated with the TLR signaling cascade',
        'TNFs bind their physiological receptors',
        'TNF receptor superfamily (TNFSF) members mediating non-canonical NF-kB pathway',
        'Regulation of TNFR1 signaling')


mut.2h <- remove.pathways(mut.2h, x1)
m2 <- reactome.barplot(mut.2h)

m2

# STM mutants 8h
mut.8h <- avg.reactome(a = 'STM_8h_over_PBS_8h',
                       b = 'SPI1_8h_over_PBS_8h',
                       c = 'SPI2_8h_over_PBS_8h')

x2 <- c('Eukaryotic Translation Elongation',
        'Signaling by Interleukins',
        'Olfactory Signaling Pathway',
        'G alpha (s) signalling events')

mut.8h <- remove.pathways(mut.8h, x2)

m8 <- reactome.barplot(mut.8h)
m8


# Serovars 2h
ser.2h <- avg.reactome(a = 'STM_2h_over_PBS_2h',
                       b = 'SE_2h_over_PBS_2h',
                       c = 'ST_2h_over_PBS_2h')

x2 <- c('Eukaryotic Translation Elongation',
        'Signaling by Interleukins',
        'Regulation of TNFR1 signaling',
        'TNF receptor superfamily (TNFSF) members mediating non-canonical NF-kB pathway',
        'Mitochondrial translation initiation',
        'Mitochondrial translation termination',
        'Mitochondrial translation elongation')

ser.2h <- remove.pathways(ser.2h, x2)

s2 <- reactome.barplot(ser.2h)
s2 



# Serovars 8h
ser.8h <- avg.reactome(a = 'STM_8h_over_PBS_8h',
                       b = 'SE_8h_over_PBS_8h',
                       c = 'ST_8h_over_PBS_8h')

x4 <- c('Toll Like Receptor 5 (TLR5) Cascade',
        'Toll Like Receptor 3 (TLR3) Cascade',
        'Toll Like Receptor 4 (TLR4) Cascade',
        'Toll Like Receptor 7/8 (TLR7/8) Cascade',
        'Toll Like Receptor 9 (TLR9) Cascade',
        'Toll Like Receptor 10 (TLR10) Cascade',
        'Toll Like Receptor 2 (TLR2) Cascade',
        'Toll Like Receptor TLR1:TLR2 Cascade',
        'Toll Like Receptor TLR6:TLR2 Cascade',
        'G alpha (s) signalling events',
        'rRNA processing in the nucleus and cytosol',
        'Striated Muscle Contraction',
        'Interleukin-1 family signaling')

ser.8h <- remove.pathways(ser.8h, x4)

s8 <- reactome.barplot(ser.8h)
s8








library(gridExtra)

png(filename = here('figures/gsea_comps1.png'),
    width = 16, height = 8, units = 'in', res = 500)

layout <- rbind(c(1, 1, 1, 2, 2, 2, 2),
                c(1, 1, 1, 2, 2, 2, 2))

grid.arrange(s2+ggtitle('Serovars 2h similar pathways'),
             s8+ggtitle('Serovars 8h similar pathways'),
             layout_matrix = layout,
             respect = F)

dev.off() 



# Save pathway lists
library(readr)
library(tidyr)

ser.2h %>% 
  spread(key = label, value = NES) %>% 
  write_csv('~/Desktop/serovars_2h_similar_pathways.csv')

ser.8h %>% 
  spread(key = label, value = NES) %>% 
  write_csv('~/Desktop/serovars_8h_similar_pathways.csv')

