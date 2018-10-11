#'---
#' title: "RNA seq workflow part 8 - GSEA REACTOME gene sets"
#' subtitle: "HIOs dual seq experiment 2"
#' author: "Ryan Berger"
#' date: "2018-08-15"
#' output: 
#'   html_document:
#'      theme: flatly
#'      highlight: tango
#'      toc: true
#'      number_sections: true
#'      toc_depth: 2	
#'      toc_float:	
#'       collapsed: false	
#'       smooth_scroll: true	  
#' ---


#' # Purpose
#' The purpose of this script is to generate a heatmap of normalized enrichment scores
#' 
#' The input needed are the differential expression files ('SE_2h_over_PBS_2h_diffexpress.csv') and the hallmark gene set [file](http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.2/h.all.v6.2.entrez.gmt) (h.all.v6.2.entrez.gmt).
#' 
#' This script was adapted from a previous version written by David Hill.


#+ 
#' # Begin script

#' ## Libraries and directories
#+ load_pkgs, error=T, message=F, warning=F
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(GSEABase)
library(here)
library(knitr)
library(magrittr)
library(org.Hs.eg.db)
library(plotly)
library(RColorBrewer)
library(ReactomePA)
library(rmarkdown)
library(stringi)


# Directory for output (in results folder)
dir.create(path =  here("results/DESeq2_human/GSEA/reactome/"), 
           recursive = TRUE)


#+ set_cache_dir, include=F
# Set a directory for all cache files to go
opts_chunk$set(cache.path = here('results/DESeq2_human/src_html_output/knitr_cache/'))

#-----------------------------------------------------------------------------

#' ## Run GSEA using REACTOME gene sets
#' This analysis is based on running gene set enrichment using differential expression files as input. These are the steps:
#' 
#' * Load differential expression file.
#' * Create vector of descending log2FoldChange values named by gene symbols.
#' * Run GSEA using the vector and `gsePathway` function.
#' * Save results as .csv file.
#' * Combine all results to single dataframe.
#' 
#' We'll write a loop to iterate these steps over each differential expression file and combine results together.



# Define list of diff expression files to run   
files <- list.files(here('results/DESeq2_human/diff_expression/'))
files


#+ loop_reactome, cache=T, warning=F
# Create empty object to put full data into
full.data <- NULL


# Iterate over all files - run reactome GSEA, save output
for(i in seq(files)){
  
  # Sample being analyzed
  sample <-  gsub('_diffexpress.csv','',files[i])
  
  # Load diff expression file, sort by descending fold change
  input <- read.csv(file = here('results/DESeq2_human/diff_expression',
                                files[i])) %>%
    arrange(-log2FoldChange)
  
  # Create vector of descending list of gene changes, named by symbol
  up.list <- input$log2FoldChange
  names(up.list) <- input$entrez
  
  # REACTOME GSEA
  gmt.gsea <- gsePathway(geneList     = up.list,
                         nPerm        = 1000,
                         minGSSize    = 10,
                         pvalueCutoff = 1,
                         verbose      = TRUE)
  
  # Convert to dataframe, add columns with sample labels
  temp <- as.data.frame(gmt.gsea) %>% 
    mutate(sample = sample, 
           label = gsub('_.*','',files[i]), 
           time = stri_extract_first_regex(files[i], '[0-9]h')) 
  
  
  # Save results
  out.dir <- here('results/DESeq2_human/GSEA/reactome/')
  write.csv(temp, 
            file = paste0(out.dir,'GSEA_reactome_',sample,'.csv'),
            row.names = F)
  
  # Combine into full dataframe
  full.data <- rbind(full.data, temp)
  
  # Monitor progress
  print(paste(i, 'of', length(files),'done:',files[i] ))
}


# Save the full data
write.csv(full.data, here('results/DESeq2_human/GSEA/reactome/GSEA_reactome_all.csv'), row.names = F)

# Load the full data
# full.data <- read.csv(here('results/DESeq2_human/GSEA/reactome/GSEA_reactome_all.csv'))


#' ## Make heatmap of GSEA REACTOME 

# Get samples for heatmap - top 50 pathways by average significance
x <- group_by(full.data, Description) %>% 
  summarise(avg.adjp = mean(p.adjust)) %>% 
  #arrange(-avg.NES) %>% 
  arrange(avg.adjp) %>% 
  head(., 50)


# Set sample order for heatmap
full.data$label <- factor(full.data$label, 
                          levels = c('SE','ST','STM','SPI1','SPI2'))


data2 <- filter(full.data, Description %in% x$Description)



# Set order for heatmap - pathways by descending average NES
NES.avg <- group_by(data2, Description) %>% 
  summarise(NES_avg = mean(NES)) %>% 
  arrange(NES_avg)

data2$Description <- factor(data2$Description, levels = NES.avg$Description)


# Round NES (remove decimals for ggplotly tooltip)
data2$NES <- round(data2$NES, 2)

# Set up palette of colors for heatmap
hm.palette <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')), space='Lab')

#+ figure, fig.height = 8.5, fig.width = 7
# Heatmap

react_heatmap <- function(input){
ggplot(input, aes(label, Description)) + 
  geom_tile(aes(fill = NES), colour = "white") + 
  scale_fill_gradientn(colors = hm.palette(100))+ 
  theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45),
        strip.text.x = element_text(size = 12, face = 'bold'),
        strip.text.y = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = 'bold'))+
  ggtitle('REACTOME gene set enrichment')+
  labs(x='', y = 'Pathway',
       subtitle = 'Top 50 pathways by avg. adjusted p')+
  facet_grid(cols = vars(time), scales = 'free')
}

p <- react_heatmap(data2)
p

#' ## Interactive heatmap
#+ interactive, fig.height = 8.5, fig.width = 7
ggplotly(p)


#' ## Save png of plot
#+ save, eval=F
png(filename = here("/img/GSEA_hallmark_heatmap.png"),
    width =7, height = 8.5, units = 'in', res = 300)
p
dev.off()



#' # Split samples
#' 

split_samples <- function(df, samples){
mut <- dplyr::filter(df, label %in% samples)
  

subset <- dplyr::group_by(mut, Description) %>% 
  dplyr::mutate(avg_adjp = mean(p.adjust))


# Get samples for heatmap - top 50 pathways by average significance
x <- group_by(mut, Description) %>% 
  summarise(avg.adjp = mean(p.adjust)) %>% 
  arrange(avg.adjp) %>% 
  head(., 50)

mut <- dplyr::filter(mut, Description %in% x$Description)



# Set order for heatmap - pathways by descending average NES
NES.avg <- group_by(mut, Description) %>% 
  summarise(NES_avg = mean(NES)) %>% 
  arrange(NES_avg)

mut$Description <- factor(mut$Description, levels = NES.avg$Description)

return(mut)
}

# STM mutants heatmap
split_samples(df = full.data, samples = c('STM','SPI1','SPI2')) %>% 
  react_heatmap()


# Serovars heatmap
split_samples(df = full.data, samples = c('STM','SE','ST')) %>% 
  react_heatmap()

#+ render, include=F
# Render source file to html 
# dir.create(here('results/DESeq2_human/src_html_output'))

# render.dir <- here('results/DESeq2_human/src_html_output/')

# render(here('src/Hs_align_src/07_gsea_hallmark.R'), output_dir = render.dir, intermediates_dir = render.dir, clean = TRUE)


# Run Reactome GSEA
STM_PBS.reactome.gsea <- gsePathway(geneList     = up_list,
                                    nPerm        = 1000,
                                    minGSSize    = 10,
                                    pvalueCutoff = 1,
                                    verbose      = TRUE)





# Heatmap for anna-lisa

full.data <- read.csv(here('results/DESeq2_human/GSEA/reactome/GSEA_reactome_all.csv'))


z <- full.data %>% 
  filter(label %in% c('STM','ST'))


# Get samples for heatmap - top 50 pathways by average significance
n <- 20

x <- group_by(z, Description) %>% 
  summarise(avg.adjp = mean(p.adjust)) %>% 
  #arrange(-avg.NES) %>% 
  arrange(avg.adjp) %>% 
  head(., 20)


# Set sample order for heatmap
full.data$label <- factor(full.data$label, 
                          levels = c('SE','ST','STM','SPI1','SPI2'))


data2 <- filter(z, Description %in% x$Description)

# Remove long description from one that uses "or"
 data2$Description <- gsub(pattern = ' or .*','',data2$Description) %>% 
   gsub(pattern = 'Interleukin', replacement = 'IL')


# Set order for heatmap - pathways by descending average NES
NES.avg <- group_by(data2, Description) %>% 
  summarise(NES_avg = mean(NES)) %>% 
  arrange(NES_avg)

data2$Description <- factor(data2$Description, levels = NES.avg$Description)



# Round NES (remove decimals for ggplotly tooltip)
data2$NES <- round(data2$NES, 2)

# Set up palette of colors for heatmap
hm.palette <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')), space='Lab')

#+ figure, fig.height = 8.5, fig.width = 7
# Heatmap

react_heatmap <- function(input){
  ggplot(input, aes(label, Description)) + 
    geom_tile(aes(fill = NES), colour = "white") + 
    scale_fill_gradientn(colors = hm.palette(100))+ 
    theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45),
          strip.text.x = element_text(size = 12, face = 'bold'),
          strip.text.y = element_text(size = 12, face = 'bold'),
          axis.title = element_text(size = 14, face = "bold"), 
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12, face = 'bold'))+
    ggtitle('REACTOME gene set enrichment')+
    labs(x='', y = 'Pathway',
         subtitle = paste('Top',n, 'pathways by avg. adjusted p'))+
    facet_grid(cols = vars(time), scales = 'free')
}

p <- react_heatmap(data2)
p


#' ## Save png of plot
#+ save, eval=F
png(filename = here("/img/GSEA_reactome_heatmap_AL.png"),
    width =7, height = 8.5, units = 'in', res = 500)
p
dev.off()

