#'---
#' title: "RNA seq workflow part 7 - GSEA hallmark gene sets"
#' subtitle: "HIOs dual seq experiment 2"
#' author: "Ryan Berger & David Hill"
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
#' The purpose of this script is to generate a heatmap of normalized enrichment scores of the 50 "hallmark pathways" in the MSigDB gene collection. According to description on the [molecular signature database site](http://software.broadinstitute.org/gsea/msigdb/collections.jsp), "Hallmark gene sets are coherently expressed signatures derived by aggregating many MSigDB gene sets to represent well-defined biological states or processes." We'll use these to get a sense of the overall global pathway changes occurring in our RNA seq samples compared to a PBS control.
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
library(readr)
library(rmarkdown)
library(stringi)


# Directory for output (in results folder)
dir.create(path =  here("results/DESeq2_human/GSEA/hallmark/"), 
           recursive = TRUE)


#+ set_cache_dir, include=F
# Set a directory for all cache files to go
opts_chunk$set(cache.path = here('results/DESeq2_human/src_html_output/knitr_cache/'))

#-----------------------------------------------------------------------------

#' ## Run GSEA using hallmark gene sets
#' This analysis is based on running gene set enrichment using differential expression files as input. These are the steps:
#' 
#' * Load differential expression file.
#' * Create vector of descending log2FoldChange values named by gene symbols.
#' * Run GSEA using the vector and hallmark gene file.
#' * Save results as .csv file.
#' * Combine all results to single dataframe
#' 
#' We'll write a loop to iterate these steps over each differential expression file and combine results together.


# load hallmark gene set
hall <- read.gmt(gmtfile = here("data/h.all.v6.2.symbols.gmt"))

# Define list of diff expression files to run   
files <- list.files(here('results/DESeq2_human/diff_expression/'))
files


#+ loop, cache=T, warning=F
# Create empty object to put full data into
full.data <- NULL


# Iterate over all files - run hallmark GSEA, save output
for(i in seq(files)){
  
  # Sample being analyzed
  sample <-  gsub('_diffexpress.csv','',files[i])

  # Load diff expression file, sort by descending fold change
 input <- read.csv(file = here('results/DESeq2_human/diff_expression',files[i])) %>%
      arrange(-log2FoldChange)
  
# Create vector of descending list of gene changes, named by symbol
   up.list <- input$log2FoldChange
   names(up.list) <- input$symbol
  
  # REACTOME GSEA
  gmt.gsea <- GSEA(geneList  = up.list,
                   nPerm        = 1000,
                   minGSSize    = 10,
                   pvalueCutoff = 1,
                   verbose      = F,
                   TERM2GENE    = hall)

  # Convert to dataframe, add columns with sample labels
  temp <- as.data.frame(gmt.gsea) %>% 
    mutate(pathway = gsub('_',' ',.$ID),
           sample = sample, 
           label = gsub('_.*','',files[i]), 
           time = stri_extract_first_regex(files[i], '[0-9]h')) 
  
  # Remove 'HALLMARK' from pathway column
  temp$pathway <- gsub('HALLMARK ', '', temp$pathway)
    
  # Save results
  out.dir <- here('results/DESeq2_human/GSEA/hallmark/')
  write.csv(temp, 
            file = paste0(out.dir,'GSEA_hallmark_',sample,'.csv'))
  
  # Combine into full dataframe
  full.data <- rbind(full.data, temp)
  
  # Monitor progress
  print(paste(i, 'of', length(files),'done:',files[i] ))
  }
  

# Save the full data
write.csv(full.data, 
          here('results/DESeq2_human/GSEA/hallmark/GSEA_hallmark_all.csv'),
          row.names = F)

# Load the full data
full.data <- read.csv(here('results/DESeq2_human/GSEA/hallmark/GSEA_hallmark_all.csv'))


#' ## Make heatmap of GSEA hallmark sets
#' Now that we have all the GSEA normalized enrichment scores for all samples in a single dataframe, we can make a heatmap to visualize changes in pathways. We'll use `ggplot` to make a heatmap (using `geom_tile()`) of the normalized enrichment scores (NES) and split the 2h and 8h samples apart vertically. The pathways will be in descending order with the most enriched pathways at the top and the least enriched at the bottom.




## Add categories for HALLMARK pathways
# Load hallmark annotation categories
cats <- read_csv(here('data/hallmark_geneset_categories.csv')) %>% 
  dplyr::select(Hallmark_name, Process_category) %>% 
  rename(pathway = Hallmark_name)

cats$Process_category <- gsub('DNA damage', 'Cellular stress', 
                              cats$Process_category) %>% 
  gsub(pattern = 'Cellular component', 'Cell comp.', x = .) %>% 
  gsub(pattern = 'Proliferation', 'Cell division', x = .)

# remove underscores to match other format
cats$pathway <- gsub(pattern = '_', replacement = ' ', 
                     x = cats$pathway)

# Check that all labels match
table(cats$pathway %in% full.data$pathway)


# Add categories to full data
full.data <- left_join(full.data, cats, by = 'pathway')

# Fix Oxygen typo
full.data$pathway <- gsub('OXIGEN', 'OXYGEN', full.data$pathway)

# Remove "PATHWAY" from reactive oxygen species label
full.data$pathway <- gsub('PATHWAY', '', full.data$pathway)




# Set sample order for heatmap
full.data$label <- factor(full.data$label, 
                          levels = c('STM','SE','ST','SPI1','SPI2'))

# Set order for heatmap - pathways by descending average NES
NES.avg <- group_by(full.data, pathway) %>% 
  summarise(NES_avg = mean(NES)) %>% 
  arrange(NES_avg)

full.data$pathway <- factor(full.data$pathway, levels = NES.avg$pathway)

str(full.data)
unique(full.data$Process_category)
# Set order of categories
full.data$Process_category <- factor(full.data$Process_category,
                                     levels = c('Immune response', 'Signaling', 
                                                'Metabolism', 'Cellular stress',
                                                'Cell division', 'DNA damage',
                                                'Development', 'Cell comp.'))


# Round NES (remove decimals for ggplotly tooltip)
full.data$NES <- round(full.data$NES, 2)

# Set up palette of colors for heatmap
hm.palette <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')), space='Lab')

# Get gene labels for interactive plot tooltip
full.data$core_enrichment <- gsub('\\/', ', ', full.data$core_enrichment)
full.data$top10genes <- sub(sprintf('^((?:[^,]*,){10}).*'), '\\1',
                            full.data$core_enrichment)


#+ figure, fig.height = 8.5, fig.width = 7
# Heatmap
hallmark_plot <- function(input){
  ggplot(input, aes(label, pathway, text = top10genes)) + 
  geom_tile(aes(fill = NES), colour = "white") + 
  scale_fill_gradientn(colors = hm.palette(100))+ 
  theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45),
        strip.text.x = element_text(size = 10, face = 'bold'),
        strip.text.y = element_text(size = 10, face = 'bold'),
        axis.title = element_text(size = 12, face = "bold"), 
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10, face = 'bold'))+
  ggtitle('GSEA hallmark gene set enrichment')+
  labs(x='', y = '')+
    facet_grid(rows = vars(Process_category), cols = vars(time),
               scale = 'free', space = 'free', labeller = label_wrap_gen(width=10))
}

p <- hallmark_plot(full.data)
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
#' The paper will be split into sections and so we'll split the heatmap into two groups of samples:
#' 
#' * STM mutants (STM, SPI-1, SPI-2)
#' * Serovars (STM, SE, ST)

#+
#' ## STM mutants GSEA hallmark

# Select samples for heatmap
muts <- dplyr::filter(full.data, label %in% c('STM', 'SPI1', 'SPI2'))

#+ fig.height = 8.5, fig.width = 7
# Heatmap
p1 <- hallmark_plot(muts)
p1


#' ## STM mutants interactive heatmap
#+ fig.height = 10, fig.width = 7
ggplotly(p1)


#' ## Serovars GSEA hallmark

#+ fig.height = 8.5, fig.width = 7
# Select samples
ser <- dplyr::filter(full.data, label %in% c('STM','SE','ST'))

p2 <- hallmark_plot(ser)
p2


#' ## Serovars interactive heatmap
#+ fig.height = 10, fig.width = 7
ggplotly(p2)


#' ## Save png of plots
#+ eval=F
png(filename = here("/img/stm_mutants/mut_GSEA_hallmark_heatmap.png"),
    width =6.5, height = 8.5, units = 'in', res = 300)
p1
dev.off()

png(filename = here("/img/serovars/ser_GSEA_hallmark_heatmap.png"),
    width =6.5, height = 8.5, units = 'in', res = 300)
p2
dev.off()

# Save ggplot objects for the plots
saveRDS(p1, file = here('img/ggplot_objects/gg_mut_hallmarkGSEA.rds'))
saveRDS(p2, file = here('img/ggplot_objects/gg_ser_hallmarkGSEA.rds'))

#+ render, include=F
# Render source file to html 
# dir.create(here('results/DESeq2_human/src_html_output'))

# render.dir <- here('results/DESeq2_human/src_html_output/')

# render(here('src/Hs_align_src/07_gsea_hallmark.R'), output_dir = render.dir, intermediates_dir = render.dir, clean = TRUE)

# Copy to dropbox
# source(here('src/Hs_align_src/XX_copy_to_dropbox.R'))
