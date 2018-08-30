# Make hallmark heatmap for Anna-Lisa
#
# Cut hallmark heatmap to only include:
# ST,STM
# Top ~20 genes by NES


library(ggplot2)
library(dplyr)
library(here)
library(RColorBrewer)


# Load data
 full.data <- read.csv(here('results/DESeq2_human/GSEA/hallmark/GSEA_hallmark_all.csv'))


#' ## Make heatmap of GSEA hallmark sets
#' Now that we have all the GSEA normalized enrichment scores for all samples in a single dataframe, we can make a heatmap to visualize changes in pathways. We'll use `ggplot` to make a heatmap (using `geom_tile()`) of the normalized enrichment scores (NES) and split the 2h and 8h samples apart vertically. The pathways will be in descending order with the most enriched pathways at the top and the least enriched at the bottom.

 
 # Filter out all samples except STM, ST
 full.data <- filter(full.data, label %in% c('ST','STM'))

# Set sample order for heatmap
full.data$label <- factor(full.data$label, 
                          levels = c('ST','STM'))

# Set order for heatmap - pathways by descending average NES
NES.avg <- group_by(full.data, pathway) %>% 
  summarise(NES_avg = mean(NES)) %>% 
  arrange(NES_avg)

full.data$pathway <- factor(full.data$pathway, levels = NES.avg$pathway)


# Round NES (remove decimals for ggplotly tooltip)
full.data$NES <- round(full.data$NES, 2)

# Set up palette of colors for heatmap
hm.palette <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')), space='Lab')




#+ figure, fig.height = 8.5, fig.width = 7
# Heatmap
hallmark_plot <- function(input){
  ggplot(input, aes(label, pathway)) + 
    geom_tile(aes(fill = NES), colour = "white") + 
    scale_fill_gradientn(colors = hm.palette(100))+ 
    theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45),
          strip.text.x = element_text(size = 12, face = 'bold'),
          strip.text.y = element_text(size = 12, face = 'bold'),
          axis.title = element_text(size = 14, face = "bold"), 
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12, face = 'bold'))+
    ggtitle('GSEA hallmark gene set enrichment')+
    labs(x='', y = 'Pathway')+
    facet_grid(cols = vars(time), scales = 'free')
}

p <- hallmark_plot(full.data)
p


#' ## Save png of plot
#+ save, eval=F
png(filename = here("/img/GSEA_hallmark_heatmap_AL.png"),
    width =7, height = 8.5, units = 'in', res = 500)
p
dev.off()
