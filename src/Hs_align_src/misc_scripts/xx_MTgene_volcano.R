#' Make a volcano plot of a list of mitochondrial genes from Basel. 

#' Run the 04_volcano_plot.R script first to get the data2 object and libraries loaded.


data2 <- read.csv(here('results/DESeq2_human/volcano_data.csv'))

library(dplyr)
library(ggplot2)
head(data2)

mtgenes <- c('MT-ATP6','MT-ATP8','MT-CO1','MT-CO2','MT-CO3','MT-CYB','MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6')

head(data2)

data3 <- filter(data2, symbol %in% mtgenes)


p <- ggplot(data = data3, aes(x=log2FoldChange, y=-log10(padj), 
                              color=colors, label=symbol))+
  geom_point(size=.85)+
  geom_vline(xintercept = 0, linetype = 'dotted')+
  geom_hline(yintercept = 0, linetype = 'dotted')+
  scale_color_manual(name = '',
                     values = c('Non significant' = 'gray',
                                'Decreasing' = 'blue',
                                'Increasing' = 'red'),
                     labels = c('Non significant',
                                'Significant decrease',
                                'Significant increase'))+
  theme_bw()+ 
  theme(panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 14, face = 'bold'),
        strip.text.y = element_text(size = 14, face = 'bold'))+
  facet_grid(rows = vars(time), cols = vars(label))

p



#' ## Save png of plot
#+ save, eval=F
png(filename = here("/img/MTgenes_volcano.png"),
    width = 900, height = 600)
print(p)
dev.off()








library(RColorBrewer)
# Set up palette of colors for heatmap
hm.palette <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')), space='Lab')

data4 <- filter(data3, label %in% c('SE','STM','ST'))
data4$label <- factor(data4$label, levels = c('STM','SE','ST'))


# Heatmap of MT genes 
p1 <- ggplot(data = data4, aes(x = label, y = symbol))+
  geom_tile(aes(fill = log2FoldChange), colour = "white")+
  scale_fill_gradientn(colors = hm.palette(11), limits = c(-1.65,1.65))+
  facet_grid(cols = vars(time))+
  ylab('')+
  xlab('') + 
  theme(plot.subtitle = element_text(vjust = 1), 
    plot.caption = element_text(vjust = 1), 
    axis.text = element_text(size = 12, face = "bold"), 
    strip.text.x = element_text(size = 12, face = 'bold'),
    legend.text = element_text(size = 12, 
        face = "bold"), legend.title = element_text(size = 12, 
        face = "bold"))

#' ## Save png of plot
#+ save, eval=F
png(filename = here("/img/MTgenes_heatmap.png"),
    width = 6, height = 5, units = 'in', res = 300)
print(p1)
dev.off()

