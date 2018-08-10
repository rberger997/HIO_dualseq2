
#'---
#' title: "RNA seq interactive volcano plots"
#' subtitle: "HIOs dual seq experiment 2"
#' author: "Ryan Berger"
#' date: "2018-08-10"
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
#' The RNA seq volcano plots that were made previously were static images but this script will convert them into interactive objects that can be used to explore the data in an html file.
#' 

#' # Interactive volcano plot
#' Hover over points to see the gene name and values. The non-significant points have been removed to increase performance and reduce lag.
#' 
#+ prep, echo=F
library(ggplot2)
library(rmarkdown)
library(knitr)
library(plotly)


# Load data
data2 <- read.csv(here('results/DESeq2_human/volcano_data.csv'))

# Filter out the non-significant points to speed up perfromance
data2 <- filter(data2, colors != 'Non significant')

#+ fig2, fig.height = 7, fig.width = 10, fig.align = 'center', echo=F
# Make plot - facet time into columns, samples into rows
p <- ggplot(data = data2, aes(x=log2FoldChange, y=-log10(padj), 
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
  facet_grid(rows = vars(time), cols = vars(label))+
  xlim(-12,12)

ggplotly(p)



#+ render, include=F
# Render source file to html 
# dir.create(here('results/DESeq2_human/src_html_output'))

# render.dir <- here('results/DESeq2_human/src_html_output/')

# render(here('src/Hs_align_src/interactive_volcano_plot.R'), output_dir = render.dir, intermediates_dir = render.dir, clean = TRUE)