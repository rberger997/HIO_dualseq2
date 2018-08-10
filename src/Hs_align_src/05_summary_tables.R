#'---
#' title: "RNA seq workflow part 5 - Summary tables"
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
#' This script is designed to summarize gene changes in a series of RNA seq differential expression analyses. The output is a table with the number of genes in the categories:
#' 
#'  * Non significant -- (adj p > 0.05)
#'  * Decreasing -- (adjp < 0.05, log2FoldChange < 0)
#'  * Increasing -- (adjp < 0.05, log2FoldChange > 0)

#+ start
#' # Begin script
 
#' ## Libraries and directories
#+ libs, message=F, warning=F
library(dplyr)
library(gtable)
library(grid)
library(gridExtra)
library(here)
library(knitr)
library(magrittr)
library(rmarkdown)
library(tidyr)



#' ## Make table

#' Load in the data from volcano plots
data2 <- read.csv(here('results/DESeq2_human/volcano_data.csv'))

head(data2)

# Set order of samples for tables
data2$label <- factor(data2$label, levels = c('SE','ST','STM','SPI1','SPI2'))
data2$colors <- factor(data2$colors, levels = c('Non significant','Decreasing','Increasing'))


# Create table for 2h gene changes
table2h <- dplyr::filter(data2, time == '2h') %>% 
  group_by(label) %>% 
  #  filter(abs(log2FoldChange) > 1) %>% 
  count(label, A = colors) %>% 
  rename(Injection = label) %>% 
  spread(key = A, value = n)

knitr::kable(table2h)

# Create table for 2h gene changes
table8h <- dplyr::filter(data2, time == '8h') %>% 
  group_by(label) %>% 
#  filter(abs(log2FoldChange) > 1) %>% 
  count(label, A = colors) %>% 
  rename(Injection = label) %>% 
  spread(key = A, value = n)

knitr::kable(table8h)



#' ## Format tables to save as images
#' The tables will be formatted and converted to .png files using the `grid` and `gridextra` packages.

# Add formatting to tables
g2 <- tableGrob(table2h, rows = NULL) %>% 
  gtable_add_grob(grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 2, b = nrow(.), l = 1, r = ncol(.)) %>% 
  gtable_add_grob(grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 1, l = 1, r = ncol(.))
grid.draw(g2)


# 8h table
g8 <- tableGrob(table8h, rows = NULL) %>% 
  gtable_add_grob(grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                  t = 2, b = nrow(.), l = 1, r = ncol(.)) %>% 
  gtable_add_grob(grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                  t = 1, l = 1, r = ncol(.))
grid.draw(g8)


# Save png of tables
png(filename = here('img/summary_table_2h.png'),
    height = 2, width = 5, units = 'in', res = 300)
grid.arrange(g2, top = 'Gene summary - 2h post injection')
dev.off()


png(filename = here('img/summary_table_8h.png'),
    height = 2, width = 5, units = 'in', res = 300)
grid.arrange(g8, top = 'Gene summary - 8h post injection')
dev.off()


#+ render, include=F
# Render source file to html 
# dir.create(here('results/DESeq2_human/src_html_output'))

# render.dir <- here('results/DESeq2_human/src_html_output/')

# render(here('src/Hs_align_src/05_summary_tables.R'), output_dir = render.dir, intermediates_dir = render.dir, clean = TRUE)