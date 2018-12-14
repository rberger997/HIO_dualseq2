# Build gene level analysis figure from ggplot objects using gridextra

library(gridExtra)
library(ggplot2)
library(here)

## -------------------------------------------------------------------------
# STM mutants figure

# Load ggplot objects
mut_volcano <- readRDS(file = here('img/ggplot_objects/gg_mut_over_stm_volcano.rds'))
spi1_2h_bar <- readRDS(file = here('img/ggplot_objects/gg_SPI1_2h_gene_barplot.rds'))
spi1_8h_bar <- readRDS(file = here('img/ggplot_objects/gg_SPI1_8h_gene_barplot.rds'))
spi2_2h_bar <- readRDS(file = here('img/ggplot_objects/gg_SPI2_2h_gene_barplot.rds'))
spi2_8h_bar <- readRDS(file = here('img/ggplot_objects/gg_SPI2_8h_gene_barplot.rds'))


mut_volcano


png(filename = here('figures/gene_level_muts_fig.png'),
    width = 8.5, height = 11, units = 'in', res = 300)

layout <- rbind(c(1, 1, 1, 1),
                c(1, 1, 1, 1),
                c(2, 2, 4, 4),
                c(3, 3, 5, 5))

grid.arrange(mut_volcano,
             spi1_2h_bar,
             spi1_8h_bar,
             spi2_2h_bar,
             spi2_8h_bar,
             layout_matrix = layout)

dev.off() 

## -------------------------------------------------------------------------
# Serovars figure

# Load ggplot objects
ser_volcano <- readRDS(file = here('img/ggplot_objects/gg_ser_over_stm_volcano.rds'))
se_2h_bar <- readRDS(file = here('img/ggplot_objects/gg_se_2h_gene_barplot.rds'))
se_8h_bar <- readRDS(file = here('img/ggplot_objects/gg_se_8h_gene_barplot.rds'))
st_2h_bar <- readRDS(file = here('img/ggplot_objects/gg_st_2h_gene_barplot.rds'))
st_8h_bar <- readRDS(file = here('img/ggplot_objects/gg_st_8h_gene_barplot.rds'))





png(filename = here('figures/gene_level_ser_fig.png'),
    width = 8.5, height = 11, units = 'in', res = 300)

layout <- rbind(c(1, 1, 1, 1),
                c(1, 1, 1, 1),
                c(2, 2, 4, 4),
                c(3, 3, 5, 5))

grid.arrange(ser_volcano,
             se_2h_bar,
             se_8h_bar,
             st_2h_bar,
             st_8h_bar,
             layout_matrix = layout)

dev.off()





