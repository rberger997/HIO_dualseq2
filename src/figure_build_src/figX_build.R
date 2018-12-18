# Build figure from ggplot objects using gridextra

library(gridExtra)
library(ggplot2)
library(here)

## -------------------------------------------------------------------------
# STM mutants figure

# Create folder for figures
dir.create(here('figures/'))


# Load ggplot objects
mut_hallmark <- readRDS(file = here('img/ggplot_objects/gg_mut_hallmarkGSEA.rds'))
mut_pca <- readRDS(file = here('img/ggplot_objects/gg_mut_pcaplot.rds'))
mut_venn_2h <- readRDS(file = here('img/ggplot_objects/gg_mut_venn_2h.rds'))
mut_venn_8h <- readRDS(file = here('img/ggplot_objects/gg_mut_venn_8h.rds'))


png(filename = here('figures/stm_muts_fig.png'),
    width = 12, height = 8, units = 'in', res = 300)

layout <- rbind(c(1, 1, 2, 2),
                c(1, 1, 2, 2),
                c(3, 4, 2, 2))

grid.arrange(mut_pca,
             mut_hallmark+ggtitle('')+ylab(''),
             mut_venn_2h,
             mut_venn_8h,
             layout_matrix = layout)

dev.off() 

## -------------------------------------------------------------------------
# Serovars figure

# Load ggplot objects
ser_hallmark <- readRDS(file = here('img/ggplot_objects/gg_ser_hallmarkGSEA.rds'))
ser_pca <- readRDS(file = here('img/ggplot_objects/gg_ser_pcaplot.rds'))
ser_venn_2h <- readRDS(file = here('img/ggplot_objects/gg_ser_venn_2h.rds'))
ser_venn_8h <- readRDS(file = here('img/ggplot_objects/gg_ser_venn_8h.rds'))


png(filename = here('figures/serovars_fig.png'),
    width = 12, height = 8, units = 'in', res = 300)

layout <- rbind(c(1, 1, 2, 2),
                c(1, 1, 2, 2),
                c(3, 4, 2, 2))


grid.arrange(ser_pca,
             ser_hallmark+ggtitle('')+ylab(''),
             #ser_volcano+theme(legend.position="none"),
             ser_venn_2h,
             ser_venn_8h,
             layout_matrix = layout)

dev.off() 





