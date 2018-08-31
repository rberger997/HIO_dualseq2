


## -----------------------------------------------------------------------------

library(gridExtra)


# Load ggplot objects
ser_hallmark <- readRDS(file = here('img/ggplot_objects/gg_ser_hallmarkGSEA.rds'))
ser_volcano <- readRDS(file = here('img/ggplot_objects/gg_ser_volcano.rds'))
ser_pca <- readRDS(file = here('img/ggplot_objects/gg_ser_pcaplot.rds'))
ser_venn_2i <- readRDS(file = here('img/ggplot_objects/gg_ser_venn_2h_incr.rds'))
ser_venn_8i <- readRDS(file = here('img/ggplot_objects/gg_ser_venn_8h_incr.rds'))
ser_venn_2d <- readRDS(file = here('img/ggplot_objects/gg_ser_venn_2h_decr.rds'))
ser_venn_8d <- readRDS(file = here('img/ggplot_objects/gg_ser_venn_8h_decr.rds'))



png(filename = here('figures/serovars_fig.png'),
    width = 12, height = 10, units = 'in', res = 300)

layout <- rbind(c(1, 1, 2, 2),
                c(1, 1, 2, 2),
                c(1, 1, 2, 2),
                c(3, 4, 2, 2),
                c(5, 6, 2, 2))

grid.arrange(ser_pca,
             ser_hallmark+ggtitle('')+ylab(''),
             #ser_volcano+theme(legend.position="none"),
             ser_venn_2i,
             ser_venn_8i,
             ser_venn_2d,
             ser_venn_8d,
             layout_matrix = layout)

dev.off() 


library(here)
library(ggplot2)
# Load ggplot objects
mut_hallmark <- readRDS(file = here('img/ggplot_objects/gg_mut_hallmarkGSEA.rds'))
mut_volcano <- readRDS(file = here('img/ggplot_objects/gg_mut_volcano.rds'))
mut_pca <- readRDS(file = here('img/ggplot_objects/gg_mut_pcaplot.rds'))
mut_venn_2i <- readRDS(file = here('img/ggplot_objects/gg_mut_venn_2h_incr.rds'))
mut_venn_8i <- readRDS(file = here('img/ggplot_objects/gg_mut_venn_8h_incr.rds'))
mut_venn_2d <- readRDS(file = here('img/ggplot_objects/gg_mut_venn_2h_decr.rds'))
mut_venn_8d <- readRDS(file = here('img/ggplot_objects/gg_mut_venn_8h_decr.rds'))

png(filename = here('figures/stm_muts_fig.png'),
    width = 12, height = 10, units = 'in', res = 300)

layout <- rbind(c(1, 1, 2, 2),
                c(1, 1, 2, 2),
                c(1, 1, 2, 2),
                c(3, 4, 2, 2),
                c(5, 6, 2, 2))

grid.arrange(mut_pca,
             mut_hallmark+ggtitle('')+ylab(''),
             #mut_volcano+theme(legend.position="none"),
             mut_venn_2i,
             mut_venn_8i,
             mut_venn_2d,
             mut_venn_8d,
             layout_matrix = layout)

dev.off() 


