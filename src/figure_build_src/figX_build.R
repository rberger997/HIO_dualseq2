


## -----------------------------------------------------------------------------

library(gridExtra)


# Load ggplot objects
ser_hallmark <- readRDS(file = here('img/ggplot_objects/gg_ser_hallmarkGSEA.rds'))
ser_volcano <- readRDS(file = here('img/ggplot_objects/gg_ser_volcano.rds'))
ser_pca <- readRDS(file = here('img/ggplot_objects/gg_ser_pcaplot.rds'))


png(filename = here('figures/testfig.png'),
    width = 12, height = 10, units = 'in', res = 300)

layout <- rbind(c(1, 1, 2, 2),
                c(3, 3, 2, 2))

grid.arrange(ser_pca,
             ser_hallmark+ggtitle('')+ylab(''),
             ser_volcano+theme(legend.position="none"),
             layout_matrix = layout)

dev.off()



