
# Making venn diagrams with the eulerr package

# Should make area proportional venn diagrams from lists or a matrix


#' ## Libraries and directories
library(dplyr)
library(eulerr)
library(here)
library(magrittr)
library(tidyr)


# Load data from volcano plots
data2 <- read.csv(here('results/DESeq2_human/volcano_data.csv'), 
                  stringsAsFactors = F)
head(data2)



# Function to make venn diagrams
make_venn <- function(df,x1,x2,x3,color,t){

  # want to compare all genes increasing
incr <- filter(df, colors == color & time == t) %>% 
  select(c(symbol,label)) %>% 
  arrange(label)


# Split into three groups by sample
a <- filter(incr, label == x1) %>% 
  .$symbol
b <- filter(incr, label == x2) %>% 
  .$symbol
c <- filter(incr, label == x3) %>% 
  .$symbol
  
# Calculate number of shared genes in each group
ABC <- length(c[c%in%a[a%in%b]])
AB <- length(a[a%in%b]) - ABC
AC <- length(a[a%in%c]) - ABC
A <- length(a) - AB -AC - ABC
BC <- length(b[b%in%c]) - ABC
BA <- length(b[b%in%a]) - ABC
B <- length(b)-BA-BC-ABC
C <- length(c)-AC-BC-ABC



eulerr_options(pointsize = 14)
options(digits = 4)
# Input in the form of a named numeric vector
fit1 <- euler(c("A" = A, "B" = B, "C" = C,
                "A&B" = AB, "A&C" = AC, "B&C" = BC,
                "A&B&C" = ABC))

test <- plot(fit1, 
             quantities = T,
             fill = c("lightblue", "lightcoral", "lemonchiffon"),
             lty = 1,
             labels = c(x1,x2,x3))
return(test)
}


# STM mutants venn diagrams
m2i <- make_venn(data2, 'STM','SPI1','SPI2','Increasing','2h')
m8i <- make_venn(data2, 'STM','SPI1','SPI2','Increasing','8h')

m2d <- make_venn(data2, 'STM','SPI1','SPI2','Decreasing','2h')
m8d <- make_venn(data2, 'STM','SPI1','SPI2','Decreasing','8h')


# Serovars venn diagrams
s2i <- make_venn(data2, 'STM','SE','ST','Increasing','2h')
s8i <- make_venn(data2, 'STM','SE','ST','Increasing','8h')

s2d <- make_venn(data2, 'STM','SE','ST','Decreasing','2h')
s8d <- make_venn(data2, 'STM','SE','ST','Decreasing','8h')


# Save venn diagram objects
saveRDS(m2i, file = here(paste0('img/ggplot_objects/gg_mut_venn_2h_incr.rds')))
saveRDS(m8i, file = here(paste0('img/ggplot_objects/gg_mut_venn_8h_incr.rds')))
saveRDS(m2d, file = here(paste0('img/ggplot_objects/gg_mut_venn_2h_decr.rds')))
saveRDS(m8d, file = here(paste0('img/ggplot_objects/gg_mut_venn_8h_decr.rds')))


saveRDS(s2i, file = here(paste0('img/ggplot_objects/gg_ser_venn_2h_incr.rds')))
saveRDS(s8i, file = here(paste0('img/ggplot_objects/gg_ser_venn_8h_incr.rds')))
saveRDS(s2d, file = here(paste0('img/ggplot_objects/gg_ser_venn_2h_decr.rds')))
saveRDS(s8d, file = here(paste0('img/ggplot_objects/gg_ser_venn_8h_decr.rds')))





eulerr::eulerr_options()





# Save as png
png(filename = here(paste0('img/venn_serovars_',i,'.png')),
    width =5, height = 5, units = 'in', res = 300)
  print(make_venn(data2, 'SE','STM','ST','Increasing',i))
dev.off()
