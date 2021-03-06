
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
make_venn <- function(df,x1,x2,x3,t, padj = 1, l2fc = 0, show_labels = T){

  # Set significance filter
  df$signif <- ifelse(df$padj < padj & abs(df$log2FoldChange) > l2fc, 'Significant', 'Non significant')
  
  # want to compare all significant genes
incr <- filter(df, signif == 'Significant' & time == t) %>% 
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


# Formatting options
eulerr_options(pointsize = 14)
options(digits = 4)


# Input in the form of a named numeric vector
fit1 <- euler(c("A" = A, "B" = B, "C" = C,
                "A&B" = AB, "A&C" = AC, "B&C" = BC,
                "A&B&C" = ABC))

diagram <- plot(fit1, 
             quantities = T,
             fill = c("lightblue", "lightcoral", "lemonchiffon"),
             lty = 1,
             labels = if(show_labels == F){
               F
             }else{
               c(x1,x2,x3)})
return(diagram)
}


# STM mutants venn diagrams 
# (leave labels off for better formatting in powerpoint)
m2 <- make_venn(data2, 'STM','SPI1','SPI2','2h', padj = 0.05, l2fc = 1,
                show_labels = F)

m8 <- make_venn(data2, 'STM','SPI1','SPI2','8h', padj = 0.05, l2fc = 1,
                show_labels = F)

# Serovars venn diagrams
# (leave labels off for better formatting in powerpoint)
s2 <- make_venn(data2, 'STM','SE','ST','2h', padj = 0.05, l2fc = 1,
                show_labels = F)
s8 <- make_venn(data2, 'STM','SE','ST','8h', padj = 0.05, l2fc = 1,
                show_labels = F)


# Save venn diagram objects
saveRDS(m2, file = here(paste0('img/ggplot_objects/gg_mut_venn_2h.rds')))
saveRDS(m8, file = here(paste0('img/ggplot_objects/gg_mut_venn_8h.rds')))

saveRDS(s2, file = here(paste0('img/ggplot_objects/gg_ser_venn_2h.rds')))
saveRDS(s8, file = here(paste0('img/ggplot_objects/gg_ser_venn_8h.rds')))


