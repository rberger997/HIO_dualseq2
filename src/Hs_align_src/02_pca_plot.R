#'---
#' title: "RNA seq workflow part 2 - Dimensionality redution and PCA plot"
#' subtitle: "HIOs dual seq experiment 2"
#' author: "Ryan Berger"
#' date: "2018-08-07"
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
#' This script is part of a workflow for RNA seq processing and analysis (modified from [Bioconductor](https://www.bioconductor.org/help/workflows/rnaseqGene/)). The purpose is to take in a DESeq2 object (dds) and return a PCA plot of all the samples.

#------------------------------------------------------------
#' # Begin script

#' ## Libraries and directories
# #+ load_pkgs, error=T
library(DESeq2)
library(ggplot2)
library(here)
library(knitr)
library(magrittr)
library(rmarkdown)
source(here("src/Hs_align_src/ggplot2-themes.R"))


#+ set_cache_dir, include=F
# Set a directory for all cache files to go
dir.create(here('results/DESeq2_human/src_html_output/knitr_cache'))
opts_chunk$set(cache.path = here('results/DESeq2_human/src_html_output/knitr_cache/'))


#' The output files will go to the `DESeq2_human` folder in the results directory.

results.dir <- here('results/DESeq2_human/')


#' ## Load the data
#' We previously saved the `dds` object in the results folder. Let's load it now:

dds <- readRDS(here('results/DESeq2_human/dds_all.rds'))


#' ## Transform data using rlog
#' Before running PCA, we'll do an rlog transformation of the data. This is a variance stabilizing transformation that "transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size" ([source](https://rdrr.io/bioc/DESeq2/man/rlog.html)). 

#+ rlog, cache=T
rld <- rlog(dds, blind = FALSE)

# Compare before and after transformation
head(assay(dds)[,1:4], 3) ; head(assay(rld)[,1:4], 3)

# Save the rlog transformed data
saveRDS(rld, file = file.path(results.dir, 'rld_all.rds'))

# Load the rlog transformed data
rld <- readRDS(here('results/DESeq2_human/rld_all.rds'))

#' ## PCA plot
#' Now we can run PCA using the `plotPCA` function from DESeq2 and use those values to make a custom plot with `ggplot`. First we'll set the order of injection to how we want it to show up in the plot.

# Set order of injection
unique(colData(rld)$Inject)

colData(rld)$Inject <- colData(rld)$Inject %>% 
  factor(., levels(.)[c(1,6,2,5,3,4)])


#+ pcaplot
# PCA plot (use just only for PC variance estimates)
pca <- plotPCA(rld, intgroup = c('code_name','Inject','hr'))

# Get PCA data
pca.df <- plotPCA(rld, intgroup = c('code_name', 'Inject', 'hr'), returnData = TRUE)


# Make plot
pca_plot <- function(input){
ggplot(data = input, aes(x = PC1, y = PC2))+ 
  geom_point(shape = 21, stroke = 1.5, 
             aes(fill = as.factor(Inject),
                 color = as.factor(hr)), 
             size = 6) +
  theme1 + 
  #scale_y_continuous(limits = c(-25,25), breaks = seq(-25, 25, by = 10)) +
  scale_fill_brewer(palette = "Set1", name = 'Injection') +
  #scale_color_brewer(palette = "Set2", name = 'Time', direction = 1) +    
  scale_color_manual(values=c("gray", "black"), name = 'Time (hr)')+
  theme(legend.position = "right") +
  geom_hline(yintercept = 0,
             size = 1, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0,
             size = 1, linetype = "dashed", color = "grey70") +
  coord_fixed(ratio = 1) +
  xlab(pca$labels$x) + #pull variance estimates from al. plotPCA call
  ylab(pca$labels$y) +
  # Move y axis
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 1, 
                                                    b = 0, l = 0))) +
  # Move x axis
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, 
                                                    b = 0, l = 0)))+
 # Shrink axis labels down
  theme(plot.caption = element_text(vjust = 1), 
    axis.title = element_text(size = 24), 
    axis.text.x = element_text(size = 18), 
    axis.text.y = element_text(size = 18), 
    plot.title = element_text(size = 30))
}

plot <- pca_plot(pca.df)

#' ## Observations
#' Looking at the PCA plot, it appears the main sources of variance are the time point of collection and the presence/absence of bacteria injected. This seems to separate samples into an interesting pattern of three clusters:
#' 
#' * The bottom left of the plot shows the PBS controls at both 2 and 8 hours clustering together and showing separation from all the HIOs injected with bacteria (except for one control that is in the upper left quadrant).
#' * The upper middle of the plot shows all of the 8 hour bacterially injected HIOs clustering together.
#'     + For the most part, it looks like 3 of 4 samples for each subtype of bacteria are grouped together while one sample shows more variance.
#' * The bottom right of the plot is where all the 2 hour HIOs injected with bacteria are clustering together.
#'
#' PC1 contains 31% of the total variance and splits samples by PBS control, 8 hour infected samples, and 2 hour infected samples. Interestingly, there appears to be greater variance from PBS control after 2 hours of exposure to bacteria than at 8 hours.
#' 
#' PC2 seems to mostly split the 8 hour infected samples from the PBS control and 2 hour infected samples.
#' 

#' ## Save plot
#' We'll save the PCA plot as a .png file in the `img` directory.

# save png of plot
#+ save, eval=F
png(filename = here("/img/pca.png"),
    width = 900, height = 500)
print(plot)
dev.off()




#' # Split samples
#' The paper will be split into sections and so we'll split the heatmap into two groups of samples:
#' 
#' * STM mutants (STM, SPI-1, SPI-2)
#' * Serovars (STM, SE, ST)
#' 



#' We previously saved the `dds` object in the results folder. Let's load it now:
dds_mut <- readRDS(here('results/DESeq2_human/dds_stm_muts.rds'))

#' ## Transform data using rlog
#+ rlog, cache=T
rld_mut <- rlog(dds_mut, blind = FALSE)

# Save the rlog transformed data
saveRDS(rld_mut, file = file.path(results.dir, 'rld_stm_muts.rds'))





#' We previously saved the `dds` object in the results folder. Let's load it now:
dds_ser <- readRDS(here('results/DESeq2_human/dds_serovars.rds'))

#' ## Transform data using rlog
#+ rlog, cache=T
rld_ser <- rlog(dds_ser, blind = FALSE)


# Compare before and after transformation
head(assay(dds1)[,1:4], 3) ; head(assay(rld1)[,1:4], 3)

# Save the rlog transformed data
saveRDS(rld_ser, file = file.path(results.dir, 'rld_serovars.rds'))




# Load the rlog transformed data
rld_mut <- readRDS(here('results/DESeq2_human/rld_stm_muts.rds'))
# Load the rlog transformed data
rld_ser <- readRDS(here('results/DESeq2_human/rld_serovars.rds'))




#' ## PCA plot
#' Now we can run PCA using the `plotPCA` function from DESeq2 and use those values to make a custom plot with `ggplot`. First we'll set the order of injection to how we want it to show up in the plot.

# Set order of injection
unique(colData(rld_mut)$Inject)

colData(rld_mut)$Inject <- colData(rld_mut)$Inject %>% 
  factor(., levels(.)[c(1,6,3,4)])



# Set order of injection
unique(colData(rld_ser)$Inject)

colData(rld_ser)$Inject <- colData(rld_ser)$Inject %>% 
  factor(., levels(.)[c(1,6,2,5)])




#+ 

pcaplot2 <- function(rld){
# PCA plot (use just only for PC variance estimates)
pca <- plotPCA(rld, intgroup = c('code_name','Inject','hr'))

# Get PCA data
pca.df <- plotPCA(rld, intgroup = c('code_name', 'Inject', 'hr'), returnData = TRUE)


# Make plot
pca_plot <- function(input){
  ggplot(data = input, aes(x = PC1, y = PC2))+ 
    geom_point(shape = 21, stroke = 1.5, 
               aes(fill = as.factor(Inject),
                   color = as.factor(hr)), 
               size = 6) +
    theme1 + 
    #scale_y_continuous(limits = c(-25,25), breaks = seq(-25, 25, by = 10)) +
    scale_fill_brewer(palette = "Set1", name = 'Injection') +
    #scale_color_brewer(palette = "Set2", name = 'Time', direction = 1) +    
    scale_color_manual(values=c("gray", "black"), name = 'Time (hr)')+
    theme(legend.position = "right") +
    geom_hline(yintercept = 0,
               size = 1, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0,
               size = 1, linetype = "dashed", color = "grey70") +
    coord_fixed(ratio = 1) +
    xlab(pca$labels$x) + #pull variance estimates from al. plotPCA call
    ylab(pca$labels$y) +
    # Move y axis
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 1, 
                                                      b = 0, l = 0))) +
    # Move x axis
    theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, 
                                                      b = 0, l = 0)))+
    # Shrink axis labels down
    theme(plot.caption = element_text(vjust = 1), 
          axis.title = element_text(size = 24), 
          axis.text.x = element_text(size = 18), 
          axis.text.y = element_text(size = 18), 
          plot.title = element_text(size = 30))
  

}
plot <- pca_plot(pca.df)
return(plot)
}

mut.plot <- pcaplot2(rld_mut)
ser.plot <- pcaplot2(rld_ser)


# save png of plot
#+ eval=F
png(filename = here("/img/stm_mutants/mut_pca.png"),
    width = 900, height = 500)
print(mut.plot)
dev.off()

png(filename = here("/img/serovars/ser_pca.png"),
    width = 900, height = 500)
print(ser.plot)
dev.off()


saveRDS(mut.plot, file = here('img/ggplot_objects/gg_mut_pcaplot.rds'))
saveRDS(ser.plot, file = here('img/ggplot_objects/gg_ser_pcaplot.rds'))


#+ render, include=F
# Render source file to html 
# dir.create(here('results/DESeq2_human/src_html_output'))

# render.dir <- here('results/DESeq2_human/src_html_output/')

# render(here('src/Hs_align_src/02_pca_plot.R'), output_dir = render.dir, intermediates_dir = render.dir, clean = TRUE)