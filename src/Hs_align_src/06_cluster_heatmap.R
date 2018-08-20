#'---
#' title: "RNA seq workflow part 6 - Cluster heatmap"
#' subtitle: "HIOs dual seq experiment 2"
#' author: "Ryan Berger"
#' date: "2018-08-13"
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
#' This script is designed to generate a summary heatmap of the top variance genes in an RNA seq experiment. The approach is to use an unbiased approach to see the genes with the most changes across all samples in the experiment. The input is DESeq counts (rlog transformed) and output is a heatmap organized by hierarchical clustering.

#+ 
#' # Begin script
#' 
#' ## Libraries and directories
library(d3heatmap)
library(dplyr)
library(genefilter)
library(here)
library(knitr)
library(magrittr)
library(pheatmap)
library(RColorBrewer)
library(rmarkdown)


# Create directories for images of split samples to go into
dir.create(here('img/stm_mutants'), recursive = T)
dir.create(here('img/serovars'), recursive = T)


#------------------------------------------------------------

#' ## Load and prep data
#' For this heatmap we want to see the top n variance genes across all samples. We also want each condition to be a single box and so we'll calculate the averages of the four replicates. There is a choice of which to do first: define the top variance genes or calculate the averages by sample. We'll calculate averages first and then use that result to pick the top variance genes to minimize single replicates that are outliers.
#' 
#' We'll iterate over all the columns in the counts matrix to calculate averages. Since they're all in order (every four columns is a single sample), we'll set a column index, use the `apply` and `mean` functions across each group of four columns, and `cbind` everything together into a new matrix.

# Load data
rld <- readRDS(here('results/DESeq2_human/rld_all.rds'))


# Counts data
counts1 <- assay(rld)

# Empty object for new matrix to go into
mat1 <- NULL


# Caculate average of sequential 4 columns
for(i in 1:12){
  # 48 samples, 4 replicates of each
  index <- seq(1, 48, by=4)
  lo <- index[i]
  hi <- index[i]+3
  
  # Get sample names for column
  label <- gsub(pattern = '-([0-9])$','',colnames(counts1)[hi])
  
  # Calculate average by sample
  temp <-  apply(counts1[,lo:hi], 1, mean)
  
  # Combine together
  mat1 <- cbind(mat1, temp)
  colnames(mat1)[i] <- label
}
head(mat1)



#' ## Heatmap of top variance genes
#' We'll select the top n variance genes and make into a heatmap using `pheatmap`. We can also make an interactive version using `d3heatmap`.



# Make function for making the heatmap
make_heatmap <- function(input, n, cutrows){

  # Define top variance genes, extract from full data
topVarGenes <- head(order(rowVars(input), decreasing = TRUE), n) 
mat <- input[topVarGenes, ]

# convert to fold over mean of all samples
mat <- mat - rowMeans(mat)  

# Get sample times to use as annotation
anno1 <- as.data.frame(gsub('.*_','',colnames(input))) %>% 
  set_colnames('Time') %>% 
  set_rownames(colnames(input))

labels <- gsub('_.*$','',rownames(anno1))

# Make heatmap
p <- pheatmap(mat,
         main = paste('Top',n,'variance genes'),
         cutree_rows = cutrows,
         cutree_cols = 3,
         treeheight_row = 25,
         treeheight_col = 25,
         annotation_col = anno1,
         annotation_names_col = F,
         labels_col = labels,
         annotation_colors = list(Time = c(`2h` = 'gray80',`8h` =  'gray30')),
         fontsize = 12,
         cellwidth = 15,
         cellheight = 15)

p
}

#+ figure, fig.height = 8, fig.width = 6
p <- make_heatmap(input = mat1, n = 25, cutrows = 3)
p


#' ## Heatmap observations
#' * The main sample clustering (vertical) splits PBS controls away from all HIOs injected with bacteria.
#' * In samples injected with bacteria, main clustering occurs by time point.
#'     + The expression profiles for these genes are more similar by time point than by the type of bacteria injected.
#'     + Matches what was seen in PCA plot where time point was the major separator of samples.
#' * The gene clustering (horizontal) shows an interesting pattern:
#'     + **Early response genes** (bottom subset): highly up-regulated at 2h, decreasing at 8h.
#'     + **Late response genes** (top subset): slightly up-regulated at 2h, higher expression at 8h. Top subset.
#'     + **Up-regulated at both times** (middle subset): up-regulated at 2h and maintaining similar expression at 8h.


#+ save
#' ## Save as png
# Save png of heatmap
png(here('img/topvar_heatmap.png'), 
    width = 5, height = 7, units = 'in', res = 300)
p
dev.off()


#' ## Interactive heatmap
#' We can use the data to make an interactive version of the heatmap using the `d3heatmap` package.

 # interactive, fig.height = 8, fig.width = 6
d3heatmap(mat,
          colors = colorRampPalette(rev(brewer.pal(n = 7, 
                                                   name = "RdYlBu")))(100),
          height = 900,
          width = 550)




#' # Split samples
#' The paper will be split into sections and so we'll split the heatmap into two groups of samples:
#' 
#' * STM mutants (STM, SPI-1, SPI-2)
#' * Serovars (STM, SE, ST)
#' 
#' Each heatmap will contain the top 25 variance genes for that set of samples.

#+ 
#' ## STM mutants heatmap

colnames(mat1)


#+ fig.height = 8, fig.width = 6
# Filter samples
mat2 <- mat1[,c(1:4,9:12)]

p1 <- make_heatmap(input = mat2, 
                   n=25, 
                   cutrows = 2)

p1

# Save png of heatmap
png(here('img/stm_mutants/mut_topvar_heatmap.png'), 
    width = 5, height = 7, units = 'in', res = 300)
p1
dev.off()



#' ## Serovars heatmap

colnames(mat1)

#+ fig.height = 8, fig.width = 6
# Filter samples
mat3 <- mat1[,c(1:8)]

p2 <- make_heatmap(input = mat3, 
                   n=25, 
                   cutrows = 2)

p2

# Save png of heatmap
png(here('img/serovars/ser_topvar_heatmap.png'), 
    width = 5, height = 7, units = 'in', res = 300)
p2
dev.off()




#+ render, include=F
# Render source file to html 
# dir.create(here('results/DESeq2_human/src_html_output'))

# render.dir <- here('results/DESeq2_human/src_html_output/')

# render(here('src/Hs_align_src/06_cluster_heatmap.R'), output_dir = render.dir, intermediates_dir = render.dir, clean = TRUE)


# Copy to dropbox
# source(here('src/Hs_align_src/XX_copy_to_dropbox.R'))



# #-----------------------------------------------------------
#+ Extra_stuff, include=F

#
#  Function for generating average by sample in counts file
#  Works as long as number of replicates are the same and columns in order by sample
#
# avg_by_sample <- function(data, n, reps){
#   samples <- n/reps
#
#   # Empty object for new matrix to go into
#   mat1 <- NULL
#
# for(i in seq(samples)){
#   # 48 samples, 4 replicates of each
#   index <- seq(1, n, by=reps)
#   lo <- index[i]
#   hi <- index[i]+(reps -1)
#
#   # Get sample names for column
#   label <- gsub(pattern = '-([0-9])$','',colnames(data)[hi])
#
#   # Calculate average by sample
#   A <-  apply(data[,lo:hi], 1, mean)
#
#   # Combine together
#   mat1 <- cbind(mat1, A)
#   colnames(mat1)[i] <- label
# }
#   return(mat1)
# }
#
# x <- assay(rld)
# new <- avg_by_sample(data = x,
#                      n = 48,
#                      reps = 4)
#