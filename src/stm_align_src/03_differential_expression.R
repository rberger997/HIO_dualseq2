#'---
#' title: "RNA seq workflow part 3 - Differential expression"
#' subtitle: "HIOs dual seq experiment 2 - STM alignment"
#' author: "Ryan Berger"
#' date: "2018-09-18"
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
#' This script is designed to calculate differential gene expression between two samples. The input is a DESeq2 `dds` object and the output is a `.csv` file for each pair of samples compared.




#+ 
#' # Begin script
## -----------------------------------------------------------------------------

#' ## Libraries and directories
library(dplyr)
library(DESeq2)
library('magrittr')
library(here)

dir.create(here('results/DESeq2_stm/diff_expression'))
results.dir <- here('results/DESeq2_stm/diff_expression/')


#' ## Load data
dds <- readRDS(file = here('results/DESeq2_stm/dds_all.rds'))


# Run DESeq
dds <- DESeq2::DESeq(dds)

# Create function for making differential expression dataframe
diff_express <- function(samp8h, samp2h){
res <- results(dds, 
               contrast = c('code_name', samp8h, samp2h))%>% 
  as.data.frame()

res$product_accession <- rownames(res) %>% 
  gsub('^.*_cds_','',.) %>% 
  gsub('_[0-9]{1,4}$','',.)



## Annotate transcripts with gene IDs/names
stm.anno <- read.csv(here('data/stm_annotation.csv'), stringsAsFactors = FALSE)
## Filter out rows with no product accession
stm.anno <-  dplyr::filter(stm.anno, !product_accession == '') %>% 
  dplyr::select(product_accession, name, symbol)


## Add annotations to results
res <- left_join(res, stm.anno, by = 'product_accession') %>% 
  arrange(padj)
res <- res[,c(7,1:6,8:9)]

## Save output
write.csv(res, file = file.path(results.dir, paste0(samp8h, "_over_",samp2h,"_diffexpress.csv")))

return(res)
}

# STM results
stm.res <- diff_express('STM_8h','STM_2h')

# SPI1 results
spi1.res <- diff_express('SPI1_8h', 'SPI1_2h')

# SPI2 results
spi2.res <- diff_express('SPI2_8h', 'SPI2_2h')





