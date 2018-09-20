#'---
#' title: "RNA seq workflow part 3 - Differential expression"
#' subtitle: "HIOs dual seq experiment 2 - ST alignment"
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
library(magrittr)
library(here)

dir.create(here('results/DESeq2_st/diff_expression'))
results.dir <- here('results/DESeq2_st/diff_expression/')


#' ## Load data
dds <- readRDS(file = here('results/DESeq2_st/dds_all.rds'))


# Run DESeq
dds <- DESeq2::DESeq(dds)

# Create function for making differential expression dataframe
res <- results(dds, 
               contrast = c('code_name', 'ST_8h', 'ST_2h'))%>% 
  as.data.frame()

res$Protein.product <- rownames(res) %>% 
  gsub('^.*_cds_','',.) %>% 
  gsub('_[0-9]{1,4}$','',.)

head(res)

## Annotate transcripts with gene IDs/names
st.anno <- read.csv(here('data/ST_annotation.csv'), 
                    stringsAsFactors = FALSE) %>% 
  dplyr::select(Protein.product, Locus, Protein.name)

nrow(res)
length(unique(res$Protein.product))

## Add annotations to results
res <- left_join(res, st.anno, by = 'Protein.product') %>% 
  arrange(padj)
res <- res[,c(7,1:6,8,9)]

res$Protein.name <- gsub('MULTISPECIES: ','', res$Protein.name)

res$Protein.name[is.na(res$Protein.name)] <- 'unknown protein'

res <- dplyr::rename(res, name = Protein.name,
                     symbol = Locus)


## Save output
write.csv(res, 
          file = file.path(results.dir, "ST_8h_over_ST_2h_diffexpress.csv"),
          row.names = F)

nrow(res)
