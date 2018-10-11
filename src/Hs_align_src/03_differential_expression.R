#'---
#' title: "RNA seq workflow part 3 - Differential expression"
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
#' This script is designed to calculate differential gene expression between two samples. The input is a DESeq2 `dds` object and the output is a `.csv` file for each pair of samples compared.

#+ 
#' # Begin script
#' 
#' ## Libraries and directories
library(AnnotationDbi)
library(DESeq2)
library(dplyr)
library(here)
library(knitr)
library(magrittr)
library(org.Hs.eg.db)
library(rmarkdown)

## Optional: install org.Hs.eg.db
#source("https://bioconductor.org/biocLite.R")
#biocLite("org.Hs.eg.db") 

# Create directory for diff expression .csv files
dir.create(here('results/DESeq2_human/diff_expression'))

#+ set_cache_dir, include=F, warning=F
# Set a directory for all cache files to go
dir.create(here('results/DESeq2_human/src_html_output/knitr_cache'))
opts_chunk$set(cache.path = here('results/DESeq2_human/src_html_output/knitr_cache/'))


#' ## Load data
#' We previously saved the `dds` object in the results folder. Let's load it now:

#+ prepdata, cache=T
dds <- readRDS(here('results/DESeq2_human/dds_all.rds'))

# Prep dds for differential expression
dds <- DESeq(dds)

#' ## Calculate differential expression
#' This experiment contains six injection conditions (PBS, STM, SE, ST, SPI1, SPI2) and two time points (2h, 8h). To start, we'll calculate the differential expression of each injection over PBS at the matching time point (e.g. STM 8h over PBS 8h, ST 2h over PBS 2h). We'll save each of these as a .csv file in the `diff_expression` folder in the results directory.

# Set up the samples for diff expression
unique(colData(dds)$code_name)

# Select all but PBS
samps <- unique(colData(dds)$code_name)[3:12] %>% 
  as.character()


#' We'll set up a loop to calculate differential expression and iterate over all the samples. We'll make sure each is matched to the proper PBS control with a logical test at the start of the loop that checks the time point of each sample.

#+ loop, cache=T
# Set up loop to calculate differential expression for all samples over PBS
for(i in seq(samps)){
  
  # Sample for i in loop
  sample <- samps[i]
  
  # Check which time point the i sample is to match PBS control
  if(grepl('2h',sample, fixed = T) == T){
    PBS <- 'PBS_2h'}
    else{
      PBS <- 'PBS_8h'
    }

# Calculate differential expression - sample over PBS
res <- results(dds, 
               contrast = c('code_name', sample, PBS))


# Add annotation - symbol and entrezID
res$symbol <- rownames(res)
# Add column for gene Entrez ID
res$entrez <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     column = 'ENTREZID',
                     keytype = 'SYMBOL',
                     multiVals = 'first')
# Add column for gene name
res$name <- mapIds(org.Hs.eg.db,
                   keys = rownames(res),
                   column = 'GENENAME',
                   keytype = 'SYMBOL',
                   multiVals = 'first')
# Make results dataframe


res.df <- as.data.frame(res) %>% 
  arrange(padj)


# Save output file
file.name <- paste0(here('results/DESeq2_human/diff_expression/'),
                    sample,'_over_',PBS,'_diffexpress.csv')
write.csv(res.df, file = file.name, row.names = F)
}



#' ## Differential expression over STM

#' Want to compare the mutants and serovars to STM so we'll calculate differential expression over STM directly and save those files in a separate folder.


# Samples to compare to STM
samps <- samps[3:10]


# Create directory for diff expression files over STM
dir.create(here('results/DESeq2_human/diff_expression_stm'))

#+ loop1, cache=T
# Set up loop to calculate differential expression for all samples over PBS
for(i in seq(samps)){
  
  # Sample for i in loop
  sample <- samps[i]
  
  # Check which time point the i sample is to match PBS control
  if(grepl('2h',sample, fixed = T) == T){
    STM <- 'STM_2h'}
  else{
    STM <- 'STM_8h'
  }
  
  # Calculate differential expression - sample over PBS
  res <- results(dds, 
                 contrast = c('code_name', sample, STM))
  
  
  # Add annotation - symbol and entrezID
  res$symbol <- rownames(res)
  # Add column for gene Entrez ID
  res$entrez <- mapIds(org.Hs.eg.db,
                       keys = rownames(res),
                       column = 'ENTREZID',
                       keytype = 'SYMBOL',
                       multiVals = 'first')
  # Add column for gene name
  res$name <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     column = 'GENENAME',
                     keytype = 'SYMBOL',
                     multiVals = 'first')
  # Make results dataframe
  
  
  res.df <- as.data.frame(res) %>% 
    arrange(padj)
  
  
  # Save output file
  file.name <- paste0(here('results/DESeq2_human/diff_expression_stm/'),
                      sample,'_over_',STM,'_diffexpress.csv')
  write.csv(res.df, file = file.name, row.names = F)
}



#' ## Ending notes
#' The differential expression outputs are now in the results directory. These will be used to make MA and volcano plots.




#+ render, include=F
# Render source file to html 
# dir.create(here('results/DESeq2_human/src_html_output'))

# render.dir <- here('results/DESeq2_human/src_html_output/')

# render(here('src/Hs_align_src/03_differential_expression.R'), output_dir = render.dir, intermediates_dir = render.dir, clean = TRUE)