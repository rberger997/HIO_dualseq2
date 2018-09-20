#'---
#' title: "RNA seq workflow part 1 - Counts to .dds"
#' subtitle: "HIOs dual seq experiment 2 - STM alignment"
#' author: "Ryan Berger"
#' date: "2018-09-17"
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
#' This script is the first in a series of RNA seq processing scripts for converting transcript counts from a kallisto alignment into gene counts used in analyses. This is modified from a workflow by [Bioconductor](https://www.bioconductor.org/help/workflows/rnaseqGene/).

#+ section2
#' # Summary
#' **Input:**
#' 
#' * RNA seq alignment results (kallisto output)
#'     + abundance.h5 
#' * Sample key for annotating files
#'     + sample_key.csv
#' * Gene reference file (tx2gene) for converting transcript counts to genes
#'     + EnsDb.Hsapiens.v75 
#'   
#' **Output:**
#' 
#' * DESeq2 dds object

#------------------------------------------------------------
#' # Begin script

#' ## Libraries
library(DESeq2)
library(dplyr)
library(EnsDb.Hsapiens.v75)
library(here)
library(knitr)
library(rmarkdown)
library(tximport)


# Optional: install packages from Bioconductor
# source("https://bioconductor.org/biocLite.R")
# biocLite("tximport")
# biocLite("rhdf5")
# biocLite("DESeq2")
# biocLite("EnsDb.Hsapiens.v75")




#' ## Set up directories
#' We'll use the `here` package to control our working directory. This will set the root directory to `'../HIO_dualseq2/'`.

here()

# Where are the input (kallisto alignment) files located?
dir <- here('results/Run_2374/STM/')


# Where should the output files go?
output_folderID <- 'DESeq2_stm'  # Name of output folder
dir.create(here('results', output_folderID))

results.dir <- here('results/DESeq2_stm/')


#+ set_cache_dir, include=F
# Set a directory for all cache files to go
dir.create(here('results/DESeq2_stm/src_html_output'))
dir.create(here('results/DESeq2_stm/src_html_output/knitr_cache'))
opts_chunk$set(cache.path = here('results/DESeq2_human/src_html_output/knitr_cache/'))



#------------------------------------------------------------
#' ## Transcript counts to gene counts

#' Need to first set up a gene reference object (tx2gene) for converting counts of ENSEMBL transcipt IDs to counts of gene IDs. We'll use the Ensembl based annotation package `EnsDb.Hsapiens.v75` to set this up.


# associate transcripts with gene IDs
# create tx2gene for STM using annotation file - product accession & name
stm.anno <- read.csv(here('data/stm_annotation.csv'), stringsAsFactors = FALSE)
## Filter out rows with no name
stm.anno <-  dplyr::filter(stm.anno, !product_accession == '')
tx2gene <- stm.anno[, c('product_accession', 'name')]
tx2gene[is.na(tx2gene$name)] <- 'unknown gene'



#' ## File import
#' Import the `kallisto` alignment output files (abundance.h5) and annotate them with a sample key. These files are located in the results directory.

# Load in sample key.
sample_key <- read.csv(here('data/sample_key.csv'))


#' This RNA seq run contains 12 samples that were from a different experiment (the first dual-seq experiment). The column `dualseq_expt` indicates which experiment each sample is from so we'll use this to filter these samples out of our analysis.

# Optional: filter out samples from the key that you don't want to analyze
sample_key <- dplyr::filter(sample_key, dualseq_expt == 2) %>% 
  dplyr::filter(code_name %in% c('STM_2h', 'STM_8h', 
                                 'SPI1_2h', 'SPI1_8h',
                                 'SPI2_2h', 'SPI2_8h'))

# Vector of samples used in the tximport
txiSamples <- as.vector(unique(sample_key$code_name))
txiSamples

# Set up path to read files into tximport
files <- file.path(dir, sample_key$filename, 'abundance.h5')
# Add sample IDs to files
names(files) <- sample_key$Description


# import the abundance.h5 files using tximport

#+ tximport, cache=T
txi <- tximport(files,
                type = "kallisto",
                tx2gene = tx2gene, 
                ignoreTxVersion = TRUE,
                txOut = TRUE)


# export abundance counts
write.csv(txi$abundance, file = file.path(results.dir, "complete_dataset_txi.csv"))


#' ## Create a DESeq2 object (dds)

#' Need to set the design of dds object to the experimental condition of each sample without replicate IDs. In our case it is the `code_name` from the `sample_key`. 

# Create DESeq dataset
dds <- DESeqDataSetFromTximport(txi, 
                                colData = sample_key,
                                design = ~code_name)
head(dds)


#' Now we'll eliminate all the genes that have zero counts across all conditions.

nrow(dds)
# Filter out rows with no counts
dds <- dds[rowSums(counts(dds)) > 1, ]
nrow(dds) 

#' This removed ~3000 genes from the dataset.

# Account for transcript length
dds <- DESeq2::estimateSizeFactors(dds)
ddscounts <- DESeq2::counts(dds, normalized = TRUE)

# Save output file:
write.csv(ddscounts, file = file.path(results.dir, "complete-dataset_DESeq2-normalized-counts.csv"))

# Save dds object
saveRDS(dds, file = file.path(results.dir, 'dds_all.rds'))

#' ## Ending notes
#' RNA seq reads are now saved in a DESeq dataset (dds object). This can be loaded for downstream analysis.




#+ render, include=F
# Render source file to html 
# dir.create(here('results/DESeq2_human/src_html_output'))

# render.dir <- here('results/DESeq2_human/src_html_output/')

# render(here('src/Hs_align_src/01_counts_to_dds.R'), output_dir = render.dir, intermediates_dir = render.dir, clean = TRUE)