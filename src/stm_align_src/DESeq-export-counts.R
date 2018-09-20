## R script to export kallisto results to matrix with DESeq2
## Ryan Berger & David R. Hill
## 5-24-2018
##
##  This script is used for importing the kallisto alignment against the Salmonella
##  typhimurium genome. The format of the input files are '.tsv'. The transcript IDs 
##  had to be string split and saved as new '_formatted.tsv' files before analyzing.
##  There is an issue with the tx2gene annotation but currently using a workaround.
## -----------------------------------------------------------------------------

## Check directory - ../STM_dual_rnaseq/src
getwd()

## load prerequisites
library(magrittr)
library(dplyr)

## Differential expression of kallisto results with DESeq2
kallisto.results.dir.hs <- "../results/Run_2286/H_sapiens"
kallisto.results.dir.stm <- "../results/Run_2286/STM"

## create directory to deposit results
results.dir <- "../results/DESeq2_STM"
dir.create(path = results.dir, recursive = TRUE)

## read in table with sample metadata
## load using the UM core provided sample submission form
samples <- readr::read_csv(file = "../data/Run_2286/Run_2286_oriordan.csv",
                           skip = 18)

## create experimental design labels from Description string
samples$group <- gsub('.{2}$', '', samples$Description)
samples$hr <- c(rep('8', 8), rep('24', 2))

## setup access to kallisto read files (human alignment)
files <- file.path(kallisto.results.dir.stm,
                   paste0(samples$Sample_Name,"_S",
                          as.numeric(rownames(samples)),
                          "_L007_R1_001.fastq"),
                   "abundance.tsv") 

## set sample names as description_rep#_seq_rep#
names(files) <- samples$Description

## String split .tsv file to give product and genomic accession ids
for(i in seq_along(files)){
  ## Split string
  x <- read.table(file = files[i], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  x[,1] <- sub('^lcl.NC_[0-9]{,6}.[0-9]_cds_', '', x[,1])
  x[,1] <- sub('_[0-9]{,4}$', '',x[,1])

  ## Save new formatted files
  name <- paste(samples$Sample_Name[i],paste('_S',i, sep = '') ,'_L007_R1_001.fastq', sep = '')
  sample.dir <- file.path('../results/Run_2286/STM',name)
  readr::write_tsv(x, path = file.path(sample.dir, 'abundance_formatted.tsv'))
}

## Switch to new formatted files
files <- gsub('abundance.tsv', 'abundance_formatted.tsv', files)

## Filter only the STM samples
samples <- samples[5:10,]
files <- files[5:10]

## check that all files are found
if (all(file.exists(files)) == FALSE) {
    print("kallisto files not found")
}else{
  print('all kallisto files are found')
}

## associate transcripts with gene IDs
## create tx2gene for STM using annotation file - product accession & name
stm.anno <- read.csv('../data/stm_annotation.csv', stringsAsFactors = FALSE)
## Filter out rows with no name
stm.anno <-  filter(stm.anno, !product_accession == '')
tx2gene <- stm.anno[, c('product_accession', 'name')]
tx2gene[is.na(tx2gene$name)] <- 'unknown gene'

## import kallisto data and generate count dataframe (dds)
## http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
library(readr)
txi <- tximport::tximport(files,
                          type = "kallisto",
                          tx2gene = tx2gene, 
                          ignoreTxVersion = TRUE,
                          txOut = TRUE)

##---------------------------------------------------
## Tximport keeps giving error: Error in summarizeToGene
## Doesn't make sense because the tx2gene is in the right format and matches
## the transcript IDs exactly.

## Workaround: 
## Set txOut = TRUE
## Still gives results with the product_accession name, can match with annotation later.
##---------------------------------------------------

## export abundance counts
write.csv(txi$abundance, file = file.path(results.dir, "complete_dataset_txi.csv"))

library(DESeq2)
dds <- DESeq2::DESeqDataSetFromTximport(txi,
                                        colData = samples,
                                        design = ~ group) # hr ~ status
## pre-filter out counts < 1
dds <- dds[rowSums(counts(dds)) > 1, ]

## write out normalized expression counts
dds <- DESeq2::estimateSizeFactors(dds)
ddscounts <- DESeq2::counts(dds, normalized = TRUE)

## write expression matrix to file
write.csv(ddscounts, file =  file.path(results.dir, "complete-dataset_DESeq2-normalized-counts.csv"))
save(dds, file = file.path(results.dir, "dds.Rdata"))

## clear working memory
rm(list = ls())
