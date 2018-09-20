## R script to generate Wald and LRT test results from DESeq2 output
## Ryan Berger & David R. Hill
## 5-24-18
## -----------------------------------------------------------------------------

## enable parallel processes
library("BiocParallel")
library('magrittr')
register(MulticoreParam(4))
results.dir <- "../results/DESeq2_STM"

## load output of DESeq-export-counts.R
load(file = file.path(results.dir, "dds.Rdata"))

## setup multifactor design
DESeq2::design(dds) <- ~ group

## Run DESeq
dds <- DESeq2::DESeq(dds)

## Apply Wald test for specific comparisons
## 'test = "Wald"' is default; specified here for clarity
res <- DESeq2::results(dds, test = "Wald",
                       contrast = c("group", "24hpi-STM", "STM")) %>% 
  as.data.frame()
res$product_accession <- rownames(res)

## Annotate transcripts with gene IDs/names
stm.anno <- read.csv('../data/stm_annotation.csv', stringsAsFactors = FALSE)
## Filter out rows with no product accession
stm.anno <-  filter(stm.anno, !product_accession == '')
stm.anno <- stm.anno[,c('product_accession','name', 'symbol')]

## Add annotations to results
res <- left_join(res, stm.anno, by = 'product_accession') %>% 
  arrange(padj)
res <- res[,c(7,1:6,8:9)]

## Save output
write.csv(res, file = file.path(results.dir, "STM24h-over-STM8h-Wald-test.csv"))

## Clean up
rm(list = ls())