## Troubleshooting tx2gene error
## Look at tx2gene format compared to input .tsv file transcript format

## Check wd - ../STM_dual_rnaseq/src
getwd()

## create tx2gene for STM using annotation file - product accession & name
#---------------------------------------------------------
stm.anno <- read.csv('../data/stm_annotation.csv', stringsAsFactors = FALSE)
## Filter out rows with no name
stm.anno <-  filter(stm.anno, !product_accession == '')
tx2gene <- stm.anno[, c('product_accession', 'name')]
tx2gene[is.na(tx2gene$name)] <- 'unknown gene'
#---------------------------------------------------------


## Load formatted .tsv file
#---------------------------------------------------------
kallisto.results.dir.stm <- "../results/Run_2286/STM"
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
                   "abundance_formatted.tsv") 
## set sample names as description_rep#_seq_rep#
names(files) <- samples$Description
## Load input .tsv file
tsv <- read.table(file = files[1], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
#---------------------------------------------------------

# Compare the formats
head(tsv)
head(tx2gene)

# Compare structures
str(tsv)
str(tx2gene)

# Check to see if the first values match
tsv$target_id[1:10] == tx2gene$product_accession[1:10]

# Check to see how many values match
table(tsv$target_id %in% tx2gene$product_accession)

# No NAs
table(is.na(tsv$target_id))
table(is.na(tx2gene))

## tx2gene needs to be a 2 column dataframe with the transcript id match in col1 and the new id in col2
## colnames don't matter, unique IDs in col2 don't matter, having perfect match for all still 
## won't fix


## Clean up
rm(list = ls())