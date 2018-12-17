
# RNA seq setup script



# Install packages
install.packages('dplyr')
install.packages('here')
install.packages('knitr')
install.packages('rmarkdown')
install.packages('tidyr')
install.packages('pheatmap')
install.packages('plotly')
install.packages('readr')
install.packages('eulerr')


# Install packages from Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("tximport")
biocLite("rhdf5")
biocLite("DESeq2")
biocLite("EnsDb.Hsapiens.v75")
biocLite("org.Hs.eg.db") 
biocLite("clusterProfiler") 
biocLite("GSEABase") 
biocLite("ReactomePA") 
