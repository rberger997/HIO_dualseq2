
# RNA seq setup script



# Install packages
install.packages('dplyr')
install.packages('here')
install.packages('knitr')
install.packages('rmarkdown')


# Install packages from Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("tximport")
biocLite("rhdf5")
biocLite("DESeq2")
biocLite("EnsDb.Hsapiens.v75")
biocLite("org.Hs.eg.db") 


