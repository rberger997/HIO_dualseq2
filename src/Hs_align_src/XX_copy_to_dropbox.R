# Copy files to dropbox folder

#' The purpose of this script is to copy files from the HIOs_dualseq2 directory into the dropbox folder. The files to be copied are:
#' 
#' * Results folder (except kallisto output)
#' * `img` folder
#' 
#' This script will be sourced at the end of each script to update the contents of the dropbox folder to match the repository.
#' 


## Move files into blog post directory ##

# Make copy of img folder and move to 'img' folder of Dropbox/DATA/ directory
file.copy(from = here('img'),
          to = '~/Dropbox/HIO Salmonella manuscript/DATA/', 
          recursive = T, overwrite = T)

# Copy results folder and move to Dropbox/DATA/ directory
file.copy(from = here('results/DESeq2_human/'),
          to = '~/Dropbox/HIO Salmonella manuscript/DATA/',
          recursive = T, overwrite = T)


# Finished message
print('Files copied to dropbox')