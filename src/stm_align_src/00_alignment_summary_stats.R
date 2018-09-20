# Create summary table of alignment results

library(readr)
library(here)
library(dplyr)

# Load stats
tab <- read_csv(here('data/alignment_stats/ALL_alignment_stats.csv')) %>% 
  rename(aligned_to = short_name)

head(tab)

tab$code_name <- tab$Description %>% 
  gsub('-[0-9]$','',.)


x <- tab %>% 
  dplyr::filter(code_name %in% c('PBS_2h', 'PBS_8h',
                                 'STM_2h', 'STM_8h',
                                 'ST_2h', 'ST_8h',
                                 'SE_2h', 'SE_8h',
                                 'SPI1_2h', 'SPI1_8h',
                                 'SPI2_2h', 'SPI2_8h')) %>% 
  group_by_(.dots=c('aligned_to','code_name')) %>% 
  summarise(avg = round(mean(n_pseudoaligned), 0)) %>% 
  tidyr::spread(key = 'aligned_to', value = avg) %>% 
  rename(Sample = code_name)

x

