#'---
#' title: "RNA seq workflow - Gene level analysis"
#' subtitle: "HIOs dual seq experiment 2"
#' author: "Ryan Berger"
#' date: "2018-11-1"
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


#+ script, include=F

# Make volcano plots of samples over STM for gene comparisons
# Use functions from volcano plot script

library(dplyr)
library(ggplot2)
library(here)
library(knitr)
library(magrittr)
library(rmarkdown)
library(stringi)

# Function to prep data
prep_df <- function(df, label, time){
  # Set label for conditions
  df$label <- label
  df$time <- time
  
  # Set up labels for colors  
  df$colors <- 'Non significant'
  df[which(df$log2FoldChange < -1 & df$padj < 0.05), 'colors'] <- 'Decreasing'
  df[which(df$log2FoldChange > 1 & df$padj < 0.05), 'colors'] <- 'Increasing'  
  df$colors <- factor(df$colors, levels = c('Non significant', 'Increasing', 'Decreasing'))
  
  return(df)
}


# Locate all the files
files <- list.files(here('results/DESeq2_human/diff_expression_stm/'))
files

# Create an empty object to put the files into
data1 <- NULL

# Iterate over all files and combine into a single dataframe
for(i in seq(files)){
  file <- files[i]
  # Extract sample label and time point
  label <- gsub('_.*$',"",file)
  time <- stri_extract_first_regex(file, '[0-9]h')
  
  # Load file
  temp <- read.csv(here(paste0('results/DESeq2_human/diff_expression_stm/',file)), header = T)
  
  temp <- prep_df(temp, label = label, time = time)
  
  data1 <- rbind(data1, temp)
  print(paste(i, 'of', length(files), 'done:', files[i]))
}


# Set order of samples for plot
data1$label <- factor(data1$label, levels = c('SE','ST','SPI1','SPI2'))



# Identify long non coding RNAs
lncrna <- dplyr::filter(data1, grepl("^RP([0-9]).*",symbol) & is.na(entrez) == T)

# Select all rows in data1 that are not in lncrna
data2 <- dplyr::anti_join(data1, lncrna)



# Make plot - facet time into rows, samples into columns

volcano_plot <- function(input, lim=5, plim = 1e-16){
  
  # Get labels for number of significant genes
  ann_decr <- input %>% 
    filter(colors == 'Decreasing') %>% 
    group_by(label, time) %>% 
    dplyr::count(label, colors = colors) %>% 
    mutate(padj = plim,
           log2FoldChange = -lim+.2)
  
  ann_incr <- input %>% 
    filter(colors == 'Increasing') %>% 
    group_by(label, time) %>% 
    dplyr::count(label, colors = colors) %>% 
    mutate(padj = plim,
           log2FoldChange = lim-.2)
  
  
  # Volcano plot
  ggplot(data = input, aes(x=log2FoldChange, y=-log10(padj), 
                           color=colors, label=symbol))+
    geom_point(size=.85)+
    geom_vline(xintercept = 0, linetype = 'dotted')+
    geom_hline(yintercept = 0, linetype = 'dotted')+
    scale_color_manual(name = '',
                       values = c('Non significant' = 'gray',
                                  'Decreasing' = 'blue',
                                  'Increasing' = 'red'),
                       labels = c('Non significant',
                                  'Significant increase',
                                  'Significant decrease'))+
    theme_bw()+ 
    theme(panel.grid.major = element_line(linetype = "blank"), 
          panel.grid.minor = element_line(linetype = "blank"), 
          axis.title = element_text(size = 14, face = "bold"), 
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 12),
          strip.text.x = element_text(size = 14, face = 'bold'),
          strip.text.y = element_text(size = 14, face = 'bold'))+
    #facet_grid(rows = vars(time), cols = vars(label))+
    geom_text(data = ann_decr,label = ann_decr$n,
              color = 'blue', fontface = 'bold', size = 5)+
    geom_text(data = ann_incr,label = ann_incr$n,
              color = 'red', fontface = 'bold', size = 5)+
    xlim(-lim,lim)+
    ylim(0, -log10(plim))+
    ggtitle(paste0(unique(input$label), '/STM ',unique(input$time),' pi'))+
    theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
  
  
}

#+ fig2, fig.height = 7, fig.width = 10, fig.align = 'center'

# Save volcano plot function
save_volcano <- function(gg, sample, time){
  # Create folder for volcano plots in img directory
  dir.create(here('img/volcano_plots'))
  
  # Save plot
  png(filename = here(paste0('img/volcano_plots/',sample,'_over_stm_', time, '_volcano_plot.png')), width = 6, height = 4, units = 'in', res = 300)
  print(gg)
  dev.off()
  
  # Save ggplot object
  saveRDS(gg, 
          file = here(paste0('img/ggplot_objects/gg_',sample,'_over_stm',time,'_volcano.rds')))
}



# Iterate through all plots
for(i in c('SPI1', 'SPI2', 'SE', 'ST')){
  for(j in c('2h', '8h')){
    # Set p limit for volcano plots
    plim <- ifelse(i %in% c('SPI1', 'SPI2'), 1e-8, 1e-15)

    # Make volcano plot
    vplot <- dplyr::filter(data2, label == i & time == j) %>% 
      arrange(colors) %>% 
      volcano_plot(., lim = 5, plim = plim)
    
    # Save volcano plot
    save_volcano(vplot, sample = i, time = j)
    
    # Update progress
    print(paste(i, j, 'done'))
  }
}



#' # STM mutants plots
mut_plot <- dplyr::filter(data2, label == c('SPI1', 'SPI2')) %>% 
  arrange(colors) %>% 
  volcano_plot(., lim = 5)+
  facet_grid(rows = vars(time), cols = vars(label))
mut_plot
save_volcano(mut_plot, sample = 'stm_muts', time = '2h')

# Save plot
png(filename = here('img/spi1_over_stm_2h_volcano_plot.png'), width = 6, height = 4, units = 'in', res = 300)
mut_plot
dev.off()

saveRDS(mut_plot, 
        file = here('img/ggplot_objects/gg_spi1_over_stm_volcano.rds'))



#' # Serovars volcano plot
ser_plot <- dplyr::filter(data2, label %in% c('SE', 'ST')) %>% 
  arrange(colors) %>% 
  volcano_plot(., lim = 5, plim = 1e-16)+
  facet_grid(rows = vars(time), cols = vars(label))
ser_plot

# Save plot
png(filename = here('img/ser_over_stm_volcano_plot.png'), width = 9, height = 5, units = 'in', res = 300)
ser_plot
dev.off()

# Save ggplot objects
saveRDS(ser_plot, 
        file = here('img/ggplot_objects/gg_ser_over_stm_volcano.rds'))




# Barplots of top gene changes

# set up 
mut_data1 <- filter(data2, label %in% c('SPI1', 'SPI2') & 
                      colors != 'Non significant')


ser_data1 <- filter(data2, label %in% c('SE', 'ST') & 
                      colors != 'Non significant')




# Function for making barplots of significant gene changes
gene_barplot <- function(df, samp, t, n, remove = NA, lim){
  # Set colors for barplot
  labs <- c('red', 'blue')
  names(labs) <- c('Increased vs. STM', 'Decreased vs. STM')
  
  
  # Filter top genes
  df %>% 
    filter(label == samp & time == t & !symbol %in% remove) %>% 
    arrange(-abs(log2FoldChange)) %>% 
    select(symbol, log2FoldChange, padj, name, label, time) %>% 
    mutate(change = as.factor(ifelse(log2FoldChange < 0, 
                                     'Decreased vs. STM', 
                                     'Increased vs. STM'))) %>% 
    head(., n) %>% 
    
    # Make barplot
    ggplot(., 
           aes(x = reorder(symbol, log2FoldChange), 
               y = log2FoldChange, fill = change))+
    geom_bar(stat = 'identity', position = 'dodge')+
    coord_flip()+
    theme(panel.grid.minor = element_line(linetype = "blank"))+
    theme_bw()+
    scale_fill_manual(values = labs, name = 'Trend')+
    labs(x = '',
         title = paste0(samp, '/STM ', t, ' pi'),
         subtitle = paste('Top', n, 'genes by fold change (padj < 0.05)'))+
    theme(plot.title = element_text(size=12, hjust=0.5, 
                                    face="bold", vjust=-2))+
    theme(plot.subtitle = element_text(size=8, hjust=0.5, vjust = -1))+
    geom_hline(yintercept = 0)+
    ylim(-lim,lim)
  
}


# Function to save barplots (png and RDS)
save_barplot <- function(p, w, h, file_label){
  # Create directory
  dir.create(here('img/gene_barplots/'), recursive = T)
  # Save image
  png(filename = here(paste0('img/gene_barplots/', file_label, '_gene_barplot.png')), width = w, height = h, units = 'in', res = 300)
  print(p)
  dev.off()
  
  # Save plot as ggplot object using saveRDS
  saveRDS(p, 
          file = here(paste0('img/ggplot_objects/gg_',file_label,'_gene_barplot.rds')))
  
}




# SE barplots 
lim <- 30
gene_barplot(ser_data1, samp = 'SE', t = '2h', n = 20, lim = lim,  # -12, 5
                   remove = c('CTB-63M22.1')) %>% 
  save_barplot(., w = 5, h = 4, file_label = 'SE_2h')

gene_barplot(ser_data1, samp = 'SE', t = '8h', n = 20, lim = lim) %>%  #-25, 5
  save_barplot(., w = 5, h = 4, file_label = 'SE_8h')


# ST barplots
gene_barplot(ser_data1, 'ST', '2h', 20, lim = lim,
                   remove = c('TFAP2B')) %>%  #-25, 30
  save_barplot(., w = 5, h = 4, file_label = 'ST_2h')

gene_barplot(ser_data1, 'ST', '8h', 20, lim = lim) %>%  #-6, 8
  save_barplot(., w = 5, h = 4, file_label = 'ST_8h')




# SPI1 barplots
lim <- 12
gene_barplot(mut_data1, 'SPI1', '2h', 20, lim = lim) %>%  #-10, 5
  save_barplot(., w = 5, h = 4, file_label = 'SPI1_2h')

gene_barplot(mut_data1, 'SPI1', '8h', 20, lim = lim, 
             remove = c('RPL3P4')) %>% # -3, 7
  save_barplot(., w = 5, h = 4, file_label = 'SPI1_8h')



# SPI2 barplots
gene_barplot(mut_data1, 'SPI2', '2h', 20, lim = lim,
             remove = c('CTB-63M22.1')) %>%  # -6, 4
  save_barplot(., w = 5, h = 4, file_label = 'SPI2_2h')

gene_barplot(mut_data1, 'SPI2', '8h', 20, lim = lim) %>%  # -3, 7.5
  save_barplot(., w = 5, h = 4, file_label = 'SPI2_8h')





#+ render, include=F
# Render source file to html 
# dir.create(here('results/DESeq2_human/src_html_output'))

# render.dir <- here('results/DESeq2_human/src_html_output/')

# render(here('src/Hs_align_src/XX_volcano_plots_stm.R'), output_dir = render.dir, intermediates_dir = render.dir, clean = TRUE)








# #+ include=F
# ## STM mutants plots with labels
# volcano_plot(mut_data)+
#   xlim(-5,5)+
#   ggrepel::geom_text_repel(data = subset(mut_data,
#                                 abs(log2FoldChange)>3 &
#                                   colors != 'Non significant' |
#                                   -log10(pvalue) > 6 &
#                                   colors != 'Non significant'),
#                   #nudge_y       = 32 - subset(nba, PTS > 25)$PTS,
#                   #segment.size  = 0.2,
#                   segment.color = "grey50")


#' #' ## STM mutants interactive plot
#' library(plotly)
#' 
#' # STM mutants interactive
#' mut_data1 <- dplyr::filter(data2, label %in% c('SPI1', 'SPI2') & colors != 'Non significant')
#' p1 <- volcano_plot(mut_data1)+
#'   xlim(-5,5)
#' ggplotly(p1)





#' #+ include=F
#' # ## Serovars plot with labels
#' # library(ggrepel)
#' # volcano_plot(ser_data)+
#' #   xlim(-5,5)+
#' #   geom_text_repel(data          = subset(ser_data, abs(log2FoldChange)>6 & colors != 'Non significant' | -log10(pvalue) > 6 & colors != 'Non significant'),
#' #                   #nudge_y       = 32 - subset(nba, PTS > 25)$PTS,
#' #                   #segment.size  = 0.2,
#' #                   segment.color = "grey50")
#' 
#' 
#' #' ## Serovars interactive plot
#' ser_data1 <- dplyr::filter(data2, label %in% c('SE', 'ST') & colors != 'Non significant')
#' p <- volcano_plot(ser_data1)+
#'   xlim(-5,5)
#' ggplotly(p)
#' 
#' p




# # Get significant gene lists
# m_sig2 <- dplyr::filter(mut_data, padj < 0.01 & abs(log2FoldChange) > 1 &
#                           baseMean > 10 & label == 'SPI1' & time == '2h')
# 
# 
# m_sig8 <- dplyr::filter(mut_data, padj < 0.01 & abs(log2FoldChange) > 1 &
#                           baseMean > 10 & label == 'SPI1' & time == '8h')
# 
# 
# 
# 
# s_sig2 <- dplyr::filter(ser_data, padj < 0.01 & abs(log2FoldChange) > 1 &
#                           baseMean > 10 & label == 'ST' & time == '2h')
# 
# s_sig8 <- dplyr::filter(ser_data, padj < 0.01 & abs(log2FoldChange) > 1 &
#                           baseMean > 10 & label == 'ST' & time == '8h')

