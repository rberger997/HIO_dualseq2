#'---
#' title: "RNA seq workflow - Volcano plots relative to STM"
#' subtitle: "HIOs dual seq experiment 2"
#' author: "Ryan Berger"
#' date: "2018-10-26"
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
  df$colors <- factor(df$colors, levels = c('Non significant', 'Decreasing', 'Increasing'))
  
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
                                  'Significant decrease',
                                  'Significant increase'))+
    theme_bw()+ 
    theme(panel.grid.major = element_line(linetype = "blank"), 
          panel.grid.minor = element_line(linetype = "blank"), 
          axis.title = element_text(size = 14, face = "bold"), 
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 12),
          strip.text.x = element_text(size = 14, face = 'bold'),
          strip.text.y = element_text(size = 14, face = 'bold'))+
    facet_grid(rows = vars(time), cols = vars(label))+
    geom_text(data = ann_decr,label = ann_decr$n,
              color = 'blue', fontface = 'bold')+
    geom_text(data = ann_incr,label = ann_incr$n,
              color = 'red', fontface = 'bold')+
    xlim(-lim,lim)
  
  
}

#+ fig2, fig.height = 7, fig.width = 10, fig.align = 'center'


#' # STM mutants plots
mut_data <- dplyr::filter(data2, label %in% c('SPI1', 'SPI2'))
volcano_plot(mut_data, lim = 5)







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



#' # Serovars volcano plot
ser_data <- dplyr::filter(data2, label %in% c('SE', 'ST'))
volcano_plot(ser_data, lim = 5, plim = 1e-16)



#+ include=F
# ## Serovars plot with labels
# library(ggrepel)
# volcano_plot(ser_data)+
#   xlim(-5,5)+
#   geom_text_repel(data          = subset(ser_data, abs(log2FoldChange)>6 & colors != 'Non significant' | -log10(pvalue) > 6 & colors != 'Non significant'),
#                   #nudge_y       = 32 - subset(nba, PTS > 25)$PTS,
#                   #segment.size  = 0.2,
#                   segment.color = "grey50")


#' ## Serovars interactive plot
ser_data1 <- dplyr::filter(data2, label %in% c('SE', 'ST') & colors != 'Non significant')
p <- volcano_plot(ser_data1)+
  xlim(-5,5)
ggplotly(p)

p



# Look at lists of significant genes by sample
head(mut_data1)

mut_data1 %>% 
  filter(label == 'SPI2') %>% 
  arrange(padj) %>% 
  select(symbol, log2FoldChange, padj, name, label, time) %>% 
  head(., 20)






head(ser_data1)

x <- ser_data1 %>% 
  filter(label == 'SE') %>% 
  arrange(padj) %>% 
  select(symbol, log2FoldChange, padj, name, label, time) %>% 
  mutate(change = as.factor(ifelse(log2FoldChange > 0, 'Higher in SE', 'Lower in SE'))) %>% 
  head(., 20)

head(x)
# SE shows most significant gene changes at 8h



ser_data1 %>% 
  filter(label == 'ST') %>% 
  arrange(padj) %>% 
  select(symbol, log2FoldChange, padj, name, label, time) %>% 
  head(., 20)

# ST shows most significant gene changes at 2h


# ST significant genes 2h
a <- ser_data1 %>% 
  filter(label == 'ST' & time == '2h') %>% 
  arrange(padj) %>% 
  select(symbol, log2FoldChange, padj, name, label, time) %>% 
  mutate(change = as.factor(ifelse(log2FoldChange > 0, 'Higher in ST', 'Lower in ST'))) %>% 
  head(., 20)


# ST significant genes 8h
b <- ser_data1 %>% 
  filter(label == 'ST' & time == '8h') %>% 
  arrange(padj) %>% 
  select(symbol, log2FoldChange, padj, name, label, time) %>% 
  mutate(change = as.factor(ifelse(log2FoldChange > 0, 'Higher in ST', 'Lower in ST'))) %>% 
  head(., 20)



# Barplot from significant genes
ggplot(a, 
       aes(x = reorder(symbol, log2FoldChange), 
           y = log2FoldChange, fill = change))+
  geom_bar(stat = 'identity', position = 'dodge')+
  coord_flip()+
  theme(panel.grid.minor = element_line(linetype = "blank"))+
  theme_bw()+
  scale_fill_manual(values = c('Higher in ST' = 'red', 'Lower in ST' = 'blue'), name = 'Trend')+
  labs(x = '')+
  #ylim(-2.5, 2.6)+
  geom_hline(yintercept = 0)+
  ggtitle('Top 20 most significant gene changes ST/STM 2h p.i.')



ggplot(x, 
       aes(x = reorder(symbol, log2FoldChange), 
           y = log2FoldChange, fill = change))+
  geom_bar(stat = 'identity', position = 'dodge')+
  coord_flip()+
  theme(panel.grid.minor = element_line(linetype = "blank"))+
  theme_bw()+
  scale_fill_manual(values = c('Higher in SE' = 'red', 'Lower in SE' = 'blue'), name = 'Trend')+
  labs(x = '')+
  #ylim(-2.5, 2.6)+
  geom_hline(yintercept = 0)+
  ggtitle('Top 20 most significant gene changes ST/STM 8h p.i.')




# Function for making barplots of significant gene changes

gene_barplot <- function(df, samp, t, n, remove = NA){
  # Set colors for barplot
  labs <- c('red', 'blue')
  names(labs) <- c('Increased vs. STM', 'Decreased vs. STM')
  
  
  # Filter top genes
  df %>% 
    filter(label == samp & time == t & !symbol %in% remove) %>% 
    arrange(padj) %>% 
    select(symbol, log2FoldChange, padj, name, label, time) %>% 
    mutate(change = as.factor(ifelse(log2FoldChange > 0, 
                                     'Increased vs. STM', 
                                     'Decreased vs. STM'))) %>% 
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
    labs(x = '')+
    geom_hline(yintercept = 0)+
    ggtitle(paste0(samp, '/STM ', t, ' p.i.'))+
    ylim(-9,9)
  
}

# Function to save barplots
save_barplot <- function(p, w, h, file_label){
  # Create directory
  dir.create(here('img/gene_barplots/'), recursive = T)
  # Save image
  png(filename = here(paste0('img/gene_barplots/', file_label, '_gene_barplot.png')), width = w, height = h, units = 'in', res = 300)
  print(p)
  dev.off()
  
  
}




# SE barplots  
p1 <- gene_barplot(ser_data1, samp = 'SE', t = '2h', n = 20, 
             remove = c('CTB-63M22.1'))
p2 <- gene_barplot(ser_data1, samp = 'SE', t = '8h', n = 20)

save_barplot(p1, w = 4, h = 4, file_label = 'SE_2h')
save_barplot(p2, w = 4, h = 4, file_label = 'SE_8h')


# ST barplots
p3 <- gene_barplot(ser_data1, 'ST', '2h', 20,
             remove = c('TFAP2B'))
p4 <- gene_barplot(ser_data1, 'ST', '8h', 20)

save_barplot(p3, w = 4, h = 4, file_label = 'ST_2h')
save_barplot(p4, w = 4, h = 4, file_label = 'ST_8h')



# SPI1 barplots
gene_barplot(mut_data1, 'SPI1', '2h', 20)
gene_barplot(mut_data1, 'SPI1', '8h', 20,
             remove = c('RPL3P4'))


# SPI2 barplots
gene_barplot(mut_data1, 'SPI2', '2h', 20,
             remove = c('CTB-63M22.1'))
gene_barplot(mut_data1, 'SPI2', '8h', 20)





# Get significant gene lists
m_sig2 <- dplyr::filter(mut_data, padj < 0.01 & abs(log2FoldChange) > 1 &
                              baseMean > 10 & label == 'SPI1' & time == '2h')


m_sig8 <- dplyr::filter(mut_data, padj < 0.01 & abs(log2FoldChange) > 1 &
                              baseMean > 10 & label == 'SPI1' & time == '8h')




s_sig2 <- dplyr::filter(ser_data, padj < 0.01 & abs(log2FoldChange) > 1 &
                              baseMean > 10 & label == 'ST' & time == '2h')

s_sig8 <- dplyr::filter(ser_data, padj < 0.01 & abs(log2FoldChange) > 1 &
                              baseMean > 10 & label == 'ST' & time == '8h')




#+ render, include=F
# Render source file to html 
# dir.create(here('results/DESeq2_human/src_html_output'))

# render.dir <- here('results/DESeq2_human/src_html_output/')

# render(here('src/Hs_align_src/XX_volcano_plots_stm.R'), output_dir = render.dir, intermediates_dir = render.dir, clean = TRUE)

