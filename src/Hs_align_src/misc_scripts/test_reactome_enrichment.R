
require(clusterProfiler)
data(gcSample)
res <- compareCluster(gcSample, fun="enrichPathway")

head(gcSample)
head(res)


class(gcSample)


here()

# Load samples
st2 <- read.csv(here('results/DESeq2_human/diff_expression_stm/ST_2h_over_STM_2h_diffexpress.csv'), stringsAsFactors = F) %>% 
  dplyr::filter(pvalue < 0.05 & abs(log2FoldChange > 1))

se2 <- read.csv(here('results/DESeq2_human/diff_expression_stm/SE_2h_over_STM_2h_diffexpress.csv'), stringsAsFactors = F) %>% 
  dplyr::filter(pvalue < 0.05 & abs(log2FoldChange > 1))


head(st2)
samp <- st2$entrez
samp <- samp[complete.cases(samp)]
samp

res.st <- compareCluster(samp, fun = 'enrichPathway')
