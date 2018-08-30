# Making venn diagrams with the eulerr package

# Should make area proportional venn diagrams from lists or a matrix



library(eulerr)
eulerr_options(pointsize = 16)
options(digits = 4)
# Input in the form of a named numeric vector
fit1 <- euler(c("A" = 25, "B" = 5, "C" = 5,
                "A&B" = 5, "A&C" = 5, "B&C" = 3,
                "A&B&C" = 3))
# Input as a matrix of logicals
set.seed(1)
mat <- cbind(
  A = sample(c(TRUE, TRUE, FALSE), 50, TRUE),
  B = sample(c(TRUE, FALSE), 50, TRUE),
  C = sample(c(TRUE, FALSE, FALSE, FALSE), 50, TRUE)
)
fit2 <- euler(mat)

plot(fit2, 
     quantities = T,
     fill = c('red','gray','blue'),
     lty = 1:3)
     #labels = c('test','A','None'))
eulerr::eulerr_options()


library(here)
# Load data from volcano plots
data2 <- read.csv(here('results/DESeq2_human/volcano_data.csv'), 
                  stringsAsFactors = F)
head(data2)


# want to compare all genes increasing
incr <- filter(data2, colors == 'Increasing' & time == '2h') %>% 
  select(c(symbol,label))
library(tidyr)

se <- filter(incr, label == 'SE')
stm <- filter(incr, label == 'STM')
st <- filter(incr, label == 'ST')


a <- se$symbol
b <- stm$symbol
c <- st$symbol



  
x <- cbind(se$symbol, stm$symbol, st$symbol)
colnames(x) <- c('SE', 'STM', 'ST')
head(x)
tail(x)
dim(x)

fit3 <- euler(c(se$symbol, stm$symbol))



a

library(VennDiagram)
z <- calculate.overlap(list(a,b,c))

v <- euler(c(A=450, B=1800, "A&B"=425)) ;plot(v)




# New venn diagram maker from http://matticklab.com/index.php?title=Weighted_Venn_diagrams_in_R

source("https://bioconductor.org/biocLite.R"); biocLite(c("RBGL","graph"))

library(devtools)
install_github("js229/Vennerable"); library(Vennerable)

w <- Venn(Sets = list(a,b,c), 
          SetNames = c('SE','STM','ST'))
plot(w,
     doWeights = T)

?Venn
vignette('Venn')
?VennThemes



library(here)
# Load data from volcano plots
data2 <- read.csv(here('results/DESeq2_human/volcano_data.csv'), 
                  stringsAsFactors = F)
head(data2)


# want to compare all genes increasing
incr <- filter(data2, colors == 'Increasing' & time == '2h') %>% 
  select(c(symbol,label))
library(tidyr)

a <- filter(incr, label == 'SE')
b <- filter(incr, label == 'STM')
c <- filter(incr, label == 'ST')


w <- Venn(Sets = list(a$symbol,b$symbol,c$symbol), 
          SetNames = c('SE','STM','ST')) 
plot(w,
     doWeights = T,
     doEuler = T,
     show = list(Faces = F, DarkMatter = F))





library(VennDiagram)
#pdf("venn_diagram.pdf")
venn.plot <- venn.diagram(list(a$symbol, b$symbol, c$symbol), NULL, fill=c("red", "green", 'blue'), alpha=c(0.5,0.5,0.5), cex = 2, cat.fontface=4, category.names=c("Condition A", "Condition B", 'C'))
grid.draw(venn.plot)
dev.off()  



