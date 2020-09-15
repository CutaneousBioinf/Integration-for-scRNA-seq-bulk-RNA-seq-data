library(data.table)
library(tidyr)
library(tidyverse)

set.seed(1000)
setwd('/Users/hayden/Desktop/project1new')
pn_matrix = fread('PN.normalized.matrix.csv')
sample = sample(2:20604,2000,replace=F)
sample = sort(sample)
spMatrix = pn_matrix[,sample, with = FALSE]
fwrite(spMatrix, file = 'sampled_raw_matrix.csv')

#adding rownames
df = read.csv('sampled_raw_matrix.csv')
rownames = pn_matrix[[1]]
row.names(df) = rownames
write.csv(df, file = "sampledrawMatrix.csv", row.names = T)

#preprocess PN cell type
df1 = read.csv('PN.celltype.csv',header = T)
View(df)
df1 = unite(df1, "cell", X, sampleID, sep = ".", remove = FALSE)
rownames(df1)<-df1[,1]
df1<-df1[,c(-1,-2)]
write.csv(df1,'PN.csv')

#sampling for cell type
pn_cellType = read.csv('PN.csv')
spCellType = pn_cellType[sample-1,]
write.csv(spCellType, file = "sampledPNcellType.csv", row.names = F)