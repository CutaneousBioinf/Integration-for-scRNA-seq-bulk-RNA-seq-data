
library(affy)
library(BiocGenerics)
library(ArrayTools)
library(xbioc)
library(tidyr)
library(MuSiC)

#scRNA
pDataPN = read.csv('sampledPNcellType.csv',header = T, row.names = 1)
pDataPN$cellTypeID = as.numeric(pDataPN$cellType)
matrixPN = read.csv('sampledrawMatrix.csv',header = T, row.names = 1)

col = colnames(matrixPN)
row = row.names(pDataPN)
sum(col == row)

set = ExpressionSet(as.matrix(matrixPN))
set = ExpressionSet(as.matrix(PNmarkers))
phenoData <- new("AnnotatedDataFrame",data=pDataPN)

scSet <- ExpressionSet(assayData=as.matrix(matrixPN),
                       phenoData=phenoData,
                       annotation="hgu95av2")
#Bulk RNA data PN
matrixNN = read.table('counts_ddsHTSeq_collapsed_NN', header = T)

bulksetNN <- ExpressionSet(assayData=as.matrix(matrixNN))



bulksetNN <- ExpressionSet(assayData=as.matrix(normalized_counts))
bulksetNN <- ExpressionSet(assayData=as.matrix(normalized_countsHi))
bulksetNN <- ExpressionSet(assayData=as.matrix(logCPM2))
bulksetNN <- ExpressionSet(assayData=as.matrix(cpm))

bulksetNN <- ExpressionSet(assayData=as.matrix(n3))
bulksetNN = ExpressionSet(assayData=as.matrix(NNmarkers))


#normalize bulkrna data
library(edgeR)
library(EnrichmentBrowser)

norm.eset = normalize(bulksetNN, data.type = 'rseq')
bulksetNN <- ExpressionSet(assay(norm.eset))


#run music
cell = unique(pDataPN$cellType)
sample_num = unique(pDataPN$sampleID)

Est.prop.GSE50244 = music_prop(bulk.eset = bulksetNN, 
                               sc.eset = scSet, clusters = 'cellType',
                               samples = 'sampleID',select.ct = cell,
                               verbose = F)

write.csv(Est.prop.GSE50244$Est.prop.weighted,'onlymarkersNormboth.csv')
write.csv(Est.prop.GSE50244$Weight.gene,'weight.csv')
View(Est.prop.GSE50244$Est.prop.weighted)


data = Est.prop.GSE50244$Est.prop.weighted
write.csv(Est.prop.GSE50244$Est.prop.weighted, 'resup.csv')
rownames = row.names(data)
##read unsampled annotation data to see the propotion
data = as.data.frame(data)
pDataPN_raw = read.csv('PN.celltype.csv',header = T)
colnames = colnames(Est.prop.GSE50244$Est.prop.weighted)
prop = rep(0, length(colnames))
for (i in 1:length(colnames)) {
  prop[i] = sum(pDataPN_raw$cellType== colnames[i])
}
prop = prop/length(pDataPN_raw$cellType)
data1 = rbind(data, prop)
write.csv(data1, 'new.csv')

#proportion plot
res = read_csv('res_withCompare.csv')
res = res %>% remove_rownames %>% column_to_rownames(var=colnames(res)[1])
res

res$category <- row.names(res)
res_2 <- melt(res, id.vars = "category")
res_2

(p <- ggplot(res_2, aes(category, value, fill = variable)) +
    geom_bar(position = "fill", stat = "identity") 
)