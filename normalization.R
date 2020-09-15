
# DESeq2
library(DESeq2)
matrixNN = read.table('counts_ddsHTSeq_collapsed_NN', header = T)

#exclude low expression data
matrixNN_hi <- matrixNN[rowSums(matrixNN)>2,]
#normalization
dim(matrixNN_hi)
condition <- factor(rep(1,84))
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = NNmarkers,
                                            colData = DataFrame(condition),  ~ 1)

dds <- DESeq(ddsFullCountTable)
normalized_counts <- counts(dds, normalized=TRUE)
#  Too few common genes!
n2 = t(apply(counts(dds,normalized=T),1,function(x){qnorm((rank(x)-(3/8))/(length(x)-2*(3/8)+1))}))
n3 = n2 - min(n2)
normalized_counts



#normalization
condition <- factor(rep(1,84))
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = matrixNN,
                                            colData = DataFrame(condition),  ~ 1)

dds <- DESeq(ddsFullCountTable)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_countsHi

# TMM and voom
library(limma)
dge <- DGEList(counts = matrixNN)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)
min(logCPM)
logCPM2 = logCPM - min(logCPM)
cpm = cpm(dge, log=F)
#https://www.bioinfo-scrounger.com/archives/115/


# calculate spearman correaltion
bulksetNN <- ExpressionSet(assayData=as.matrix(matrixNN))
bulksetNN <- ExpressionSet(assayData=as.matrix(normalized_counts))
bulksetNN <- ExpressionSet(assayData=as.matrix(normalized_countsHi))
bulksetNN <- ExpressionSet(assayData=as.matrix(logCPM))
bulksetNN <- ExpressionSet(assayData=as.matrix(cpm))
Est.prop.GSE50244 = music_prop(bulk.eset = bulksetNN, 
                               sc.eset = scSet, clusters = 'cellType',
                               samples = 'sampleID',select.ct = cell,
                               verbose = F)
datade = Est.prop.GSE50244$Est.prop.weighted
dataraw = Est.prop.GSE50244$Est.prop.weighted
datatm = Est.prop.GSE50244$Est.prop.weighted
datav = Est.prop.GSE50244$Est.prop.weighted
dim(data)
dim(data1)
corrTMM = rep(0, 84)
corrRaw = rep(0, 84)
corrVoom = rep(0, 84)
corrDEseq2 = rep(0, 84)
for (i in 1:84) {
  corrVoom[i] = as.numeric(cor.test(as.numeric(datav[i,]), prop,method = "spearman")$estimate)
}
for (i in 1:84) {
  corrDEseq2[i] = as.numeric(cor.test(as.numeric(data[i,]), prop,method = "spearman")$estimate)
}
for (i in 1:84) {
  corrTMM[i] = as.numeric(cor.test(as.numeric(datatm[i,]), prop,method = "spearman")$estimate)
}
for (i in 1:84) {
  corrRaw[i] = as.numeric(cor.test(as.numeric(Est.prop.GSE50244$Est.prop.weighted[i,]), prop,method = "spearman")$estimate)
}
output = as.data.frame(corrTMM)
output$TMM = corrTMM
output$rawcount = corrRaw
output$DEseq2 = corrDEseq2
output$voom = corrVoom
# Distribution of the above gene candidates (add the KRT14) 
# in bulk RNA-seq raw count / different normalized counts
# Raw
par(mfrow=c(3,2))
plot(density(as.matrix(matrixNN['FLG',])), main = 'FLG')
plot(density(as.matrix(matrixNN['LOR',])), main = 'LOR')
plot(density(as.matrix(matrixNN['FLG2',])), main = 'FLG2')
plot(density(as.matrix(matrixNN['CDSN',])), main = 'CDSN')
plot(density(as.matrix(matrixNN['ARG1',])), main = 'ARG1')
plot(density(as.matrix(matrixNN['SLURP1',])), main = 'SLURP1')
par(mfrow=c(3,2))
plot(density(as.matrix(matrixNN['PTGS1',])), main = 'PTGS1')
plot(density(as.matrix(matrixNN['KRT10',])), main = 'KRT10')
plot(density(as.matrix(matrixNN['KRT6B',])), main = 'KRT6B')
plot(density(as.matrix(matrixNN['KRT1',])), main = 'KRT1')
plot(density(as.matrix(matrixNN['KRT14',])), main = 'KRT14')
par(mfrow=c(2,2))
plot(density(as.matrix(matrixNN['KRT14',])), main = 'KRT14')
plot(density(as.matrix(matrixNN['MT2A',])), main = 'MT2A')
plot(density(as.matrix(matrixNN['ASS1',])), main = 'ASS1')
plot(density(as.matrix(matrixNN['WNT10A',])), main = 'WNT10A')

# Normalized
par(mfrow=c(3,2))
matrixNN = normalized_counts
matrixNN = cpm
matrixNN = logCPM
plot(density(as.matrix(matrixNN['FLG',])), main = 'FLG')
plot(density(as.matrix(matrixNN['LOR',])), main = 'LOR')
plot(density(as.matrix(matrixNN['FLG2',])), main = 'FLG2')
plot(density(as.matrix(matrixNN['CDSN',])), main = 'CDSN')
plot(density(as.matrix(matrixNN['ARG1',])), main = 'ARG1')
plot(density(as.matrix(matrixNN['SLURP1',])), main = 'SLURP1')
par(mfrow=c(3,2))
plot(density(as.matrix(matrixNN['PTGS1',])), main = 'PTGS1')
plot(density(as.matrix(matrixNN['KRT10',])), main = 'KRT10')
plot(density(as.matrix(matrixNN['KRT6B',])), main = 'KRT6B')
plot(density(as.matrix(matrixNN['KRT1',])), main = 'KRT1')
plot(density(as.matrix(matrixNN['KRT14',])), main = 'KRT14')
par(mfrow=c(2,2))
plot(density(as.matrix(matrixNN['KRT14',])), main = 'KRT14')
plot(density(as.matrix(matrixNN['MT2A',])), main = 'MT2A')
plot(density(as.matrix(matrixNN['ASS1',])), main = 'ASS1')
plot(density(as.matrix(matrixNN['WNT10A',])), main = 'WNT10A')


# check how many FLG +ve cells in scRNA-seq by checking if they have >1 or 2 scaled values

matrixPN = read.csv('sampledrawMatrix.csv',header = T, row.names = 1)
sum(as.numeric(matrixPN['KRT14',] > 1))
length(matrixPN['FLG',])


# gene markers
df = read.csv('markers.csv')
markers = as.vector(as.matrix(df))
markers = na.omit(markers)
NNmarkers = matrixNN[markers,]
PNmarkers = matrixPN[markers,]
NNmarkers = na.omit(NNmarkers)
PNmarkers = na.omit(PNmarkers)
