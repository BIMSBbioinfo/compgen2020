counts_file <- system.file('extdata/rna-seq/SRP049988.raw_counts.tsv', 
                           package = 'compGenomRData')
colData_file <- system.file('extdata/rna-seq/SRP049988.colData.tsv', 
                            package = 'compGenomRData')

counts <- read.table(counts_file)
colData <- read.table(colData_file, header = T, 
                      sep = '\t', stringsAsFactors = TRUE)
# simplify condition descriptions
colData$source_name <- ifelse(colData$group == 'CASE', 
                              'EHF_overexpression', 'Empty_Vector')

library(DESeq2)
# remove the 'width' column from the counts matrix
countData <- as.matrix(subset(counts, select = c(-width)))
# set up a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, 
                              design = ~ 1)
dds <- estimateSizeFactors(dds)
counts_normalized <- counts(dds, normalized = TRUE)

library(EDASeq)
# create a seqExpressionSet object using EDASeq package 
set <- newSeqExpressionSet(counts = countData,
                           phenoData = colData)
par(mfrow = c(1,2))
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=as.numeric(colData$group))
plotPCA(set, col = as.numeric(colData$group), adj = 0.5, 
        ylim = c(-0.7, 0.5), xlim = c(-0.5, 0.5))

par(mfrow = c(1,2))
plotRLE(counts_normalized, outline=FALSE, ylim=c(-4, 4), col=as.numeric(colData$group))
plotPCA(counts_normalized, col=as.numeric(colData$group), adj = 0.5, 
        ylim = c(-0.3, 1), xlim = c(-0.5, 0.5))


library(RUVSeq)

differences <- makeGroups(colData$group)
## looking for two different sources of unwanted variation (k = 2)
## use information from all genes in the expression object
par(mfrow = c(2, 2))
for(k in 1:4) {
  set_s <- RUVs(set, unique(rownames(set)), 
                k=k, differences) #all genes
  plotPCA(set_s, col=as.numeric(colData$group), 
          cex = 0.9, adj = 0.5, 
          main = paste0('with RUVs, k = ',k), 
          ylim = c(-1, 1), xlim = c(-0.6, 0.6))
}

# choose k = 2
set_s <- RUVs(set, unique(rownames(set)), k=2, differences) #


par(mfrow = c(1,2))
plotRLE(set, outline=FALSE, ylim=c(-4, 4), 
        col=as.numeric(colData$group), 
        main = 'without RUVs')
plotRLE(set_s, outline=FALSE, ylim=c(-4, 4), 
        col=as.numeric(colData$group),
        main = 'with RUVs')

par(mfrow = c(1,2))
plotPCA(set, col=as.numeric(colData$group), 
        main = 'without RUVs', adj = 0.5,
        ylim = c(-0.75, 0.75), xlim = c(-0.75, 0.75))
plotPCA(set_s, col=as.numeric(colData$group), 
        main = 'with RUVs', adj = 0.5,
        ylim = c(-0.75, 0.75), xlim = c(-0.75, 0.75))


library(EDASeq)
library(pheatmap)
# extract normalized counts that are cleared from unwanted variation using RUVs
normCountData <- normCounts(set_s)
plot_heatmap(normCountData, colData)


# How to integrate estimated sources of variation in DESeq2 analysis 

library(DESeq2)
colData_updated <- cbind(colData, 
                         pData(set_s)[,3:4])

#set up DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData_updated, 
                              design = ~ W_1 + W_2 + group)
# filter for low count genes
dds <- dds[rowSums(DESeq2::counts(dds)) > 10]
dds <- DESeq(dds)



