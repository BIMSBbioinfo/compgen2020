counts_file <- system.file("extdata/rna-seq/SRP029880.raw_counts.tsv",
                           package = "compGenomRData")
coldata_file <- system.file("extdata/rna-seq/SRP029880.colData.tsv",
                            package = "compGenomRData")

counts <- as.matrix(read.table(counts_file, header = T, sep = '\t'))
colData <- read.table(coldata_file, header = T, sep = '\t', 
                      stringsAsFactors = TRUE)

counts <- subset(counts, select = -width)

#compute the standard deviation of each gene across samples
V <- apply(counts, 1, sd)
#sort the results by std. deviation in decreasing order 
#and select the top 100 genes 
selectedGenes <- names(V[order(V, decreasing = T)][1:100])

require(pheatmap)
pheatmap(counts[selectedGenes,], show_rownames = FALSE)

pheatmap(counts[selectedGenes,], scale = 'row', show_rownames = FALSE)

pheatmap(counts[selectedGenes,], scale = 'row', 
         show_rownames = FALSE, 
         annotation_col = colData, 
         cutree_cols = 2, cutree_rows = 2)


library(ggplot2)
require(ggfortify)
#transpose the matrix 
M <- t(counts[selectedGenes,])
# transform the counts to log2 scale 
M <- log2(M + 1)
#compute PCA 
pcaResults <- prcomp(M)

#plot PCA results making use of ggplot2's autoplot function
#ggfortify is needed to let ggplot2 know about PCA data structure. 
autoplot(pcaResults, data = colData, colour = 'group') +
  theme_bw(base_size = 14)


## Correlation plots
library(corrplot)
correlationMatrix <- cor(counts)
corrplot(correlationMatrix, order = 'hclust', 
         addrect = 2, addCoef.col = 'white', 
         number.cex = 0.7) 

corrplot(correlationMatrix, order = 'hclust', 
         addrect = 2, addCoef.col = 'white', 
         number.cex = 0.5, type = 'upper')

# split the clusters into two based on the clustering similarity 
pheatmap(correlationMatrix,  
         annotation_col = colData, 
         cutree_cols = 2, cutree_rows = 2, display_numbers = T)


