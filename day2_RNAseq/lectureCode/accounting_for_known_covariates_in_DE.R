counts_file <- system.file('extdata/rna-seq/SRP021193.raw_counts.tsv', 
                           package = 'compGenomRData')
colData_file <- system.file('extdata/rna-seq/SRP021193.colData.tsv', 
                            package = 'compGenomRData')

counts <- read.table(counts_file)
colData <- read.table(colData_file, header = T, sep = '\t', 
                      stringsAsFactors = TRUE)

plot_pca <- function(M, colData, by = NULL) {
  selected <- names(sort(apply(M, 1, sd), decreasing = T)[1:500])
  M <- t(M[selected,])
  M <- log(M+1)
  pca <- prcomp(M)
  pca.df <- as.data.frame(pca$x)
  pca.df <- merge(pca.df, colData, by = 'row.names')
  #compute variable explained by each principle component
  var_exp <- round(diag(cov(pca$x))/sum(diag(cov(pca$x))) * 100, 1)
  ggplot(pca.df, aes(x = PC1, y = PC2)) + 
    geom_point(aes_string(color = by), size = 5) +
    theme_bw(base_size = 14) + 
    labs(x = paste0('PC1 (',var_exp[['PC1']],'%)'), 
         y = paste0('PC2 (',var_exp[['PC2']],'%)')) +
    scale_color_brewer(type = 'qual', palette = 2)
}

plot_heatmap <- function(M, colData) {
  selected <- names(sort(apply(M, 1, sd), decreasing = T)[1:500])
  pheatmap(M[selected,], scale = 'row', show_rownames = F,
           annotation_col = colData)
}

library(DESeq2)
# remove the 'width' column from the counts matrix
countData <- as.matrix(subset(counts, select = c(-width)))
# set up a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, 
                              design = ~ 1)
dds <- estimateSizeFactors(dds)
counts_normalized <- counts(dds, normalized = TRUE)

plot_heatmap(counts_normalized, colData)
plot_pca(counts_normalized, colData, by = 'LibrarySelection')

dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, 
                              design = ~ LibrarySelection + group)
dds <- DESeq(dds)


