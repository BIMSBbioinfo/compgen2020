library(DESeq2)
counts_file <- system.file("extdata/rna-seq/SRP029880.raw_counts.tsv",
                           package = "compGenomRData")
coldata_file <- system.file("extdata/rna-seq/SRP029880.colData.tsv",
                            package = "compGenomRData")
counts <- as.matrix(read.table(counts_file, header = T, sep = '\t'))
colData <- read.table(coldata_file, header = T, sep = '\t', stringsAsFactors = TRUE)
countData <- as.matrix(subset(counts, select = c(-width)))

dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, 
                              design = ~ group)
dds <- dds[ rowSums(DESeq2::counts(dds)) > 1, ]
dds <- DESeq(dds)
DEresults = results(dds, contrast = c("group", 'CASE', 'CTRL'))
DEresults <- DEresults[order(DEresults$pvalue),]

#remove genes with NA values 
DE <- DEresults[!is.na(DEresults$padj),]
#select genes with adjusted p-values below 0.1
DE <- DE[DE$padj < 0.1,]
#select genes with absolute log2 fold change above 1 (two-fold change)
DE <- DE[abs(DE$log2FoldChange) > 1,]
#calculate enriched GO terms
library(gprofiler2)
goResults <- gprofiler2::gost(query = rownames(DE)[1:500], 
                       organism = 'hsapiens', sources = 'GO:BP', evcodes = T)

p <- gprofiler2::gostplot(goResults, interactive = F, capped = F)

library(ggplot2)
library(ggrepel)
p + geom_text_repel(aes(label = ifelse(-log10(p_value) > 20, term_name, '')))

# GSEA analysis 
GO <- goResults$result
GO <- GO[order(GO$p_value),]
#restrict the terms that have at most 100 genes
GO <- GO[GO$term_size < 100,]
# use the top term from this table to create a gene set 
geneSet1 <- unlist(strsplit(GO[1,]$intersection, ','))

#Define another gene set by just randomly selecting 25 genes from the counts
#table get normalized counts from DESeq2 results
normalizedCounts <- DESeq2::counts(dds, normalized = TRUE)
geneSet2 <- sample(rownames(normalizedCounts), 25)

geneSets <- list('top_GO_term' = geneSet1,
                 'random_set' = geneSet2)

library(gage)
#use the normalized counts to carry out a GSEA. 
gseaResults <- gage(exprs = log2(normalizedCounts+1), 
                    ref = 6:10, 
                    samp = 1:5,
                    gsets = geneSets, compare = 'as.group')
library(pheatmap)
# get the expression data for the gene set of interest
M <- normalizedCounts[rownames(normalizedCounts) %in% geneSet1, ]
# log transform the counts for visualization scaling by row helps visualizing
# relative change of expression of a gene in multiple conditions
pheatmap(log2(M+1), 
         annotation_col = colData, 
         show_rownames = TRUE, 
         fontsize_row = 8,
         scale = 'row', 
         cutree_cols = 2)

hallmark_genesets <- do.call(c, lapply(readLines('../data/msigdb.hallmark_genesets.v7.1.gmt'), 
                            function(x) {
                              x <- unlist(strsplit(x, "\t"))
                              l <- list(x[3:length(x)])
                              names(l) <- x[1]
                              return(l)
                              }))

gsea_hallmarks <- gage(exprs = log2(normalizedCounts+1), 
                       ref = 6:10, samp = 1:5,  
                       gsets = hallmark_genesets, 
                       compare = 'as.group')

M <- normalizedCounts[rownames(normalizedCounts) %in% hallmark_genesets$HALLMARK_COAGULATION, ]
# log transform the counts for visualization scaling by row helps visualizing
# relative change of expression of a gene in multiple conditions
pheatmap(log2(M+1), 
         annotation_col = colData, 
         show_rownames = TRUE, 
         fontsize_row = 3,
         scale = 'row', 
         cutree_cols = 2)

df <- data.frame('case' = rowMeans(log(M[,1:5]+1)), 
                'ctrl' = rowMeans(log(M[,6:10]+1)))

df$log2FoldChange <- DEresults[rownames(df),]$log2FoldChange
df$gene <- rownames(df)

ggplot(df, aes(x = case, y = ctrl)) + geom_point(aes(color = log2FoldChange)) + 
  coord_fixed() +
  geom_abline(intercept = 1) + 
  theme_bw(base_size = 16) +
  geom_text_repel(aes(label = ifelse(log2FoldChange > 5, gene, '')))
  









