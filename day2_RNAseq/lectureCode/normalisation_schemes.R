counts_file <- system.file("extdata/rna-seq/SRP029880.raw_counts.tsv",
                           package = "compGenomRData")
coldata_file <- system.file("extdata/rna-seq/SRP029880.colData.tsv",
                            package = "compGenomRData")

counts <- as.matrix(read.table(counts_file, header = T, sep = '\t'))
colData <- read.table(coldata_file, header = T, sep = '\t', 
                      stringsAsFactors = TRUE)

boxplot(log10(counts[,c(1:3,8:10)]+1))

# Computing CPM (counts per million reads)
cpm <- apply(subset(counts, select = c(-width)), 2, 
             function(x) x/sum(x) * 10^6)

## Computing RPKM  
# create a vector of gene lengths 
geneLengths <- as.vector(subset(counts, select = c(width)))

# compute rpkm (reads per kilobase per million)
rpkm <- apply(X = subset(counts, select = c(-width)),
              MARGIN = 2, 
              FUN = function(x) {
                x <- x / sum(x) * 10^6
                x <- x / geneLengths * 10^3
              })
colSums(rpkm)


## Computing TPM (transcripts per million)
tpm <- apply(X = subset(counts, select = c(-width)),
              MARGIN = 2, 
              FUN = function(x) {
                x <- x / geneLengths * 10^3
                x <- x / sum(x) * 10^6
              })
colSums(tpm)


## DESeq2 normalisation

dds <- DESeq2::DESeqDataSetFromMatrix(countData = subset(counts, select = -width),
                                      colData = colData,
                                      design = ~1)
dds <- DESeq2::estimateSizeFactors(dds)
counts_deseq <- counts(dds, normalized = TRUE)

# statquest video for deseq normalisation: https://www.youtube.com/watch?v=UFB993xufUU


