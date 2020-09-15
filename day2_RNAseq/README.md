# compgen2020

code and exercises for the computational genomics 2020 course

## Required packages

```
- Bioconductor and CRAN packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install('DESeq2', 'EDASeq', 'RUVSeq','corrplot', 
                      'gprofiler2', 'gage', 'ggfortify', 'ggplot2', 
                      'knitr', 'pheatmap', 'rmarkdown', 
                      'ExpressionAtlas', 'devtools')

- Data package from the book [at](http://compgenomr.github.io/book/)
devtools::install_github("compgenomr/compGenomRData")

```
