library(ExpressionAtlas)

datasets <- ExpressionAtlas::searchAtlasExperiments(properties = 'lung cancer',
                                                    species = 'mouse')
download <- getAtlasData(experimentAccessions = 'E-GEOD-59831')

exp <- download$`E-GEOD-59831`$rnaseq

counts <- assay(exp, 'counts')

cd <- colData(exp)

