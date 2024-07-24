# Single-Cell-Proteomics-data-processing-and-analysis

library("scp")

library("ggplot2")

library("dplyr")

# Read in SCP data

## The first step is to read in the PSM quantification table  

## a dummy MaxQuant dataset (Tyanova, Temu, and Cox (2016)).

data("mqScpData")

# sample annotation

data("sampleAnnotation")

table(sampleAnnotation$SampleType)

# Using readSCP, we combine both tables in a QFeatures object formatted as described above.




               scp <- readSCP(assayData = mqScpData,
               colData = sampleAnnotation,
               
               runCol = "Raw.file",
               
               removeEmptyCols = TRUE)
