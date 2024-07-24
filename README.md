# Single-Cell-Proteomics-data-processing-and-analysis

               
               library("scp")
               
               library("ggplot2")

               library("dplyr")


# Read in SCP data

## The first step is to read in the PSM quantification table  

## a dummy MaxQuant dataset (Tyanova, Temu, and Cox (2016)).
               
               data("mqScpData")

               data("sampleAnnotation")

                table(sampleAnnotation$SampleType)



table(sampleAnnotation$SampleType)

# Using readSCP, we combine both tables in a QFeatures object formatted as described above.


               scp <- readSCP(assayData = mqScpData,
               
               colData = sampleAnnotation,
               
               runCol = "Raw.file",
               
               removeEmptyCols = TRUE)

               scp


An instance of class QFeatures containing 4 assays:
 
 [1] 190222S_LCA9_X_FP94BM: SingleCellExperiment with 395 rows and 11 columns 
 
 [2] 190321S_LCA10_X_FP97_blank_01: SingleCellExperiment with 109 rows and 11 columns 
 
 [3] 190321S_LCA10_X_FP97AG: SingleCellExperiment with 487 rows and 11 columns 
 
 [4] 190914S_LCB3_X_16plex_Set_21: SingleCellExperiment with 370 rows and 16 columns 

                
               plot(scp)

  ![image](https://github.com/user-attachments/assets/ecfa0b0b-2f8c-486c-9e87-5091fd0b0268)

