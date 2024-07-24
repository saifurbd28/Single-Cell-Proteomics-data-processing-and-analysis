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

  # Clean missing data

               scp <- zeroIsNA(scp, i = 1:4)


 # Filtering Option-1: Filter PSMs - A common steps in SCP is to filter out low-confidence PSMs

               
               rowDataNames(scp)
               
CharacterList of length 4

[["190222S_LCA9_X_FP94BM"]] uid Sequence Length Modifications ... exclude residual participated

[["190321S_LCA10_X_FP97_blank_01"]] uid Sequence Length Modifications ... exclude residual participated

[["190321S_LCA10_X_FP97AG"]] uid Sequence Length Modifications ... exclude residual participated

[["190914S_LCB3_X_16plex_Set_21"]] uid Sequence Length Modifications ... exclude residual participated

# Filter features based on feature annotations:

## 1. Remove PSMs that are matched to contaminants

## 2. Remove PSMs that are matched to the decoy database

## 3. Keep PSMs that exhibit a high PIF (parental ion fraction), indicative of the purity of a spectrum


               
               scp <- filterFeatures(scp,
                      ~ Reverse != "+" &
                        Potential.contaminant != "+" &
                        !is.na(PIF) & PIF > 0.8)

# Filtering option-2: Filter assays based on detected features


                      dims(scp)
                      keepAssay <- dims(scp)[1, ] > 150 # assays that have sufficient PSMs (> 150 rows)
                      scp <- scp[, , keepAssay] # subset the scp object for the assays that meet the criterion

  # Filtering option-3: Filter features based on SCP metrics


                      table(colData(scp)[, "SampleType"])
                      scp <- computeSCR(scp,
                                        i = 1:3,
                                        colvar = "SampleType",
                                        carrierPattern = "Carrier",
                                        samplePattern = "Macrophage|Monocyte",
                                        sampleFUN = "mean",
                                        rowDataName = "MeanSCR")

                    rbindRowData(scp, i = 1:3) |>
                      data.frame() |>
                      ggplot(aes(x = MeanSCR)) +
                      geom_histogram() +
                      geom_vline(xintercept = c(1/200, 0.1),
                                 lty = c(2, 1)) +
                      scale_x_log10()


                  scp <- filterFeatures(scp,
                                        ~ !is.na(MeanSCR) &
                                        MeanSCR < 0.1)
  ![image](https://github.com/user-attachments/assets/fa46b8de-b74f-4738-b874-43d93a033c30)

