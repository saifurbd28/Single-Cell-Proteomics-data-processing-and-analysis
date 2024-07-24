# Single-Cell-Proteomics-data-processing-and-analysis

      if (!requireNamespace("BiocManager"))
      install.packages("BiocManager")
      BiocManager::install("scp")
      BiocManager::install(version = "devel")
      BiocManager::install("remotes")
      BiocManager::install("UCLouvain-CBIO/scp")


               
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

# Filter features to control for FDR
                  scp <- pep2qvalue(scp,
                                    i = 1:3,
                                    PEP = "dart_PEP",
                                    groupBy = "Leading.razor.protein",
                                    rowDataName = "qvalue_proteins")


                  scp <- filterFeatures(scp,
                                        ~ qvalue_proteins < 0.01)
# Process the PSM data

                  scp <- divideByReference(scp,
                                           i = 1:3,
                                           colvar = "SampleType",
                                           samplePattern = ".",
                                           refPattern = "Reference")
# Aggregate PSM data to peptide data
                scp

                scp <- aggregateFeaturesOverAssays(scp,
                                   i = 1:3,
                                   fcol = "Modified.sequence",
                                   name = paste0("peptides_", names(scp)),
                                   fun = matrixStats::colMedians, na.rm = TRUE)

                
                scp

                plot(scp)

![image](https://github.com/user-attachments/assets/402536c3-98b9-4de7-9c66-162bd378d5f0)

# Join the SCoPE2 sets in one assay
              scp <- joinAssays(scp,
                                i = 4:6,
                                name = "peptides")

              plot(scp)
![image](https://github.com/user-attachments/assets/ae03cf75-1d8d-4d5b-91f7-235d580b681a)

# Filter single-cells 
            scp
            scp <- scp[, scp$SampleType %in% c("Blank", "Macrophage", "Monocyte"), ] 
### subsetting removes unwanted samples from all assays
            scp

# Filter based on the median relative intensity

          medians <- colMedians(assay(scp[["peptides"]]), na.rm = TRUE)
          scp$MedianRI <- medians

          colData(scp) |>
            data.frame() |>
            ggplot() +
            aes(x = MedianRI,
                y = SampleType,
                fill = SampleType) +
            geom_boxplot() +
            scale_x_log10()
![image](https://github.com/user-attachments/assets/abf348de-cb3b-424e-ae0f-01b5f9eb6c74)

# Process the peptide data

## Normalization

### Divide columns by median
            scp <- sweep(scp,
                         i = "peptides",
                         MARGIN = 2,
                         FUN = "/",
                         STATS = colMedians(assay(scp[["peptides"]]), na.rm = TRUE),
                         name = "peptides_norm_col")

            plot(scp)

![image](https://github.com/user-attachments/assets/f86decc1-dbf2-4d63-bcd4-c639eb634318)

# Remove peptides with high missing rate
          scp <- filterNA(scp,
                i = "peptides_norm",
                pNA = 0.99)
# Log-transformation
          scp <- logTransform(scp,
                    base = 2,
                    i = "peptides_norm",
                    name = "peptides_log")

# Aggregate peptide data to protein data
          scp <- aggregateFeatures(scp,
                         i = "peptides_log",
                         name = "proteins",
                         fcol = "Leading.razor.protein",
                         fun = matrixStats::colMedians, na.rm = TRUE)

          scp

# Process the protein data

## Normalization
        ## Center columns with median
        scp <- sweep(scp, i = "proteins",
             MARGIN = 2,
             FUN = "-",
             STATS = colMedians(assay(scp[["proteins"]]),
                                na.rm = TRUE),
             name = "proteins_norm_col")
        ## Center rows with mean
        scp <- sweep(scp, i = "proteins_norm_col",
             MARGIN = 1,
             FUN = "-",
             STATS = rowMeans(assay(scp[["proteins_norm_col"]]),
                              na.rm = TRUE),
             name = "proteins_norm")

# Imputation missing values
        scp[["proteins_norm"]] |>
            assay() |>
            is.na() |>
            mean()


        scp <- impute(scp,
              i = "proteins_norm",
              name = "proteins_imptd",
              method = "knn",
              k = 3, rowmax = 1, colmax= 1,
              maxp = Inf, rng.seed = 1234)

      scp[["proteins_imptd"]] |>
              assay() |>
              is.na() |>
              mean()

# Dimension reduction
      library(scater)

      scp[["proteins_batchC"]] <- runPCA(scp[["proteins_batchC"]],
                                   ncomponents = 5,
                                   ntop = Inf,
                                   scale = TRUE,
                                   exprs_values = 1,
                                   name = "PCA")


      plotReducedDim(scp[["proteins_batchC"]],
               dimred = "PCA",
               colour_by = "SampleType",
               point_alpha = 1)

![plot_PCA-1](https://github.com/user-attachments/assets/7ea03d7c-ab00-4697-ba3f-dd8c5c5c9a79)

                                           
