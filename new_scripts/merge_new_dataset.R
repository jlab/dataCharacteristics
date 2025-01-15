new_csv_file <- "/home/ubuntu/dataCharacteristics/new_data/dataCharacteristics_jlab_microbiome.csv"
data_type_name <- "JLAB Microbiome"
previous_base <- "/home/ubuntu/dataCharacteristics/shinyApps/datasets_results_clean_renamed.csv"
new_data_base <- "/home/ubuntu/dataCharacteristics/new_data/datasets_results_appended.csv"

new_csv_file <- "/home/ubuntu/dataCharacteristics/new_data/dataCharacteristics_marbel.csv"
data_type_name <- "Marbel"
previous_base <- new_data_base
new_data_base <- "/home/ubuntu/dataCharacteristics/new_data/datasets_results_microbiome_marbel.csv"

new_csv_file <- "/home/ubuntu/dataCharacteristics/new_data/dataCharacteristics_metatranscriptomics.csv"
data_type_name <- "Metatranscriptomic"
previous_base <- new_data_base
new_data_base <- "/home/ubuntu/dataCharacteristics/new_data/datasets_results_microbiome_marbel_metatranscriptomics.csv"



data <- read.csv(new_csv_file, check.names = FALSE)

data$prctPC1.wNAs.log2 <- 100 * data$prctPC1.wNAs.log2
data$prctPC2.wNAs.log2 <- 100 * data$prctPC2.wNAs.log2



data$dataType <- data$dataTypeSubgroups <- data_type_name

cols_to_zero_point_zero_one <- c("Corr(Mean vs. % NA) (Samples) (p-Value)", 
                  "Corr(Mean vs. % NA) (Analytes) (p-Value)", 
                  "Bimodality of sample correlations", 
                  "prctnDistinctValues")
data[cols_to_zero_point_zero_one] <- 0.01

Oldnames <- c("datasetID", "dataType", "dataTypeSubgroups", "nSamples", 
              "nAnalytes", "minRowNaPercentage.log2", 
              "maxRowNaPercentage.log2", "minColNaPercentage.log2", 
              "maxColNaPercentage.log2", 
              "percNATotal.log2", "percOfRowsWithNAs.log2", "percOfColsWithNAs.log2", 
              "corSampleMeanNA.log2", 
              # "corSampleMeanNAPval", 
              "corAnalyteMeanNA.log2", 
              # "corAnalyteMeanNAPval", 
              "mean.wNAs.log2", "median.wNAs.log2", "min.wNAs.log2", "max.wNAs.log2", "medianSampleVariance.wNAs.log2", 
              "medianAnalyteVariance.wNAs.log2", 
              "variance.wNAs.log2", "kurtosis.wNAs.log2", "skewness.wNAs.log2", "prctPC1.wNAs.log2", "prctPC2.wNAs.log2", 
              # "bimodalityColCorr", 
              "linearCoefPoly2Row.wNAs.log2", "quadraticCoefPoly2Row.wNAs.log2", "coefHclustRows.wNAs.log2", 
              "intensityNAProb50.sd.wNAs.log2", "intensityNAProb90.sd.wNAs.log2", 
              "intensityNAProbnSamplesWithProbValue.wNAs.log2", 
              "prctnDistinctValues")

Newnames <- c("Dataset ID", "Data type", "Data type subgroups", 
              "# Samples", "# Analytes", "min(% NA in analytes)", 
              "max(% NA in analytes)", "min(% NA in samples)", 
              "max(% NA in samples)", 
              "% NA", "% Analytes with NAs", "% Samples with NAs", 
              "Corr(Mean vs. % NA) (Samples)", 
              # "Corr(Mean vs. % NA) (Samples) (p-Value)", 
              "Corr(Mean vs. % NA) (Analytes)", 
              # "Corr(Mean vs. % NA) (Analytes) (p-Value)", 
              "Mean", "Median", "Min", "Max", "median(Variance of samples)", 
              "median(Variance of analytes)", 
              "Variance", "Kurtosis", "Skewness", "% Var. explained by PC1",
              "% Var. explained by PC2", 
              # "Bimodality of sample correlations", 
              "Lin. coef. of Poly2(Means vs. Vars) (Analytes)", 
              "Quadr. coef. of Poly2(Means vs. Vars) (Analytes)", 
              "Agglom. coef. hierarch. analyte clustering", 
              "sd(Intensity w/ prob(NA) = 50% for sample)", 
              "sd(Intensity w/ prob(NA) = 90% for sample)", 
              "# Samples w/ intensityNAProb50.sd or intensityNAProb90.sd", 
              "% Distinct values")

renameTable <- data.frame(Oldnames = Oldnames,
                          Newnames = Newnames)

data <- data %>% dplyr::rename_with(~ Newnames[which(Oldnames == .x)], 
                                    .cols = Oldnames)
data <- data %>% dplyr::mutate(`|Skewness|` = abs(Skewness))

base_data <- read.csv(previous_base, check.names = FALSE)

filtered_data <- data[colnames(data) %in% colnames(base_data)]

data.allDatasets <- rbind(base_data, filtered_data)

write.csv(data.allDatasets, new_data_base, row.names = FALSE)
