setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(plyr)
library(dplyr)
library(ggplot2)
library(GGally)
library(ggforce)
library(stringr)
theme_set(theme_bw())

dataset <- ldply(list.files("20230614_results", pattern = ".csv", full.names = TRUE), read.csv, header=TRUE)
#dataset <- plyr::ldply(list.files(pattern = ".csv", full.names = TRUE), read.csv, header=TRUE)
dataset <- dataset[dataset$nSamples>0,]
dataset <- dataset[dataset$variance>0,]
names(dataset) <- gsub(x = names(dataset), pattern = "\\.log2|\\.wNAs", replacement = "")  
#names(dataset)[names(dataset) == 'corColR'] <- 'corSampleMeanNA'
#names(dataset)[names(dataset) == 'corRowR'] <- 'corAnalyteMeanNA'

dataset <- dataset %>% 
  mutate(
    # dataType2 = case_when(grepl("\\^iBAQ", datasetID) ~ paste0(dataType, "_iBAQ"),
    #                            grepl("\\^LFQ", datasetID) ~ paste0(dataType, "_LFQ"),
    #                            grepl("\\^Intensity", datasetID) ~ paste0(dataType, "_Intensity"),
    #                            grepl("LC", datasetID, ignore.case = TRUE) & dataType == "metabolomics_MS" ~ paste0(dataType, "_LC"),
    #                            grepl("GC", datasetID, ignore.case = TRUE) & dataType == "metabolomics_MS" ~ paste0(dataType, "_GC"),
    #                            TRUE ~ dataType),
         # dataType3 = case_when(grepl("LC", datasetID, ignore.case = TRUE) & 
         #                         grepl("pos", datasetID, ignore.case = TRUE) & 
         #                         dataType == "metabolomics_MS" ~ paste0(dataType, "_LC_pos"),
         #                       grepl("LC", datasetID, ignore.case = TRUE) & 
         #                         grepl("neg", datasetID, ignore.case = TRUE) & 
         #                         dataType == "metabolomics_MS" ~ paste0(dataType, "_LC_neg"),
         #                       grepl("GC", datasetID, ignore.case = TRUE) & 
         #                         grepl("pos", datasetID, ignore.case = TRUE) & 
         #                         dataType == "metabolomics_MS" ~ paste0(dataType, "_GC_pos"),
         #                       grepl("GC", datasetID, ignore.case = TRUE) & 
         #                         grepl("neg", datasetID, ignore.case = TRUE) & 
         #                         dataType == "metabolomics_MS" ~ paste0(dataType, "_GC_neg"),
         #                       TRUE ~ dataType4),
         # dataType4 = case_when(grepl("\\^iBAQ", datasetID) ~ paste0(dataType, "_iBAQ"),
         #                       grepl("\\^LFQ", datasetID) ~ paste0(dataType, "_LFQ"),
         #                       grepl("\\^Intensity", datasetID) ~ paste0(dataType, "_Intensity"),
         #                       TRUE ~ dataType),
         dataType = case_when(grepl("\\^iBAQ", datasetID) ~ paste0(dataType, "_iBAQ"),
                               grepl("\\^LFQ", datasetID) ~ paste0(dataType, "_LFQ"),
                               grepl("\\^Intensity", datasetID) ~ paste0(dataType, "_Intensity"),
                               dataType == "RNAseq_raw_undecorated" ~ "RNAseq_raw",
                               TRUE ~ dataType) # ,
         # dataTypeMetabolomicsCombined = case_when(grepl("\\^iBAQ", datasetID) ~ paste0(dataType, "_iBAQ"),
         #                                          grepl("\\^LFQ", datasetID) ~ paste0(dataType, "_LFQ"),
         #                                          grepl("\\^Intensity", datasetID) ~ paste0(dataType, "_Intensity"),
         #                                          grepl("metabolomics", dataType, ignore.case = TRUE) ~ "metabolomics",
         #                                          TRUE ~ dataType)
         )


# Remove outliers: 
# E-GEOD-152766.aggregated_filtered_counts.mtx, (also E-GEOD-152766.aggregated_filtered_normalised_counts.mtx to make
# numbers of normalized and unnormalized sc datasets even.)
# MTBLS407_m_MTBLS407_21022_metabolite_profiling_mass_spectrometry_v2_maf.tsv,
# MTBLS407_m_MTBLS407_21548_metabolite_profiling_mass_spectrometry_v2_maf.tsv
data <- dataset %>% dplyr::filter(!(datasetID %in% c("E-GEOD-152766.aggregated_filtered_counts.mtx", 
                                                     "E-GEOD-152766.aggregated_filtered_normalised_counts.mtx",
                                                   "MTBLS407_m_MTBLS407_21022_metabolite_profiling_mass_spectrometry_v2_maf.tsv",
                                                   "MTBLS407_m_MTBLS407_21548_metabolite_profiling_mass_spectrometry_v2_maf.tsv")))


# selectedDataTypeCol <- "dataType"
# data <- dataset %>% dplyr::group_by(!!sym(selectedDataTypeCol)) %>% filter(n() > 5) %>% ungroup
data <- data %>% dplyr::group_by(dataType) %>% filter(n() > 5) %>% ungroup


write.csv(data, "datasets_results.csv", row.names = FALSE)
write.csv(data.frame(table(data[, "dataType"])), "numberOfDatasets.csv", row.names = FALSE)
######

# data <- dataset %>% dplyr::select(-nNegativeNumbers)
# Get %NA per column grouped by data type
naPrctPerCol <- data %>%
  dplyr::group_by(dataType) %>%
  summarise_all(~sum(is.na(.))/n()*100)

# Remove "bimodalityRowCorr" because of too many NAs
# data <- dataset %>% dplyr::select(-c(nNegativeNumbers, bimodalityRowCorr, bimodalityColCorr))
data$bimodalityRowCorr <- NULL

data <- data %>%
  mutate_if(is.integer, as.numeric) %>% 
  mutate("prctnDistinctValues" = nDistinctValues/((nSamples * nAnalytes) * ((100 - percNATotal)/100)) * 100 ) %>%
  dplyr::select(-c(nDistinctValues, nNegativeNumbers))



seedCols <- c("bimodalityRowCorrSeed", "bimodalityColCorrSeed", 
                                                          "coefHclustRowsSeed", "intensityNAProbSeed"#,
                                                          #"minRowNaPercentage", "maxRowNaPercentage", 
                                                          #"minColNaPercentage", "maxColNaPercentage", 
                                                          #"percOfRowsWithNAs", "percOfColsWithNAs"
)

data <- data[,!(colnames(data) %in% seedCols)]

# data <- data %>%
#   mutate_if(is.integer, as.numeric) %>% 
#   mutate("prctminRowNonNaNumber" = (nSamples - minRowNonNaNumber)/nSamples * 100 ) %>%
#   mutate("prctmaxRowNonNaNumber" = (nSamples - maxRowNonNaNumber)/nSamples * 100 ) 
# # prctminRowNonNaNumber is the same as maxRowNaPercentage, prctmaxRowNonNaNumber is the same as minRowNaPercentage
data <- data[,!(colnames(data) %in% c("minRowNonNaNumber", "maxRowNonNaNumber"))] 

write.csv(data, "datasets_results_clean.csv", row.names = FALSE)
#############################################

# Rename Data types
OldDataTypeNames <- c("metabolomics_MS", "metabolomics_NMR", "microarray", "microbiome", 
                      "proteomics_expressionatlas_iBAQ", "proteomics_expressionatlas_Intensity", 
                      "proteomics_pride_LFQ", "proteomics_pride_Intensity", "proteomics_pride_iBAQ", 
                      "RNAseq_fpkms_median", "RNAseq_raw", "RNAseq_tpms_median", "sc_normalized", 
                      "sc_unnormalized", "scProteomics")

NewDataTypeNames <- c("Metabolomics (MS)", "Metabolomics (NMR)", "Microarray", "Microbiome", 
              "Proteomics (iBAQ, Expression Atlas)", "Proteomics (Intensity, Expression Atlas)", 
              "Proteomics (LFQ, PRIDE)", "Proteomics (Intensity, PRIDE)", "Proteomics (iBAQ, PRIDE)", 
              "RNA-seq (FPKM)", "RNA-seq (raw)", "RNA-seq (TPM)", "scRNA-seq (normalized)", 
              "scRNA-seq (unnormalized)", "scProteomics")

renameDataTypeTable <- data.frame(OldDataTypeNames = OldDataTypeNames,
                                  NewDataTypeNames = NewDataTypeNames)

data$dataType <- renameDataTypeTable$NewDataTypeNames[
  match(data$dataType, renameDataTypeTable$OldDataTypeNames)]

# Rename Variables
Oldnames <- c("datasetID", "dataType", "nSamples", "nAnalytes", "minRowNaPercentage", 
              "maxRowNaPercentage", "minColNaPercentage", "maxColNaPercentage", 
              "percNATotal", "percOfRowsWithNAs", "percOfColsWithNAs", "corSampleMeanNA", 
              "corSampleMeanNAPval", "corAnalyteMeanNA", "corAnalyteMeanNAPval", 
              "mean", "median", "min", "max", "medianSampleVariance", "medianAnalyteVariance", 
              "variance", "kurtosis", "skewness", "prctPC1", "prctPC2", "bimodalityColCorr", 
              "linearCoefPoly2Row", "quadraticCoefPoly2Row", "coefHclustRows", 
              "intensityNAProb50.sd", "intensityNAProb90.sd", "intensityNAProbnSamplesWithProbValue", 
              "prctnDistinctValues")

Newnames <- c("Dataset ID", "Data type", "# Samples", "# Analytes", "min(% NA in analytes)", 
              "max(% NA in analytes)", "min(% NA in samples)", "max(% NA in samples)", 
              "% NA", "% Analytes with NAs", "% Samples with NAs", 
              "Corr(Mean vs. % NA) (Samples)", "Corr(Mean vs. % NA) (Samples) (p-Value)", 
              "Corr(Mean vs. % NA) (Analytes)", "Corr(Mean vs. % NA) (Analytes) (p-Value)", 
              "Mean", "Median", "Min", "Max", "median(Variance of samples)", "median(Variance of analytes)", 
              "Variance", "Kurtosis", "Skewness", "% Var. explained by PC1", "% Var. explained by PC2", 
              "Bimodality of sample correlations", 
              "Lin. coef. of Poly2(Means vs. Vars) (Analytes)", 
              "Quadr. coef. of Poly2(Means vs. Vars) (Analytes)", 
              "Agglom. coef. hierarch. analyte clustering", 
              "sd(Intensity w/ prob(NA) = 50% for sample)", "sd(Intensity w/ prob(NA) = 90% for sample)", 
              "# Samples w/ intensityNAProb50.sd or intensityNAProb90.sd", 
              "% Dinstinct values")

renameTable <- data.frame(Oldnames = Oldnames,
                          Newnames = Newnames)

data <- data %>% dplyr::rename_with(~ Newnames[which(Oldnames == .x)], .cols = Oldnames)


data$`Data type` <- factor(data$`Data type`, levels = c( "Metabolomics (NMR)", "Metabolomics (MS)", 
                                                   "Proteomics (iBAQ, PRIDE)", "Proteomics (Intensity, PRIDE)", "Proteomics (LFQ, PRIDE)", 
                                                   "Proteomics (iBAQ, Expression Atlas)", "Proteomics (Intensity, Expression Atlas)",
                                                  "scProteomics",
                                                  "Microarray", 
                                                  "RNA-seq (raw)", "RNA-seq (FPKM)", "RNA-seq (TPM)", 
                                                  "scRNA-seq (unnormalized)",
                                                  "scRNA-seq (normalized)", 
                                                  "Microbiome"))

#############################################
boxplotCols <- setdiff(unique(c("Dataset ID", "Data type", "# Samples", "# Analytes", "min(% NA in analytes)", 
                                  "max(% NA in analytes)", "min(% NA in samples)", "max(% NA in samples)", 
                                  "% NA", "% Analytes with NAs", "% Samples with NAs", "Mean", 
                                  "Median", "Min", "Max", "median(Variance of samples)", "median(Variance of analytes)", 
                                  "Variance", "Kurtosis", "Skewness", "% Var. explained by PC1", 
                                  "% Var. explained by PC2", "Lin. coef. of Poly2(Means vs. Vars) (Analytes)", 
                                  "Quadr. coef. of Poly2(Means vs. Vars) (Analytes)", "Agglom. coef. hierarch. analyte clustering", 
                                  "% Dinstinct values", "Corr(Mean vs. % NA) (Samples)", "Corr(Mean vs. % NA) (Analytes)",
                                "sd(Intensity w/ prob(NA) = 50% for sample)", 
                                "sd(Intensity w/ prob(NA) = 90% for sample)")), c("Data type", "Dataset ID"))

# boxplotCols <-  c("nSamples", "nAnalytes", 
#                    "medianSampleVariance", "medianAnalyteVariance",
#                    "median", "skewness", "prctPC1", "prctPC2",
#                    "corSampleMeanNA", "corAnalyteMeanNA", "percNATotal")


data2 <- data[, c("Data type", 
                  boxplotCols)]

# To be log2-transformed:
toBeLog2Transformed <- c("# Samples", "# Analytes", 
                         "min(% NA in analytes)", "max(% NA in analytes)", 
                         "% Analytes with NAs", "% Samples with NAs",
                         "median(Variance of samples)", "median(Variance of analytes)", 
                         "Variance", "Kurtosis", "Skewness",
                         "Lin. coef. of Poly2(Means vs. Vars) (Analytes)", 
                         "Quadr. coef. of Poly2(Means vs. Vars) (Analytes)",
                         "sd(Intensity w/ prob(NA) = 50% for sample)", 
                         "sd(Intensity w/ prob(NA) = 90% for sample)")


colsWithNegativeNumbers <- colnames(data2[, sapply(data2, FUN = function(x) any(x <= 0, na.rm = TRUE))])
# [1] "min(% NA in analytes)"                            "max(% NA in analytes)"                           
# [3] "min(% NA in samples)"                             "max(% NA in samples)"                            
# [5] "% NA"                                             "% Analytes with NAs"                             
# [7] "% Samples with NAs"                               "Mean"                                            
# [9] "Median"                                           "Min"                                             
# [11] "Max"                                              "median(Variance of samples)"                     
# [13] "median(Variance of analytes)"                     "Kurtosis"                                        
# [15] "Skewness"                                         "Lin. coef. of Poly2(Means vs. Vars) (Analytes)"  
# [17] "Quadr. coef. of Poly2(Means vs. Vars) (Analytes)" "Corr(Mean vs. % NA) (Samples)"                   
# [19] "Corr(Mean vs. % NA) (Analytes)"                   "sd(Intensity w/ prob(NA) = 50% for sample)"      
# [21] "sd(Intensity w/ prob(NA) = 90% for sample)"      

toBeLog2Transformed <- setdiff(toBeLog2Transformed, colsWithNegativeNumbers)
# "# Samples"  "# Analytes" "Variance"  

log2transform <- function(df, variable){
  df[[paste0("log2(", variable, ")")]] <- log2(df[[variable]])
  df[[variable]] <- NULL
  df
}

for (var in toBeLog2Transformed){
  # print(var)
  data2 <- log2transform(df = data2, variable = var)
}

write.csv(data2, "data2.csv", row.names = FALSE)
data2[sapply(data2, is.infinite)] <- NA

dput(colnames(data2))


naRelatedCols <- c("Corr(Mean vs. % NA) (Samples)", 
                   "Corr(Mean vs. % NA) (Samples) (p-Value)", 
                   "Corr(Mean vs. % NA) (Analytes)", 
                   "Corr(Mean vs. % NA) (Analytes) (p-Value)", 
                   "sd(Intensity w/ prob(NA) = 50% for sample)", 
                   "sd(Intensity w/ prob(NA) = 90% for sample)",
                   "# Samples w/ intensityNAProb50.sd or intensityNAProb90.sd",
                   "Bimodality of sample correlations"
)

data2.complete <- data2[,!(colnames(data2) %in% naRelatedCols)]

data2.complete <- data2.complete[complete.cases(data2.complete), ]

#############################################

data2.long <- reshape2::melt(data2)

neworder <- c("Data type", 
              "log2(# Analytes)", "log2(# Samples)", 
              "Mean", "Median", "Min", "Max", 
              "median(Variance of samples)", "median(Variance of analytes)", "log2(Variance)",
              "Kurtosis", "Skewness", "% Dinstinct values",
              "% NA",
              "min(% NA in samples)", "max(% NA in samples)",  
              "min(% NA in analytes)", "max(% NA in analytes)", 
              "% Analytes with NAs", "% Samples with NAs", 
              "sd(Intensity w/ prob(NA) = 50% for sample)", 
              "sd(Intensity w/ prob(NA) = 90% for sample)",
              "Corr(Mean vs. % NA) (Samples)", "Corr(Mean vs. % NA) (Analytes)", 
              "% Var. explained by PC1", "% Var. explained by PC2", 
              "Agglom. coef. hierarch. analyte clustering",
              "Lin. coef. of Poly2(Means vs. Vars) (Analytes)", 
              "Quadr. coef. of Poly2(Means vs. Vars) (Analytes)"
)

data2.long <- dplyr::arrange(dplyr::mutate(data2.long,
                                           variable=factor(variable,levels=neworder)), variable)


numberOfDatasetsIncluded <- data2.long %>% 
  na.omit() %>%
  group_by(`Data type`, variable) %>%  tally() %>%
  ungroup()
colnames(numberOfDatasetsIncluded) <- c("Data type", "variable", "count")
write.csv(numberOfDatasetsIncluded, "numberOfDatasetsIncluded.csv", row.names = FALSE)

plotBoxplots <- function(data2.long, fileNameAddition = "") {
  # For each data characteristic: Median of the median of each data type
  medianValues <- data2.long %>%
    group_by(`Data type`, variable) %>%  
    summarise(medianValue = median(value, na.rm=T)) %>%
    group_by(variable) %>%  
    summarise(medianValue = median(medianValue, na.rm=T))
  write.csv(medianValues, paste0("medianValues", fileNameAddition,".csv"), row.names = FALSE)
  
  ggplot.charact <- ggplot(data2.long, aes(forcats::fct_rev(`Data type`), value)) +
    # geom_boxplot(aes(fill = `Data type`), alpha=0.5, outlier.size=0.5) +
    geom_violin(aes(fill = `Data type`), alpha=0.5, scale = "width") +
    geom_boxplot(aes(fill = `Data type`), width=0.5, alpha=0.25, outlier.size=0.5) +
    coord_flip() +
    xlab("") +
    ylab("") +
    geom_hline(aes(yintercept = medianValue), medianValues, colour = 'red') +
    facet_wrap( ~ variable, scales = "free_x", ncol=7, strip.position = "bottom") +
    ggplot2::theme_bw() +
    #theme_minimal() + 
    theme(panel.spacing.y=unit(1.5, "lines"),
          legend.title = element_blank(), axis.text.y = element_text(hjust=0, face = "bold"), 
          # legend.position="bottom",
          legend.position = "none",
          legend.justification = "left", legend.direction = "horizontal",
          strip.text = element_text(face="bold", size = 6),
          strip.placement = "outside",                      # Place facet labels outside x axis labels.
          strip.background = element_blank(),  # Make facet label background white.
          axis.title = element_blank()) +     
    guides(fill = guide_legend(reverse = TRUE))
  ggsave(file=paste0("boxplots", fileNameAddition, ".pdf"), ggplot.charact, height=12, width=18)
}

plotBoxplots(data2.long, fileNameAddition = "_allDataTypes")
# plotBoxplots(data2.long %>% filter(dataType != "metabolomics_MS"), fileNameAddition = "_metabolomics_MSRemoved")

plotBoxplots(data2.long %>% na.omit() %>% dplyr::group_by(variable) %>% dplyr::mutate(value = rank(value)), 
             fileNameAddition = "_allDataTypes_ranks")

#############################################
colorForCorValues <- function(data, mapping, method="spearman", use="pairwise", ...){
  
  # grab data
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # calculate correlation
  corr <- cor(x, y, method=method, use=use)
  
  # calculate colour based on correlation value
  # Here I have set a correlation of minus one to blue, 
  # zero to white, and one to red 
  # Change this to suit: possibly extend to add as an argument of `my_fn`
  colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
  
  ggally_cor(data = data, mapping = mapping, 
             method = method,
             colour = "black", ...) + 
    theme_void() +
    theme(panel.background = element_rect(fill=fill))
}



plotPairsPlotForTypes <- function(df, groupColName = "Data type", colsForCorr = c(), corrMethod = "spearman", 
                                  smooth = FALSE, width = 19, height = 19, addStr = "") {
  for (group in unlist(unique(df[, groupColName]))){
    print(group)
    df.group <- df[df[, groupColName] == group,]
    df.group <- df.group[, colsForCorr]
    colnames(df.group) <- stringr::str_wrap(colnames(df.group), width = 15)
    pdf(paste0("ggpairs_", gsub(" ", "_", group), addStr, ".pdf"), width = width, height = height)
    if (smooth) {
      ggpairsObj <- ggpairs(df.group, 
                            lower = list(continuous = wrap("smooth", alpha = 0.1, size = 0.1)),
                            upper = list(continuous = wrap(colorForCorValues, method = corrMethod, size = 2.5),
                                         continuous = "smooth") #,
                            # mapping=aes(# color = get(groupColName),
                            #   # fill= get(groupColName),
                            #   alpha=0.3)
      )
    } else {
      ggpairsObj <- ggpairs(df.group, 
                            # lower = list(continuous = wrap("smooth", alpha = 0.1, size = 0.1)),
                            lower = list(continuous = wrap("points", alpha = 0.1, size = 0.1)),
                            upper = list(continuous = wrap(colorForCorValues, method = corrMethod, size = 2.5),
                                         continuous = "smooth") #,
                            # mapping=aes(# color = get(groupColName),
                            #   # fill= get(groupColName),
                            #   alpha=0.3)
      )
    }

    print(ggpairsObj + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
      ggtitle(group) +
      # theme_bw() +
      theme(strip.placement = "outside", text = element_text(size = 6))
    )
    
    dev.off()
  }
}

colsForCorr <- setdiff(colnames(data2), "Data type")
corrMethod <- "spearman"
plotPairsPlotForTypes(df = data2, groupColName = "Data type", colsForCorr = colsForCorr, 
                      corrMethod = corrMethod, addStr = paste0("_allcols_", corrMethod)) 

#############################################
# Interesting features based on correlation plot of all columns:
# "Lin. coef. of Poly2(Means vs. Vars) (Analytes)", 
# "Quadr. coef. of Poly2(Means vs. Vars) (Analytes)", 
# "sd(Intensity w/ prob(NA) = 50% for sample)", "sd(Intensity w/ prob(NA) = 90% for sample)", 
# "max(% NA in analytes)", "log2(# Samples)",  "% Dinstinct values",
# "Agglom. coef. hierarch. analyte clustering", "log2(# Analytes)",
# "Kurtosis", "Skewness",

# One of: "median(Variance of samples)", "log2(Variance)"
# One of: "Mean", "Median", "Min", "Max"
# One of : "min(% NA in samples)", "max(% NA in samples)", "% NA", "% Analytes with NAs", "% Samples with NAs", 

selForCorr <- c("% Dinstinct values", 
                "log2(# Analytes)", "log2(# Samples)", 
                "Mean", "log2(Variance)", 
                "% NA", "max(% NA in analytes)", "Corr(Mean vs. % NA) (Samples)", "Corr(Mean vs. % NA) (Analytes)",
                "Skewness", "Kurtosis", 
                "% Var. explained by PC1", 
                "% Var. explained by PC2",  
                "Lin. coef. of Poly2(Means vs. Vars) (Analytes)", "Quadr. coef. of Poly2(Means vs. Vars) (Analytes)", 
                "Agglom. coef. hierarch. analyte clustering", 
                "sd(Intensity w/ prob(NA) = 50% for sample)", "sd(Intensity w/ prob(NA) = 90% for sample)"
                )

# selForCorr <- selForCorr[order(match(selForCorr, colnames(data2)))]

corrMethod <- "spearman"
plotPairsPlotForTypes(df = data2, groupColName = "Data type", colsForCorr = selForCorr, 
                      corrMethod = corrMethod, width = 14, height = 14, addStr = paste0("_selectedCols_", corrMethod)) 


# c("Data type", "min(% NA in analytes)",
#   "median(Variance of analytes)",  "% Var. explained by PC1", 
#   "% Var. explained by PC2",  
#  "Corr(Mean vs. % NA) (Samples)", "Corr(Mean vs. % NA) (Analytes)", 
#   )
#############################################

# margPlot <- ggplot(data, aes(x = corSampleMeanNA, y = corAnalyteMeanNA, colour = get(selectedDataTypeCol))) +
# geom_point(aes(fill = get(selectedDataTypeCol)),  size = 0.8, alpha = 0.5) +
#   theme_minimal() +
#   theme(legend.title=element_blank())

margPlot <- ggplot(data, aes(x = `Corr(Mean vs. % NA) (Samples)`, 
                             y = `Corr(Mean vs. % NA) (Analytes)`, colour = `Data type`)) +
  geom_point(aes(fill = `Data type`),  size = 0.8, alpha = 0.5) +
  theme_minimal() +
  theme(legend.title=element_blank())

pdf("marginPlot.pdf", width = 12, height = 10)
ggExtra::ggMarginal(margPlot, groupFill = TRUE, groupColour = TRUE)
dev.off()

#############################################

plotPCABiplot <- function(df, groups= c(), alpha = 0.5, 
                             pcaMethod = "nipals",
                             coordRatio = 0.5, 
                             facetZoom = TRUE, 
                             xlimLower = NA, xlimUpper = NA,
                             ylimLower = NA, ylimUpper = NA,
                             PCchoices = 1:2) {
  # See https://stackoverflow.com/a/49788251
  #   # devtools::install_github("vqv/ggbiplot")

  library(ggbiplot)
  
  iris_dummy<-df
  iris_dummy[is.na(iris_dummy)]<-7777 #swap out your NAs with a dummy number so prcomp will run
  pca.obj <- prcomp(iris_dummy, center=TRUE, scale.=TRUE)
  
  # scale: One of "UV" (unit variance a=a/\sigma_{a}), 
  # "vector" (vector normalisation b=b/|b|), 
  # "pareto" (sqrt UV) or "none" 
  # to indicate which scaling should be used to scale the matrix with aa variables and b samples. 
  pca.obj2 <- pcaMethods::pca(df, method=pcaMethod, nPcs=4, center=TRUE
                              , scale = "uv"
  )
  
  pca.obj$x<-pca.obj2@scores 
  pca.obj$rotation<-pca.obj2@loadings 
  pca.obj$sdev<-pca.obj2@sDev
  pca.obj$center<-pca.obj2@center
  pca.obj$scale<-pca.obj2@scale
  
  P2 <- ggbiplot::ggbiplot(pca.obj,
                           choices = PCchoices,
                           obs.scale = 1, 
                           var.scale=1,
                           ellipse=T,
                           circle=F,
                           varname.size=3,
                           var.axes=T,
                           groups=groups, 
                           alpha=0)  +
    scale_color_discrete(name = '') +  
    coord_fixed(ratio = coordRatio)
  
  if (facetZoom) {
    P2 <- P2 + ggforce::facet_zoom(xlim = c(xlimLower, xlimUpper), ylim = c(ylimLower, ylimUpper))
  } else {
    if (!all(is.na(c(xlimLower, xlimUpper))))
      P2 <- P2 + xlim(c(xlimLower, xlimUpper))
    
    if (!all(is.na(c(ylimLower, ylimUpper))))
      P2 <- P2 + ylim(c(ylimLower, ylimUpper))
  }
  

  P2 <- P2 + theme(legend.direction ='horizontal', 
                   legend.position = 'bottom')
  
  P2$layers <- c(geom_point(aes(colour=groups), cex=1, alpha = alpha), P2$layers)
  
  P2
}

plotPCABiplots <- function(df, groupColName = "", addStr = "", pcaMethod = "nipals") {
    
  pdf(file = paste0("biplot_", pcaMethod, "_facetZoom_PC1vs2", addStr, ".pdf"), width = 12, height = 10)
  print(plotPCABiplot(df = df %>% dplyr::select(-!!groupColName), 
                         groups= df[[groupColName]],
                         alpha = 0.3,
                         pcaMethod = pcaMethod,
                         facetZoom = TRUE,
                         PCchoices = 1:2,
                         xlimLower = -4, xlimUpper = 5,
                         ylimLower = -5, ylimUpper = 8))
  dev.off()
  
  pdf(file = paste0("biplot_", pcaMethod, "_facetZoom_PC1vs3", addStr, ".pdf"), width = 12, height = 10)
  print(plotPCABiplot(df = df %>% dplyr::select(-!!groupColName), 
                         groups= df[[groupColName]],
                         alpha = 0.3,
                         pcaMethod = pcaMethod,
                         facetZoom = TRUE,
                         PCchoices = c(1, 3),
                         xlimLower = -5, xlimUpper = 5,
                         ylimLower = -4, ylimUpper = 8))
  dev.off()
  
  pdf(file = paste0("biplot_", pcaMethod, "_facetZoom_PC2vs3", addStr, ".pdf"), width = 12, height = 10)
  print(plotPCABiplot(df = df %>% dplyr::select(-!!groupColName), 
                         groups= df[[groupColName]],
                         alpha = 0.3,
                         pcaMethod = pcaMethod,
                         facetZoom = TRUE,
                         PCchoices = c(2, 3),
                         xlimLower = -5, xlimUpper = 7,
                         ylimLower = -4, ylimUpper = 8
                         ))
  dev.off()
  
  ## Remove outlier "E-GEOD-152766.aggregated_filtered_counts.mtx"
  # df <- df[row.names(df) != "E-GEOD-152766.aggregated_filtered_counts.mtx", ]
  pdf(file = paste0("biplot_", pcaMethod, "_", addStr, ".pdf"), width = 12, height = 10)
  print(plotPCABiplot(df = df %>% dplyr::select(-!!groupColName), 
                         groups= df[[groupColName]],
                         alpha = 0.3,
                         pcaMethod = pcaMethod,
                         facetZoom = FALSE))
  dev.off()
  
  pca <- pcaMethods::pca(df %>% dplyr::select(-!!groupColName), method=pcaMethod, nPcs=4, center=TRUE
                         , scale = "uv")
  write.csv(pca@loadings, paste0("loadings_",  pcaMethod, "_", addStr,".csv"))
  dat <- merge(pcaMethods::scores(pca), df, by=0)
  
  library(GGally)
  pdf(paste0("ggpairs_", pcaMethod, "_", addStr, ".pdf"), width = 12, height = 10)
  print(GGally::ggpairs(dat, columns = 2:5, ggplot2::aes(colour=get(groupColName)),
                lower = list(continuous = wrap("smooth", alpha = 0.3, size = 1), 
                             combo = wrap("dot_no_facet", alpha = 0.4)),
                upper=list(continuous = wrap("cor", method = "spearman", size = 3)),
                mapping=aes(color = get(groupColName),
                            fill= get(groupColName), 
                            alpha=0.5)) +
          ggtitle(groupColName) +
          theme_bw())
  dev.off()
  
  library(plotly)
  library(htmlwidgets)
  fig <- plot_ly(dat, x = ~PC1, y = ~PC2, z = ~PC3, 
                 color = ~as.factor(dat[[groupColName]]), 
                 type="scatter3d", mode="markers",
                 #colors = c('#636EFA','#EF553B') , 
                 marker = list(size = 2)
                 # , alpha = 0.75
  ) #%>%
  # add_markers(size = 5, marker=list(sizeref=8, sizemode="area"))
  fig <- fig %>%
    layout(
      title = "3D PCA",
      scene = list(bgcolor = "#e5ecf6"
      )
    )
  
  htmlwidgets::saveWidget(fig, paste0("plotly_", pcaMethod, "_", addStr,".html"), selfcontained = F, libdir = "lib")
}

plotPCABiplots(df = data2.complete, groupColName = "Data type", addStr = "", pcaMethod = "svd") # "svd" or "nipals"
plotPCABiplots(df = data2, groupColName = "Data type", addStr = "", pcaMethod = "nipals") # "svd" or "nipals"


# ####################
# y <- as.factor(data.complete$`Data type`)
# levels(y)
# x <- data.complete %>% select(-c("Data type", "Dataset ID")) 
# 
# for (var in toBeLog2Transformed){
#   # print(var)
#   x <- log2transform(df = x, variable = var)
# }
# 
# x <- x %>% as.matrix()
# # Fit parallel cumulative logit model
# fit1 <- ordinalNet::ordinalNet(x, y, family="cumulative", link="logit",
#                                parallelTerms=TRUE, nonparallelTerms=FALSE)
# summary(fit1)
# coefs <- data.frame(coef(fit1))
# coefsOver005 <- row.names(coefs)[which(abs(coefs$coef.fit1.) > 0.05)]
# 
# # colsForCorr <- c(coefsOver005[!grepl("(Intercept)", coefsOver005)],
# #                  "Corr(Mean vs. % NA) (Samples)",
# #                  "Corr(Mean vs. % NA) (Analytes)"
# #                  )

#####################


session <- sessionInfo()
sink(paste0("visualizations_sessionInfo.txt"))
print(session)
sink()
