#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(shiny)

library(scales)
library(ggforce)

library(data.table)
library(pcaMethods)
library(Matrix)
library(dplyr)
library(ggplot2)
ggplot2::theme_set(ggplot2::theme_bw())
library(cluster)
library(amap)
library(mclust)
library(mlr3misc)
library(biglm)
library(reshape2)
library(speedglm)
library(plotly)
library(tibble)
# library(gridExtra)
library(shinydashboardPlus)
library(ggpubr)
library(umap)

# devtools::install_github("vqv/ggbiplot")
library(ggbiplot)

library(shinycssloaders)
options(spinner.color = "#0275D8", spinner.color.background = "#ffffff", 
        spinner.size = 2)
################################################################################

source("/home/ubuntu/dataCharacteristics/shinyApps/inialize_helper.R")

`%notin%` <- Negate(`%in%`)

customColors <- c("Metabolomics (NMR)" = "#4BADF1", 
                  "Metabolomics (MS)" = "#0033CC",
                  "Lipidomics (MS)" = "#000000",
                  
                  "Proteomics (LFQ, PRIDE)" = "#800000",
                  "Proteomics (Intensity, PRIDE)" =  "#DC143C",
                  "Proteomics (iBAQ, PRIDE)" = "#DE4B7E",
                  "Proteomics (Intensity, Expression Atlas)" = "#FF6600",
                  "Proteomics (iBAQ, Expression Atlas)" = "#FFC0CB",
                  
                  "RNA-seq (raw)" = "#2E5F72",
                  "RNA-seq (FPKM)" = "#467741",
                  "RNA-seq (TPM)" = "#32CD32",
                  "Microarray" = "#5BB3B1",
                  
                  "scProteomics" = "#806FC4",
                  "scRNA-seq (unnormalized)" = "#DDA0DD",
                  "scRNA-seq (normalized)" = "#BA30B5",
                  "Microbiome" = "#6911D3",
                  "JLAB Microbiome" = "#99C59E",
                  "Marbel" = "#f73939",
                  "Metatranscriptomic" = "#fff100"
)

removeEmptyRowsAndColumns <- function(mtx, zerosToNA = FALSE){
  if (zerosToNA) mtx[mtx == 0] <- NA
  mtx[mtx %in% c("NaN", "N/A", "<-->")] <- NA
  mtx <- subset(mtx, rowSums(is.na(mtx) | mtx == 0) != ncol(mtx))
  if (!is.vector(mtx)) mtx <- mtx[, colSums(is.na(mtx) | mtx == 0) != nrow(mtx)]
  mtx
}

################################################################################
# DATA CHARACTERISTICS

calc_kurtosis <- function(data){
  n <- length(data[!is.na(data)])
  data <- data - mean(data, na.rm = TRUE)
  r <- n * sum(data^4, na.rm = TRUE) / (sum(data^2, na.rm = TRUE)^2)
  return(r * (1 - 1/n)^2 - 3)
}

calc_variance <- function(data){
  # var(c(data), na.rm=TRUE) # was replaced because for single cell data large 
  # vectors can become a problem
  1/(sum(!is.na(data)) - 1)*sum((data - mean(data, na.rm = TRUE))^2, 
                                na.rm = TRUE) 
}

calc_sd <- function(data){
  sqrt(calc_variance(data))
}

calc_skewness <- function(data) {
  data <- data[!is.na(data)]
  # return(sum((data - mean(data))^3)/(length(data) * sd(data)^3))
  return(sum((data - mean(data))^3)/(length(data) * calc_sd(data)^3))
}

calc_RMS <- function(data) sqrt(mean(data^2, na.rm = TRUE))

applyFunctionWithSeed <- function(functionName, seed = 123,  ...){
  
  oldseed <- mlr3misc::get_seed()
  
  if (is.null(seed)) {
    seed <- mlr3misc::get_seed()
  }
  
  set.seed(seed)
  res <- functionName(...)
  
  .Random.seed <- oldseed
  
  return(list(res = res, seed = seed))
}

# row means 
get_rowMeans <- function(mtx, ...) return(apply(mtx, 1, mean, na.rm = TRUE))

# row standard deviations
get_rowSd <- function(mtx, ...) {
  return(apply(mtx, 1, calc_sd))
}

# Poly2 (features)
# linear coefficient of 2nd order polynomial fit with x = row means,
# y = row variances
get_LinearCoefPoly2XRowMeansYRowVars <- function(mtx, ...){
  return(unname(coef(biglm::biglm(y~x+I(x^2), data = data.frame(
    y = get_rowSd(mtx)^2,
    x = get_rowMeans(mtx))))[2]))
}

# quadratic coefficient of 2rd order polynomial fit with x = row means, 
# y = row variances
get_QuadraticCoefPoly2XRowMeansYRowVars <- function(mtx, ...){ 
  return(unname(coef(biglm::biglm(y~x+I(x^2), data = data.frame(
    y = get_rowSd(mtx)^2,
    x = get_rowMeans(mtx))))[3]))
}

# linear and quadratic coefficient of 2rd order polynomial fit with
# x = row means, y = row variances
get_CoefPoly2XRowMeansYRowVars <- function(mtx, ...){
  
  rowSd <- rowMean <- linearCoef <- quadraticCoef <- NA
  try({rowSd <- get_rowSd(mtx)^2})
  try({rowMean <- get_rowMeans(mtx)})
  
  try({ 
    if (!(length(rowSd) == 1 && is.na(rowSd)) & 
        !(length(rowMean) == 1 && is.na(rowMean))) {
      df <- data.frame(
        y = rowSd,
        x = rowMean)
      model <- biglm::biglm(y~x+I(x^2), data = df)
      coefs <- coef(model)
      linearCoef <- unname(coefs[2])
      quadraticCoef <- unname(coefs[3])
    }
  })
  return(list(linearCoef = linearCoef, quadraticCoef = quadraticCoef))
}

# Coef.hclust (features)
get_coefHclustRowsWithFewestNAs <- function(mtx, naToZero = FALSE, 
                                            nMax = 500, ...) {
  mtx <- mtx[order(rowSums(is.na(mtx)), -rowMeans(mtx, na.rm = TRUE)), ]
  
  mtx <- mtx[order(rowMeans(mtx, na.rm = TRUE)), ]
  mtx <- mtx[order(rowSums(is.na(mtx))), ]
  if (naToZero) mtx[is.na(mtx)] <- 0
  randomCols <- applyFunctionWithSeed(sample, x = 1:ncol(mtx), 
                                      size = min(nMax, ncol(mtx)))
  # max. 500 features and samples
  mtx <- mtx[1:min(500, nrow(mtx)), randomCols$res]
  mtx <- removeEmptyRowsAndColumns(mtx, zerosToNA = FALSE)
  return(list(res = cluster::coef.hclust(amap::hcluster(mtx)), 
              seed = randomCols$seed))
}

calculateIntensityNAProbability5090 <- function(mtx, nmaxSamples = 200) {
  if (nrow(mtx) > 1000)
    mtx <- mtx[sample(1:nrow(mtx), 1000),] # downsampling in order to decrease RAM usage
  
  seedUsed <- NA
  if (ncol(mtx) > nmaxSamples) { # random subset of features are selected
    randomCols <- applyFunctionWithSeed(sample, x = 1:ncol(mtx), 
                                        size = min(nmaxSamples, ncol(mtx)))
    seedUsed <- randomCols$seed
    mtx <- mtx[, randomCols$res]
    mtx <- removeEmptyRowsAndColumns(mtx, zerosToNA = TRUE)
  }
  
  data.long <- reshape2::melt(mtx)
  mtx <- NULL
  gc()
  colnames(data.long) <- c("Feature", "Sample", "Value")
  data.long$isNA <- as.integer(is.na(data.long$Value))
  
  # Group by mean using dplyr
  featureMean.df <- data.long %>% dplyr::group_by(Feature) %>% 
    dplyr::summarise(mean = mean(Value, na.rm = TRUE))
  
  imputed <- data.long$Value
  
  data.long <- data.long %>% dplyr::left_join(featureMean.df, by = 'Feature')
  imputed[is.na(imputed)] <- data.long$mean[is.na(imputed)]
  data.long$imputed <- imputed
  
  
  x5090.lst <- lapply(unique(data.long$Sample), function(sample){
    data.long.Sample <- data.long %>% dplyr::filter(Sample == sample)
    x5090 <- NA
    glmModelSample <- NA
    if (any(data.long.Sample$isNA == 1)) {
      glmModelSample <- speedglm::speedglm(isNA ~ imputed, 
                                           family = binomial(link = 'logit'), 
                                           data = data.long.Sample, 
                                           fitted = TRUE)
      
      if (unname(coefficients(glmModelSample)[2]) < 0) {
        intensities <- seq(min(data.long$imputed) - 20, max(data.long$imputed),
                           length.out = 400)
        x5090 <- approx(predict(glmModelSample, 
                                newdata = data.frame(imputed = intensities),
                                type = "response"), intensities, c(0.5, 0.9))$y
      }
    }
    x5090
  })
  
  data.long <- NULL
  gc()
  
  x5090.df <- data.frame(do.call(rbind, x5090.lst))
  colnames(x5090.df) <- c("IntensityNAp50", "IntensityNAp90")
  x5090.sds <- apply(x5090.df, 2, calc_sd)
  
  nSamplesWithProbValue <- max(c(sum(!is.na(x5090.df$IntensityNAp50)), 
                                 sum(!is.na(x5090.df$IntensityNAp90))))
  
  list(IntensityNAp50 = x5090.sds[["IntensityNAp50"]], 
       IntensityNAp90 = x5090.sds[["IntensityNAp90"]],
       nSamplesWithProbValue = nSamplesWithProbValue,
       seed = seedUsed
  )
}


getCharacteristicsHelper <- function(mtx){
  
  if (is.vector(mtx)) {
    nFeatures <- 1
  } else {
    nFeatures <- nrow(mtx) # number of proteins with no NAs
  }
  
  mean <- mean(mtx, na.rm = TRUE)
  median <- median(mtx, na.rm = TRUE)
  min <- min(mtx, na.rm = TRUE)
  max <- max(mtx, na.rm = TRUE)
  
  medianSampleVariance <- medianAnalyteVariance <- skewness <- 
    kurtosis <- variance <-
    prctPC1 <- prctPC2 <- 
    linearCoefPoly2Row <- quadraticCoefPoly2Row <- 
    coefHclustRows <- 
    intensityNAProb50.sd <- intensityNAProb90.sd <- 
    intensityNAProbnSamplesWithProbValue <- NA
  
  try({medianSampleVariance <- median(apply(mtx, 2, calc_variance), 
                                      na.rm = TRUE)})
  try({medianAnalyteVariance <- median(unname(apply(mtx, 1, calc_variance)), 
                                       na.rm = TRUE)})
  
  try({skewness <- calc_skewness(mtx)})
  try({kurtosis <- calc_kurtosis(mtx)})
  
  try({variance <- calc_variance(mtx)})
  
  # Poly2 (features)
  try({
    coefs <- get_CoefPoly2XRowMeansYRowVars(mtx)
    linearCoefPoly2Row <- coefs[["linearCoef"]]
    quadraticCoefPoly2Row <- coefs[["quadraticCoef"]]
  })
  
  # Coef.hclust (features)
  try({
    coefHclustRowsRes <- get_coefHclustRowsWithFewestNAs(mtx, naToZero = TRUE)
    coefHclustRows <- coefHclustRowsRes$res
  })
  
  try({
    x5090.sds <- calculateIntensityNAProbability5090(mtx)
    intensityNAProb50.sd <- x5090.sds[["IntensityNAp50"]]
    intensityNAProb90.sd <- x5090.sds[["IntensityNAp90"]]
    intensityNAProbnSamplesWithProbValue <- x5090.sds[["nSamplesWithProbValue"]]
  })
  
  mtx <- mtx %>% t()
  # Remove zero variance columns 
  mtx <- mtx[, which(apply(mtx, 2, calc_variance) != 0)] 
  
  if (!is.vector(mtx)) {
    try({
      if (nrow(mtx) > 1000) # samples
        mtx <- mtx[sample(1:nrow(mtx), 1000),] # downsampling in order to decrease RAM usage
      if (ncol(mtx) > 10000) # features
        mtx <- mtx[, sample(1:ncol(mtx), 10000)] # downsampling in order to decrease RAM usage
      
      pca <- pcaMethods::pca(mtx, method = "nipals", center = TRUE, 
                             maxSteps = 5000)
      prctPC1 <- pca@R2[1]
      prctPC2 <- pca@R2[2]
    })
  }
  
  mtx <- NULL
  gc()
  resultvec <- c(
    mean = mean,
    median = median,
    min = min,
    max = max,
    medianSampleVariance = medianSampleVariance,
    medianAnalyteVariance = medianAnalyteVariance,
    variance = variance, 
    kurtosis = kurtosis,
    skewness = skewness,
    prctPC1 = prctPC1, 
    prctPC2 = prctPC2,
    linearCoefPoly2Row = linearCoefPoly2Row,
    quadraticCoefPoly2Row = quadraticCoefPoly2Row,
    coefHclustRows = coefHclustRows,
    intensityNAProb50.sd = intensityNAProb50.sd,
    intensityNAProb90.sd = intensityNAProb90.sd,
    intensityNAProbnSamplesWithProbValue = intensityNAProbnSamplesWithProbValue
  )
  
  return(resultvec)
}

getNaFeatures <- function(mtx) {
  colNaPercentage <- colMeans(is.na(mtx))*100
  rowNaPercentage <- rowMeans(is.na(mtx))*100
  
  colMeans2 <- colMeans(mtx, na.rm = TRUE)
  rowMeans2 <- rowMeans(mtx, na.rm = TRUE)
  
  corSampleMeanNA <- corAnalyteMeanNA <- NA
  corCoefType <- "spearman"
  try({
    cortestCol <- cor.test(colNaPercentage, colMeans2, method = corCoefType)
    corSampleMeanNA <- unname(cortestCol$estimate)
  })
  
  try({
    cortestRow <- cor.test(rowNaPercentage, rowMeans2, method = corCoefType)
    corAnalyteMeanNA <- unname(cortestRow$estimate)
  })
  
  percNATotal <- mean(is.na(mtx)) * 100
  percOfRowsWithNAs <- sum(apply(mtx, 1, anyNA))/nrow(mtx) * 100
  percOfColsWithNAs <- sum(apply(mtx, 2, anyNA))/ncol(mtx) * 100
  
  # mtx <- NULL
  # gc()
  c(
    minRowNaPercentage = min(rowNaPercentage),            
    maxRowNaPercentage = max(rowNaPercentage),
    minColNaPercentage = min(colNaPercentage),
    maxColNaPercentage = max(colNaPercentage),
    percNATotal = percNATotal,
    percOfRowsWithNAs = percOfRowsWithNAs,
    percOfColsWithNAs = percOfColsWithNAs,
    corSampleMeanNA = corSampleMeanNA,
    corAnalyteMeanNA = corAnalyteMeanNA
  )
}


getDataCharacteristics <- function(mtx) {
  
  mtx[mtx == 0] <- NA
  mtx[mtx == Inf] <- NA
  mtx <- mtx[, colSums(is.na(mtx)) != nrow(mtx)]
  
  nSamples <- ncol(mtx)
  nAnalytes <- nrow(mtx)
  
  nNegativeNumbers <- sum(mtx < 0, na.rm = TRUE)
  mtx <- log2(mtx)
  mtx[mtx == "NaN"] <- NA
  
  nDistinctValues <- length(unique(c(mtx[!is.na(mtx)])))
  naFeatures <- getNaFeatures(mtx)
  gc()
  
  characts <- getCharacteristicsHelper(mtx)
  mtx <- NULL
  gc()
  charact.log <- c(naFeatures, characts, nDistinctValues = nDistinctValues, 
                   nNegativeNumbers = nNegativeNumbers)
  
  res <- c(nSamples = nSamples, 
           nAnalytes = nAnalytes,
           # charact.noLog, 
           charact.log)
  res
  
}

logTransform <- function(df, variable, logBase = c("log2", "log1p")){
  df[[paste0(logBase, "(", variable, ")")]] <- get(logBase)(df[[variable]])
  df[[variable]] <- NULL
  df
}

plotPCABiplotNewDataset <- function(df, groups= c(), alpha = 0.5, 
                                    pcaMethod = "nipals",
                                    coordRatio = 0.5, 
                                    facetZoom = TRUE, 
                                    xlimLower = NA, xlimUpper = NA,
                                    ylimLower = NA, ylimUpper = NA,
                                    PCchoices = 1:2,
                                    ellipse = TRUE) {
  # See https://stackoverflow.com/a/49788251
  
  iris_dummy <- df
  #swap out your NAs with a dummy number so prcomp will run
  iris_dummy[is.na(iris_dummy)] <- 7777 
  pca.obj <- prcomp(iris_dummy, center = TRUE, scale.=TRUE)
  
  iris_dummy <- NULL
  gc()
  pca.obj2 <- pcaMethods::pca(df, method = pcaMethod, nPcs = 4, center = TRUE, 
                              scale = "uv")
  df <- NULL
  gc()
  pca.obj$x <- pca.obj2@scores 
  pca.obj$rotation <- pca.obj2@loadings 
  pca.obj$sdev <- pca.obj2@sDev
  pca.obj$center <- pca.obj2@center
  pca.obj$scale <- pca.obj2@scale
  
  pca.obj2 <- NULL
  gc()
  cond <- pca.obj$x[which(groups == "newDataset"),]
  
  P2 <- ggbiplot::ggbiplot(pca.obj,
                           choices = PCchoices,
                           obs.scale = 1, 
                           var.scale = 1,
                           ellipse = ellipse,
                           circle = F,
                           varname.size = 3,
                           var.axes = T,
                           groups = groups, 
                           alpha = 0)  +
    coord_fixed(ratio = coordRatio)
  
  pca.obj <- NULL
  gc()
  if (facetZoom) {
    P2 <- P2 + ggforce::facet_zoom(xlim = c(xlimLower, xlimUpper), 
                                   ylim = c(ylimLower, ylimUpper))
  } else {
    if (!all(is.na(c(xlimLower, xlimUpper))))
      P2 <- P2 + xlim(c(xlimLower, xlimUpper))
    
    if (!all(is.na(c(ylimLower, ylimUpper))))
      P2 <- P2 + ylim(c(ylimLower, ylimUpper))
  }
  
  P2 <- P2 + theme(legend.direction = 'horizontal', 
                   legend.position = 'bottom')
  
  P2$layers <- c(geom_point(aes(colour = groups), cex = 1, alpha = alpha), 
                 P2$layers)
  
  my_colors <- scales::hue_pal()(length(unique(groups)) - 1)
  
  P2 <- P2 + 
    scale_color_manual(name = '',
                       values = my_colors,
                       limits = setdiff(unique(groups), "newDataset")
    ) +
    geom_point(aes(x = cond[PCchoices[1]], y = cond[PCchoices[2]]), 
               col = "red", size = 3) 
  
  P2
}

plot3DPCA <- function(df, groupColName = "", addStr = "", 
                      pcaMethod = "nipals") {
  
  pca <- pcaMethods::pca(df %>% dplyr::select(-!!groupColName), 
                         method = pcaMethod, nPcs = 4, center = TRUE
                         , scale = "uv")
  dat <- merge(pcaMethods::scores(pca), df, by = 0)
  df <- NULL
  gc()
  # library(htmlwidgets)
  dat.woNewDataset <- dat %>% dplyr::filter(Row.names != "newDataset")
  newDataset <- dat[dat$Row.names == "newDataset",]
  dat <- NULL
  gc()
  fig <- plotly::plot_ly(dat.woNewDataset, x = ~PC1, y = ~PC2, z = ~PC3,
                         color = ~as.factor(dat.woNewDataset[[groupColName]]),
                         colors = customColors,
                         type = "scatter3d", mode = "markers",
                         marker = list(size = 3, opacity = 0.7),
                         hovertext = paste("Dataset ID :", 
                                           dat.woNewDataset$Row.names)
                         # , alpha = 0.75
  ) 
  
  fig <- fig %>%
    plotly::layout(
      hoverlabel = list(namelength = -1),
      legend = list(itemsizing = "constant")
    )
  
  
  
  fig
}

integrateNewDataset <- function(mtx,
                                nNegativeNumbers) {
  # data <- getDataCharacteristics(mtx)
  
  ##############################
  # getDataCharacteristics()
  # mtx[mtx == 0] <- NA
  # mtx[mtx == Inf] <- NA
  # mtx <- mtx[, colSums(is.na(mtx)) != nrow(mtx)]
  
  nSamples <- ncol(mtx)
  nAnalytes <- nrow(mtx)
  
  # nNegativeNumbers <- sum(mtx < 0, na.rm = TRUE)
  # mtx <- log2(mtx)
  # mtx[mtx == "NaN"] <- NA
  
  nDistinctValues <- length(unique(c(mtx[!is.na(mtx)])))
  naFeatures <- getNaFeatures(mtx)
  gc()
  
  characts <- getCharacteristicsHelper(mtx)
  gc()
  charact.log <- c(naFeatures, characts, nDistinctValues = nDistinctValues, 
                   nNegativeNumbers = nNegativeNumbers)
  
  data <- c(nSamples = nSamples, 
           nAnalytes = nAnalytes,
           # charact.noLog, 
           charact.log)
  
  ##############################
  
  # mtx <- NULL
  # gc()
  data <- data.frame(as.list(data))
  
  data$prctPC1 <- 100 * data$prctPC1
  data$prctPC2 <- 100 * data$prctPC2
  
  data <- data %>%
    dplyr::mutate_if(is.integer, as.numeric) %>% 
    dplyr::mutate("prctnDistinctValues" = 
                    nDistinctValues/((nSamples * nAnalytes) * 
                                       ((100 - percNATotal)/100)) * 100 ) %>%
    dplyr::select(-c(nDistinctValues, nNegativeNumbers))
  
  data$datasetID <- data$dataType <- data$dataTypeSubgroups <- "newDataset"
  
  Oldnames <- c("datasetID", "dataType", "dataTypeSubgroups", "nSamples", 
                "nAnalytes", "minRowNaPercentage", 
                "maxRowNaPercentage", "minColNaPercentage", 
                "maxColNaPercentage", 
                "percNATotal", "percOfRowsWithNAs", "percOfColsWithNAs", 
                "corSampleMeanNA", 
                # "corSampleMeanNAPval", 
                "corAnalyteMeanNA", 
                # "corAnalyteMeanNAPval", 
                "mean", "median", "min", "max", "medianSampleVariance", 
                "medianAnalyteVariance", 
                "variance", "kurtosis", "skewness", "prctPC1", "prctPC2", 
                # "bimodalityColCorr", 
                "linearCoefPoly2Row", "quadraticCoefPoly2Row", "coefHclustRows", 
                "intensityNAProb50.sd", "intensityNAProb90.sd", 
                "intensityNAProbnSamplesWithProbValue", 
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
  
 # data.allDatasets <-  read.csv("datasets_results_clean_renamed.csv",
#                                check.names = FALSE)
  data.allDatasets <-  read.csv("/home/ubuntu/dataCharacteristics/new_data/datasets_results_microbiome_marbel_metatranscriptomics.csv", check.names = FALSE)
  data.allDatasets$`Corr(Mean vs. % NA) (Samples) (p-Value)` <- 
  data.allDatasets$`Corr(Mean vs. % NA) (Analytes) (p-Value)` <-
  data.allDatasets$`Bimodality of sample correlations` <- NULL

  
  
  levels <- c('Metabolomics (NMR)', 
              'Metabolomics (MS)',
              
              "Lipidomics (NMR)",
              "Lipidomics (MS)",
              
              'Proteomics (LFQ, PRIDE)',
              'Proteomics (Intensity, PRIDE)',
              'Proteomics (iBAQ, PRIDE)',
              'Proteomics (Intensity, Expression Atlas)',
              'Proteomics (iBAQ, Expression Atlas)', 
              
              'RNA-seq (raw)',
              'RNA-seq (FPKM)',
              'RNA-seq (TPM)',
              'Microarray',
              
              'scProteomics',
              'scRNA-seq (unnormalized)',
              'scRNA-seq (normalized)',
              'Microbiome',
              'JLAB Microbiome',
              'Marbel',
              'Metatranscriptomic'
  )
  
  data.allDatasets[, "Data type"] <- factor(data.allDatasets[, "Data type"], levels = levels)
  
  data.allDatasets <- rbind(data.allDatasets, data)
  
  data <- NULL
  gc()
  data.allDatasets
}

plotCorrelation <- function(mtx.corr, plotTitle = "", plotSubtitle = "", 
                            ySampleMeanMin = NULL, ySampleMeanMax = NULL) {
  colMeans2 <- colMeans(mtx.corr, na.rm = TRUE)
  gc()
  colNaPercentage <- colMeans(is.na(mtx.corr)) * 100
  
  #mtx.corr <- NULL
  #gc()
  corRes <- cor.test(colMeans2, colNaPercentage, method = "spearman")
  corrPlot <- ggplot(data.frame(colNaPercentage, colMeans2), 
                     aes(colNaPercentage, colMeans2)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm") +
    theme_bw() +
    labs(x = "%NA in sample", y = "Sample mean") +
    ggtitle(plotTitle, subtitle = plotSubtitle) +
    xlim(0, 100) 
  
  if (!is.null(ySampleMeanMin) & !is.null(ySampleMeanMax))
    corrPlot <- corrPlot + ylim(c(ySampleMeanMin, ySampleMeanMax))
  
  if (!is.na(corRes$p.value) & corRes$p.value < 0.05) 
    corrPlot <- corrPlot + ggpubr::stat_cor(method = "spearman", label.x = 3)
  
  corrPlot
}


getUMAPNewDataset <- function(df, groupColName = "") {
  
  set.seed(142)
  groupVec <- df[[groupColName]]
  rownames.df <- row.names(df)
  
  umap_fit <- df %>% dplyr::mutate(ID = row_number()) %>% 
    dplyr::select(-!!groupColName) %>% dplyr::select_if(~ !any(is.na(.))) %>%
    remove_rownames() %>% column_to_rownames("ID") %>%
    scale() %>%
    umap::umap(n_components = 3)
  
  df <- NULL
  gc()
  
  umap_df <- umap_fit$layout %>%
    as.data.frame() %>%
    dplyr::rename(UMAP1 = "V1",
                  UMAP2 = "V2",
                  UMAP3 = "V3") %>%
    dplyr::mutate(!!groupColName := !!groupVec) 
  row.names(umap_df) <- rownames.df
  
  newDataset <- umap_df[row.names(umap_df) == "newDataset",]
  dat.woNewDataset <- umap_df %>% 
    dplyr::filter(!!as.name(groupColName) != "newDataset")
  
  highlightTypes <- c("Metatranscriptomic", "Marbel", "JLAB Microbiome")
  
  # Add a column to indicate highlighted types
  dat.woNewDataset <- dat.woNewDataset %>%
    mutate(highlight = ifelse(`Data type`  %in% highlightTypes, "mark", "not_mark"))
  
  # Define custom marker properties for highlighted and normal points
  markerProps <- list(
    size = ~ifelse(highlight == "mark", 6, 3),  # Larger size for highlighted
    opacity = ~ifelse(highlight == "mark", 1, 0.5),  # Higher opacity for highlighted
    symbol = ~ifelse(highlight == "mark", "circle", "circle-open")  # Different symbol
  )  
  
  fig <- plotly::plot_ly(dat.woNewDataset, x = ~UMAP1, y = ~UMAP2, z = ~UMAP3,
                         color = ~as.factor(dat.woNewDataset[[groupColName]]),
                         colors = customColors,
                         type = "scatter3d", mode = "markers",
                         marker = markerProps,
                         hovertext = paste("Dataset ID :", 
                                           row.names(dat.woNewDataset))
  ) 
  
  fig <- fig %>%
    plotly::layout(
      hoverlabel = list(namelength = -1),
      legend = list(itemsizing = "constant")
    )
  
  fig
}

generatePlots <- function(mtx, output, omicsTypes){
  
  mtx[mtx == 0] <- NA
  mtx[mtx == Inf] <- NA
  mtx <- mtx[, colSums(is.na(mtx)) != nrow(mtx)]
  
  allDataTypeLevels <- c("Data type", "Data type subgroups")
  selectedDataTypeLevel <- "Data type"
  
  nNegativeNumbers <- sum(mtx < 0, na.rm = TRUE)
  mtx <- log2(mtx)
  mtx[mtx == "NaN"] <- NA
  
  correlationplot <- plotCorrelation(mtx)
  gc()

  data <- integrateNewDataset(mtx, nNegativeNumbers = nNegativeNumbers) %>% 
    dplyr::select(-setdiff(!!allDataTypeLevels, !!selectedDataTypeLevel)) %>% 
    dplyr::rename("Data type" = !!selectedDataTypeLevel)
  mtx <- NULL
  gc()
  
  data <- data[
    !(data$`Data type` %in% c("Metabolomics (MS, Other)",
                              "Lipidomics (MS, Other)",
                              "Proteomics (LFQ, PRIDE, Multiple)", 
                              "Proteomics (LFQ, PRIDE, Other)", 
                              "Proteomics (Intensity, PRIDE, Multiple)",
                              "Proteomics (Intensity, PRIDE, Other)",
                              "Proteomics (iBAQ, PRIDE, Multiple)", 
                              "Proteomics (iBAQ, PRIDE, Other)"
    )),]
  
  #TODO: add metatranscriptomics when need be
  data <- data %>% dplyr::group_by(`Data type`) %>% 
    filter(n() > 5 | `Data type` == "newDataset" |  `Data type` == "Marbel" | `Data type` == "Metatranscriptomic") %>% ungroup
  data$`Data type` <- droplevels(data$`Data type`)

  ##############################################################################
  boxplotCols <- setdiff(
    unique(c("Dataset ID", "Data type", "# Samples", "# Analytes", 
             "min(% NA in analytes)", 
             "max(% NA in analytes)", "min(% NA in samples)", 
             "max(% NA in samples)", 
             "% NA", "% Analytes with NAs", "% Samples with NAs", "Mean", 
             "Median", "Min", "Max", "median(Variance of samples)", 
             "median(Variance of analytes)", 
             "Variance", "Kurtosis", "Skewness", "|Skewness|", 
             "% Var. explained by PC1", 
             "% Var. explained by PC2", 
             "Lin. coef. of Poly2(Means vs. Vars) (Analytes)", 
             "Quadr. coef. of Poly2(Means vs. Vars) (Analytes)", 
             "Agglom. coef. hierarch. analyte clustering", 
             "% Distinct values", "Corr(Mean vs. % NA) (Samples)", 
             "Corr(Mean vs. % NA) (Analytes)",
             "sd(Intensity w/ prob(NA) = 50% for sample)", 
             "sd(Intensity w/ prob(NA) = 90% for sample)")), 
    c("Data type", "Dataset ID"))
  
  data <- data %>% column_to_rownames("Dataset ID") %>% 
    dplyr::select(c("Data type", !!boxplotCols))
  
  # To be log2-transformed:
  toBeLog2Transformed <- c("# Samples", "# Analytes", 
                           "min(% NA in analytes)", "max(% NA in analytes)", 
                           "% Analytes with NAs", "% Samples with NAs",
                           "median(Variance of samples)", 
                           "median(Variance of analytes)", 
                           "Variance", "Kurtosis", "|Skewness|", "Skewness",
                           "Lin. coef. of Poly2(Means vs. Vars) (Analytes)", 
                           "Quadr. coef. of Poly2(Means vs. Vars) (Analytes)",
                           "sd(Intensity w/ prob(NA) = 50% for sample)", 
                           "sd(Intensity w/ prob(NA) = 90% for sample)")
  
  colsWithNegativeNumbers <- colnames(
    data[, sapply(data, FUN = function(x) any(x <= 0, na.rm = TRUE))])
  
  toBeLog2Transformed <- setdiff(toBeLog2Transformed, colsWithNegativeNumbers)  
  
  for (var in toBeLog2Transformed) {
    data <- logTransform(df = data, variable = var, logBase = "log2")
  }
  
  data[sapply(data, is.infinite)] <- NA
  
  groupColName <- "Data type"
  plotlyPCA <- plot3DPCA(df = data, groupColName = groupColName)
  
  newDatasetValues <- data %>% filter(`Data type` == "newDataset") %>% 
    t() %>% as.data.frame()
  
  df.newDatasetOnly <- data %>% filter(`Data type` == "newDataset") %>% 
    dplyr::select(-`Data type`)
  df.newDatasetOnly.vec <- as.vector(unlist(df.newDatasetOnly))
  names(df.newDatasetOnly.vec) <- colnames(df.newDatasetOnly)
  df.newDatasetOnly.vec <- paste(names(df.newDatasetOnly.vec), 
                                 round(df.newDatasetOnly.vec, 3), 
                                 sep = ": ", collapse = "<br/>")
  
  pcaMethod <- "nipals"
  
  plotlyUMAP <- getUMAPNewDataset(data, groupColName = groupColName)
  selectedDataType <- omicsTypes # "Proteomics (LFQ, PRIDE)"
  boxplot.df <- data %>% filter(`Data type` == !!selectedDataType)
  
  data <- NULL
  gc()
  
  boxplot.df.long <- reshape2::melt(boxplot.df)
  lev <- c("log2(# Analytes)", "log2(# Samples)", 
           "Mean", "Median", "Min", "Max", 
           "log2(Variance)", "log2(median(Variance of samples))", 
           "log2(median(Variance of analytes))",
           "Kurtosis", "Skewness", "log2(|Skewness|)", "% Distinct values",
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
  
  newDatasetValues <- 
    tibble::rownames_to_column(newDatasetValues, "variable") %>% 
    filter(variable != "Data type")
  colnames(newDatasetValues) <- c("variable", "value")
  newDatasetValues$value <- as.numeric(newDatasetValues$value)
  
  quantile.df <- data.frame(t(sapply(boxplot.df[, 2:ncol(boxplot.df)],
                                     function(x) quantile(x[!is.na(x)], 
                                                          c(.05, .95)) )))
  quantile.df$variable <- row.names(quantile.df)
  newDatasetValues <- dplyr::left_join(newDatasetValues, quantile.df, 
                                       by = "variable")
  newDatasetValues$included <- ifelse(
    (newDatasetValues$value >= newDatasetValues$X5.) &
      (newDatasetValues$value <= newDatasetValues$X95.), 'green', 'red')
  
  gg.boxplot <- ggplot(boxplot.df.long, aes(x = value, y = factor(1))) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.5, alpha = 0.25, outlier.size = 0.5) +
    geom_rect(data = newDatasetValues, 
              aes(fill = factor(included, levels = c('red', 'green'))), 
              xmin = -Inf, xmax = Inf,
              ymin = -Inf, ymax = Inf, alpha = 0.3) +
    geom_vline(data = newDatasetValues, aes(xintercept = as.numeric(value)), 
               colour = 'blue') +
    facet_wrap(. ~ factor(variable, levels=lev), ncol = 4, scales = "free_x") +
    ggplot2::theme_bw() +
    theme(axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none")
  
  gc()
  
  list(plotlyUMAP = plotlyUMAP,
       gg.boxplot = gg.boxplot, 
       plotlyPCA = plotlyPCA,
       correlationplot = correlationplot,
       df.newDatasetOnly.vec = df.newDatasetOnly.vec)
}

ui <- fluidPage(
  title = "Analysis",
  sidebarLayout(
    sidebarPanel(
      title = "Inputs",
      fileInput("tsv_input", 
                "Select .tsv file to import (Format: samples in columns and analytes in rows, max. 100 MB)",
                accept = ".tsv"),
      checkboxInput("header", "Header", TRUE),
      checkboxInput("firstColumnRownames", 
                    "First column contains row names", TRUE),
      h6("Example file from Fröhlich et al. https://doi.org/10.1038/s41467-022-30094-0:"),
      downloadButton("downloadData", "Download example file"),
      
      hr(),
      
      radioButtons("omicsTypes", "Omics types",
                   choices =  c('Metabolomics (NMR)', 
                                'Metabolomics (MS)',
                                'Lipidomics (MS)',
                                'Proteomics (LFQ, PRIDE)',
                                'Proteomics (Intensity, PRIDE)',
                                'Proteomics (iBAQ, PRIDE)',
                                'Proteomics (Intensity, Expression Atlas)',
                                'Proteomics (iBAQ, Expression Atlas)', 
                                'RNA-seq (raw)',
                                'RNA-seq (FPKM)',
                                'RNA-seq (TPM)',
                                'Microarray',
                                'scProteomics',
                                'scRNA-seq (unnormalized)',
                                'scRNA-seq (normalized)',
                                'Microbiome', 
                                'JLAB Microbiome',
                                'Marbel',
                                'Metatranscriptomic'),
                   selected = "Proteomics (LFQ, PRIDE)"),
      br(),
      actionButton("run_button", "Run Analysis", icon = icon("play")),
      h5("Please be patient, analysis can take a while.")
    ),
    mainPanel(
      
      h2("Data characteristics:"),
      withSpinner(htmlOutput(outputId = "dataCharacteristics"), type = 5),
      hr(),
      h2("How close is the provided dataset to the selected omics type?"),
      h3("If the dataset lies between the 5th and 95th percentile, the box is colored green, else red."),
      withSpinner(plotOutput(outputId = "boxplot", 
                             width = "1100px",height = "500px"),
                  type = 5),
      hr(),
      h2("PCA"),
      h3("Omics types can be switched on and off by clicking on the respective symbol in the legend."),
      withSpinner(plotlyOutput('plotlyPCA', width = "1100px", height = "500px"), 
                  type = 5),
      hr(),
      h2("UMAP"),
      h3("Omics types can be switched on and off by clicking on the respective symbol in the legend."),
      withSpinner(plotlyOutput('plotlyUMAP', width = "1100px", height = "500px"), 
                  type = 5),
      hr(),
      h2("Sample mean vs. %NA for provided dataset"),
      withSpinner(plotOutput(outputId = "correlationplot", 
                             width = "1100px", height = "500px"), type = 5)
      
    )
  )
)


server <- function(input, output){
  
  options(shiny.maxRequestSize = 100*1024^2) 
  
  data_input <- reactive({
    default_file_path <- "/home/ubuntu/Data/init/testiecountie.tsv"
    
    # Check if a user uploads a file; otherwise, use the default path
    file_path <- if (!is.null(input$tsv_input)) {
      input$tsv_input
    } else {
      default_file_path
    }
    
    mtx <- data.table::fread(file_path,
                             header = TRUE,
                             data.table = FALSE,
                             check.names = FALSE,
                             sep = '\t')
    
    if (input$firstColumnRownames) mtx <- mtx[, 2:ncol(mtx)]
    mtx <- as.matrix(mtx)
    mtx <- removeEmptyRowsAndColumns(mtx, zerosToNA = TRUE)
    
    validate(
      need(is.numeric(mtx), "Please provide numeric dataset")
    )
    mtx
  })
  
  omicsTypes <- eventReactive(input$run_button, input$omicsTypes)
  
  plots_initial <- reactive({
    init_data = init_database()
    generateInitialPCAandUMAP(init_data)
  })
  
  plots <- reactive({
    generatePlots(data_input(), output, omicsTypes = "Default")
  })
  
  gc() 
  
  output$plotlyPCA <- renderPlotly({
    req(plots())
    plots()[["plotlyPCA"]]
  })
  
  output$plotlyUMAP <- renderPlotly({
    req(plots())
    plots()[["plotlyUMAP"]]
  })
  
  
  output$correlationplot <- renderPlot(plots()[["correlationplot"]])
  output$boxplot <- renderPlot(plots()[["gg.boxplot"]])
  output$dataCharacteristics <- renderUI(HTML(
    plots()[["df.newDatasetOnly.vec"]]))
  
  output$downloadData <- downloadHandler(
    filename = function() {
      "DIANN_DIANN_AI_GPF_example.tsv"
    },
    content = function(file) {
      write.table(read.table("DIANN_DIANN_AI_GPF_example.tsv", 
                             check.names = FALSE, sep = '\t', header = TRUE), 
                  file, 
                  sep = '\t', 
                  row.names = FALSE)
    }
  )
  
}

shinyApp(ui = ui, server = server)