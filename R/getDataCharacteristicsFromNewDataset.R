setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##########################################################################################
library(data.table)
library(pcaMethods)
library(Matrix)
library(dplyr)
# library(foreach)
# library(BimodalIndex)
library(ggplot2)
theme_set(theme_bw())
library(cluster)
library(amap)
library(mclust)
library(mlr3misc)
library(biglm)
library(reshape2)
library(speedglm)

`%notin%` <- Negate(`%in%`)

removeEmptyRowsAndColumns <- function(mtx, zerosToNA = FALSE){
  if (zerosToNA) mtx[mtx == 0] <- NA
  mtx[mtx %in% c("NaN", "N/A", "<-->")] <- NA
  # remove rows with only NAs
  # mtx <- mtx[rowSums(is.na(mtx) | mtx == 0) != ncol(mtx), ]
  mtx <- subset(mtx, rowSums(is.na(mtx) | mtx == 0) != ncol(mtx))
  if (!is.vector(mtx)) mtx <- mtx[, colSums(is.na(mtx) | mtx == 0) != nrow(mtx)]
  mtx
}

#################################################################################
# DATA CHARACTERISTICS

calc_kurtosis <- function(data){
  n <- length(data[!is.na(data)])
  data <- data - mean(data, na.rm=TRUE)
  r <- n * sum(data^4, na.rm=TRUE) / (sum(data^2, na.rm=TRUE)^2)
  return(r * (1 - 1/n)^2 - 3)
}

calc_variance <- function(data){
  # var(c(data), na.rm=TRUE) # was replaced because for single cell data large vectors can become a problem
  1/(sum(!is.na(data))-1)*sum((data-mean(data, na.rm = TRUE))^2, na.rm = TRUE) 
}

calc_sd <- function(data){
  sqrt(calc_variance(data))
}

calc_skewness <- function (data){
  data <- data[!is.na(data)]
  # return(sum((data - mean(data))^3)/(length(data) * sd(data)^3))
  return(sum((data - mean(data))^3)/(length(data) * calc_sd(data)^3))
}

calc_RMS <- function(data) sqrt(mean(data^2, na.rm=TRUE))

applyFunctionWithSeed <- function(functionName, seed = 123,  ...){
  
  #oldseed <- .Random.seed
  oldseed <- mlr3misc::get_seed()
  
  if (is.null(seed)) {
    # seed <- .Random.seed
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
  # return(apply(mtx,1,sd,na.rm=T))
  return(apply(mtx, 1, calc_sd))
}

# Poly2 (features)
# linear coefficient of 2nd order polynomial fit with x = row means, y = row variances
get_LinearCoefPoly2XRowMeansYRowVars <- function(mtx, ...){
  return(unname(coef(biglm::biglm(y~x+I(x^2), data=data.frame(
    y=get_rowSd(mtx)^2,
    x=get_rowMeans(mtx))))[2]))
}

# quadratic coefficient of 2rd order polynomial fit with x = row means, y = row variances
get_QuadraticCoefPoly2XRowMeansYRowVars <- function(mtx, ...){ 
  return(unname(coef(biglm::biglm(y~x+I(x^2), data=data.frame(
    y=get_rowSd(mtx)^2,
    x=get_rowMeans(mtx))))[3]))
}

# linear and quadratic coefficient of 2rd order polynomial fit with x = row means, y = row variances
get_CoefPoly2XRowMeansYRowVars <- function(mtx, ...){
  
  rowSd <- rowMean <- linearCoef <- quadraticCoef <- NA
  try({rowSd <- get_rowSd(mtx)^2})
  try({rowMean <- get_rowMeans(mtx)})
  
  try({ 
    if (!(length(rowSd) == 1 && is.na(rowSd)) & 
        !(length(rowMean) == 1 && is.na(rowMean))) {
      df <- data.frame(
        y=rowSd,
        x=rowMean)
      model <- biglm::biglm(y~x+I(x^2), data=df)
      coefs <- coef(model)
      linearCoef <- unname(coefs[2])
      quadraticCoef <- unname(coefs[3])
    }
  })
  return(list(linearCoef = linearCoef, quadraticCoef = quadraticCoef))
}

# Coef.hclust (features)
get_coefHclustRowsWithFewestNAs <- function(mtx, naToZero = FALSE, nMax = 500, ...) {
  mtx <- mtx[order(rowSums(is.na(mtx)), -rowMeans(mtx, na.rm = TRUE)), ]
  
  mtx <- mtx[order(rowMeans(mtx, na.rm = TRUE)), ]
  mtx <- mtx[order(rowSums(is.na(mtx))), ]
  if (naToZero) mtx[is.na(mtx)] <- 0
  randomCols <- applyFunctionWithSeed(sample, x = 1:ncol(mtx), size= min(nMax, ncol(mtx)))
  mtx <- mtx[1:min(500, nrow(mtx)), randomCols$res] # max. 500 features and samples
  mtx <- removeEmptyRowsAndColumns(mtx, zerosToNA = FALSE)
  return(list(res = cluster::coef.hclust(amap::hcluster(mtx)), seed = randomCols$seed))
}

calculateIntensityNAProbability5090 <- function(mtx, nmaxSamples = 200) {
  
  seedUsed <- NA
  if (ncol(mtx) > nmaxSamples){ # random subset of features are selected
    randomCols <- applyFunctionWithSeed(sample, x = 1:ncol(mtx), size= min(nmaxSamples, ncol(mtx)))
    seedUsed <- randomCols$seed
    mtx <- mtx[, randomCols$res]
    mtx <- removeEmptyRowsAndColumns(mtx, zerosToNA = TRUE)
  }
  
  data.long <- reshape2::melt(mtx)
  colnames(data.long) <- c("Feature", "Sample", "Value")
  data.long$isNA <- as.integer(is.na(data.long$Value))
  
  # Group by mean using dplyr
  featureMean.df <- data.long %>% dplyr::group_by(Feature) %>% 
    dplyr::summarise(mean=mean(Value, na.rm=TRUE))
  
  imputed <- data.long$Value
  
  data.long <- data.long %>% dplyr::left_join(featureMean.df, by='Feature')
  imputed[is.na(imputed)] <- data.long$mean[is.na(imputed)]
  data.long$imputed <- imputed
  
  
  x5090.lst <- lapply(unique(data.long$Sample), function(sample){
    data.long.Sample <- data.long %>% dplyr::filter(Sample == sample)
    x5090 <- NA
    glmModelSample <- NA
    if (any(data.long.Sample$isNA == 1)){
      glmModelSample <- speedglm::speedglm(isNA ~ imputed, 
                                           family=binomial(link='logit'), 
                                           data = data.long.Sample, fitted = TRUE)
      
      if (unname(coefficients(glmModelSample)[2]) < 0){
        intensities <- seq(min(data.long$imputed)-20, max(data.long$imputed),
                           length.out = 400)
        x5090 <- approx(predict(glmModelSample, 
                                newdata=data.frame(imputed=intensities),
                                type="response"), intensities, c(0.5, 0.9))$y
      }
    }
    x5090
  })
  
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
  #KS.SignProp <- kolSmirTestSignProp(mtx)
  
  if (is.vector(mtx)){
    nFeatures <- 1
  } else {
    nFeatures <- nrow(mtx) # number of proteins with no NAs
  }
  
  mean <- mean(mtx, na.rm = TRUE)
  median <- median(mtx, na.rm = TRUE)
  min <- min(mtx, na.rm = TRUE)
  max <- max(mtx, na.rm = TRUE)
  
  medianSampleVariance <- medianAnalyteVariance <- skewness <- kurtosis <- variance <-
    prctPC1 <- prctPC2 <- 
    linearCoefPoly2Row <- quadraticCoefPoly2Row <- 
    coefHclustRows <- 
    intensityNAProb50.sd <- intensityNAProb90.sd <- intensityNAProbnSamplesWithProbValue <- NA
  
  try({medianSampleVariance <- median(apply(mtx, 2, calc_variance), na.rm = TRUE)})
  try({medianAnalyteVariance <- median(unname(apply(mtx, 1, calc_variance)), na.rm = TRUE)})
  
  try({skewness <- calc_skewness(mtx)})
  try({kurtosis <- calc_kurtosis(mtx)})
  
  try({variance <- calc_variance(mtx)})
  
  # Poly2 (features)
  # try(linearCoefPoly2Row <- get_LinearCoefPoly2XRowMeansYRowVars(mtx))
  # try(quadraticCoefPoly2Row <- get_QuadraticCoefPoly2XRowMeansYRowVars(mtx))
  try({
    coefs <- get_CoefPoly2XRowMeansYRowVars(mtx)
    linearCoefPoly2Row <- coefs[["linearCoef"]]
    quadraticCoefPoly2Row <- coefs[["quadraticCoef"]]
  })

  # Coef.hclust (features)
  try({
    coefHclustRowsRes <- get_coefHclustRowsWithFewestNAs(mtx, naToZero = TRUE)
    coefHclustRows <- coefHclustRowsRes$res
    # coefHclustRowsSeed <- coefHclustRowsRes$seed
  })

  try({
    x5090.sds <- calculateIntensityNAProbability5090(mtx)
    intensityNAProb50.sd <- x5090.sds[["IntensityNAp50"]]
    intensityNAProb90.sd <- x5090.sds[["IntensityNAp90"]]
    # intensityNAProbSeed <- x5090.sds[["seed"]]
    intensityNAProbnSamplesWithProbValue <- x5090.sds[["nSamplesWithProbValue"]]
  })
  
  
  mtx <- mtx %>% t()
  mtx <- mtx[ , which(apply(mtx, 2, calc_variance) != 0)] # Remove zero variance columns 
  
  if (!is.vector(mtx)){
    try({
      pca <- pcaMethods::pca(mtx, method="nipals", center = TRUE, maxSteps=5000)
      prctPC1 <- pca@R2[1]
      prctPC2 <- pca@R2[2]
    })
  }
  
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
  rowNonNaNumber <- rowSums(!is.na(mtx))
  
  colMeans <- colMeans(mtx, na.rm = TRUE)
  rowMeans <- rowMeans(mtx, na.rm = TRUE)
  
  corSampleMeanNA <- corAnalyteMeanNA <- NA
  corCoefType <- "spearman"
  try({
    cortestCol <- cor.test(colNaPercentage, colMeans, method = corCoefType)
    corSampleMeanNA <- unname(cortestCol$estimate)
  })
  
  try({
    cortestRow <- cor.test(rowNaPercentage, rowMeans, method = corCoefType)
    corAnalyteMeanNA <- unname(cortestRow$estimate)
  })
  
  c(
    minRowNaPercentage = min(rowNaPercentage),            
    maxRowNaPercentage = max(rowNaPercentage),
    minColNaPercentage = min(colNaPercentage),
    maxColNaPercentage = max(colNaPercentage),
    percNATotal = mean(is.na(mtx)) * 100,
    percOfRowsWithNAs = sum(apply(mtx, 1, anyNA))/nrow(mtx) * 100,
    percOfColsWithNAs = sum(apply(mtx, 2, anyNA))/ncol(mtx) * 100,
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
    
    characts <- getCharacteristicsHelper(mtx)
    
    charact.log <- c(naFeatures, characts, nDistinctValues = nDistinctValues, nNegativeNumbers = nNegativeNumbers)
    
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
  
  cond <- pca.obj$x[which(groups == "newDataset"),]
  
  P2 <- ggbiplot::ggbiplot(pca.obj,
                           choices = PCchoices,
                           obs.scale = 1, 
                           var.scale=1,
                           ellipse=ellipse,
                           circle=F,
                           varname.size=3,
                           var.axes=T,
                           groups=groups, 
                           alpha=0)  +
    # scale_color_discrete(name = '') +  
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
  
  my_colors <- scales::hue_pal()(length(unique(groups))-1)
  
  P2 <- P2 + 
    scale_color_manual(name = '',
      values = my_colors,
      limits = setdiff(unique(groups), "newDataset")
    ) +
    geom_point(aes(x=cond[PCchoices[1]], y=cond[PCchoices[2]]), col="red", size=3) 
  
  P2
}

plotPCABiplotsNewDataset <- function(df, groupColName = "", addStr = "", pcaMethod = "nipals") {
  
  # pdf(file = paste0("biplot_", pcaMethod, "_facetZoom_PC1vs2_", addStr, ".pdf"), width = 12, height = 10)
  # print(plotPCABiplot(df = df %>% dplyr::select(-!!groupColName), 
  #                     groups= df[[groupColName]],
  #                     alpha = 0.3,
  #                     pcaMethod = pcaMethod,
  #                     facetZoom = TRUE,
  #                     PCchoices = 1:2,
  #                     xlimLower = -4, xlimUpper = 5,
  #                     ylimLower = -5, ylimUpper = 8))
  # dev.off()
  # 
  # pdf(file = paste0("biplot_", pcaMethod, "_facetZoom_PC1vs3_", addStr, ".pdf"), width = 12, height = 10)
  # print(plotPCABiplot(df = df %>% dplyr::select(-!!groupColName), 
  #                     groups= df[[groupColName]],
  #                     alpha = 0.3,
  #                     pcaMethod = pcaMethod,
  #                     facetZoom = TRUE,
  #                     PCchoices = c(1, 3),
  #                     xlimLower = -5, xlimUpper = 5,
  #                     ylimLower = -4, ylimUpper = 8))
  # dev.off()
  # 
  # pdf(file = paste0("biplot_", pcaMethod, "_facetZoom_PC2vs3_", addStr, ".pdf"), width = 12, height = 10)
  # print(plotPCABiplot(df = df %>% dplyr::select(-!!groupColName), 
  #                     groups= df[[groupColName]],
  #                     alpha = 0.3,
  #                     pcaMethod = pcaMethod,
  #                     facetZoom = TRUE,
  #                     PCchoices = c(2, 3),
  #                     xlimLower = -5, xlimUpper = 7,
  #                     ylimLower = -4, ylimUpper = 8
  # ))
  # dev.off()
  
  ## Remove outlier "E-GEOD-152766.aggregated_filtered_counts.mtx"
  # df <- df[row.names(df) != "E-GEOD-152766.aggregated_filtered_counts.mtx", ]
  pdf(file = paste0("biplot_", pcaMethod, "_", addStr, ".pdf"), width = 12, height = 10)
  print(plotPCABiplotNewDataset(df = df %>% dplyr::select(-!!groupColName), 
                      groups= df[[groupColName]],
                      alpha = 0.3,
                      pcaMethod = pcaMethod,
                      coordRatio = 1/3,
                      facetZoom = FALSE))
  dev.off()
  
  pca <- pcaMethods::pca(df %>% dplyr::select(-!!groupColName), method=pcaMethod, nPcs=4, center=TRUE
                         , scale = "uv")
  write.csv(pca@loadings, paste0("loadings_",  pcaMethod, "_", addStr,".csv"))
  dat <- merge(pcaMethods::scores(pca), df, by=0)
  # 
  # # library(GGally)
  # # pdf(paste0("ggpairs_", pcaMethod, "_", addStr, ".pdf"), width = 12, height = 10)
  # # print(GGally::ggpairs(dat, columns = 2:5, ggplot2::aes(colour=get(groupColName)),
  # #                       lower = list(continuous = wrap("smooth", alpha = 0.3, size = 1), 
  # #                                    combo = wrap("dot_no_facet", alpha = 0.4)),
  # #                       upper=list(continuous = wrap("cor", method = "spearman", size = 3)),
  # #                       mapping=aes(color = get(groupColName),
  # #                                   fill= get(groupColName), 
  # #                                   alpha=0.5)) +
  # #         ggtitle(groupColName) +
  # #         theme_bw())
  # # dev.off()
  # 
  # library(plotly)
  # library(htmlwidgets)
  # fig <- plot_ly(dat, x = ~PC1, y = ~PC2, z = ~PC3, 
  #                color = ~as.factor(dat[[groupColName]]), 
  #                type="scatter3d", mode="markers",
  #                #colors = c('#636EFA','#EF553B') , 
  #                marker = list(size = 2)
  #                # , alpha = 0.75
  # ) #%>%
  # # add_markers(size = 5, marker=list(sizeref=8, sizemode="area"))
  # fig <- fig %>%
  #   layout(
  #     title = "3D PCA",
  #     scene = list(bgcolor = "#e5ecf6"
  #     )
  #   )
  # 
  # htmlwidgets::saveWidget(fig, paste0("plotly_", pcaMethod, "_", addStr,".html"), selfcontained = F, libdir = "lib")
}

################################################################################

mtx <- readRDS("diaWorkflowResults_allDilutions.rds")
mtx <- 2^as.matrix(mtx[["DIANN_DIANN_AI_GPF"]])

# if (nrow(mtx) > 9 & ncol(mtx) > 4) {}
mtx <- removeEmptyRowsAndColumns(mtx, zerosToNA = TRUE)

data <- getDataCharacteristics(mtx)
data <- data.frame(as.list(data))


data$prctPC1 <- 100 * data$prctPC1
data$prctPC2 <- 100 * data$prctPC2

data <- data %>%
  mutate_if(is.integer, as.numeric) %>% 
  mutate("prctnDistinctValues" = nDistinctValues/((nSamples * nAnalytes) * ((100 - percNATotal)/100)) * 100 ) %>%
  dplyr::select(-c(nDistinctValues, nNegativeNumbers))

data$datasetID <- data$dataType <- data$dataTypeSubgroups <- "newDataset"
# seedCols <- c("coefHclustRowsSeed", "intensityNAProbSeed")
# data <- data[,!(colnames(data) %in% seedCols)]

Oldnames <- c("datasetID", "dataType", "dataTypeSubgroups", "nSamples", "nAnalytes", "minRowNaPercentage", 
              "maxRowNaPercentage", "minColNaPercentage", "maxColNaPercentage", 
              "percNATotal", "percOfRowsWithNAs", "percOfColsWithNAs", "corSampleMeanNA", 
              # "corSampleMeanNAPval", 
              "corAnalyteMeanNA", 
              # "corAnalyteMeanNAPval", 
              "mean", "median", "min", "max", "medianSampleVariance", "medianAnalyteVariance", 
              "variance", "kurtosis", "skewness", "prctPC1", "prctPC2", 
              # "bimodalityColCorr", 
              "linearCoefPoly2Row", "quadraticCoefPoly2Row", "coefHclustRows", 
              "intensityNAProb50.sd", "intensityNAProb90.sd", "intensityNAProbnSamplesWithProbValue", 
              "prctnDistinctValues")

Newnames <- c("Dataset ID", "Data type", "Data type subgroups", "# Samples", "# Analytes", "min(% NA in analytes)", 
              "max(% NA in analytes)", "min(% NA in samples)", "max(% NA in samples)", 
              "% NA", "% Analytes with NAs", "% Samples with NAs", 
              "Corr(Mean vs. % NA) (Samples)", 
              # "Corr(Mean vs. % NA) (Samples) (p-Value)", 
              "Corr(Mean vs. % NA) (Analytes)", 
              # "Corr(Mean vs. % NA) (Analytes) (p-Value)", 
              "Mean", "Median", "Min", "Max", "median(Variance of samples)", "median(Variance of analytes)", 
              "Variance", "Kurtosis", "Skewness", "% Var. explained by PC1", "% Var. explained by PC2", 
              # "Bimodality of sample correlations", 
              "Lin. coef. of Poly2(Means vs. Vars) (Analytes)", 
              "Quadr. coef. of Poly2(Means vs. Vars) (Analytes)", 
              "Agglom. coef. hierarch. analyte clustering", 
              "sd(Intensity w/ prob(NA) = 50% for sample)", "sd(Intensity w/ prob(NA) = 90% for sample)", 
              "# Samples w/ intensityNAProb50.sd or intensityNAProb90.sd", 
              "% Distinct values")

renameTable <- data.frame(Oldnames = Oldnames,
                          Newnames = Newnames)

data <- data %>% dplyr::rename_with(~ Newnames[which(Oldnames == .x)], .cols = Oldnames)
data <- data %>% dplyr::mutate(`|Skewness|` = abs(Skewness))



data.allDatasets <-  read.csv("datasets_results_clean_renamed.csv", check.names = FALSE)
data.allDatasets$`Corr(Mean vs. % NA) (Samples) (p-Value)` <- 
  data.allDatasets$`Corr(Mean vs. % NA) (Analytes) (p-Value)` <-
  data.allDatasets$`Bimodality of sample correlations` <- NULL

data.allDatasets <- rbind(data.allDatasets, data)
################################################################################


allDataTypeLevels <- c("Data type", "Data type subgroups")
# selectedDataTypeLevel <- "Data type subgroups"
data.copy <- data.allDatasets
for (selectedDataTypeLevel in allDataTypeLevels) {
  
  data <- data.copy %>% dplyr::select(-setdiff(!!allDataTypeLevels, !!selectedDataTypeLevel)) %>% dplyr::rename("Data type" = !!selectedDataTypeLevel)
  data <- data[!(data$`Data type` %in% c("Metabolomics (Undefined-MS)",
                                         "Metabolomics (Other ionization-MS)",
                                         "Proteomics (iBAQ, PRIDE, Undefined)", 
                                         "Proteomics (Intensity, PRIDE, Undefined)", 
                                         "Proteomics (LFQ, PRIDE, Undefined)")),]
  
  data <- data %>% dplyr::group_by(`Data type`) %>% filter(n() > 5 | `Data type` == "newDataset") %>% ungroup
  
  ################################################################################
  boxplotCols <- setdiff(unique(c("Dataset ID", "Data type", "# Samples", "# Analytes", "min(% NA in analytes)", 
                                  "max(% NA in analytes)", "min(% NA in samples)", "max(% NA in samples)", 
                                  "% NA", "% Analytes with NAs", "% Samples with NAs", "Mean", 
                                  "Median", "Min", "Max", "median(Variance of samples)", "median(Variance of analytes)", 
                                  "Variance", "Kurtosis", "Skewness", "|Skewness|", "% Var. explained by PC1", 
                                  "% Var. explained by PC2", "Lin. coef. of Poly2(Means vs. Vars) (Analytes)", 
                                  "Quadr. coef. of Poly2(Means vs. Vars) (Analytes)", "Agglom. coef. hierarch. analyte clustering", 
                                  "% Distinct values", "Corr(Mean vs. % NA) (Samples)", "Corr(Mean vs. % NA) (Analytes)",
                                  "sd(Intensity w/ prob(NA) = 50% for sample)", 
                                  "sd(Intensity w/ prob(NA) = 90% for sample)")), c("Data type", "Dataset ID"))
  
  data2 <- data[, c("Data type", 
                    boxplotCols)]
  
  # To be log2-transformed:
  toBeLog2Transformed <- c("# Samples", "# Analytes", 
                           "min(% NA in analytes)", "max(% NA in analytes)", 
                           "% Analytes with NAs", "% Samples with NAs",
                           "median(Variance of samples)", "median(Variance of analytes)", 
                           "Variance", "Kurtosis", "|Skewness|", "Skewness",
                           "Lin. coef. of Poly2(Means vs. Vars) (Analytes)", 
                           "Quadr. coef. of Poly2(Means vs. Vars) (Analytes)",
                           "sd(Intensity w/ prob(NA) = 50% for sample)", 
                           "sd(Intensity w/ prob(NA) = 90% for sample)")
  
  colsWithNegativeNumbers <- colnames(data2[, sapply(data2, FUN = function(x) any(x <= 0, na.rm = TRUE))])
  
  toBeLog2Transformed <- setdiff(toBeLog2Transformed, colsWithNegativeNumbers)  
  
  for (var in toBeLog2Transformed){
    # print(var)
    data2 <- logTransform(df = data2, variable = var, logBase = "log2")
  }
  
  #write.csv(data2, paste0("data2_", gsub(" ", "_", selectedDataTypeLevel), ".csv"), row.names = FALSE)
  data2[sapply(data2, is.infinite)] <- NA
  

  
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
  
  


  # One of: "median(Variance of samples)", "log2(Variance)"
  # One of: "Mean", "Median", "Min", "Max"
  # One of : "min(% NA in samples)", "max(% NA in samples)", "% NA", "% Analytes with NAs", "% Samples with NAs", 
  
  


  
  if (selectedDataTypeLevel == "Data type") {
    # plotPCABiplotsNewDataset(df = data2.complete, groupColName = "Data type", addStr = paste0(gsub(" ", "_", selectedDataTypeLevel), "_newDataset"), pcaMethod = "svd") # "svd" or "nipals"
    plotPCABiplotsNewDataset(df = data2, groupColName = "Data type", addStr = paste0(gsub(" ", "_", selectedDataTypeLevel), "_newDataset"), pcaMethod = "nipals") # "svd" or "nipals"
    
    
    df = data2
    groupColName = "Data type"
    pcaMethod = "nipals"
    
    plotPCABiplotNewDataset(df = df %>% dplyr::select(-!!groupColName), 
                            groups= df[[groupColName]],
                            alpha = 0.3,
                            pcaMethod = pcaMethod,
                            coordRatio = 1/3,
                            facetZoom = FALSE)
    
    groups= df[[groupColName]]
    
    pca <- pcaMethods::pca(df %>% dplyr::select(-!!groupColName), method=pcaMethod, nPcs=4, center=TRUE
                           , scale = "uv")
    
    
    loadings <- data.frame(pca@loadings)[, 1:2]
    loadings.abs <- abs(loadings)
    loadings.abs <- loadings.abs/colSums(loadings.abs)
    loadings.abs <- loadings.abs %>% mutate(sumPC1PC2 = PC1 + PC2)
    loadings.abs$variable <- row.names(loadings.abs)
    loadings.abs.sorted <- loadings.abs %>% arrange(-sumPC1PC2)
    
    scores <- data.frame(pca@scores)[,1:2]
    cond <- scores[which(groups == "newDataset"),]
    
    
    dataTypes <- setdiff(unique(df$`Data type`), "newDataset")
    densVals <- c()
    for (type in dataTypes) {
      
      scores.sub <- scores[which(groups == type),]
      
      x <- scores.sub$PC1
      y <- scores.sub$PC2
      
      library(MASS)
      dens <- kde2d(x,y)
      
      # create a new data frame of that 2d density grid
      # (needs checking that I haven't stuffed up the order here of z?)
      gr <- data.frame(with(dens, expand.grid(x,y)), as.vector(dens$z))
      names(gr) <- c("xgr", "ygr", "zgr")
      gr$zgr <- gr$zgr/sum(gr$zgr)
      
      # Fit a model
      mod <- loess(log2(zgr)~xgr*ygr, data=gr, span = 0.1
                   , control=loess.control(surface="direct")
      )
      
      plot <- TRUE
      if (plot) {
        all.data <- expand.grid(xgr = seq(from = min(gr$xgr), to = max(gr$xgr), length.out = 500),
                                ygr = seq(from = min(gr$ygr), to = max(gr$ygr), length.out = 500))
        
        # all.data <- expand.grid(xgr = seq(from = min(scores$PC1), to = max(scores$PC1), length.out = 500),
        #                         ygr = seq(from = min(scores$PC2), to = max(scores$PC2), length.out = 500))
        predicted <- predict(mod, newdata =  all.data)
        all.data = data.frame(W = 2^as.vector(predicted),
                              all.data)
        
        gg <- ggplot(all.data, aes(x = xgr, y = ygr, z = W)) +
          stat_contour(geom = "polygon", aes(fill = after_stat(level)) ) +
          geom_tile(aes(fill = W)) +
          stat_contour(bins = 10) +
          xlab("PC1") +
          ylab("PC2") +
          xlim(c(min(scores$PC1), max(scores$PC1))) +
          ylim(c(min(scores$PC2), max(scores$PC2))) +
          guides(fill = guide_colorbar(title = "Probability")) +
          ggtitle(type) +
          theme_bw() +
          geom_point(aes(x=cond$PC1, y=cond$PC2), col="red", size=3) 
        
        ggsave(file=paste0("dataType_2dDensity_", type, ".pdf"), gg)
      }
      
      # Apply the model to the original data to estimate density at that point
      pointdens <- predict(mod, newdata=data.frame(xgr=cond$PC1, ygr=cond$PC2))
      pointdens <- 2^pointdens
      # pointdens <- ifelse(is.na(pointdens), 0, pointdens)
      print(paste0(type, ": ", pointdens))
      densVals <- c(densVals, pointdens)
    }
    names(densVals) <- dataTypes
    densVals <- sort(densVals, decreasing = TRUE) 
    densVals
    # densVals <- ifelse(densVals<0, 0, densVals)
    top3 <- names(densVals[1:3])
    
    boxplot.df <- data2 %>% filter(`Data type` %in% top3)
    boxplot.df.long <- reshape2::melt(boxplot.df)
    
    newDatasetValues <- data2 %>% filter(`Data type` == "newDataset") %>% t() %>% as.data.frame()
    newDatasetValues <- tibble::rownames_to_column(newDatasetValues, "variable") %>% 
      filter(variable != "Data type")
    colnames(newDatasetValues) <- c("variable", "value")
    newDatasetValues$value <- as.numeric(newDatasetValues$value)
    
    
    boxplot.df.long$`Data type` <- factor(boxplot.df.long$`Data type`, levels = rev(top3))
    ggplot(boxplot.df.long, aes(`Data type`, value)) +
      geom_violin(aes(fill = `Data type`), alpha=0.5) +
      geom_boxplot(aes(fill = `Data type`), width=0.5, alpha=0.25, outlier.size=0.5) +
      coord_flip() +
      xlab("") +
      ylab("") +
      geom_hline(aes(yintercept = value), newDatasetValues, colour = 'red') +
      facet_wrap( ~ factor(variable, levels = loadings.abs.sorted$variable), scales = "free_x", ncol=6, strip.position = "bottom") +
      ggplot2::theme_bw() +
      ggtitle("Top 3") +
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
  } 
  
  
  # selForCorr <- c("% Distinct values", 
  #                 "log2(# Analytes)", "log2(# Samples)", 
  #                 "Mean", "log2(Variance)", 
  #                 "% NA", "max(% NA in analytes)", "Corr(Mean vs. % NA) (Samples)", "Corr(Mean vs. % NA) (Analytes)",
  #                 "Skewness", "log2(|Skewness|)", "Kurtosis", 
  #                 "% Var. explained by PC1", 
  #                 "% Var. explained by PC2",  
  #                 "Lin. coef. of Poly2(Means vs. Vars) (Analytes)", "Quadr. coef. of Poly2(Means vs. Vars) (Analytes)", 
  #                 "Agglom. coef. hierarch. analyte clustering", 
  #                 "sd(Intensity w/ prob(NA) = 50% for sample)", "sd(Intensity w/ prob(NA) = 90% for sample)"
  # )
  # 
  # data3 <- data2[, c("Data type", selForCorr)]
  # 
  # data3$`Data type` <- as.factor(data3$`Data type`)
  # library("partykit")
  # ct <- ctree(`Data type` ~ ., data = data3 %>% filter(`Data type` != "newDataset") %>% select("Data type", "% Distinct values", 
  #                                                                                              "log2(# Analytes)", "log2(# Samples)", 
  #                                                                                              "Mean", "log2(Variance)", 
  #                                                                                              "% NA"),
  #             control = ctree_control(alpha=0.0000000000000001))
  # plot(ct, tp_args = list(text = TRUE))
  
}

################################################################################
################################################################################
session <- sessionInfo()
sink(paste0("newDataset_getDataCharacteristics_sessionInfo_ ", dataType, ".txt"))
print(session)
sink()