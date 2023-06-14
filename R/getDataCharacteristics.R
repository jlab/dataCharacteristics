
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# All datasets of Metabolights project MTBLS1541
# (m_MTBLS1541_AF_NMR_metabolite_profiling_v2_maf.tsv , m_MTBLS1541_Plasma_NMR_metabolite_profiling_v2_maf.tsv , m_MTBLS1541_Urine_NMR_metabolite_profiling_v2_maf.tsv)
# are excluded because it was unclear what the quantitative data corresponds to

# m_MTBLS1497_Serum_NMR_metabolite_profiling_v2_maf.tsv and m_MTBLS1497_Urine_NMR_metabolite_profiling_v2_maf.tsv have been changed manually because of 
# meta data that had to be removed from numeric data

##########################################################################################
library(data.table)
library(pcaMethods)
library(Matrix)
library(dplyr)
# library(foreach)
# library(BimodalIndex)
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

# # Kolmogorov-Smirnov Test, 
# # Returns percentage of significant samples when KS test was applied on single samples compared to all samples combined
# kolSmirTestSignProp <- function(mtx) {
#   pvals.mtx <- apply(mtx, 2, function(x) stats::ks.test(x, mtx)$p.value)
#   cnt <- sum(pvals.mtx<0.05)
#   signProportion <- cnt/length(pvals.mtx)
#   return(signProportion)
# }

# calc_ functions are from radiomics package
calc_energy <- function(data){
  #TODO: Add dim check for 2D vs 3D
  return(sum(as.numeric(data)*as.numeric(data), na.rm=TRUE))
}

#' @describeIn first_order_features Entropy
#' @param base The base for which the logarithm is calculate
#' @param nbins The number of bins the histogram is discretized into
calc_entropy <- function(data, base=2, nbins=length(unique(c(data)))){
  # Break data into a hist
  im_range <- range(data, na.rm=TRUE)
  cuts <- table(cut(data, seq(im_range[1], im_range[2], by=diff(im_range)/nbins), include.lowest=TRUE))/length(data[!is.na(data)])
  
  #Logs cannot take 0 values, so let = 0 if no value
  entropy_vals <- vapply(cuts, function(data) ifelse(data != 0, data*logb(data, base=base), 0), FUN.VALUE = 1)
  return(-1*sum(entropy_vals))
}

calc_kurtosis <- function(data){
  n <- length(data[!is.na(data)])
  data <- data - mean(data, na.rm=TRUE)
  r <- n * sum(data^4, na.rm=TRUE) / (sum(data^2, na.rm=TRUE)^2)
  return(r * (1 - 1/n)^2 - 3)
}

calc_meanDeviation <- function(data){
  #scale <- 1/prod(dim(data))
  scale <- 1/length(data)
  mu <- mean(data, na.rm=TRUE)
  return(scale * sum(abs(data - mu), na.rm=TRUE))
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

calc_uniformity <- function(data, nbins=length(unique(c(data)))){
  # Break data into a hist
  data <- data[!is.na(data)]
  im_range <- range(data, na.rm=TRUE)
  cuts <- table(cut(data, seq(im_range[1], im_range[2], by=diff(im_range)/nbins), include.lowest=TRUE))/length(data)
  function_vals <- vapply(cuts, function(data) data^2, FUN.VALUE = 1)
  return(sum(function_vals))
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

# Pairwise correlation of rows
# For more than nmaxFeature feature a subset of nmaxFeatu random features is selected to speed up runtime. 
get_rowCorr <- function(mtx, nmaxFeature=100, corMethod = "spearman", ...){
  res <- seedUsed <- NA
  
  if (nrow(mtx)>nmaxFeature){ # random subset of features are selected
    randomRows <- applyFunctionWithSeed(sample, x = 1:nrow(mtx), size= min(nmaxFeature, nrow(mtx)))
    seedUsed <- randomRows$seed
    mtx <- removeEmptyRowsAndColumns(mtx[randomRows$res,], zerosToNA = TRUE)
  } 
  
  res <- cor(t(mtx), method = corMethod, use = "pairwise.complete.obs")

  return(list(res = res, seed = seedUsed))  
}

# pairwise correlation of columns
get_colCorr <- function(mtx, nmaxSamples=100, corMethod = "spearman", ...){
  res <- seedUsed <- NA
  
  if (ncol(mtx) > nmaxSamples){ # random subset of features are selected
    randomCols <- applyFunctionWithSeed(sample, x = 1:ncol(mtx), size= min(nmaxSamples, ncol(mtx)))
    seedUsed <- randomCols$seed
    mtx <- removeEmptyRowsAndColumns(mtx[, randomCols$res], zerosToNA = TRUE)
  } 
  res <- cor(mtx, method = corMethod, use = "pairwise.complete.obs")
  return(list(res = res, seed = seedUsed))  
}

# Adjusted from BimodalIndex::bimodalIndex such that itmax can befined
bimodalIndexWithIterDef <- function (dataset, verbose = TRUE, itmax = 10000, ...) {
  bim <- matrix(NA, nrow = nrow(dataset), ncol = 6)
  if (verbose) 
    cat("1 ")
  for (i in 1:nrow(dataset)) {
    if (verbose && 0 == i%%100) 
      cat(".")
    if (verbose && 0 == i%%1000) 
      cat(paste("\n", 1 + i/1000, " ", sep = ""))
    x <- as.vector(as.matrix(dataset[i, ]))
    if (any(is.na(x))) 
      next
    mc <- mclust::Mclust(x, G = 2, modelNames = "E", verbose = FALSE, control = mclust::emControl(itmax = itmax))
    sigma <- sqrt(mc$parameters$variance$sigmasq)
    delta <- abs(diff(mc$parameters$mean))/sigma
    pi <- mc$parameters$pro[1]
    bi <- delta * sqrt(pi * (1 - pi))
    bim[i, ] <- c(mc$parameters$mean, sigma = sigma, delta = delta, 
                  pi = pi, bim = bi)
  }
  if (verbose) 
    cat("\n")
  dimnames(bim) <- list(rownames(dataset), c("mu1", "mu2", 
                                             "sigma", "delta", "pi", "BI"))
  bim <- as.data.frame(bim)
  bim
}

# bimodalIndex 
# bimodality of row correlations
get_bimodalityRowCorr <- function(mtx, naToZero = FALSE, ...) {
  if (naToZero) mtx[is.na(mtx)] <- 0
  corrRes <- get_rowCorr(mtx)
  res <- bimodalIndexWithIterDef(matrix(corrRes$res, nrow=1), verbose=F)$BI
  return(list(res = res, seed = corrRes$seed))
}

# bimodality of column correlations
#get_bimodalityColCorr <- function(mtx, ...) return(BimodalIndex::bimodalIndex(matrix(get_colCorr(mtx),nrow=1),verbose=F)$BI)
get_bimodalityColCorr <- function(mtx, naToZero = FALSE, ...) {
  if (naToZero) mtx[is.na(mtx)] <- 0
  corrRes <- get_colCorr(mtx)
  res <- bimodalIndexWithIterDef(matrix(corrRes$res, nrow=1), verbose=F)$BI
  return(list(res = res, seed = corrRes$seed))
}


# Poly2 (features)
# linear coefficient of 2nd order polynomial fit with x = row means, y = row variances
get_LinearCoefPoly2XRowMeansYRowVars <- function(mtx, ...){
  # return(unname(lm(y~x+I(x^2), data=data.frame(
  #   y=get_rowSd(mtx)^2,
  #   x=get_rowMeans(mtx)))$coefficients[2]))
  return(unname(coef(biglm::biglm(y~x+I(x^2), data=data.frame(
    y=get_rowSd(mtx)^2,
    x=get_rowMeans(mtx))))[2]))
}

# quadratic coefficient of 2rd order polynomial fit with x = row means, y = row variances
get_QuadraticCoefPoly2XRowMeansYRowVars <- function(mtx, ...){ 
  # return(unname(lm(y~x+I(x^2), data=data.frame(
  #   y=get_rowSd(mtx)^2,
  #   x=get_rowMeans(mtx)))$coefficients[3]))
  return(unname(coef(biglm::biglm(y~x+I(x^2), data=data.frame(
    y=get_rowSd(mtx)^2,
    x=get_rowMeans(mtx))))[3]))
}

# linear and quadratic coefficient of 2rd order polynomial fit with x = row means, y = row variances
get_CoefPoly2XRowMeansYRowVars <- function(mtx, ...){ 
  # coefs <- lm(y~x+I(x^2), data=data.frame(
  #   y=get_rowSd(mtx)^2,
  #   x=get_rowMeans(mtx)))$coefficients
  
  # model <- biglm::biglm(y~x+I(x^2), data=data.frame(
  #   y=get_rowSd(mtx)^2,
  #   x=get_rowMeans(mtx)))
  
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
get_coefHclustRows <- function(mtx, naToZero = FALSE, ...) {
  if (naToZero) mtx[is.na(mtx)] <- 0
  mtx.tmp <- mtx[order(rowSums(is.na(mtx))), ]
  # mtx.tmp <- mtx.tmp[1:min(500,dim(mtx)[1]),1:min(500,dim(mtx)[2])] # max. 500 features and samples
  randomCols <- applyFunctionWithSeed(sample, x = 1:ncol(mtx), size= min(500, ncol(mtx)))
  seed <- randomCols$seed
  randomRows <- applyFunctionWithSeed(sample, x = 1:nrow(mtx), size= min(500, nrow(mtx)), seed = seed)
  mtx.tmp <- mtx.tmp[randomRows$res, randomCols$res] # max. 500 features and samples
  
  mtx.tmp <- removeEmptyRowsAndColumns(mtx.tmp, zerosToNA = FALSE)
  return(list(res = cluster::coef.hclust(amap::hcluster(mtx.tmp)), seed = seed))
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
    glmModelSample <- speedglm::speedglm(isNA ~ imputed, 
                                         family=binomial(link='logit'), 
                                         data = data.long.Sample, fitted = TRUE)
    
    intensities <- seq(-20, max(data.long$imputed), length.out=200)
    approx(predict(glmModelSample, 
                   newdata=data.frame(imputed=intensities),
                   type="response"), intensities, c(0.5, 0.9))$y
  })
  
  x5090.lst2 <- list()
  for (sample in unique(data.long$Sample)){
    x5090 <- NULL
    data.long.Sample <- data.long %>% dplyr::filter(Sample == sample)
    glmModelSample <- speedglm::speedglm(isNA ~ imputed,
                                         family=binomial(link='logit'),
                                         data = data.long.Sample, fitted = TRUE)

    intensities <- seq(-20, max(data.long$imputed), length.out=200)
    x5090 <- approx(predict(glmModelSample, newdata=data.frame(imputed=intensities),type="response"),intensities, c(0.5, 0.9))$y
    x5090.lst2 <- append(x5090.lst2, list(x5090))
  }
  
  x5090.df <- do.call(rbind, x5090.lst)
  colnames(x5090.df) <- c("IntensityNAp50", "IntensityNAp90")
  x5090.sds <- apply(x5090.df, 2, calc_sd)
  
  list(IntensityNAp50 = x5090.sds[["IntensityNAp50"]], 
       IntensityNAp90 = x5090.sds[["IntensityNAp90"]],
       seed = seedUsed
  )
}


getCharacteristicsHelper <- function(mtx, fast = TRUE){
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
    prctPC1 <- prctPC2 <- bimodalityRowCorr <- bimodalityRowCorrSeed <- bimodalityColCorr <- bimodalityColCorrSeed <- 
    linearCoefPoly2Row <- quadraticCoefPoly2Row <- 
    coefHclustRows <- coefHclustRowsSeed <- 
    intensityNAProbability50.sd <- intensityNAProbability90.sd <- intensityNAProbabilitySeed <- NA
  
  try({medianSampleVariance <- median(apply(mtx, 2, calc_variance), na.rm = TRUE)})
  try({medianAnalyteVariance <- median(unname(apply(mtx, 1, calc_variance)), na.rm = TRUE)})
  
  try({skewness <- calc_skewness(mtx)})
  try({kurtosis <- calc_kurtosis(mtx)})
  
  if (!fast){
    entropy <- calc_entropy(mtx) #
    meanDeviation <- calc_meanDeviation(mtx)
    uniformity <- calc_uniformity(mtx) #
    RMS <- calc_RMS(mtx)
  }
  
  try({variance <- calc_variance(mtx)})
  
  # group.size <- ncol(mtx)/2
  
  # var.groups.ratio <- median(matrixStats::rowVars(mtx[, 1:group.size], na.rm = TRUE)/matrixStats::rowVars(mtx[, (group.size+1):ncol(mtx)], na.rm = TRUE), na.rm = TRUE)
  
  # bimodalIndex 
  try({
    bimodalityRowCorrRes <- get_bimodalityRowCorr(mtx, naToZero = TRUE)
    bimodalityRowCorr <- bimodalityRowCorrRes$res
    bimodalityRowCorrSeed <- bimodalityRowCorrRes$seed
  })
  
  gc()
  
  try({
    bimodalityColCorrRes <- get_bimodalityColCorr(mtx, naToZero = TRUE)
    bimodalityColCorr <- bimodalityColCorrRes$res
    bimodalityColCorrSeed <- bimodalityColCorrRes$seed
  })
  
  gc()
  
  # Poly2 (features)
  # try(linearCoefPoly2Row <- get_LinearCoefPoly2XRowMeansYRowVars(mtx))
  # try(quadraticCoefPoly2Row <- get_QuadraticCoefPoly2XRowMeansYRowVars(mtx))
  try({
    coefs <- get_CoefPoly2XRowMeansYRowVars(mtx)
    linearCoefPoly2Row <- coefs[["linearCoef"]]
    quadraticCoefPoly2Row <- coefs[["quadraticCoef"]]
  })

  gc()
  
  # Coef.hclust (features)
  try({
    coefHclustRowsRes <- get_coefHclustRows(mtx, naToZero = TRUE)
    coefHclustRows <- coefHclustRowsRes$res
    coefHclustRowsSeed <- coefHclustRowsRes$seed
  })
  
  gc()
  
  
  try({
    x5090.sds <- calculateIntensityNAProbability5090(mtx)
    intensityNAProbability50.sd <- x5090.sds[["IntensityNAp50"]]
    intensityNAProbability90.sd <- x5090.sds[["IntensityNAp90"]]
    intensityNAProbabilitySeed <- x5090.sds[["seed"]]
  })
  
  gc()
  
  
  mtx <- mtx %>% t()
  mtx <- mtx[ , which(apply(mtx, 2, calc_variance) != 0)] # Remove zero variance columns 
  
  if (!is.vector(mtx)){
    try({
      pca <- pcaMethods::pca(mtx, method="nipals", center = TRUE, maxSteps=5000)
      prctPC1 <- pca@R2[1]
      prctPC2 <- pca@R2[2]
    })
  }
  
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
    bimodalityRowCorr = bimodalityRowCorr,
    bimodalityRowCorrSeed = bimodalityRowCorrSeed,
    bimodalityColCorr = bimodalityColCorr,
    bimodalityColCorrSeed = bimodalityColCorrSeed,
    linearCoefPoly2Row = linearCoefPoly2Row,
    quadraticCoefPoly2Row = quadraticCoefPoly2Row,
    coefHclustRows = coefHclustRows,
    coefHclustRowsSeed = coefHclustRowsSeed,
    intensityNAProbability50.sd = intensityNAProbability50.sd,
    intensityNAProbability90.sd = intensityNAProbability90.sd,
    intensityNAProbabilitySeed = intensityNAProbabilitySeed
  )
  
  if (!fast){
    resultvec <- c(resultvec, 
                   entropy = entropy,
                   meanDeviation = meanDeviation,
                   uniformity = uniformity,
                   RMS = RMS)
    
  }
  
  return(resultvec)
}

getNaFeatures <- function(mtx) {
  colNaPercentage <- colMeans(is.na(mtx))*100
  rowNaPercentage <- rowMeans(is.na(mtx))*100
  rowNonNaNumber <- rowSums(!is.na(mtx))
  
  colMeans <- colMeans(mtx, na.rm = TRUE)
  rowMeans <- rowMeans(mtx, na.rm = TRUE)
  
  corSampleMeanNAPval <- corAnalyteMeanNAPval <- corSampleMeanNA <- corAnalyteMeanNA <- NA
  corCoefType <- "spearman"
  try({
    cortestCol <- cor.test(colNaPercentage, colMeans, method = corCoefType)
    corSampleMeanNAPval <- cortestCol$p.value
    corSampleMeanNA <- unname(cortestCol$estimate)
  })
  
  try({
    cortestRow <- cor.test(rowNaPercentage, rowMeans, method = corCoefType)
    corAnalyteMeanNAPval <- cortestRow$p.value
    corAnalyteMeanNA <- unname(cortestRow$estimate)
  })
  
  c(
    minRowNonNaNumber = min(rowNonNaNumber),
    maxRowNonNaNumber = max(rowNonNaNumber),
    minRowNaPercentage = min(rowNaPercentage),            
    maxRowNaPercentage = max(rowNaPercentage),
    minColNaPercentage = min(colNaPercentage),
    maxColNaPercentage = max(colNaPercentage),
    percNATotal = mean(is.na(mtx)) * 100,
    percOfRowsWithNAs = sum(apply(mtx, 1, anyNA))/nrow(mtx) * 100,
    corSampleMeanNA = corSampleMeanNA,
    corSampleMeanNAPval = corSampleMeanNAPval,
    corAnalyteMeanNA = corAnalyteMeanNA,
    corAnalyteMeanNAPval = corAnalyteMeanNAPval
  )
}

getDataCharacteristicsLogNoLog <- function(mtx, takeLog2 = FALSE, fast = TRUE) {
  
  nNegativeNumbers <- sum(mtx < 0, na.rm = TRUE)
  if (takeLog2) mtx <- log2(mtx)
  mtx[mtx == "NaN"] <- NA
  
  nDistinctValues <- length(unique(c(mtx[!is.na(mtx)])))
  naFeatures <- getNaFeatures(mtx)
  # if (!is.vector(mtx)) mtx <- mtx[, colSums(is.na(mtx) | mtx == 0) != nrow(mtx)]
  
  
  # Positive correlation of pride datasets only due to LFQ?
  # Compare correlation when taking "Intensity " and "LFQ " columns for PXD020490/proteinGroups.txt
  #KS.SignProp <- kolSmirTestSignProp(as.matrix(df))
  
  characts.wNAs <- getCharacteristicsHelper(mtx, fast = fast)
  names(characts.wNAs) <- paste0(names(characts.wNAs), ".wNAs")
  
  #   # Remove rows with NAs
  #   mtx <- mtx[rowSums(is.na(mtx)) == 0, ]
  #   
  #   if (is.vector(mtx)){
  #     nFeatures.woNAs <- 1
  #   } else {
  #     nFeatures.woNAs <- nrow(mtx) # number of proteins with no NAs
  #   }
  #   
  #   runCalculation <- TRUE
  #   if (is.vector(mtx) & length(mtx)>1){
  #     runCalculation <- FALSE
  #   } else {
  #     if (!is.vector(mtx) &  nFeatures.woNAs == 0){
  #       runCalculation <- FALSE
  #     } 
  #   }
  #   
  #   if (runCalculation) {
  #     characts.woNAs <- getCharacteristicsHelper(mtx, withNAs=FALSE)
  #   } else {
  #     characts.woNAs <- c(
  #       mean = NA,
  #       median = NA,
  #       min = NA,
  #       max = NA,
  #       medianSampleVariance = NA,
  #       medianAnalyteVariance = NA,
  #       entropy = NA,
  #       kurtosis = NA, 
  #       meanDeviation = NA,
  #       skewness = NA,
  #       uniformity = NA,
  #       variance = NA, 
  #       RMS = NA,
  #       prctPC1 = NA, 
  #       prctPC2 = NA)
  #   }
  #   
  #   names(characts.woNAs) <- paste0(names(characts.woNAs), ".woNAs")
  
  # c(naFeatures, characts.wNAs, characts.woNAs)
  c(naFeatures, characts.wNAs, nDistinctValues = nDistinctValues, nNegativeNumbers = nNegativeNumbers)
}

getDataCharacteristics <- function(mtx, datasetID="test", dataType="test", ignoreRDS=FALSE) {
  resultName <- paste0(dataType,"__",datasetID)
  
  dir.create("Results",showWarnings = F)
  rdsName <- paste0("Results/",resultName,".RDS")

  ## Load if already calcuated  
  if(file.exists(rdsName) && !ignoreRDS){
    lst <- readRDS(file=rdsName)
  
    
  ## Calculate and save if already calcuated:
  }else{
    mtx[mtx == 0] <- NA
    mtx[mtx == Inf] <- NA
    mtx <- mtx[, colSums(is.na(mtx)) != nrow(mtx)]
    
    nSamples <- ncol(mtx)
    nAnalytes <- nrow(mtx)
    
    # charact.noLog <- getDataCharacteristicsLogNoLog(mtx, takeLog2 = FALSE)
    # names(charact.noLog) <- paste0(names(charact.noLog), ".noLog2")
    charact.log <- getDataCharacteristicsLogNoLog(mtx, takeLog2 = TRUE)
    names(charact.log) <- paste0(names(charact.log), ".log2")
    
    lst <- list(c(list(datasetID = datasetID, dataType = gsub("^\\./", "", dataType)), 
                  c(nSamples = nSamples, 
                    nAnalytes = nAnalytes,
                    # charact.noLog, 
                    charact.log)))
    saveRDS(lst,file=rdsName)

    session <- sessionInfo()
    sink(sub(".RDS",".log",rdsName))
    print(session)
    sink()
  }
  
  return(lst)
}
################################################################################


readInMetabolightsFiles <- function(filePath, zerosToNA = FALSE) {
  dat <- read.csv(filePath, check.names = FALSE, sep = "\t")
  
  # dat <- data.frame(data.table::fread(filePath, check.names = FALSE, sep = "\t", 
  #                                     integer64 = "character"), 
  #                   check.names = FALSE)
  
  remove <- c(
    "database_identifier",
    "chemical_formula",
    "smiles",
    "inchi",
    "metabolite_identification",
    "mass_to_charge",
    "fragmentation",
    "modifications",
    "charge",
    "retention_time",
    "retention time",
    "taxid",
    "species",
    "database",
    "database_version",
    "uri",
    "search_engine",
    "search_engine_score",
    "smallmolecule_abundance_sub",
    "smallmolecule_abu0ance_sub",
    "smallmolecule_abundance_stdev_sub",
    "smallmolecule_abu0ance_stdev_sub",
    "smallmolecule_abundance_std_error_sub",
    "smallmolecule_abu0ance_std_error_sub",
    "pattern ID",
    "pattern_ID",
    "chemical_shift",
    "chemical_shift F2 (ppm)",
    "chemical_shift F1 (ppm)",
    "multiplicity",
    "reliability",
    "Reliability",
    "Functional group",
    "Molecule source",
    "KEGG_database_identifier",
    "row_id", 
    "pubchem_first_synonym", 
    "final_external_id", 
    "pubchem_cid", 
    "pubchem_cid_ik", 
    "csid_ik", 
    "final_smiles", 
    "final_inchi", 
    "final_inchi_key", 
    "search_type", 
    "ID", 
    "NAME", 
    "IUPAC_NAME", 
    "DATABASE_ACCESSION", 
    "pubchem_smiles", 
    "pubchem_inchi", 
    "pubchem_inchi_key", 
    "pubchem_formula", 
    "unichem_id", 
    "dime_id",
    "FDR",
    "p value",
    "Adjusted p value"
  )
  
  metabolite_identification <- dat$metabolite_identification
  dat <- dat[, colnames(dat) %notin% remove]
  mtx <- as.matrix(dat)
  
  try({
    if(sum(metabolite_identification != "") == nrow(dat)) 
      row.names(mtx) <- make.names(metabolite_identification, unique=TRUE)
  })
  
  dat <- NULL
  
  mtx[mtx == "BLQ"] <- 0
  
  mtx <- removeEmptyRowsAndColumns(mtx, zerosToNA = zerosToNA)
  
  mtx <- apply(mtx, 2, gsub, pattern=",", replacement=".")
  class(mtx) <- "numeric"
  
  mtx <- removeEmptyRowsAndColumns(mtx, zerosToNA = zerosToNA)
  mtx
}

readInAllMetabolightsFiles <- function(dataTypePath, lst = list(), zerosToNA = FALSE) {
  dataTypeFilePaths <- list.files(dataTypePath, full.names = TRUE)
  for (dataTypeFilePath in dataTypeFilePaths){
    files <- list.files(dataTypeFilePath, recursive = TRUE, full.names = TRUE)
    # read in all files in subfolder, name dataset according to folder name plus filename
    for (file in files) {
      print(file)
      try({
        mtx <- readInMetabolightsFiles(filePath = file, zerosToNA = zerosToNA)
        # if (!is.vector(mtx)) {
        if (!is.vector(mtx) & ((sum(mtx < 0, na.rm = TRUE) / length(mtx)) < 0.01)) {
          if (nrow(mtx) > 9 & ncol(mtx) > 4) lst <- append(lst, getDataCharacteristics(mtx=mtx, 
                                                                                       datasetID=gsub(" ", "_", paste0(basename(dataTypeFilePath), "_", basename(file))), 
                                                                                       dataType=dataTypePath))
        }
      })
    }
  }
  lst
}

readInFile <- function(filePath, rowLabelCol, colsToRemove = c(), zerosToNA = FALSE, alternativeRowLabelCol = "") {
  dat <- read.csv(filePath, check.names = FALSE, sep = "\t")
  
  # dat <- data.frame(data.table::fread(filePath, check.names = FALSE, sep = "\t",
  #                                     integer64 = "numeric"), 
  #                   check.names = FALSE)
  
  if (is.character(rowLabelCol)){
    geneId <- dat[[rowLabelCol]]
  } else {
    if (rowLabelCol %notin% colnames(dat) & alternativeRowLabelCol != "") rowLabelCol <- alternativeRowLabelCol
    geneId <- dat[, rowLabelCol]
    colsToRemove <- c(colsToRemove, colnames(dat)[rowLabelCol])
  }
  
  dat <- data.frame(dat[, setdiff(colnames(dat), colsToRemove)])
  mtx <- as.matrix(dat)
  dat <- NULL
  
  row.names(mtx) <- make.names(unlist(geneId), unique=TRUE)
  mtx <- removeEmptyRowsAndColumns(mtx, zerosToNA = zerosToNA)
  mtx
}

readInAllDataTypeFiles <- function(dataTypePath, rowLabelCol, colsToRemove, zerosToNA = FALSE, lst = list(), alternativeRowLabelCol = "") {
  dataTypeFilePaths <- list.files(dataTypePath, full.names = TRUE)
  for (dataTypeFilePath in dataTypeFilePaths){
    print(dataTypeFilePath)
    try({
      mtx <- readInFile(filePath = dataTypeFilePath, 
                        rowLabelCol = rowLabelCol, 
                        colsToRemove = colsToRemove,
                        zerosToNA = zerosToNA,
                        alternativeRowLabelCol = alternativeRowLabelCol)
      # Only keep datasets with two dimensions
      if (!is.vector(mtx)) {
        if (nrow(mtx) > 9 & ncol(mtx) > 4)  lst <- append(lst, getDataCharacteristics(mtx=mtx, datasetID=gsub(" ", "_", basename(dataTypeFilePath)), dataType=dataTypePath))
      }
    })
  } 
  lst
}

readInMaxQuantFiles <- function (filePath, quantColPattern = c("^LFQ ", "^iBAQ ", "^Intensity "), zerosToNA = FALSE) {
  dat <- read.csv(filePath, allowEscapes = TRUE, check.names = FALSE,
                  sep = "\t")
  
  # mtx <- as.matrix(dat[, grepl("^LFQ ", names(dat))])
  
  if (quantColPattern == "^iBAQ "){
    # Because of PXD003834/proteinGroups_IF_151005_Volcanoplot_iBAQ_modSEENOTE_151022.txt
    cols.selected <- grepl(quantColPattern, names(dat)) & 
      !grepl("^iBAQ p value ", names(dat)) &
      !grepl("^iBAQ log2 ", names(dat))
  } else {
    cols.selected <- grepl(quantColPattern, names(dat))
  }
  mtx <- as.matrix(dat[, cols.selected])
  
  if (zerosToNA) mtx[mtx == 0] <- NA
  mtx[mtx == "NaN"] <- NA
  mtx <- `dimnames<-`(`dim<-`(as.numeric(mtx), dim(mtx)), dimnames(mtx))
  allColNA <- as.vector(apply(mtx, 1, function(r) {
    return(all(is.na(r)|r==0))
  }))
  message(paste("Number of proteins with empty entries:", length(which(allColNA))))
  if (!("Protein IDs" %in% colnames(dat))) {
    colnames(dat)[1] <- "Protein IDs"
  }
  featureAnnotations <- data.frame(proteinName = dat[, "Protein IDs"])
  annotations <- c("Fasta headers", "Q-value", "Reverse", "Peptides",
                   "Reverse", paste(c("Potential contaminant", "Contaminant"),
                                    collapse = "|"), "Only identified by site")
  fieldnames <- c("proteinDescription", "idScore", "isDecoy",
                  "nbPeptides", "isFiltered", "isPotential.contaminant",
                  "isIdentified.by.site")
  
  if (identical(grep(annotations[6], colnames(dat), value = TRUE), character(0))){
    bool1 <- logical(length = 0)
  } else {
    bool1 <- dat[, grep(annotations[6], colnames(dat), value = TRUE)] ==
      "+" & !is.na(dat[, grep(annotations[6], colnames(dat),
                              value = TRUE)])
  }
  
  bool2 <- dat[["Only identified by site"]] == "+" & !is.na(dat[["Only identified by site"]])
  bool3 <- dat[["Reverse"]] == "+" & !is.na(dat[["Reverse"]])
  ixs <- bool1 | bool2 | bool3
  if (identical(ixs, logical(0))){
    ixs <- rep(FALSE, nrow(dat))
  }
  
  for (i in seq_len(length(fieldnames))) {
    if (length(grep(annotations[i], colnames(dat))) > 0) {
      if (length(grep(fieldnames[i], c("isDecoy", "isPotential.contaminant",
                                       "Only identified by site"))) > 0) {
        featureAnnotations[[fieldnames[i]]] <- dat[,
                                                   grep(annotations[i], colnames(dat), value = TRUE)] ==
          "+"
      }
      else {
        featureAnnotations[[fieldnames[i]]] <- dat[,
                                                   annotations[i]]
      }
    }
  }
  row.names(mtx) <- dat[, "Protein IDs"]
  mtx <- subset(mtx, rowSums(is.na(mtx) | mtx == 0) != ncol(mtx))
  if (!is.vector(mtx)) mtx <- mtx[, colSums(is.na(mtx) | mtx == 0) != nrow(mtx)]
  return(mtx)
}


getDataCharacteristicsForDataType <- function(dataType) {
  print(dataType)
  
  lst <- list()
  path <- "./"
  dataTypePath <- paste0(path, dataType)
  if (dataType == "microarray"){
    lst <- readInAllDataTypeFiles(dataTypePath = dataTypePath, 
                                  rowLabelCol = "Gene ID", 
                                  colsToRemove = c("Gene ID", "Gene Name", "DesignElementAccession", "Design Element"),
                                  zerosToNA = TRUE,
                                  lst = lst)
  } else if (dataType == "RNAseq_raw"){
    lst <- readInAllDataTypeFiles(dataTypePath = dataTypePath, 
                                  rowLabelCol = "Gene ID", 
                                  colsToRemove = c("Gene ID", "Gene Name"),
                                  zerosToNA = TRUE,
                                  lst = lst)
  } else if (dataType == "RNAseq_raw_undecorated"){
    lst <- readInAllDataTypeFiles(dataTypePath = dataTypePath, 
                                  rowLabelCol = "Gene ID", 
                                  colsToRemove = c("Gene ID", "Gene"),
                                  zerosToNA = TRUE,
                                  lst = lst,
                                  alternativeRowLabelCol = "Gene")
  } else if (dataType == "RNAseq_transcripts_raw_undecorated"){
    lst <- readInAllDataTypeFiles(dataTypePath = dataTypePath, 
                                  rowLabelCol = "Transcript ID", 
                                  colsToRemove = c("Transcript ID"),
                                  zerosToNA = TRUE,
                                  lst = lst)
  } else if (dataType == "RNAseq_transcripts_tpms"){
    lst <- readInAllDataTypeFiles(dataTypePath = dataTypePath, 
                                  rowLabelCol = "GeneID", 
                                  colsToRemove = c("Gene ID", "Gene Name", "GeneID"),
                                  zerosToNA = TRUE,
                                  lst = lst)
  } else if (dataType %in% c("RNAseq_fpkms_median", "RNAseq_tpms_median", "microbiome")){
    lst <- readInAllDataTypeFiles(dataTypePath = dataTypePath, 
                                  rowLabelCol = 1, 
                                  colsToRemove = c(),
                                  zerosToNA = TRUE,
                                  lst = lst)
  } else if (dataType == "proteomics_expressionatlas"){
    dataTypeFilePaths <- list.files(dataTypePath, full.names = TRUE)
    for (dataTypeFilePath in dataTypeFilePaths){
      print(dataTypeFilePath)
      for (quantColPattern in c("^LFQ ", "^iBAQ ", "^Intensity ")){
        try({
          mtx <- readInMaxQuantFiles(filePath = dataTypeFilePath, 
                                     quantColPattern = quantColPattern,
                                     zerosToNA = TRUE)
          if (!is.vector(mtx)) {
            if (nrow(mtx) > 9 & ncol(mtx) > 4) {
              lst <- append(lst, 
                            getDataCharacteristics(
                              mtx=mtx, 
                              datasetID=gsub(" ", "_", paste0(basename(dataTypeFilePath), "_", quantColPattern)), 
                              dataType=dataTypePath))
            }
          }
        })
      }
    }
  } else if (dataType == "proteomics_pride"){
    dataTypeFilePaths <- list.files(dataTypePath, full.names = TRUE)
    for (dataTypeFilePath in dataTypeFilePaths){
      files <- list.files(dataTypeFilePath, recursive = TRUE, full.names = TRUE)
      # read in all files in subfolder, name dataset according to folder name plus filename
      for (file in files) {
        if (grepl("\\.txt$", file, ignore.case=TRUE)){ 
          print(file)
          for (quantColPattern in c("^LFQ ", "^iBAQ ", "^Intensity ")){
            try({
              mtx <- readInMaxQuantFiles(filePath = file,
                                         quantColPattern = quantColPattern,
                                         zerosToNA = TRUE)
              
              if (!is.vector(mtx)) {
                if (nrow(mtx) > 9 & ncol(mtx) > 4)  {
                  lst <- append(lst, 
                                getDataCharacteristics(
                                  mtx=mtx, 
                                  datasetID=gsub(" ", "_", paste0(basename(dataTypeFilePath), "_", basename(file), "_", quantColPattern)), 
                                  dataType=dataTypePath))
                }
              }
            })
          }
        }
      }
    }
  } else if (dataType %in% c("metabolomics_MS", "metabolomics_NMR")){
    lst <- readInAllMetabolightsFiles(dataTypePath, lst = lst,  zerosToNA = TRUE)
  } else if (dataType %in% c("sc_normalized", "sc_unnormalized")){
    dataTypeFilePaths <- list.files(dataTypePath, full.names = TRUE)
    for (dataTypeFilePath in dataTypeFilePaths){
      datasetID <- gsub(" ", "_", basename(dataTypeFilePath))
      rdsName <- paste0("Results/", dataTypePath, "__", datasetID, ".RDS")
      
      ## Load if already calcuated  
      if(!file.exists(rdsName)){
        print(dataTypeFilePath)
        ## Open MatrixMarket file
        # library(Matrix)
        mtx <- Matrix::readMM(dataTypeFilePath)
        mtx <- as.matrix(mtx)
        mtx <- removeEmptyRowsAndColumns(mtx, zerosToNA = TRUE)
        # Skip datasets if they contain more than 1% negative numbers
        if (!is.vector(mtx)) {
          if (nrow(mtx) > 9 & ncol(mtx) > 4) lst <- append(lst, getDataCharacteristics(mtx=mtx, datasetID=datasetID, dataType=dataTypePath))
        }
        mtx <- NULL
        
        gc()
      }
    }
  } else if (dataType == "scProteomics"){
    dataTypeFilePaths <- list.files(dataTypePath, full.names = TRUE)
    for (dataTypeFilePath in dataTypeFilePaths){
      print(dataTypeFilePath)
      mtx <- read.csv(dataTypeFilePath, row.names = 1)
      mtx <- as.matrix(mtx)
      mtx <- removeEmptyRowsAndColumns(mtx, zerosToNA = TRUE)
      if (!is.vector(mtx)) {
        if (nrow(mtx) > 9 & ncol(mtx) > 4) lst <- append(lst, getDataCharacteristics(mtx=mtx, datasetID=gsub(" ", "_", basename(dataTypeFilePath)), dataType=dataTypePath))
      }
    }
  } 
  
  lst.df <- rbindlist(lst, fill = TRUE)
  write.csv(lst.df, paste0("dataCharacteristics_", dataType, ".csv"), row.names = FALSE)
  
  session <- sessionInfo()
  sink(paste0("getDataCharacteristics_sessionInfo_ ", dataType, ".txt"))
  print(session)
  sink()
}

dataTypes <- c(
  "scProteomics",
  "metabolomics_NMR", "metabolomics_MS", 
  "proteomics_expressionatlas", "proteomics_pride",
  "microbiome",
  "RNAseq_fpkms_median", "RNAseq_tpms_median",
  "RNAseq_raw", "RNAseq_raw_undecorated", 
  # "RNAseq_transcripts_tpms",
  # "RNAseq_transcripts_raw_undecorated",
  "microarray",
  "sc_normalized",
  "sc_unnormalized"
)

# path <- "exampleFiles/"
#path <- "./"

# lst <- list()
for (dataType in dataTypes){
  getDataCharacteristicsForDataType(dataType)
}

# cl <- parallel::makeCluster(7, outfile="")
# doParallel::registerDoParallel(cl)
# 
# foreach(i = seq_along(dataTypes)) %dopar% {
#   getDataCharacteristicsForDataType(dataTypes[[i]])
# }
# parallel::stopCluster(cl)