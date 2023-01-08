
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##########################################################################################
# library(dplyr)
# library(matrixTests)
# library(qvalue)
# library(DescTools)
# library(matrixcalc)
# library(psych) 
# library(stats)
# library(foreach)
# library(rlist)

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

calc_skewness <- function (data){
  data <- data[!is.na(data)]
  return(sum((data - mean(data))^3)/(length(data) * sd(data)^3))
}

calc_uniformity <- function(data, nbins=length(unique(c(data)))){
  # Break data into a hist
  data <- data[!is.na(data)]
  im_range <- range(data, na.rm=TRUE)
  cuts <- table(cut(data, seq(im_range[1], im_range[2], by=diff(im_range)/nbins), include.lowest=TRUE))/length(data)
  function_vals <- vapply(cuts, function(data) data^2, FUN.VALUE = 1)
  return(sum(function_vals))
}

calc_variance <- function(data) var(c(data), na.rm=TRUE)
calc_RMS <- function(data) sqrt(mean(data^2, na.rm=TRUE))

getCharacteristicsHelper <- function(mtx, withNAs=TRUE){
  #KS.SignProp <- kolSmirTestSignProp(mtx)
  
  entropy <- calc_entropy(mtx)
  kurtosis <- calc_kurtosis(mtx)
  meanDeviation <- calc_meanDeviation(mtx)
  skewness <- calc_skewness(mtx)
  uniformity <- calc_uniformity(mtx)
  variance <- calc_variance(mtx)
  RMS <- calc_RMS(mtx)
  
  # group.size <- ncol(mtx)/2
  
  # var.groups.ratio <- median(matrixStats::rowVars(mtx[, 1:group.size], na.rm = TRUE)/matrixStats::rowVars(mtx[, (group.size+1):ncol(mtx)], na.rm = TRUE), na.rm = TRUE)
  
  if (withNAs){
    resultvec <- c(#KS.SignProp = KS.SignProp,
                   entropy = entropy,
                   kurtosis = kurtosis, 
                   meanDeviation = meanDeviation,
                   skewness = skewness,
                   uniformity = uniformity,
                   variance = variance, 
                   RMS = RMS#,
                   #var.groups.ratio = var.groups.ratio
    )
  } else {
    t.mtx <- t(mtx)
    t.mtx <- t.mtx[ , which(apply(t.mtx, 2, var) != 0)] # Remove zero variance columns 
    pca <- stats::prcomp(t.mtx, scale.=T)
    eigs <- pca$sdev^2
    prctPC1 <- eigs[1]/sum(eigs)
    prctPC2 <- eigs[2]/sum(eigs)
    
    elongation <- sqrt(eigs[2] / eigs[1]) # elongation
    flatness <- sqrt(eigs[length(eigs)]/eigs[1]) # flatness
    resultvec <- c(#KS.SignProp = KS.SignProp,
                   entropy = entropy,
                   kurtosis = kurtosis, 
                   meanDeviation = meanDeviation,
                   skewness = skewness,
                   uniformity = uniformity,
                   variance = variance, 
                   RMS = RMS,
                   #var.groups.ratio = var.groups.ratio,
                   prctPC1 = prctPC1, 
                   prctPC2 = prctPC2,
                   elongation = elongation,
                   flatness = flatness)
  }
  return(resultvec)
}

getDataCharacteristics <- function(mtx, datasetID="test", dataType="test") {
  mtx[mtx == 0] <- NA
  df <- data.frame(mtx)
  
  colMeans <- colMeans(mtx, na.rm = TRUE)
  #colMedians <- colMedians(mtx, na.rm = TRUE)
  colNaPercentage <- colMeans(is.na(mtx))*100
  rowNaPercentage <- rowMeans(is.na(mtx))*100
  minRowNonNaNumber <- min(rowSums(!is.na(mtx)))
  maxRowNonNaNumber <- max(rowSums(!is.na(mtx)))
  # numberMissingMatrix <- sum(is.na(mtx))
  minRowNaPercentage <- min(rowNaPercentage)            
  maxRowNaPercentage <- max(rowNaPercentage)
  minColNaPercentage <- min(colNaPercentage)
  maxColNaPercentage <- max(colNaPercentage)
  
  # Positive correlation of pride datasets only due to LFQ?
  # Compare correlation when taking "Intensity " and "LFQ " columns for PXD020490/proteinGroups.txt
  corCoefType <- "spearman"
  cortest <- cor.test(colNaPercentage, colMeans, method = corCoefType)
  corPval <- cortest$p.value
  corR <- cortest$estimate
  
  medianSampleVariance <- median(apply(df, 2, var,  na.rm=TRUE), na.rm = TRUE)
  medianProteinVariance <- median(unname(apply(df, 1, var,  na.rm=TRUE)), na.rm = TRUE)
  #KS.SignProp <- kolSmirTestSignProp(as.matrix(df))
  percNATotal <- mean(is.na(df)) * 100 
  percOfRowsWithNAs <- sum(apply(df, 1, anyNA))/nrow(df) * 100
  characts.wNAs <- getCharacteristicsHelper(mtx, withNAs=TRUE)
  names(characts.wNAs) <- paste0(names(characts.wNAs), ".wNAs")
  
  nSamples <- ncol(mtx)
  nFeatures.wNAs <- nrow(mtx)
  mtx <- mtx[rowSums(is.na(mtx)) == 0, ]
  nFeatures.woNAs <- nrow(mtx) # number of proteins with no NAs
  
  if (nFeatures.woNAs == 0){
    characts.woNAs <- c(# KS.SignProp = NA,
                        entropy = NA,
                        kurtosis = NA, 
                        meanDeviation = NA,
                        skewness = NA,
                        uniformity = NA,
                        variance = NA, 
                        RMS = NA)
  } else {
    characts.woNAs <- getCharacteristicsHelper(mtx, withNAs=FALSE)
  }

  names(characts.woNAs) <- paste0(names(characts.woNAs), ".woNAs")
  list(c(list(datasetID = datasetID, dataType = dataType), 
         c(nSamples = nSamples, nFeatures.wNAs = nFeatures.wNAs, nFeatures.woNAs = nFeatures.woNAs, 
           minRowNonNaNumber = minRowNonNaNumber, maxRowNonNaNumber = maxRowNonNaNumber,
           minRowNaPercentage = minRowNaPercentage, maxRowNaPercentage = maxRowNaPercentage, 
           minColNaPercentage = minColNaPercentage, maxColNaPercentage = maxColNaPercentage,
           corPval = corPval, corR = corR,
           medianSampleVariance = medianSampleVariance, medianProteinVariance = medianProteinVariance, 
           percNATotal = percNATotal, percOfRowsWithNAs = percOfRowsWithNAs, 
           characts.wNAs, characts.woNAs)))
}
################################################################################

`%notin%` <- Negate(`%in%`)

readInMetabolightsFiles <- function(filePath, zerosToNA = FALSE) {
  dat <- read.csv(filePath, allowEscapes = TRUE, check.names = FALSE, 
                  sep = "\t")
  
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
    "multiplicity",
    "reliability",
    "Functional group",
    "Molecule source"
    )
  
  dat2 <- dat[,colnames(dat) %notin% remove]

  mtx <- as.matrix(dat2)
  
  if(sum(dat$metabolite_identification != "") == nrow(dat)) 
    row.names(mtx) <- make.names(dat$metabolite_identification, unique=TRUE)
  
  if (zerosToNA) mtx[mtx == 0] <- NA
  mtx[mtx == "NaN"] <- NA
  # remove rows with only NAs
  # mtx <- mtx[rowSums(is.na(mtx) | mtx == 0) != ncol(mtx), ]
  mtx <- subset(mtx, rowSums(is.na(mtx) | mtx == 0) != ncol(mtx))
  if (!is.vector(mtx)) mtx <- mtx[, colSums(is.na(mtx) | mtx == 0) != nrow(mtx)]
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
        if (nrow(mtx) != 0) lst <- append(lst, getDataCharacteristics(mtx=mtx, 
                                                  datasetID=gsub(" ", "_", paste0(basename(dataTypeFilePath), "_", basename(file))), 
                                                  dataType=dataTypePath))
      })
    }
  }
  lst
}

readInFile <- function(filePath, rowLabelCol, colsToRemove = c(), zerosToNA = FALSE) {
  dat <- read.csv(filePath, check.names = FALSE, sep = "\t")
  if (is.character(rowLabelCol)){
    geneId <- dat[[rowLabelCol]]
  } else {
    geneId <- dat[, rowLabelCol]
    colsToRemove <- c(colsToRemove, colnames(dat)[rowLabelCol])
  }
    
  dat <- data.frame(dat[,setdiff(colnames(dat), colsToRemove)])
  mtx <- as.matrix(dat)
  row.names(mtx) <- make.names(geneId, unique=TRUE)
  
  if (zerosToNA) mtx[mtx == 0] <- NA
  mtx[mtx == "NaN"] <- NA
  # remove rows with only NAs
  # mtx <- mtx[rowSums(is.na(mtx) | mtx == 0) != ncol(mtx), ]
  mtx <- subset(mtx, rowSums(is.na(mtx) | mtx == 0) != ncol(mtx))
  if (!is.vector(mtx)) mtx <- mtx[, colSums(is.na(mtx) | mtx == 0) != nrow(mtx)]
  mtx
}

readInAllDataTypeFiles <- function(dataTypePath, rowLabelCol, colsToRemove, zerosToNA = FALSE, lst = list()) {
  dataTypeFilePaths <- list.files(dataTypePath, full.names = TRUE)
  for (dataTypeFilePath in dataTypeFilePaths){
    print(dataTypeFilePath)
    try({
      mtx <- readInFile(filePath = dataTypeFilePath, 
                        rowLabelCol = rowLabelCol, 
                        colsToRemove = colsToRemove,
                        zerosToNA = zerosToNA)
      # Only keep datasets with two dimensions
      if (!is.vector(mtx)) {
        lst <- append(lst, getDataCharacteristics(mtx=mtx, datasetID=gsub(" ", "_", basename(dataTypeFilePath)), dataType=dataTypePath))
      }
    })
  } 
  lst
}

readInMaxQuantFiles <- function (filePath, quantColPattern = c("^LFQ ", "^iBAQ ", "^Intensity "), zerosToNA = FALSE) {
  dat <- read.csv(filePath, allowEscapes = TRUE, check.names = FALSE,
                  sep = "\t")
  
  # mtx <- as.matrix(dat[, grepl("^LFQ ", names(dat))])
  mtx <- as.matrix(dat[, grepl(quantColPattern, names(dat))])
  
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

dataTypes <- c("metabolomics_MS", "metabolomics_NMR", "microarray",
               "proteomics_expressionatlas", "proteomics_pride",
               "RNAseq_fpkms_median",
               "RNAseq_raw", "RNAseq_raw_undecorated", "RNAseq_tpms_median",
               "RNAseq_transcripts_raw_undecorated", "RNAseq_transcripts_tpms",
               "sc_normalized", "sc_unnormalized", "microbiome")

path <- "exampleFiles/"

lst <- list()
for (dataType in dataTypes){
  print(dataType)
  dataTypePath <- paste0(path, dataType)
  if (dataType == "microarray"){
    lst <- readInAllDataTypeFiles(dataTypePath = dataTypePath, 
                           rowLabelCol = "Gene ID", 
                           colsToRemove = c("Gene ID", "Gene Name", "DesignElementAccession"),
                           lst = lst)
  } else if (dataType == "RNAseq_raw"){
    lst <- readInAllDataTypeFiles(dataTypePath = dataTypePath, 
                           rowLabelCol = "Gene ID", 
                           colsToRemove = c("Gene ID", "Gene Name"),
                           lst = lst)
  } else if (dataType == "RNAseq_raw_undecorated"){
    lst <- readInAllDataTypeFiles(dataTypePath = dataTypePath, 
                           rowLabelCol = "Gene ID", 
                           colsToRemove = c("Gene ID"),
                           lst = lst)
  } else if (dataType == "RNAseq_transcripts_raw_undecorated"){
    lst <- readInAllDataTypeFiles(dataTypePath = dataTypePath, 
                           rowLabelCol = "Transcript ID", 
                           colsToRemove = c("Transcript ID"),
                           lst = lst)
  } else if (dataType == "RNAseq_transcripts_tpms"){
    lst <- readInAllDataTypeFiles(dataTypePath = dataTypePath, 
                           rowLabelCol = "GeneID", 
                           colsToRemove = c("Gene ID", "Gene Name", "GeneID"),
                           lst = lst)
  } else if (dataType %in% c("RNAseq_fpkms_median", "RNAseq_tpms_median", "microbiome")){
    lst <- readInAllDataTypeFiles(dataTypePath = dataTypePath, 
                           rowLabelCol = 1, 
                           colsToRemove = c(),
                           lst = lst)
  } else if (dataType == "proteomics_expressionatlas"){
    dataTypeFilePaths <- list.files(dataTypePath, full.names = TRUE)
    for (dataTypeFilePath in dataTypeFilePaths){
      print(dataTypeFilePath)
      for (quantColPattern in c("^LFQ ", "^iBAQ ", "^Intensity ")){
        try({
          mtx <- readInMaxQuantFiles(filePath = dataTypeFilePath, 
                                     quantColPattern = quantColPattern,
                                     zerosToNA = FALSE)
          lst <- append(lst, 
                        getDataCharacteristics(
                          mtx=mtx, 
                          datasetID=gsub(" ", "_", paste0(basename(dataTypeFilePath), "_", quantColPattern)), 
                          dataType=dataTypePath))
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
                                         zerosToNA = FALSE)
              lst <- append(lst, 
                            getDataCharacteristics(
                              mtx=mtx, 
                              datasetID=gsub(" ", "_", paste0(basename(dataTypeFilePath), "_", basename(file), "_", quantColPattern)), 
                              dataType=dataTypePath))
            })
          }
        }
      }
    }
  } else if (dataType %in% c("metabolomics_MS", "metabolomics_NMR")){
    lst <- readInAllMetabolightsFiles(dataTypePath, lst = lst)
  } else if (dataType %in% c("sc_normalized", "sc_unnormalized")){
    dataTypeFilePaths <- list.files(dataTypePath, full.names = TRUE)
    for (dataTypeFilePath in dataTypeFilePaths){
      print(dataTypeFilePath)
      ## Open MatrixMarket file
      # library(Matrix)
      mtx <- Matrix::readMM(dataTypeFilePath)
      mtx <- as.matrix(mtx)
      lst <- append(lst, getDataCharacteristics(mtx=mtx, datasetID=gsub(" ", "_", basename(dataTypeFilePath)), dataType=dataTypePath))
    }
  } 
}

library(data.table)
lst.df <- rbindlist(lst, fill = TRUE)

write.csv(lst.df, "dataCharacteristics.csv", row.names = FALSE)
