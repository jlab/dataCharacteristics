


##################################################


`%notin%` <- Negate(`%in%`)

colMedians_drop0 <- function (dgCMat) {
  nnz_per_col <- diff(dgCMat@p)
  ColInd <- rep.int(1:ncol(dgCMat), nnz_per_col)
  sapply(split(dgCMat@x, ColInd), median)
}

plotCorrelations <- function(colStats, colNaPercentage, nrow, ncol, modus, pxd_id, corCoefType){
  statsString <- colStats[[1]]
  colStat <- colStats[[2]]
  
  tryCatch({
    corrPlot<- ggscatter(data.frame(colNaPercentage, colStat), x = "colNaPercentage", y = "colStat",
                         title=paste0(pxd_id, ", N(Samples)=", ncol, ", N(Features)=", nrow), 
                         xlab = "NA Percentage", 
                         ylab=statsString,
                         alpha = 0.5,
                         add = "reg.line",                                 # Add regression line
                         conf.int = TRUE,                                  # Add confidence interval
                         add.params = list(color = "blue",
                                           fill = "lightgray")) +
      stat_cor(method = corCoefType, label.x = 3) + # Add correlation coefficient
      xlim(0, 100) +
      theme(plot.title = element_text(margin = margin(10, 0, 10, 0)))
    
    ggsave(file=paste0("naMeanCorrelation_", modus, "_", statsString, "_", pxd_id, ".pdf"), plot = ggMarginal(corrPlot, type = "density"), height = 5, width = 6)
  },error = function(err){
    print(paste0("ERROR: ", pxd_id, ": ", err))
  })
}

naMeanCorrelation.pre <- function(index, pxd_ids, modus, path) {
  pxd_id <- pxd_ids[index]
  print(pxd_id)
  tryCatch({
    if (modus == "pride"){
      out <- readInMaxQuant(pxd_id, source.path = path)
      #if (pxd_id %notin% c("PXD004179", "PXD005477", "PXD009342")){
      # filter for potential contaminants and identified only by site features
      out <- out[!rowData(out)[["ixs"]],] 
      
      # extract data and feature annotation
      mtx <- assays(out)[["data"]]
      colnames(mtx) <- gsub("LFQ intensity","",colnames(mtx))
    } else if (modus == "sc"){
      
      # # Here, because it is a dgCMatrix, zeros can't be converted to NAs
      ds <- pxd_id[[2]]
      pxd_id <- pxd_id[[1]]
      # dgCMatrix
      mtx <- assay(ds)
      # remove rows with only 0s
      mtx <- mtx[rowSums(mtx)>0,]
      mtx <- as(mtx, "dgCMatrix")
    } else if (modus == "metabolights"){
      mtx <- getMetabolightsMatrix(pxd_id, source.path = path)
      mtx[mtx == 0] <- NA
      mtx[mtx == "NaN"] <- NA
      # Remove empty columns
      mtx <- mtx[, colSums( is.na(mtx) ) < nrow(mtx)]
    } 
    
    if (!is.null(ncol(mtx))){
      if (ncol(mtx) > 2){    
        # for proteomics data with NAs
        if (modus %in% c("pride", "metabolights", "mw")){
          ncol <- ncol(mtx)
          nrow <- nrow(mtx)
          colMeans <- colMeans(mtx, na.rm = TRUE)
          #colMedians <- colMedians(mtx, na.rm = TRUE)
          colNaPercentage <- colMeans(is.na(mtx))*100
          rowNaPercentage <- rowMeans(is.na(mtx))*100
          minRowNonNaNumber <- min(rowSums(!is.na(mtx)))
          maxRowNonNaNumber <- max(rowSums(!is.na(mtx)))
          numberMissingMatrix <- sum(is.na(mtx))
          
          #library(fBasics)
          #stats <- basicStats(mtx)
          
          # for scRNAseq data for which zeros correspiónd to missing values
        } else if (modus == "sc"){
          ncol <- ncol(mtx)
          nrow <- nrow(mtx)
          rowSumsNotZero <- Matrix::rowSums(mtx!=0) # mtx!=0 faster than mtx==0 due to sparsity of matrix
          rowSums <- rep(ncol, nrow)-rowSumsNotZero
          colNaPercentage <- (rep(nrow, ncol)-Matrix::colSums(mtx!=0))/nrow(mtx)*100
          rowNaPercentage <- (rep(ncol, nrow)-rowSumsNotZero)/nrow(mtx)*100
          minRowNonNaNumber <- min(rowSumsNotZero) 
          maxRowNonNaNumber <- max(rowSumsNotZero)
          numberMissingMatrix <- sum(rowSums) # as sum(mtx == 0) but faster
          colMeans <- colSums(mtx)/colSums(!!mtx)
          #colMedians <- colMedians_drop0(mtx)
          
          
          #library(fBasics)
          mtx <- as.matrix(mtx)
          mtx[mtx == 0] <- NA
          #stats <- basicStats(mtx)
        }
        
        minRowNaPercentage <- min(rowNaPercentage)            
        maxRowNaPercentage <- max(rowNaPercentage)
        minColNaPercentage <- min(colNaPercentage)
        maxColNaPercentage <- max(colNaPercentage)
        
        print(numberMissingMatrix)
        
        # stats <- data.frame(t(stats))
        # col1stQuartiles <- stats$X1..Quartile
        # col3rdQuartiles <- stats$X3..Quartile
        # colIQRs <- col3rdQuartiles - col1stQuartiles
        # colMeans <- stats$Mean
        # colMedians <- stats$Median
        # colSEMeans <- stats$SE.Mean
        # colLCLMeans <- stats$LCL.Mean
        # colUCLMeans <- stats$UCL.Mean
        # colVariances <- stats$Variance
        # colSkewness <- stats$Skewness
        # colKurtosis <- stats$Kurtosis
        
        # statsList <- list(
        #     list("1stQuartile", col1stQuartiles),
        #     list("3rdQuartile", col3rdQuartiles),
        #     list("IQR", colIQRs),
        #     list("Mean", colMeans),
        #     list("Median", colMedians),
        #     list("SEMean", colSEMeans),
        #     list("LCLMean", colLCLMeans),
        #     list("UCLMean", colUCLMeans),
        #     list("Variance", colVariances),
        #     list("Skewness", colSkewness),
        #     list("Kurtosis", colKurtosis))
        
        # mind 1 feature, das bei Samples x>50% (ausser bei sc) --> maxRowNaPercentage > 50
        # datensatz mit mind 10 samples --> ncol>9
        # mind 1 sample, bei dem mind 10% Nas --> maxColNaPercentage > 10
        
        #print(pxd_id)
        print(ncol)
        print(maxColNaPercentage)
        print(maxRowNaPercentage)
        
        corCoefType <- "spearman"
        
        cortest <- cor.test(colNaPercentage, colMeans, method = corCoefType)
        pval <- cortest$p.value
        cor <- cortest$estimate
        
        #if (ncol>9  &  maxColNaPercentage > 10 & (maxRowNaPercentage > 50 | modus == "sc")){
        if (numberMissingMatrix > 0){  
          corrPlot<- ggscatter(data.frame(colNaPercentage, colMeans), x = "colNaPercentage", y = "colMeans",
                               title=paste0(pxd_id, ", N(Samples)=", ncol, ", N(Features)=", nrow), 
                               xlab = "NA Percentage", 
                               ylab= "Mean",
                               add = "reg.line",                                 # Add regression line
                               conf.int = TRUE,                                  # Add confidence interval
                               add.params = list(color = "blue",
                                                 fill = "lightgray")) +
            stat_cor(method = corCoefType, label.x = 3) + # Add correlation coefficient
            xlim(0, 100) +
            theme(plot.title = element_text(margin = margin(10, 0, 10, 0)))
          
          ggsave(file=paste0("naMeanCorrelation_", modus, "_colMeans_", pxd_id, ".pdf"), plot = ggMarginal(corrPlot, type = "density"), height = 5, width = 6)
          
          
          data.chars <- getDataCharacteristics(as.data.frame(mtx))
          return(list(df=data.frame(
            ID = pxd_id,
            #sampleSize = ncol,
            #featureNumber = nrow,
            R = cor,
            p = pval,
            numberMissingMatrix = numberMissingMatrix,
            minRowNaPercentage = minRowNaPercentage,
            maxRowNaPercentage = maxRowNaPercentage,
            minRowNonNaNumber = minRowNonNaNumber,
            maxRowNonNaNumber = maxRowNonNaNumber,
            minColNaPercentage = minColNaPercentage,
            maxColNaPercentage = maxColNaPercentage,
            as.list(data.chars),
            stringsAsFactors=FALSE)))
        }
      }
    }
    
  },
  error = function(err){
    print(paste0("ERROR: ", pxd_id, ": ", err))
  })
}
########################################

setwd(('/Users/admin/Desktop/PhD/2022_OmicsCorrelationDataCharacteristics'))


datasources <- c("pride", "sc", "mw", "metabolights")

for (datasource in datasources){

  cl <- makeCluster(detectCores()/2)
  #sink('/Users/admin/Desktop/PhD/2022_OmicsCorrelationDataCharacteristics/analysis-output.txt')
  registerDoParallel(cl)
  naMeanCorrelation.pre.list <- foreach(i = 1:length(pxd_ids), .packages=c('SummarizedExperiment',
                                                                           'scRNAseq', 
                                                                           'lipidr',
                                                                           'jsonlite')) %dopar% {
                                                                             naMeanCorrelation.pre(index=i, pxd_ids = pxd_ids, modus=datasource, path = path)
                                                                           }
  stopCluster(cl)
  
  cor.df <- do.call( rbind, naMeanCorrelation.pre.list )[,1]
  cor.df <- do.call(rbind, cor.df)
  
  write.csv(cor.df, file=paste0("naMeanCorrelation_", datasource, ".csv"), row.names = FALSE)
  write.xlsx(cor.df,paste0("naMeanCorrelation_", datasource, ".xlsx"),
             col.names = TRUE, row.names = FALSE, keepNA = TRUE,
             na.string="NA")
}


####################
facToNum <- function(vec){
  vec <- as.numeric(as.character(vec))
  vec
}

plotPairs <- function(datasource){
  print(datasource)
  path <- paste0("/Users/Eva/Desktop/PhD/pride_2013to2018/20200518_NAs_mean_correlation_allDatasources_spearman/", datasource, "/naMeanCorrelation.csv")
  df <- read.csv(path)
  if (any(grepl("ERROR", df$ID))) {
    df <- df[- grep("ERROR", df$ID),] 
  }
  
  df[,-1] <- sapply(df[,-1], facToNum)
  df$ID <- as.character(df$ID)
  
  # max. Anzahl von NA >x% (der Samples)
  medianMaxRowNaPercentage <-  median(facToNum(df$maxRowNaPercentage))
  
  # mind >y Samples 
  medianMaxRowNonNaNumber <- median(facToNum(df$maxRowNonNaNumber))
  
  # >z% (der Features) --> result value x has to be seen as 100-x
  medianMinColNaPercentage <- median(facToNum(df$minColNaPercentage))
  
  pdf(file=paste0("pairs_", datasource, ".pdf"), height = 10, width = 14)
  pairs.panels(df[,-1], 
               main = paste0(datasource," medianMaxRowNaPercentage: ", round(medianMaxRowNaPercentage, 3), 
                             " medianMaxRowNonNaNumber: ", round(medianMaxRowNonNaNumber, 3), 
                             "  medianMinColNaPercentage: ", round(medianMinColNaPercentage, 3)),
               method = "pearson", # correlation method
               hist.col = "#00AFBB",
               density = TRUE,  # show density plots
               ellipses = TRUE # show correlation ellipses
  )
  dev.off()
  
  return(list(df=data.frame(
    datasource = datasource,
    medianSampleSize = median(facToNum(df$sampleSize)),
    medianFeatureNumber = median(facToNum(df$featureNumber)),
    medianR = median(facToNum(df$R), na.rm = TRUE),
    medianP = median(facToNum(df$p), na.rm = TRUE),
    medianNumberMissingMatrix = median(facToNum(df$numberMissingMatrix)),
    medianMinRowNaPercentage = median(facToNum(df$minRowNaPercentage)),
    medianMaxRowNaPercentage = medianMaxRowNaPercentage,
    medianMinRowNonNaNumber = median(facToNum(df$minRowNonNaNumber)),
    medianMaxRowNonNaNumber = medianMaxRowNonNaNumber,
    medianMinColNaPercentage = median(facToNum(df$minColNaPercentage)),
    medianMaxColNaPercentage =  median(facToNum(df$maxColNaPercentage)),
    stringsAsFactors=FALSE)))
  
}
pairs.list <- lapply(datasources, plotPairs)
pairs.df <- do.call(rbind, pairs.list )[,1]
pairs.df <- do.call(rbind, pairs.df)

# write.csv(pairs.df, file="pairs.csv", row.names = FALSE)
# library("openxlsx")
# write.xlsx(pairs.df, "pairs.xlsx",
#     col.names = TRUE, row.names = FALSE, keepNA = TRUE,
#     na.string="NA")

getNoOfDatasets <- function(datasource){
  print(datasource)
  path <- paste0("/Users/Eva/Desktop/PhD/pride_2013to2018/20200518_NAs_mean_correlation_allDatasources_spearman/", datasource, "/naMeanCorrelation.csv")
  df <- read.csv(path)
  if (any(grepl("ERROR", df$ID))) {
    df <- df[- grep("ERROR", df$ID),] 
  }
  df[,-1] <- sapply(df[,-1], facToNum)
  df$ID <- as.character(df$ID)
  print(nrow(df))
}

lapply(datasources, getNoOfDatasets)



