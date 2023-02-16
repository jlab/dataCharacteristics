#  From 5-number summary (min, lower quartile, median, higher quartile, max) get
# the third value corresponding to median value
getThirdNumber <- function(str1){
  as.numeric(unlist(strsplit(str1, ",")))[3]
}

writeMedianRNAseqFpkmsTpmsFiles <- function(filePath, rowLabelCol, colsToRemove, zerosToNA = FALSE) {
  print(filePath)
  dat <- read.csv(filePath, check.names = FALSE, sep = "\t")
  geneId <- dat[[rowLabelCol]]
  dat <- data.frame(dat[,setdiff(colnames(dat), colsToRemove)])
  row.names(dat) <- geneId
  
  dat <- subset(dat, rowSums(dat == "0,0,0,0,0") != ncol(dat))
  
  # dat <- data.frame(dat[rowSums(dat == "0,0,0,0,0") != ncol(dat), ])
  geneId <- row.names(dat)
  dat.lst <- as.list(dat)
  
  # Take third value (= median of 5-number summary) of the comma separated numbers in each cell
  # The same seems to have been done to get the values in the FPKM file on the 
  # Downloads tab of the ExpressionAtlas website 
  # (see e.g. https://www.ebi.ac.uk/gxa/experiments/E-MTAB-2812/Downloads)
  dat.lst <- rapply(
    dat.lst, 
    function(x) lapply(x, getThirdNumber), 
    how = "replace")
  
  mtx <- t(matrix(unlist(dat.lst), nrow=length(dat.lst),byrow=TRUE))
  colnames(mtx) <- colnames(dat)
  row.names(mtx) <- make.names(geneId, unique=TRUE)
  
  if (zerosToNA) mtx[mtx == 0] <- NA
  mtx[mtx == "NaN"] <- NA
  
  mtx <- subset(mtx, rowSums(is.na(mtx)) != ncol(mtx))
  
  fileName <- paste0(sub('\\.tsv$', '', basename(filePath)), "_median.tsv")
  write.table(mtx, file = fileName, sep = "\t", col.names=NA)
}

writeMedianFiles <- function(dataTypePath, rowLabelCol, colsToRemove) {
  dataTypeFilePaths <- list.files(dataTypePath, full.names = TRUE)
  dataTypeFilePaths <- dataTypeFilePaths[!grepl("\\.txt$", dataTypeFilePaths)]
  
  for (dataTypeFilePath in dataTypeFilePaths){
    writeMedianRNAseqFpkmsTpmsFiles(filePath = dataTypeFilePath,
                                           rowLabelCol = rowLabelCol,
                                           colsToRemove = colsToRemove,
                                           zerosToNA = FALSE)
  }
}

writeMedianFiles(dataTypePath = "RNAseq_fpkms", rowLabelCol = "GeneID", colsToRemove = c("GeneID", "Gene Name")) 
writeMedianFiles(dataTypePath = "RNAseq_tpms", rowLabelCol = "GeneID", colsToRemove = c("GeneID", "Gene Name")) 

session <- sessionInfo()
sink("writeFpkmsTpmsMedianFiles_sessionInfo.txt")
print(session)
sink()
