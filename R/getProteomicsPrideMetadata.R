
library(dplyr)
library(rpx)

urls <- readLines("proteomics_pride_urls.txt")

datasetPrideIDs <- unique(sub(".*/(.*?)/[^/]*$", "\\1", urls))
prideMetaLst <- list()
for (datasetPrideID in datasetPrideIDs) {
  
  try({
    px <- rpx::PXDataset(datasetPrideID)
    prideMetaLst <- append(prideMetaLst, 
                           list(c(
                             prideID = datasetPrideID,
                             tax = paste(px@px_tax, collapse = ";"),
                             organismParts = paste(px@px_metadata[["organismParts"]][["name"]], collapse = ";"), 
                             experimentTypes = paste(px@px_metadata[["experimentTypes"]][["name"]], collapse = ";"),
                             instruments = paste(px@px_metadata[["instruments"]][["name"]], collapse = ";")
                           )))
  })
}

prideMeta.df <- do.call(rbind.data.frame, prideMetaLst) 
colnames(prideMeta.df) <- c("prideID", "tax", "organismParts", "experimentTypes", "instruments")
uniqueInstruments <- sort(unique( unlist(strsplit(prideMeta.df$instruments, ";"))))

write.csv(prideMeta.df, "proteomicsPride_metadata.csv", row.names = FALSE)

session <- sessionInfo()
sink(paste0("proteomicsPrideMetadata_sessionInfo.txt"))
print(session)
sink()