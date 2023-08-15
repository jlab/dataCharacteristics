
library("jsonlite")

jsonContent <- jsonlite::fromJSON("https://www.ebi.ac.uk/metabolights/ws/studies/technology")

write.csv(jsonContent, "metabolomics_metadata.csv", row.names = FALSE)

session <- sessionInfo()
sink(paste0("metabolomicsMetadata_sessionInfo.txt"))
print(session)
sink()