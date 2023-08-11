library(dplyr)
library("jsonlite")

urls <- readLines("microarray_urls.txt")

chipIDs <- sub(".*/", "", urls)
chipIDs <- unique(sub(".*_(.*?)\\-normalized-expressions\\.tsv", "\\1", chipIDs))

#chipIDs <- paste0("A-", sub("microarray_", "", dataset %>% filter(dataType == "microarray") %>% select(dataType2) %>% unique() %>% as.vector() %>% unlist() %>% unname()))
modulos <- formatC(as.numeric(sub(".+-(.*)", "\\1", chipIDs)) %% 1000, width=3, flag="0")
fourLetterCodes <- sub("^(.*)-.*$", "\\1", chipIDs)
# jsons <- system(paste0('lftp ', "http://ftp.ebi.ac.uk/biostudies/nfs/A-AFFY-/",' <<<\'find| grep "\\\\.json$"; exit;\';'),intern=T);

ftpPath <- "http://ftp.ebi.ac.uk/biostudies/nfs/"
jsonPaths <- paste0(ftpPath, fourLetterCodes, "-/", modulos, "/", chipIDs, "/", chipIDs, ".json")

jsonLst <- list()
for (jsonPath in jsonPaths) {
  print(jsonPath)
  jsonContent <- jsonlite::fromJSON(jsonPath)
  jsonAttributes <- jsonContent[["section"]][["attributes"]]
  jsonAttributesVec <- with(jsonAttributes, setNames(value, name))
  
  Organism <- ifelse("Organism" %in% names(jsonAttributesVec), 
                     jsonAttributesVec[["Organism"]],
                     NA)
  chipID <- gsub(".*/(.*).json", "\\1", jsonPath)
  
  jsonLst <- append(jsonLst, list(c(chipID = chipID,
                                    Title = jsonAttributesVec[["Title"]], 
                                    Description = jsonAttributesVec[["Description"]],
                                    Organism = Organism
  )))
}

json.df <- do.call(rbind.data.frame, jsonLst)
colnames(json.df) <- c("chipID", "Title", "Description", "Organism")

json.df <- json.df %>% mutate(
  Manufacturer = case_when(
    grepl("Affymetrix", Title) ~ "Affymetrix",
    grepl("Agilent", Title) ~ "Agilent",
    grepl("Illumina", Title) ~ "Illumina"),
  Organism2 = case_when(
    grepl("hugenefl", tolower(Title)) | grepl("HT_HG-U133_Plus_PM", Title) | grepl("human", tolower(Title)) ~ "Human",
    grepl("chicken", tolower(Title)) ~ "Chicken",
    grepl("rice", tolower(Title)) ~ "Rice",
    grepl("bovine", tolower(Title)) ~ "Bovine",
    grepl("mouse", tolower(Title)) ~ "Mouse",
    grepl("poplar", tolower(Title)) ~ "Poplar",
    grepl("rhesus macaque", tolower(Title)) ~ "Rhesus Macaque",
    grepl("canine", tolower(Title)) ~ "Canine",
    grepl("rat", tolower(Title)) ~ "Rat",
    grepl("drosophila", tolower(Title)) ~ "Drosophila",
    grepl("arabidopsis", tolower(Title)) ~ "Arabidopsis",
    grepl("murine", tolower(Title)) ~ "Murine",
    grepl("barley", tolower(Title)) ~ "Barley",
    grepl("zebrafish", tolower(Title)) ~ "Zebrafish",
    grepl("yeast", tolower(Title)) ~ "Yeast",
    grepl("c. elegans", tolower(Title)) ~ "Caenorhabditis elegans",
    grepl("porcine", tolower(Title)) ~ "Porcine",
    grepl("maize", tolower(Title)) ~ "Maize",
    grepl("grape", tolower(Title)) ~ "Grape"
  )
)

write.csv(json.df, "microarray_metadata.csv", row.names = FALSE)

session <- sessionInfo()
sink(paste0("microarrayMetadata_sessionInfo.txt"))
print(session)
sink()
