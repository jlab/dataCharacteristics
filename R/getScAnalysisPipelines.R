urls <- readLines("sc_normalized_urls.txt")

pattern <- ".*/(.*?)/[^/]+$"  # Regular expression pattern to capture substring between last two slashes
projectIds <- sub(pattern, "\\1", urls)
alvinLst <- list()
for (projectId in projectIds){
  filePath <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments/", projectId, "/", projectId, ".software.tsv")
  fileStr <- paste(readLines(filePath), collapse="\n")
  processedWAlevin <- grepl("Alevin", fileStr)
  alvinLst <- append(alvinLst, list(c(projectId = projectId, processedWAlevin = processedWAlevin)))
} 

df <- data.frame(t(sapply(alvinLst,c)))
# Alevin --> Droplet technologies, else --> SMART-like technologies
df$technology <- "SMART-like"
df$technology[df$processedWAlevin == TRUE] <- "Droplet-based"

write.csv(df, file = "scAnalysisPipelines.csv", row.names = FALSE)

session <- sessionInfo()
sink(paste0("scAnalysisPipelines_sessionInfo.txt"))
print(session)
sink()