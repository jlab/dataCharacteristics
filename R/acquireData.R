setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(RCurl)
library(XML)

writeURL <- function(contentPath, outFileName) {
  print(contentPath)
  cat(contentPath,
      file = outFileName,
      append = TRUE,
      sep = "\n")
}

writeURLsToFileHelper <- function(projectIds, dataType = c("proteomics", "metabolomics", "transcriptomics", "sc")) {
  for (projectID in projectIds){
    print(projectID)
    if (dataType == "sc"){
      # Filtered because on the EBI website in the download section there are only filtered data available (see e.g. https://www.ebi.ac.uk/gxa/sc/experiments/E-GEOD-99058/downloads)
      cat(paste0("ftp://", projectID, "/", basename(projectID), ".aggregated_filtered_counts.mtx"), 
        file = "sc_unnormalized_urls.txt",
        append = TRUE, 
        sep = "\n")
      cat(paste0("ftp://", projectID, "/", basename(projectID), ".aggregated_filtered_normalised_counts.mtx"), 
          file = "sc_normalized_urls.txt",
          append = TRUE, 
          sep = "\n")
    } else {
      folderpath <- paste0("http://", projectID, "/")
      try({
        folderContent <- getHTMLLinks(folderpath)
      
        # print(folderContent)
        for (content in folderContent){
          # expType <- sub("^[^-]*-([^-]*).*", "\\1", content)
          
          contentPath <- paste0("ftp://", projectID, "/", content)
          
          if (dataType == "proteomics_expressionatlas") {
            if (grepl("-proteinGroups.txt$", content, ignore.case = TRUE))
              writeURL(contentPath, outFileName="proteomics_expressionatlas_urls.txt")
          } else if (dataType == "proteomics_pride") {
            if (grepl("-proteinGroups.txt$", content, ignore.case = TRUE))
              writeURL(contentPath, outFileName="proteomics_pride_urls.txt")
            
          } else if (dataType == "transcriptomics"){
            if (grepl("-raw-counts.tsv$", content, ignore.case = TRUE))
              writeURL(contentPath, outFileName="RNAseq_raw_urls.txt")

            # RNAseq projects with -raw-counts.tsv files 
            # are different from the RNAseq projects that include all of the following files
            if (content == paste0(basename(projectID), "-raw-counts.tsv.undecorated"))
              writeURL(contentPath, outFileName="RNAseq_raw_undecorated_urls.txt")
            
            if (content == paste0(basename(projectID), "-transcripts-raw-counts.tsv.undecorated"))
              writeURL(contentPath, outFileName="RNAseq_transcripts_raw_undecorated_urls.txt")
            
            if (content == paste0(basename(projectID), "-fpkms.tsv"))
              writeURL(contentPath, outFileName="RNAseq_fpkms_urls.txt")
            
            if (content == paste0(basename(projectID), "-tpms.tsv"))
              writeURL(contentPath, outFileName="RNAseq_tpms_urls.txt")
            
            if (content == paste0(basename(projectID), "-transcripts-tpms.tsv"))
              writeURL(contentPath, outFileName="RNAseq_transcripts_tpms_urls.txt")
            
            if (grepl("-normalized-expressions.tsv$", content, ignore.case = TRUE))
              writeURL(contentPath, outFileName="microarray_urls.txt")

          } else if (dataType == "metabolomics"){
            if (grepl("^m_.*\\.tsv$", content, ignore.case = TRUE) & 
                !grepl("NMR", content, ignore.case = TRUE))
              writeURL(contentPath, outFileName="metabolomics_MS_urls.txt")
            
            if (grepl("^m_.*\\.tsv$", content, ignore.case = TRUE) & 
                grepl("NMR", content, ignore.case = TRUE))
              writeURL(contentPath, outFileName="metabolomics_NMR_urls.txt")
          } 
        }
      })
    }
  }
}

getProjectIds <- function(url) {
  projectIds <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
  projectIds <- paste(url, strsplit(projectIds, "\r*\n")[[1]], sep = "")
  projectIds <- projectIds[grepl("\\d$", projectIds)]
  projectIds
}

# !!!Install command line tool ncftp on your system to run this function!!!
writePrideProteinGroupsToFile <- function() {
  urlMain <- 'ftp://ftp.pride.ebi.ac.uk/pride/data/archive/'
  years <- system(paste0('ncftpls -1 ', urlMain), intern=TRUE)
  
  for (year in years){
    print(year)
    url <- paste0(urlMain, year, "/")
    prideFiles <- system(paste0('ncftpls -1R ', url), intern=TRUE)
    #prideFiles <- system(paste0("lftp -c 'open  ", url,"  && find -l && exit'"), intern = TRUE)
    #prideFiles <- strsplit(prideFiles, split= "\\./")
    #prideFiles <- prideFiles[lengths(prideFiles) == 2]
    #prideFiles <-  sapply(prideFiles, "[[", 2)
    prideFiles <- prideFiles[grepl("proteinGroups", prideFiles, ignore.case = TRUE)]
    prideFiles <- sub("^\\./", "", prideFiles)
    # prideFiles <- paste0(url, prideFiles)
    for (prideFile in prideFiles){
      writeURL(paste0(url, prideFile), outFileName="proteomics_pride_urls.txt")
    }
  }
}

writeURLsToFile <- function() {
  
  url <- 'ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments/'
  projectIds <- getProjectIds(url)
  writeURLsToFileHelper(projectIds = gsub("ftp://", "", projectIds), dataType = "sc") 
  
  url <- 'ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/'
  projectIds <- getProjectIds(url)
  writeURLsToFileHelper(projectIds = gsub("ftp://", "", projectIds), dataType = "metabolomics") 
  
  url <- 'ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/'
  projectIds <- getProjectIds(url)
  # projectIds <- sapply(strsplit(projectIds, "/\\s*"), tail, 1)
  projectIds.prot <- projectIds[grepl("-PROT-", projectIds, ignore.case = TRUE)]
  projectIds.geneexpr <- projectIds[!grepl("-PROT-", projectIds)]
  writeURLsToFileHelper(projectIds = gsub("ftp://", "", projectIds.prot), dataType = "proteomics_expressionatlas") 
  writeURLsToFileHelper(projectIds = gsub("ftp://", "", projectIds.geneexpr), dataType = "transcriptomics") 
}

writePrideProteinGroupsToFile()
writeURLsToFile()

session <- sessionInfo()
sink("acquireURLs_sessionInfo.txt")
print(session)
sink()

# Then in shell: wget -i sc_normalized_urls.txt
#system("wget -i sc_normalized_urls.txt")
# or wget -r -nH --cut-dirs=5 -i metabolomics_urls.txt for metabolomics datasets
