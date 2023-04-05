#BiocManager::install("scpdata")
library("scpdata")
eh <- ExperimentHub()
query(eh, "scpdata")
info <- scpdata()
# snapshotDate(): 2022-10-31

for (dataset in row.names(info)){
  print(paste0("Dataset: ", dataset))
  scp <- eh[[dataset]]
  content <- names(scp@assayLinks@listData)
  
  for (el in content[startsWith(content, "proteins")]){
    df <- assay(scp[[el]])
    if (min(dim(df)) > 0)
      write.csv(df, file = paste0("scpdata_prot_", dataset, "_", el, ".csv"))
  }

  #write.csv(assay(scp[["proteins"]]), file = paste0("scpdata_prot_", dataset, ".csv"))
}

session <- sessionInfo()
sink("scpdata_sessionInfo.txt")
print(session)
sink()


# "NA" --> NA
# 0 --> NA

# Definitely not taken log:   
# all EH7295 (williams2020_lfq), EH7085 (schoof2021), 
# EH6011 (liang2020_hela), all EH3909 (zhu2019EL), EH3908 (cong2020AC), 
# EH3905 (zhu2018NC_hela), EH3903 (dou2019_boosting), EH3902 (dou2019_mouse), EH3901 (dou2019_lysates),
# EH7713 (brunner2022) 

# Positive but small
# all EH7296 (williams2020_tmt)

# Negative numbers
# EH3900 (specht2019v3), EH3899 (specht2019v2), all EH7711 (leduc2022), EH7712 (derks2022)

# Specht2019: The [protein] data was again normalized by subtracting the column then the row medians (log2 scale).
# Leduc2022: Single-cell data was normalized relative to the mean within each protein and log2 transformed and imputed values were marked as NA
# derks2022: log2 in code

# williams2020_tmt: For TMT-based quantification, corrected reporter ion intensities were extracted. Reporter ion intensities from single cells were normalized to the reference channel containing 0.2 ng of a peptide mixture from the three AML cell lines at PSM level using MaxQuant (v. 1.6.12.0). (26) To minimize the batch effect from multiple TMT experiments, the relative abundances from 19 TMT plexes were (log 2)-transformed and the data matrices from all of the TMT experiments were combined after sequentially centering the column and row values to their median values. A criterion of >70% valid values and at least two identified peptides per protein were required to generate the “quantifiable” protein list. 
