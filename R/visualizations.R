setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(plyr)
library(dplyr)
library(ggplot2)
library(GGally)
library(ggforce)
library(stringr)
library(tibble)
library(ggbiplot)
theme_set(theme_bw())

library(plotly)
library(htmlwidgets)
################################################################################
################################################################################
# FUNCTIONS

translateTermsInColumns <- function(oldAndID.df, translation.df, 
                                    oldCol = "old", idCol = "ID", 
                                    multipleHandling = "Multiple") {
  replaced <- strsplit(oldAndID.df[,oldCol], ";")
  for (row in seq(nrow(translation.df))) {
    replaced <- rapply(replaced, function(x){
      gsub(paste0("^", translation.df[row, ]$old, "$"),
           translation.df[row, ]$new, x)
    }, how = "list")
  }
  
  replaced2 <- rapply(replaced, function(x){ifelse(length(unique(x)) > 1,
                                                   multipleHandling, unique(x))
  }, how = "list")
  
  replaced2.df <- data.frame(ID = oldAndID.df[,idCol], 
                             old = oldAndID.df[,oldCol], 
                             new = unlist(replaced2))
  replaced2.df
}

logTransform <- function(df, variable, logBase = c("log2", "log1p")){
  df[[paste0(logBase, "(", variable, ")")]] <- get(logBase)(df[[variable]])
  df[[variable]] <- NULL
  df
}

plotBoxplots <- function(data2.long, fileNameAddition = "", height = 12, 
                         width = 18, varsToRank = c(), order = NULL) {
  if (length(varsToRank) > 0) {
    data2.long <- data2.long %>% na.omit() %>% 
      dplyr::group_by(variable) %>%
      dplyr::mutate(value = ifelse(variable %in% varsToRank, 
                                   rank(value), value)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(variable = ifelse(
        variable %in% varsToRank, paste0("rank(", variable, ")"), 
        paste0(variable)))
    
    for (vars in varsToRank) {
      order[order == vars] <- paste0("rank(", vars, ")")
    }
    
    data2.long <- dplyr::arrange(
      dplyr::mutate(data2.long, variable = factor(variable, levels = order)), 
      variable)
  }
  # For each data characteristic: Median of the median of each data type
  medianValues <- data2.long %>%
    group_by(`Data type`, variable) %>%  
    summarise(medianValue = median(value, na.rm = T)) %>%
    group_by(variable) %>%  
    summarise(medianValue = median(medianValue, na.rm = T))
  write.csv(medianValues, paste0("medianValues", fileNameAddition,".csv"), 
            row.names = FALSE)
  
  ggplot.charact <- ggplot(
    data2.long, aes(forcats::fct_rev(`Data type`), value)) +
    # geom_boxplot(aes(fill = `Data type`), alpha=0.5, outlier.size=0.5) +
    geom_violin(aes(fill = `Data type`), alpha = 0.5, scale = "width") +
    geom_boxplot(aes(fill = `Data type`), width = 0.5, alpha = 0.25, 
                 outlier.size = 0.5) +
    coord_flip() +
    xlab("") +
    ylab("") +
    geom_hline(aes(yintercept = medianValue), medianValues, colour = 'red') +
    facet_wrap( ~ variable, scales = "free_x", ncol = 6, 
                strip.position = "bottom", 
                labeller = ggplot2::label_wrap_gen(width = 30)) +
    ggplot2::theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.spacing.y = unit(0.5, "lines"),
          legend.title = element_blank(), 
      axis.text.y = element_text(hjust = 0, face = "bold"), 
          # legend.position="bottom",
          legend.position = "none",
          legend.justification = "left", legend.direction = "horizontal",
          strip.text = element_text(face = "bold",
                                    vjust = 1#, size = 6
                                    ),
          strip.placement = "outside", # Place facet labels outside x axis labels.
          strip.background = element_blank(),  # Make facet label background white.
          axis.title = element_blank()) +     
    guides(fill = ggplot2::guide_legend(reverse = TRUE))
  ggsave(file = paste0("boxplots", fileNameAddition, ".pdf"), 
         ggplot.charact, height = height, width = width)
}

colorForCorValues <- function(data, mapping, method="spearman", 
                              use="pairwise", ...){
  
  # grab data
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # calculate correlation
  corr <- cor(x, y, method = method, use = use)
  
  # calculate colour based on correlation value
  # Here I have set a correlation of minus one to blue, 
  # zero to white, and one to red 
  # Change this to suit: possibly extend to add as an argument of `my_fn`
  colFn <- colorRampPalette(c("blue", "white", "red"), interpolate = 'spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length = 100))]
  
  ggally_cor(data = data, mapping = mapping, 
             method = method,
             colour = "black", ...) + 
    theme_void() +
    theme(panel.background = element_rect(fill = fill))
}

plotPairsPlotForTypes <- function(df, groupColName = "Data type", 
                                  colsForCorr = c(), corrMethod = "spearman", 
                                  smooth = FALSE, width = 19, height = 19, 
                                  addStr = "") {
  for (group in unlist(unique(df[, groupColName]))) {
    print(group)
    df.group <- df[df[, groupColName] == group,]
    df.group <- df.group[, colsForCorr]
    colnames(df.group) <- stringr::str_wrap(colnames(df.group), width = 15)
    pdf(paste0("ggpairs_", gsub(" ", "_", group), addStr, ".pdf"), 
        width = width, height = height)
    if (smooth) {
      ggpairsObj <- ggpairs(
        df.group, 
        lower = list(continuous = wrap("smooth", alpha = 0.1, size = 0.1)),
        upper = list(continuous = wrap(colorForCorValues, 
                                       method = corrMethod, size = 2.5),
                     continuous = "smooth")
      )
    } else {
      ggpairsObj <- ggpairs(
        df.group, 
        # lower = list(continuous = wrap("smooth", alpha = 0.1, size = 0.1)),
        lower = list(continuous = wrap("points", alpha = 0.1, size = 0.1)),
        upper = list(continuous = wrap(colorForCorValues, 
                                       method = corrMethod, size = 2.5),
                     continuous = "smooth")
      )
    }
    
    print(ggpairsObj + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
            ggtitle(group) +
            theme(strip.placement = "outside", text = element_text(size = 6))
    )
    
    dev.off()
  }
}

plotPCABiplot <- function(df, groups= c(), alpha = 0.5, 
                          pcaMethod = "nipals",
                          coordRatio = 0.5, 
                          facetZoom = TRUE, 
                          xlimLower = NA, xlimUpper = NA,
                          ylimLower = NA, ylimUpper = NA,
                          PCchoices = 1:2,
                          ellipse = TRUE,
                          addStr = "") {
  # See https://stackoverflow.com/a/49788251
  #   # devtools::install_github("vqv/ggbiplot")
  # library(ggbiplot)
  
  iris_dummy <- df
  #swap out your NAs with a dummy number so prcomp will run
  iris_dummy[is.na(iris_dummy)] <- 7777
  pca.obj <- prcomp(iris_dummy, center = TRUE, scale.=TRUE)
  
  # scale: One of "UV" (unit variance a=a/\sigma_{a}), 
  # "vector" (vector normalisation b=b/|b|), 
  # "pareto" (sqrt UV) or "none" 
  # to indicate which scaling should be used to scale the matrix with 
  # aa variables and b samples. 
  pca.obj2 <- pcaMethods::pca(df, method = pcaMethod, nPcs = 4, center = TRUE
                              , scale = "uv"
  )
  
  pca.obj$x <- pca.obj2@scores 
  pca.obj$rotation <- pca.obj2@loadings 
  pca.obj$sdev <- pca.obj2@sDev
  pca.obj$center <- pca.obj2@center
  pca.obj$scale <- pca.obj2@scale
  
  PCA.data <- data.frame(pca.obj[["x"]])
  write.csv(PCA.data, file = paste0("biplot_",  pcaMethod, "_", addStr,".csv"))
  
  P2 <- ggbiplot::ggbiplot(pca.obj,
                           choices = PCchoices,
                           obs.scale = 1, 
                           var.scale = 1,
                           ellipse = ellipse,
                           circle = FALSE,
                           varname.size = 2.5,
                           var.axes = T,
                           groups = groups, 
                           alpha = 0)  +
    scale_color_discrete(name = '') +  
    coord_fixed(ratio = coordRatio) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(alpha = 1, size = 3)))
  
  if (all(PCchoices == 1:2)) P2 <- P2 + 
    ggplot2::scale_y_continuous(limits = c(-5, 10.3))
  
  # For nipals:
  # row.names(PCA.data)[which(PCA.data$PC2<(-5) | PCA.data$PC2>10.3)] 
  # # [1] "MTBLS2801_m_MTBLS2801_LC-MS_negative_reverse-phase_metabolite_profiling_v2_maf.tsv"
  # # [2] "MTBLS2801_m_MTBLS2801_LC-MS_positive_reverse-phase_metabolite_profiling_v2_maf.tsv"
  # # [3] "MTBLS4186_m_MTBLS4186_LC-MS_positive_reverse-phase_metabolite_profiling_v2_maf.tsv"
  
  # For svd:
  # row.names(PCA.data)[which(PCA.data$PC2<(-5) | PCA.data$PC2>10.3)] 
  # [1] "MTBLS1183_m_clinic_metabolite_profiling_mass_spectrometry-1_v2_maf.tsv"            
  # [2] "MTBLS2801_m_MTBLS2801_LC-MS_negative_reverse-phase_metabolite_profiling_v2_maf.tsv"
  # [3] "MTBLS2801_m_MTBLS2801_LC-MS_positive_reverse-phase_metabolite_profiling_v2_maf.tsv"
  # [4] "MTBLS4186_m_MTBLS4186_LC-MS_negative_reverse-phase_metabolite_profiling_v2_maf.tsv"
  # [5] "MTBLS4186_m_MTBLS4186_LC-MS_positive_reverse-phase_metabolite_profiling_v2_maf.tsv"
  # [6] "E-MTAB-7869.aggregated_filtered_counts.mtx"   
  
  if (facetZoom) {
    P2 <- P2 + ggforce::facet_zoom(xlim = c(xlimLower, xlimUpper), 
                                   ylim = c(ylimLower, ylimUpper))
  } else {
    if (!all(is.na(c(xlimLower, xlimUpper))))
      P2 <- P2 + xlim(c(xlimLower, xlimUpper))
    
    if (!all(is.na(c(ylimLower, ylimUpper))))
      P2 <- P2 + ylim(c(ylimLower, ylimUpper))
  }
  
  
  P2 <- P2 + theme(legend.direction = 'horizontal', 
                    legend.position = 'bottom')
  
  P2$layers <- c(geom_point(aes(colour = groups), cex = 1, alpha = alpha), 
                 P2$layers)
  
  P2
}

plotPCABiplots <- function(df, groupColName = "", addStr = "", 
                           pcaMethod = "nipals") {
  
  pdf(file = paste0("biplot_", pcaMethod, "_", addStr, ".pdf"),
      width = 11, height = 9)
  print(plotPCABiplot(df = df %>% dplyr::select(-!!groupColName), 
                      groups = df[[groupColName]],
                      alpha = 0.3,
                      pcaMethod = pcaMethod,
                      coordRatio = 1/2,
                      facetZoom = FALSE,
                      addStr = addStr))
  dev.off()
  
  pca <- pcaMethods::pca(df %>% dplyr::select(-!!groupColName), 
                         method = pcaMethod, nPcs = 4, center = TRUE
                         , scale = "uv")
  write.csv(pca@loadings, paste0("loadings_",  pcaMethod, "_", addStr,".csv"))
  dat <- merge(pcaMethods::scores(pca), df, by = 0)
  
  pdf(paste0("ggpairs_", pcaMethod, "_", addStr, ".pdf"), 
      width = 12, height = 10)
  print(GGally::ggpairs(dat, 
                        columns = 2:5, ggplot2::aes(colour = get(groupColName)),
                        lower = list(
                          continuous = wrap("smooth", alpha = 0.3, size = 1),
                          combo = wrap("dot_no_facet", alpha = 0.4)),
                        upper = list(continuous = wrap("cor", 
                                                     method = "spearman", 
                                                     size = 3)),
                        mapping = aes(color = get(groupColName),
                                    fill = get(groupColName), 
                                    alpha = 0.5)) +
          ggtitle(groupColName) +
          theme_bw())
  dev.off()
  
  fig <- plotly::plot_ly(dat, x = ~PC1, y = ~PC2, z = ~PC3, 
                 color = ~as.factor(dat[[groupColName]]), 
                 type = "scatter3d", mode = "markers",
                 #colors = c('#636EFA','#EF553B') , 
                 marker = list(size = 2)
                 # , alpha = 0.75
  ) #%>%
  # add_markers(size = 5, marker=list(sizeref=8, sizemode="area"))
  fig <- fig %>%
    plotly::layout(
      title = "3D PCA",
      scene = list(bgcolor = "#e5ecf6"
      )
    )
  
  htmlwidgets::saveWidget(fig, 
                          paste0("plotly_", pcaMethod, "_", addStr,".html"), 
                          selfcontained = F, libdir = "lib")
}

plotUMAPplots <- function(df, groupColName = "", addStr = "", alpha = 0.3) {
  set.seed(142)
  umap_fit <- df %>% dplyr::mutate(ID = row_number())  %>%
    dplyr::select(-!!groupColName) %>% 
    dplyr::select_if(~ !any(is.na(.))) %>%
    remove_rownames() %>% column_to_rownames("ID") %>%
    scale() %>%
    umap::umap(n_components = 3)
  
  groupVec <- df[[groupColName]]
  umap_df <- umap_fit$layout %>%
    as.data.frame() %>%
    dplyr::rename(UMAP1 = "V1",
                  UMAP2 = "V2",
                  UMAP3 = "V3") %>%
    dplyr::mutate(!!groupColName := !!groupVec)
  row.names(umap_df) <- row.names(df)
  
  write.csv(umap_df, file = paste0("UMAP_", addStr,".csv"))
  
  pdf(file = paste0("umap_Dim1vsDim2_", addStr, ".pdf"), width = 11, height = 9)
  print(
    ggplot2::ggplot(umap_df, ggplot2::aes(
      x = UMAP1, y = UMAP2, color = get(groupColName))) +
      ggplot2::geom_point(size = 0.8, alpha = alpha)  + 
      ggplot2::theme(legend.direction = 'horizontal', 
                     legend.position = 'bottom') +
      ggplot2::scale_color_discrete(name = '') +
      ggplot2::guides(colour = ggplot2::guide_legend(
        override.aes = list(alpha = 1, size = 3)))
  )
  dev.off()
  
  pdf(file = paste0("umap_Dim2vsDim3_", addStr, ".pdf"), width = 11, height = 9)
  print(
    ggplot2::ggplot(umap_df, ggplot2::aes(x = UMAP2, y = UMAP3, 
                                          color = get(groupColName))) +
      ggplot2::geom_point(size = 0.8, alpha = alpha)  + 
      ggplot2::theme(legend.direction = 'horizontal', 
                     legend.position = 'bottom') +
      ggplot2::scale_color_discrete(name = '') +
      ggplot2::guides(colour = ggplot2::guide_legend(
        override.aes = list(alpha = 1, size = 3)))
  )
  dev.off()
}

################################################################################
################################################################################

dataset <- ldply(list.files("20230614_results", pattern = ".csv", 
                            full.names = TRUE), read.csv, header = TRUE)
dataset <- dataset[dataset$nSamples > 0,]
dataset <- dataset[dataset$variance > 0,]
dataset$prctPC1.wNAs.log2 <- 100 * dataset$prctPC1.wNAs.log2
dataset$prctPC2.wNAs.log2 <- 100 * dataset$prctPC2.wNAs.log2

names(dataset) <- gsub(x = names(dataset), pattern = "\\.log2|\\.wNAs", 
                       replacement = "")  
scPipeline.df <- read.csv("scAnalysisPipelines.csv")

dataset <- dataset %>% 
  dplyr::mutate(
    dataType = case_when(
      grepl("\\^iBAQ", datasetID) ~ paste0(dataType, "_iBAQ"),
      grepl("\\^LFQ", datasetID) ~ paste0(dataType, "_LFQ"),
      grepl("\\^Intensity", datasetID) ~ paste0(dataType, "_Intensity"),
      datasetID == "MTBLS242_m_mtbls242_v2_maf.tsv" ~ "metabolomics_NMR", 
      dataType == "RNAseq_raw_undecorated" ~ "RNAseq_raw",
      TRUE ~ dataType),
    dataTypeSubgroups = case_when(
      # For LC-MS: Add negative and positive?
      grepl("GC", datasetID, ignore.case = TRUE) & 
        dataType == "metabolomics_MS" ~ paste0(dataType, "_GC"),
      grepl("LC|HILIC|RP", datasetID, ignore.case = TRUE) & 
        dataType == "metabolomics_MS" ~ paste0(dataType, "_LC"),
      grepl("FIA", datasetID, ignore.case = TRUE) & 
        dataType == "metabolomics_MS" ~ paste0(dataType, "_FIA"),
      grepl("MALDI|MSImaging|DESI", datasetID, ignore.case = TRUE) & 
        dataType == "metabolomics_MS" ~ paste0(dataType, "_otherIonization"),
      grepl("_DI-", datasetID, ignore.case = TRUE) & 
        dataType == "metabolomics_MS" ~ paste0(dataType, "_DI"),
      dataType == "metabolomics_MS" ~ paste0(dataType, "_Undefined"),
      dataType %in% c(
        "sc_normalized", "sc_unnormalized") ~ paste0(
          dataType, "_", 
          scPipeline.df$technology[match(gsub("\\..*", "", datasetID), 
                                         scPipeline.df$projectId)]),
      dataType == "microarray" ~ paste0(
        dataType, "_", sub(".*-([A-Z]+-\\d+).*", "\\1", datasetID)),
      dataType == "microbiome" ~ paste0(
        dataType, "_", ifelse(grepl("16s", datasetID, 
                                    ignore.case = TRUE), "16S", "WGS")),
      TRUE ~ dataType)
  ) 

json.df <- read.csv("microarray_metadata.csv")

manufacturer.df <- json.df %>% 
  dplyr::select("chipID", "Manufacturer")  %>%  
  dplyr::mutate(dataTypeSubgroups = paste0(
    "microarray_", sub("A-", "", chipID))) %>% 
  dplyr::select(-chipID)

dataset <- dplyr::left_join(dataset, manufacturer.df, 
                            by = "dataTypeSubgroups") %>%
  dplyr::mutate(dataTypeSubgroups = case_when(
    grepl("microarray_", dataTypeSubgroups) ~ paste0(
      "microarray_", Manufacturer), TRUE ~ dataTypeSubgroups)
  ) %>% dplyr::select(-Manufacturer)

prideMeta.df <- read.csv("proteomicsPride_metadata.csv")

includedPrideIDs <- sort(
  unique(sub("^(.*?)_.*$", "\\1", 
             dataset[grepl("_pride_", dataset$dataType),]$datasetID)))
prideMeta.df.forIncludedPrideIDs <- prideMeta.df[
  prideMeta.df$prideID %in% includedPrideIDs, ]
setdiff(includedPrideIDs, prideMeta.df.forIncludedPrideIDs$prideID)
# [1] "PXD021464" "PXD022721" "PXD024748" "PXD024937"
# Add information manually for PXDs that could  not be found automatically
prideMeta.df.forIncludedPrideIDs <- rbind(
  prideMeta.df.forIncludedPrideIDs, 
  c(prideID = "PXD021464", 
    tax = "Mus musculus (mouse)", 
    organismParts = "Hippocampal cell line",
    experimentTypes = "Shotgun proteomics",
    instruments = "LTQ Orbitrap Velos"),
  c(prideID = "PXD022721", 
    tax = "Homo sapiens (human)", 
    organismParts = "Cell culture",
    experimentTypes = "Shotgun proteomics",
    instruments = "LTQ Orbitrap Elite;Orbitrap Fusion;Q Exactive HF"),
  c(prideID = "PXD024748", 
    tax = "Homo sapiens (human)", 
    organismParts = "Substantia nigra",
    experimentTypes = "Shotgun proteomics",
    instruments = "LTQ Orbitrap Velos"),
  c(prideID = "PXD024937", 
    tax = "Homo sapiens (human)", 
    organismParts = "Brain",
    experimentTypes = "Shotgun proteomics",
    instruments = "Q Exactive")
)
prideMeta.df.forIncludedPrideIDs$instruments <- gsub(
  ";MaxQuant", "", prideMeta.df.forIncludedPrideIDs$instruments)
write.csv(data.frame(instrument = sort(
  unique(unlist(strsplit(prideMeta.df.forIncludedPrideIDs$instruments, ";"))))), 
          "proteomicsPride_instruments.csv",
          row.names = FALSE)

prideMetaInstrumentsManufacturerAdded.df <- read.csv(
  "proteomicsPride_instruments_manufacturers.csv")

prideTranslation.df <- translateTermsInColumns(
  oldAndID.df = prideMeta.df.forIncludedPrideIDs, 
  translation.df = prideMetaInstrumentsManufacturerAdded.df %>%
    dplyr::rename("old" = "instrument", "new" = "manufacturer"), 
  oldCol = "instruments", idCol = "prideID", 
  multipleHandling = "Multiple")

colnames(prideTranslation.df) <- c("prideID", "instrument", "manufacturer")

prideMeta.df.forIncludedPrideIDs <- dplyr::left_join(
  prideMeta.df.forIncludedPrideIDs, prideTranslation.df, by = join_by(prideID))

dataset <- dataset %>% 
  dplyr::mutate(dataTypeSubgroups = case_when(
    grepl("_pride_", dataType) ~ 
      paste0(dataTypeSubgroups, "_", 
             prideMeta.df.forIncludedPrideIDs$manufacturer[
               match(sub("^(.*?)_.*$", "\\1", datasetID), 
                     prideMeta.df.forIncludedPrideIDs$prideID)]),
    TRUE ~ dataTypeSubgroups)
  )

metabolomicsMeta.df <- read.csv("metabolomics_metadata.csv")
metabolomicsDatasets <- dataset %>% 
  filter(dataType %in% c("metabolomics_MS", "metabolomics_NMR")) %>%
  dplyr::select(datasetID, dataType, dataTypeSubgroups) %>%
  dplyr::mutate(accession = sub("^(.*?)_.*$", "\\1", datasetID)) %>%
  dplyr::left_join(metabolomicsMeta.df, by = join_by(accession == accession))

replaced <- strsplit(metabolomicsDatasets$technology, ";")

replaced2 <- rapply(replaced, function(x){
  x[grepl("LC|UPLD", x)] <- "LC"
  x[grepl("GC", x)] <- "GC"
  x[grepl("FIA", x)] <- "FIA"
  x[grepl("DI-", x)] <- "DI"
  x[grepl("CE-", x)] <- "CE"
  x[grepl("NMR|MRS", x)] <- "NMR"
  x[grepl("MALDI|DART", x)] <- "otherIonization"
  x[grepl("Insufficient data supplied", x)] <- NA
  x <- sort(unique(x))
  x
}, how = "list") 

replaced2 <- lapply(replaced2, function(x) paste0(x, collapse = ";"))

metabolomicsTranslation.df <- cbind(
  metabolomicsDatasets %>% dplyr::select(- dataType, accession), 
  technology2 = unlist(replaced2))

metabolomicsTranslation.df[grepl(
  ";", metabolomicsTranslation.df$technology2) | 
    metabolomicsTranslation.df$technology2 == "", ]$technology2 <- NA
metabolomicsTranslation.df <- metabolomicsTranslation.df %>% 
  dplyr::mutate(
    dataTypeSubgroupsNew = case_when(
      dataTypeSubgroups == "metabolomics_MS_Undefined" & 
        !is.na(technology2) ~ paste0("metabolomics_MS_", technology2), 
      TRUE ~ dataTypeSubgroups)
)

# MTBLS440 --> LC
# MTBLS728 --> LC
metabolomicsTranslation.df[
  metabolomicsTranslation.df$accession %in% 
    c("MTBLS440", "MTBLS728"),]$dataTypeSubgroupsNew <- "metabolomics_MS_LC"
dataset <- dataset %>%
  dplyr::mutate(
    dataTypeSubgroups = case_when(
      grepl("metabolomics_", dataType) ~ 
        metabolomicsTranslation.df$dataTypeSubgroupsNew[match(
          datasetID, metabolomicsTranslation.df$datasetID)],
      TRUE ~ dataTypeSubgroups)
  )

# Remove outliers: 
# E-GEOD-152766.aggregated_filtered_counts.mtx, (also E-GEOD-152766.aggregated_filtered_normalised_counts.mtx to make
# numbers of normalized and unnormalized sc datasets even.)
# MTBLS407_m_MTBLS407_21022_metabolite_profiling_mass_spectrometry_v2_maf.tsv,
# MTBLS407_m_MTBLS407_21548_metabolite_profiling_mass_spectrometry_v2_maf.tsv
data <- dataset %>% 
  dplyr::filter(!(
    datasetID %in% 
      c("E-GEOD-152766.aggregated_filtered_counts.mtx", 
        "E-GEOD-152766.aggregated_filtered_normalised_counts.mtx",
        "MTBLS407_m_MTBLS407_21022_metabolite_profiling_mass_spectrometry_v2_maf.tsv",
        "MTBLS407_m_MTBLS407_21548_metabolite_profiling_mass_spectrometry_v2_maf.tsv")))


# # Check for duplicates --> e.g. because a dataset was re-uploaded under a different 
# # file name and the old and new file were both present on the FTP server, 
# # or if dataset was reuploaded under a different project ID
data.duplicateRowsOnly <- data[which(
  duplicated(data %>% dplyr::select(-datasetID)) | 
    duplicated(data %>% dplyr::select(-datasetID), 
               fromLast = TRUE)), ] %>% 
  arrange(nAnalytes, minRowNonNaNumber, mean)
data.duplicateRowsOnly2 <- data[which(
  duplicated(data %>% dplyr::select(-c(datasetID, dataType, dataTypeSubgroups))) | 
    duplicated(data %>% dplyr::select(-c(datasetID, dataType, dataTypeSubgroups)), 
               fromLast = TRUE)), ] %>% 
  arrange(nAnalytes, minRowNonNaNumber, mean)

dataDiff <- data[data$datasetID %in% setdiff(
  data.duplicateRowsOnly2$datasetID, data.duplicateRowsOnly$datasetID), ]

# Remove as they are datasets present also in scProteomics category:
# "PXD006847_CulturedCells_proteinGroups.txt_^Intensity_",
# "PXD006847_CulturedCells_proteinGroups.txt_^iBAQ_",
# "PXD006847_CulturedCells_proteinGroups.txt_^Intensity_"
# "PXD021882_proteinGroups.txt_^LFQ_"
# "PXD021882_proteinGroups.txt_^Intensity_"

# Remove as it has duplicate in different category
# "MTBLS103_m_ibanez_02_metabolite_profiling_mass_spectrometry_v2_maf.tsv"

data <- data[!(
  data$datasetID %in% 
    c("PXD006847_CulturedCells_proteinGroups.txt_^LFQ_",
      "PXD006847_CulturedCells_proteinGroups.txt_^iBAQ_",
      "PXD006847_CulturedCells_proteinGroups.txt_^Intensity_",
      "PXD021882_proteinGroups.txt_^LFQ_",
      "PXD021882_proteinGroups.txt_^Intensity_",
      "MTBLS103_m_ibanez_02_metabolite_profiling_mass_spectrometry_v2_maf.tsv")),]

##  Duplicate rows removed
data <- data[-which(duplicated(data %>% dplyr::select(-datasetID))), ]
write.csv(data, "datasets_results.csv", row.names = FALSE)

################################################################################
################################################################################

# Remove datasets with negative values --> potentially due to log already taken
data <- data[data$nNegativeNumbers == 0,]

# Remove datasets with medianAnalyteVariance or medianSampleVariance of zero or NA
data <- data[!is.na(data$medianAnalyteVariance) & 
               (data$medianAnalyteVariance > 0),]
data <- data[!is.na(data$medianSampleVariance) & 
               (data$medianSampleVariance > 0),]

# Remove "bimodalityRowCorr" because of too many NAs
# data <- dataset %>% dplyr::select(-c(nNegativeNumbers, bimodalityRowCorr, bimodalityColCorr))
data$bimodalityRowCorr <- NULL

data <- data %>%
  dplyr::mutate_if(is.integer, as.numeric) %>% 
  dplyr::mutate("prctnDistinctValues" = nDistinctValues/(
    (nSamples * nAnalytes) * ((100 - percNATotal)/100)) * 100 ) %>%
  dplyr::select(-c(nDistinctValues, nNegativeNumbers))

seedCols <- c("bimodalityRowCorrSeed", "bimodalityColCorrSeed", 
              "coefHclustRowsSeed", "intensityNAProbSeed"
)

data <- data[,!(colnames(data) %in% seedCols)]

# data <- data %>%
#   mutate_if(is.integer, as.numeric) %>% 
#   mutate("prctminRowNonNaNumber" = (nSamples - minRowNonNaNumber)/nSamples * 100 ) %>%
#   mutate("prctmaxRowNonNaNumber" = (nSamples - maxRowNonNaNumber)/nSamples * 100 ) 
# # prctminRowNonNaNumber is the same as maxRowNaPercentage, prctmaxRowNonNaNumber is the same as minRowNaPercentage
data <- data[,!(colnames(data) %in% 
                  c("minRowNonNaNumber", "maxRowNonNaNumber"))] 

write.csv(data, "datasets_results_clean.csv", row.names = FALSE)
################################################################################

# Rename Data types
for (dataTypeLevel in c("dataType", "dataTypeSubgroups")) {
  if (dataTypeLevel == "dataType") {
    OldDataTypeNames <- c("metabolomics_MS", "metabolomics_NMR", 
                          "microarray", "microbiome", 
                          "proteomics_expressionatlas_iBAQ", 
                          "proteomics_expressionatlas_Intensity", 
                          "proteomics_pride_LFQ", "proteomics_pride_Intensity",
                          "proteomics_pride_iBAQ", 
                          "RNAseq_fpkms_median", "RNAseq_raw", 
                          "RNAseq_tpms_median", "sc_normalized", 
                          "sc_unnormalized", "scProteomics")
    NewDataTypeNames <- c("Metabolomics (MS)", "Metabolomics (NMR)", 
                          "Microarray", "Microbiome", 
                          "Proteomics (iBAQ, Expression Atlas)", 
                          "Proteomics (Intensity, Expression Atlas)", 
                          "Proteomics (LFQ, PRIDE)", 
                          "Proteomics (Intensity, PRIDE)", 
                          "Proteomics (iBAQ, PRIDE)", 
                          "RNA-seq (FPKM)", "RNA-seq (raw)", "RNA-seq (TPM)", 
                          "scRNA-seq (normalized)", 
                          "scRNA-seq (unnormalized)", "scProteomics")
    levels <- c( "Metabolomics (NMR)", "Metabolomics (MS)", 
                 "Proteomics (iBAQ, Expression Atlas)", 
                 "Proteomics (Intensity, Expression Atlas)",
                 "Proteomics (iBAQ, PRIDE)", 
                 "Proteomics (Intensity, PRIDE)", "Proteomics (LFQ, PRIDE)", 
                 "scProteomics",
                 "Microarray", 
                 "RNA-seq (raw)", "RNA-seq (FPKM)", "RNA-seq (TPM)", 
                 "scRNA-seq (unnormalized)",
                 "scRNA-seq (normalized)", 
                 "Microbiome")
  } else if (dataTypeLevel == "dataTypeSubgroups") {
    OldDataTypeNames <- c(
      "metabolomics_MS_LC", "metabolomics_MS_GC", "metabolomics_MS_Undefined", 
      "metabolomics_MS_DI", "metabolomics_MS_FIA", "metabolomics_NMR", 
      "metabolomics_MS_CE", "metabolomics_MS_otherIonization", 
      "microarray_Affymetrix", 
      "microarray_Illumina", "microarray_Agilent", "microbiome_16S", 
      "microbiome_WGS", "proteomics_expressionatlas_iBAQ", 
      "proteomics_expressionatlas_Intensity", 
      "proteomics_expressionatlas_LFQ", "proteomics_pride_LFQ_Thermo", 
      "proteomics_pride_Intensity_Thermo", "proteomics_pride_iBAQ_Thermo", 
      "proteomics_pride_LFQ_Multiple", "proteomics_pride_iBAQ_Multiple", 
      "proteomics_pride_Intensity_Multiple", 
      "proteomics_pride_Intensity_Bruker", 
      "proteomics_pride_LFQ_Agilent", "proteomics_pride_Intensity_Agilent", 
      "proteomics_pride_LFQ_Bruker", "proteomics_pride_LFQ_SCIEX", 
      "proteomics_pride_Intensity_SCIEX", "proteomics_pride_iBAQ_SCIEX", 
      "proteomics_pride_iBAQ_Bruker", "proteomics_pride_iBAQ_Agilent", 
      "RNAseq_fpkms_median", 
      "RNAseq_raw", 
      "RNAseq_tpms_median", 
      "sc_normalized_SMART-like",
      "sc_normalized_Droplet-based", "sc_unnormalized_SMART-like", 
      "sc_unnormalized_Droplet-based", "scProteomics")
    NewDataTypeNames <- c(
      "Metabolomics (LC-MS)", "Metabolomics (GC-MS)", 
      "Metabolomics (Undefined-MS)", 
      "Metabolomics (DI-MS)", "Metabolomics (FIA-MS)", "Metabolomics (NMR)", 
      "Metabolomics (CE-MS)", "Metabolomics (Other ionization-MS)", 
      "Microarray (Affymetrix)", 
      "Microarray (Illumina)", "Microarray (Agilent)", "Microbiome (16S)", 
      "Microbiome (WGS)", "Proteomics (iBAQ, Expression Atlas)", 
      "Proteomics (Intensity, Expression Atlas)", 
      "Proteomics (LFQ, Expression Atlas)", "Proteomics (LFQ, PRIDE, Thermo)", 
      "Proteomics (Intensity, PRIDE, Thermo)", 
      "Proteomics (iBAQ, PRIDE, Thermo)", 
      "Proteomics (LFQ, PRIDE, Undefined)", 
      "Proteomics (iBAQ, PRIDE, Undefined)", 
      "Proteomics (Intensity, PRIDE, Undefined)", 
      "Proteomics (Intensity, PRIDE, Bruker)", 
      "Proteomics (LFQ, PRIDE, Agilent)", 
      "Proteomics (Intensity, PRIDE, Agilent)", 
      "Proteomics (LFQ, PRIDE, Bruker)", "Proteomics (LFQ, PRIDE, SCIEX)", 
      "Proteomics (Intensity, PRIDE, SCIEX)", "Proteomics (iBAQ, PRIDE, SCIEX)", 
      "Proteomics (iBAQ, PRIDE, Bruker)", "Proteomics (iBAQ, PRIDE, Agilent)", 
      "RNA-seq (FPKM)", 
      "RNA-seq (raw)", 
      "RNA-seq (TPM)", 
      "scRNA-seq (SMART-like, normalized)",
      "scRNA-seq (Droplet-based, normalized)", 
      "scRNA-seq (SMART-like, unnormalized)",
      "scRNA-seq (Droplet-based, unnormalized)", "scProteomics")
    
    levels <- c(
      "Metabolomics (NMR)", 
      "Metabolomics (GC-MS)", "Metabolomics (LC-MS)", "Metabolomics (DI-MS)", 
      "Metabolomics (FIA-MS)", 
      "Metabolomics (CE-MS)", "Metabolomics (Other ionization-MS)", 
      "Metabolomics (Undefined-MS)", 
      
      "Proteomics (iBAQ, Expression Atlas)", 
      "Proteomics (Intensity, Expression Atlas)", 
      "Proteomics (LFQ, Expression Atlas)", 
      
      "Proteomics (iBAQ, PRIDE, Agilent)", "Proteomics (iBAQ, PRIDE, Thermo)", 
      "Proteomics (iBAQ, PRIDE, Bruker)", 
      "Proteomics (iBAQ, PRIDE, SCIEX)", "Proteomics (iBAQ, PRIDE, Undefined)", 
      
      "Proteomics (Intensity, PRIDE, Agilent)", 
      "Proteomics (Intensity, PRIDE, Thermo)", 
      "Proteomics (Intensity, PRIDE, Bruker)", 
      "Proteomics (Intensity, PRIDE, SCIEX)", 
      "Proteomics (Intensity, PRIDE, Undefined)", 
      
      "Proteomics (LFQ, PRIDE, Agilent)", "Proteomics (LFQ, PRIDE, Thermo)", 
      "Proteomics (LFQ, PRIDE, Bruker)", 
      "Proteomics (LFQ, PRIDE, SCIEX)", "Proteomics (LFQ, PRIDE, Undefined)",
      
      "scProteomics",
      
      "Microarray (Affymetrix)", "Microarray (Illumina)", 
      "Microarray (Agilent)", 
      
      "RNA-seq (raw)", 
      "RNA-seq (FPKM)", 
      "RNA-seq (TPM)", 
      
      "scRNA-seq (SMART-like, unnormalized)",
      "scRNA-seq (Droplet-based, unnormalized)", 
      "scRNA-seq (SMART-like, normalized)", 
      "scRNA-seq (Droplet-based, normalized)",
      
      "Microbiome (16S)", "Microbiome (WGS)"
    )
  }
  
  renameDataTypeTable <- data.frame(OldDataTypeNames = OldDataTypeNames,
                                    NewDataTypeNames = NewDataTypeNames)
  
  data[, dataTypeLevel] <- renameDataTypeTable$NewDataTypeNames[
    match(data[, dataTypeLevel], renameDataTypeTable$OldDataTypeNames)]
  
  data[, dataTypeLevel] <- factor(data[, dataTypeLevel], levels = levels)
}  

# Rename Variables
Oldnames <- c("datasetID", "dataType", "dataTypeSubgroups", "nSamples", 
              "nAnalytes", "minRowNaPercentage", 
              "maxRowNaPercentage", "minColNaPercentage", "maxColNaPercentage", 
              "percNATotal", "percOfRowsWithNAs", "percOfColsWithNAs", 
              "corSampleMeanNA", 
              "corSampleMeanNAPval", "corAnalyteMeanNA", "corAnalyteMeanNAPval", 
              "mean", "median", "min", "max", "medianSampleVariance", 
              "medianAnalyteVariance", 
              "variance", "kurtosis", "skewness", "prctPC1", "prctPC2", 
              "bimodalityColCorr", 
              "linearCoefPoly2Row", "quadraticCoefPoly2Row", "coefHclustRows", 
              "intensityNAProb50.sd", "intensityNAProb90.sd", 
              "intensityNAProbnSamplesWithProbValue", 
              "prctnDistinctValues")

Newnames <- c("Dataset ID", "Data type", "Data type subgroups", 
              "# Samples", "# Analytes", "min(% NA in analytes)", 
              "max(% NA in analytes)", "min(% NA in samples)", 
              "max(% NA in samples)", 
              "% NA", "% Analytes with NAs", "% Samples with NAs", 
              "Corr(Mean vs. % NA) (Samples)", 
              "Corr(Mean vs. % NA) (Samples) (p-Value)", 
              "Corr(Mean vs. % NA) (Analytes)", 
              "Corr(Mean vs. % NA) (Analytes) (p-Value)", 
              "Mean", "Median", "Min", "Max", "median(Variance of samples)", 
              "median(Variance of analytes)", 
              "Variance", "Kurtosis", "Skewness", 
              "% Var. explained by PC1", "% Var. explained by PC2", 
              "Bimodality of sample correlations", 
              "Lin. coef. of Poly2(Means vs. Vars) (Analytes)", 
              "Quadr. coef. of Poly2(Means vs. Vars) (Analytes)", 
              "Agglom. coef. hierarch. analyte clustering", 
              "sd(Intensity w/ prob(NA) = 50% for sample)", 
              "sd(Intensity w/ prob(NA) = 90% for sample)", 
              "# Samples w/ intensityNAProb50.sd or intensityNAProb90.sd", 
              "% Distinct values")

renameTable <- data.frame(Oldnames = Oldnames,
                          Newnames = Newnames)

data <- data %>% dplyr::rename_with(~ Newnames[which(Oldnames == .x)], 
                                    .cols = Oldnames)
data <- data %>% dplyr::mutate(`|Skewness|` = abs(Skewness))

write.csv(data, "datasets_results_clean_renamed.csv", row.names = FALSE)

write.csv(data.frame(table(data[, "Data type" ])), "numberOfDatasets.csv", 
          row.names = FALSE)
write.csv(data.frame(table(data[, "Data type subgroups" ])), 
          "numberOfDatasets_subgroups.csv", row.names = FALSE)
################################################################################

allDataTypeLevels <- c("Data type", "Data type subgroups")
# selectedDataTypeLevel <- "Data type subgroups"
data.copy <- data
for (selectedDataTypeLevel in allDataTypeLevels) {
  
  data <- data.copy %>% 
    dplyr::select(-setdiff(
      !!allDataTypeLevels, !!selectedDataTypeLevel)) %>% 
    dplyr::rename("Data type" = !!selectedDataTypeLevel)
  data <- data[
    !(data$`Data type` %in% c("Metabolomics (Undefined-MS)",
                              "Metabolomics (Other ionization-MS)",
                              "Proteomics (iBAQ, PRIDE, Undefined)", 
                              "Proteomics (Intensity, PRIDE, Undefined)", 
                              "Proteomics (LFQ, PRIDE, Undefined)")),]
  
  data <- data %>% dplyr::group_by(`Data type`) %>% filter(n() > 5) %>% ungroup
  
  ################################################################################
  boxplotCols <- setdiff(
    unique(c("Dataset ID", "Data type", "# Samples", 
             "# Analytes", "min(% NA in analytes)", 
             "max(% NA in analytes)", "min(% NA in samples)", 
             "max(% NA in samples)", 
             "% NA", "% Analytes with NAs", "% Samples with NAs", "Mean", 
             "Median", "Min", "Max", "median(Variance of samples)", 
             "median(Variance of analytes)", 
             "Variance", "Kurtosis", "Skewness", "|Skewness|", 
             "% Var. explained by PC1", 
             "% Var. explained by PC2", 
             "Lin. coef. of Poly2(Means vs. Vars) (Analytes)", 
             "Quadr. coef. of Poly2(Means vs. Vars) (Analytes)", 
             "Agglom. coef. hierarch. analyte clustering", 
             "% Distinct values", "Corr(Mean vs. % NA) (Samples)", 
             "Corr(Mean vs. % NA) (Analytes)",
             "sd(Intensity w/ prob(NA) = 50% for sample)", 
             "sd(Intensity w/ prob(NA) = 90% for sample)")),
    c("Data type", "Dataset ID"))
  
  data2 <- data %>% tibble::column_to_rownames("Dataset ID") %>% 
    dplyr::select(c("Data type", !!boxplotCols))
  
  # To be log2-transformed:
  toBeLog2Transformed <- c("# Samples", "# Analytes", 
                           "min(% NA in analytes)", "max(% NA in analytes)", 
                           "% Analytes with NAs", "% Samples with NAs",
                           "median(Variance of samples)", 
                           "median(Variance of analytes)", 
                           "Variance", "Kurtosis", "|Skewness|", "Skewness",
                           "Lin. coef. of Poly2(Means vs. Vars) (Analytes)", 
                           "Quadr. coef. of Poly2(Means vs. Vars) (Analytes)",
                           "sd(Intensity w/ prob(NA) = 50% for sample)", 
                           "sd(Intensity w/ prob(NA) = 90% for sample)")
  
  colsWithNegativeNumbers <- colnames(data2[, sapply(
    data2, FUN = function(x) any(x <= 0, na.rm = TRUE))])
  # [1] "min(% NA in analytes)"                            "max(% NA in analytes)"                           
  # [3] "min(% NA in samples)"                             "max(% NA in samples)"                            
  # [5] "% NA"                                             "% Analytes with NAs"                             
  # [7] "% Samples with NAs"                               "Mean"                                            
  # [9] "Median"                                           "Min"                                             
  # [11] "Max"                                              "median(Variance of samples)"                     
  # [13] "median(Variance of analytes)"                     "Kurtosis"                                        
  # [15] "Skewness"                                         "Lin. coef. of Poly2(Means vs. Vars) (Analytes)"  
  # [17] "Quadr. coef. of Poly2(Means vs. Vars) (Analytes)" "Corr(Mean vs. % NA) (Samples)"                   
  # [19] "Corr(Mean vs. % NA) (Analytes)"                   "sd(Intensity w/ prob(NA) = 50% for sample)"      
  # [21] "sd(Intensity w/ prob(NA) = 90% for sample)"      
  
  toBeLog2Transformed <- setdiff(toBeLog2Transformed, colsWithNegativeNumbers)
  # "# Samples"  "# Analytes" "Variance"  

  for (var in toBeLog2Transformed){
    # print(var)
    data2 <- logTransform(df = data2, variable = var, logBase = "log2")
  }
  
  write.csv(data2, 
            paste0("data2_", gsub(" ", "_", selectedDataTypeLevel), ".csv"), 
            row.names = FALSE)
  data2[sapply(data2, is.infinite)] <- NA

  naRelatedCols <- c(
    "Corr(Mean vs. % NA) (Samples)", 
    "Corr(Mean vs. % NA) (Samples) (p-Value)", 
    "Corr(Mean vs. % NA) (Analytes)", 
    "Corr(Mean vs. % NA) (Analytes) (p-Value)", 
    "sd(Intensity w/ prob(NA) = 50% for sample)", 
    "sd(Intensity w/ prob(NA) = 90% for sample)",
    "# Samples w/ intensityNAProb50.sd or intensityNAProb90.sd",
    "Bimodality of sample correlations"
  )
  
  data2.complete <- data2[,!(colnames(data2) %in% naRelatedCols)]
  data2.complete <- data2.complete[complete.cases(data2.complete), ]
  
  #############################################
  
  data2.long <- reshape2::melt(data2)
  
  neworder <- c("Data type", 
                "log2(# Analytes)", "log2(# Samples)", 
                "Mean", "Median", "Min", "Max", 
                "log2(Variance)", "log2(median(Variance of samples))", 
                "log2(median(Variance of analytes))",
                "Kurtosis", "Skewness", "log2(|Skewness|)", "% Distinct values",
                "% NA",
                "min(% NA in samples)", "max(% NA in samples)",  
                "min(% NA in analytes)", "max(% NA in analytes)", 
                "% Analytes with NAs", "% Samples with NAs", 
                "sd(Intensity w/ prob(NA) = 50% for sample)", 
                "sd(Intensity w/ prob(NA) = 90% for sample)",
                "Corr(Mean vs. % NA) (Samples)", 
                "Corr(Mean vs. % NA) (Analytes)", 
                "% Var. explained by PC1", "% Var. explained by PC2", 
                "Agglom. coef. hierarch. analyte clustering",
                "Lin. coef. of Poly2(Means vs. Vars) (Analytes)", 
                "Quadr. coef. of Poly2(Means vs. Vars) (Analytes)"
  )
  
  data2.long <- dplyr::arrange(
    dplyr::mutate(data2.long, variable = factor(variable, levels = neworder)), 
    variable)
  
  numberOfDatasetsIncluded <- data2.long %>% 
    na.omit() %>%
    group_by(`Data type`, variable) %>%  tally() %>%
    ungroup()
  colnames(numberOfDatasetsIncluded) <- c("Data type", "variable", "count")
  write.csv(numberOfDatasetsIncluded, paste0("numberOfDatasetsIncluded_", 
                                             selectedDataTypeLevel, ".csv"), 
            row.names = FALSE)
  
  if (selectedDataTypeLevel == "Data type") {
    height <- 14# 12
    width <- 15# 18
  } else {
    height <- 21
    width <- 15
  }
  
  plotBoxplots(data2.long, 
               fileNameAddition = paste0(
                 "_allDataTypes_", 
                 gsub(" ", "_", selectedDataTypeLevel)), 
               height = height, width = width)
  
  plotBoxplots(data2.long, 
               fileNameAddition = paste0("_allDataTypes_ranks_", 
                                         gsub(" ", "_", selectedDataTypeLevel)), 
               height = height, width = width,
               varsToRank = unique(data2.long$variable),
               order = neworder)
  
  plotBoxplots(data2.long, 
               fileNameAddition = paste0(
                 "_allDataTypes_selectedRanks_", 
                 gsub(" ", "_", selectedDataTypeLevel)), 
               height = height, width = width, 
               varsToRank = c(
                 "Kurtosis",
                 "Skewness",
                 "Lin. coef. of Poly2(Means vs. Vars) (Analytes)",
                 "Quadr. coef. of Poly2(Means vs. Vars) (Analytes)"),
               order = neworder)
  
  #############################################
  
  # One of: "median(Variance of samples)", "log2(Variance)"
  # One of: "Mean", "Median", "Min", "Max"
  # One of : "min(% NA in samples)", "max(% NA in samples)", "% NA", "% Analytes with NAs", "% Samples with NAs", 
  
  
  selForCorr <- c("% Distinct values", 
                  "log2(# Analytes)", "log2(# Samples)", 
                  "Mean", "log2(Variance)", 
                  "% NA", "max(% NA in analytes)", 
                  "Corr(Mean vs. % NA) (Samples)", 
                  "Corr(Mean vs. % NA) (Analytes)",
                  "Skewness", "log2(|Skewness|)", "Kurtosis", 
                  "% Var. explained by PC1", 
                  "% Var. explained by PC2",  
                  "Lin. coef. of Poly2(Means vs. Vars) (Analytes)", 
                  "Quadr. coef. of Poly2(Means vs. Vars) (Analytes)", 
                  "Agglom. coef. hierarch. analyte clustering", 
                  "sd(Intensity w/ prob(NA) = 50% for sample)", 
                  "sd(Intensity w/ prob(NA) = 90% for sample)"
  )
  
  corrMethod <- "spearman"
  plotPairsPlotForTypes(df = data2, groupColName = "Data type", 
                        colsForCorr = selForCorr, 
                        corrMethod = corrMethod, width = 14, height = 14, 
                        addStr = paste0("_selectedCols_", corrMethod, "_", 
                                        gsub(" ", "_", selectedDataTypeLevel))) 
  
  #############################################
  
  margPlot <- ggplot(data, aes(x = `Corr(Mean vs. % NA) (Samples)`, 
                               y = `Corr(Mean vs. % NA) (Analytes)`, 
                               colour = `Data type`)) +
    geom_point(aes(fill = `Data type`), size = 0.8, alpha = 0.5) +
    theme_minimal() +
    theme(legend.title=element_blank())
  
  pdf(paste0("marginPlot_", gsub(" ", "_", selectedDataTypeLevel), ".pdf"), 
      width = 12, height = 10)
  print(ggExtra::ggMarginal(margPlot, groupFill = TRUE, groupColour = TRUE))
  dev.off()
  
  #############################################
  
  if (selectedDataTypeLevel == "Data type") {
    plotPCABiplots(df = data2.complete, groupColName = "Data type", 
                   addStr = gsub(" ", "_", selectedDataTypeLevel), 
                   pcaMethod = "svd") # "svd" or "nipals"
    plotPCABiplots(df = data2, groupColName = "Data type", 
                   addStr = gsub(" ", "_", selectedDataTypeLevel), 
                   pcaMethod = "nipals") # "svd" or "nipals"
    
    plotUMAPplots(df = data2, groupColName = "Data type", 
                  addStr = gsub(" ", "_", selectedDataTypeLevel))
  } 
}

session <- sessionInfo()
sink(paste0("visualizations_sessionInfo.txt"))
print(session)
sink()
