setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(plyr)
library(dplyr)
library(ggplot2)
library(GGally)
library(ggforce)
library(stringr)
library(tibble)

# devtools::install_github("vqv/ggbiplot")
library(ggbiplot)
theme_set(theme_bw())

library(patchwork)
library(plotly)
library(htmlwidgets)

library(pcaMethods)
################################################################################
################################################################################
# FUNCTIONS

logTransform <- function(df, variable, logBase = c("log2", "log1p")){
  df[[paste0(logBase, "(", variable, ")")]] <- get(logBase)(df[[variable]])
  df[[variable]] <- NULL
  df
}

plotBoxplots <- function(data2.long, fileNameAddition = "", height = 12, 
                         width = 18, varsToRank = c(), order = NULL,
                         customColors = c()) {
  
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
    geom_violin(aes(fill = `Data type`), alpha = 0.9, scale = "width") +
    geom_boxplot(aes(fill = `Data type`), width = 0.5, alpha = 0.9, 
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
  
  if (length(customColors) > 0) {
    ggplot.charact <- ggplot.charact +
      ggplot2::scale_fill_manual(values = customColors)
  }
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
                          addStr = "",
                          groupColName = "",
                          customColors = customColors,
                          axisLabelSize = 10,
                          axisTitleSize = 12) {
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
                           alpha = 0,
                           varname.adjust = 1.25)  +
    scale_color_discrete(name = '') +  
    coord_fixed(ratio = coordRatio) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(alpha = 1, size = 3)))
  
    if (groupColName == "Data type") {
      P2 <- P2 + ggplot2::scale_colour_manual(values = customColors) +
        ggplot2::theme(legend.title = element_blank())
    }
  
  #if (all(PCchoices == 1:2)) P2 <- P2 + 
  #  ggplot2::scale_y_continuous(limits = c(-5, 10.3))
  
  # For nipals:
  # row.names(PCA.data)[which(PCA.data$PC2<(-5) | PCA.data$PC2>10.3)] 
  # [1] "MTBLS2801_m_MTBLS2801_LC-MS_negative_reverse-phase_metabolite_profiling_v2_maf.tsv"
  # [2] "MTBLS2801_m_MTBLS2801_LC-MS_positive_reverse-phase_metabolite_profiling_v2_maf.tsv"
  # [3] "MTBLS4186_m_MTBLS4186_LC-MS_positive_reverse-phase_metabolite_profiling_v2_maf.tsv"
  
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
                   legend.position = 'bottom',
                   axis.text = element_text(size = axisLabelSize),
                   axis.title = element_text(size = axisTitleSize))
  
  P2$layers <- c(geom_point(aes(colour = groups), cex = 1, alpha = alpha), 
                 P2$layers)
  
  P2
}

plotPCABiplots <- function(df, groupColName = "", addStr = "", 
                           pcaMethod = "nipals",
                          coordRatio = 0.5,
                          xlimLower = NA, xlimUpper = NA,
                          ylimLower = NA, ylimUpper = NA,
                          plotWidth = 12, plotheight = 9,
                          axisLabelSize = 10,
                          axisTitleSize = 12) {
  
  customColors <- c("Metabolomics (NMR)" = "#4BADF1", 
                    "Metabolomics (MS)" = "#0033CC",
                    "Lipidomics (MS)" = "#000000",
                    
                    "Proteomics (LFQ, PRIDE)" = "#800000",
                    "Proteomics (Intensity, PRIDE)" =  "#DC143C",
                    "Proteomics (iBAQ, PRIDE)" = "#DE4B7E",
                    "Proteomics (Intensity, Expression Atlas)" = "#FF6600",
                    "Proteomics (iBAQ, Expression Atlas)" = "#FFC0CB",
                    
                    "RNA-seq (raw)" = "#2E5F72",
                    "RNA-seq (FPKM)" = "#467741",
                    "RNA-seq (TPM)" = "#32CD32",
                    "Microarray" = "#5BB3B1",
                    
                    "scProteomics" = "#806FC4",
                    "scRNA-seq (unnormalized)" = "#DDA0DD",
                    "scRNA-seq (normalized)" = "#BA30B5",
                    "Microbiome" = "#6911D3"
  )
  
  pdf(file = paste0("biplot_", pcaMethod, "_", addStr, ".pdf"),
      width = plotWidth, height = plotheight)
  print(plotPCABiplot(df = df %>% dplyr::select(-!!groupColName), 
                      groups = df[[groupColName]],
                      alpha = 0.3,
                      pcaMethod = pcaMethod,
                      coordRatio = coordRatio, 
                      xlimLower = xlimLower, 
                      xlimUpper = xlimUpper,
                      ylimLower = ylimLower, 
                      ylimUpper = ylimUpper,
                      facetZoom = FALSE,
                      addStr = addStr,
                      groupColName = groupColName,
                      customColors = customColors,
                      axisLabelSize = axisLabelSize,
                      axisTitleSize = axisTitleSize))
  dev.off()
  
  pca <- pcaMethods::pca(df %>% dplyr::select(-!!groupColName), 
                         method = pcaMethod, nPcs = 4, center = TRUE
                         , scale = "uv")
  write.csv(pca@loadings, paste0("loadings_",  pcaMethod, "_", addStr,".csv"))
}


plotUMAPplots <- function(df, groupColName = "", addStr = "", pointAlpha = 1,
                          pointSize = 0.8#, 
                          # pointStroke = 0.3
                          ) {
  
  df <- df[, apply(df, 2, function(x) {length(unique(x)) > 1}) ]
  
  customColors <- c("Metabolomics (NMR)" = "#4BADF1", 
                    "Metabolomics (MS)" = "#0033CC",
                    "Lipidomics (MS)" = "#000000",
                    
                    "Proteomics (LFQ, PRIDE)" = "#800000",
                    "Proteomics (Intensity, PRIDE)" =  "#DC143C",
                    "Proteomics (iBAQ, PRIDE)" = "#DE4B7E",
                    "Proteomics (Intensity, Expression Atlas)" = "#FF6600",
                    "Proteomics (iBAQ, Expression Atlas)" = "#FFC0CB",
                    
                    "RNA-seq (raw)" = "#2E5F72",
                    "RNA-seq (FPKM)" = "#467741",
                    "RNA-seq (TPM)" = "#32CD32",
                    "Microarray" = "#5BB3B1",
                    
                    "scProteomics" = "#806FC4",
                    "scRNA-seq (unnormalized)" = "#DDA0DD",
                    "scRNA-seq (normalized)" = "#BA30B5",
                    "Microbiome" = "#6911D3"
  )
  
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
  
  gg12 <- ggplot2::ggplot(umap_df, ggplot2::aes(
    x = UMAP1, y = UMAP2, color = get(groupColName))) +
    ggplot2::geom_point(size = pointSize, alpha = pointAlpha#, 
                        #shape = 21, aes(stroke = pointStroke)
                        )  + 
    ggplot2::theme(legend.direction = 'horizontal', 
                   legend.position = 'bottom') +
    ggplot2::scale_color_discrete(name = '') +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(alpha = 1, size = 3)))
  
  
  umapDataTypes <- lapply(
    levels(groupVec),
    function(x, umap_df, groupColName) {
      ggplot2::ggplot(umap_df[umap_df[groupColName] == as.character(x), ], ggplot2::aes(
        x = UMAP1, y = UMAP2)) +
        ggplot2::geom_point(size = 0.7, alpha = 0.25)  + 
        ggplot2::theme(legend.direction = 'horizontal', 
                       legend.position = 'bottom') +
        ggplot2::scale_color_discrete(name = '') +
        ggplot2::guides(colour = ggplot2::guide_legend(
          override.aes = list(alpha = 1, size = 3))) +
        xlim(min(umap_df$UMAP1), max(umap_df$UMAP1)) + 
        ylim(min(umap_df$UMAP2), max(umap_df$UMAP2)) +
        ggplot2::ggtitle(as.character(x))
    }, umap_df = umap_df, groupColName = groupColName)
  
  
  umapDataTypes <- list(
    umapDataTypes[[1]], umapDataTypes[[2]], umapDataTypes[[3]],
    patchwork::plot_spacer(), patchwork::plot_spacer(), 
    umapDataTypes[[4]], umapDataTypes[[5]], umapDataTypes[[6]],
    umapDataTypes[[7]], umapDataTypes[[8]],
    umapDataTypes[[9]], umapDataTypes[[10]], umapDataTypes[[11]], umapDataTypes[[12]],
    patchwork::plot_spacer(),
    umapDataTypes[[13]], umapDataTypes[[14]], umapDataTypes[[15]], umapDataTypes[[16]])
  ggsave(file = paste0("umap_Dim1vsDim2_", addStr, "_dataTypes.pdf"), 
         patchwork::wrap_plots(umapDataTypes, ncol = 5), 
         width = 20, height = 12)
  
  
  gg23 <- ggplot2::ggplot(umap_df, ggplot2::aes(x = UMAP2, y = UMAP3, 
                                        color = get(groupColName))) +
    ggplot2::geom_point(size = pointSize, alpha = pointAlpha#, 
                        #shape = 21, aes(stroke = pointStroke)
                        )  + 
    ggplot2::theme(legend.direction = 'horizontal', 
                   legend.position = 'bottom') +
    ggplot2::scale_color_discrete(name = '') +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(alpha = 1, size = 3)))
  
  if (groupColName == "Data type") {
    gg12 <- gg12 + ggplot2::scale_colour_manual(values = customColors) +
      ggplot2::theme(legend.title = element_blank())
    gg23 <- gg23 + ggplot2::scale_colour_manual(values = customColors) +
      ggplot2::theme(legend.title = element_blank())
  }
  
  pdf(file = paste0("umap_Dim1vsDim2_", addStr, ".pdf"), width = 12, height = 9)
  print(gg12)
  dev.off()
  
  pdf(file = paste0("umap_Dim2vsDim3_", addStr, ".pdf"), width = 12, height = 9)
  print(gg23)
  dev.off()
  
  gg12
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

proteomics.replaced <- strsplit(prideMeta.df.forIncludedPrideIDs$instruments, ";")

proteomics.replaced2 <- rapply(proteomics.replaced, function(x){
  x[grepl("Q Exactive", x)] <- "QExactive"
  x[grepl("Orbitrap", x)] <- "Orbitrap"
  x[grepl("timsTOF", x)] <- "timsTOF"
  x[grepl("TripleTOF", x)] <- "TripleTOF"
  x <- sort(unique(x))
  x
}, how = "list") 

replaced2 <- lapply(proteomics.replaced2, function(x) paste0(x, collapse = ";"))

proteomicsReplaced.df <- data.frame(technology = unlist(replaced2))
  
proteomicsReplaced.df$technology[grepl(";", proteomicsReplaced.df$technology)] <- 'Multiple'  

proteomicsReplaced.df$technology[
  proteomicsReplaced.df$technology %in% 
    names(which(table(proteomicsReplaced.df$technology) < 6))] <- 'Other'  

proteomicsReplaced.df <- data.frame(prideID = prideMeta.df.forIncludedPrideIDs$prideID, technology = proteomicsReplaced.df$technology)

dataset <- dataset %>% 
  dplyr::mutate(dataTypeSubgroups = case_when(
    grepl("_pride_", dataType) ~ 
      paste0(dataTypeSubgroups, "_", 
             proteomicsReplaced.df$technology[
               match(sub("^(.*?)_.*$", "\\1", datasetID), 
                     proteomicsReplaced.df$prideID)]),
    TRUE ~ dataTypeSubgroups)
  )

metabolomicsMeta.df <- read.csv("metabolomics_metadata.csv")
metabolomicsDatasets <- dataset %>% 
  filter(dataType %in% c("metabolomics_MS", "metabolomics_NMR")) %>%
  dplyr::select(datasetID, dataType, dataTypeSubgroups) %>%
  dplyr::mutate(accession = sub("^(.*?)_.*$", "\\1", datasetID)) %>%
  dplyr::left_join(metabolomicsMeta.df, by = join_by(accession == accession))

replaced <- strsplit(metabolomicsDatasets$technology, ";")

replacedTechnologies <- rapply(replaced, function(x){
  x[grepl("-(qTOF|TOF)-", x, ignore.case = TRUE)] <- "TOF"
  x[grepl("-(LTQ)-", x, ignore.case = TRUE)] <- "LTQ"
  x[grepl("-(TQ)-", x, ignore.case = TRUE)] <- "TQ"
  x[grepl("-(Q)-", x, ignore.case = TRUE)] <- "Q"
  x[grepl("-(FT-ICR)-", x, ignore.case = TRUE)] <- "FTICR"  
  x[grepl("UPLC-MS", x, ignore.case = TRUE)] <- NA  
  x[grepl("Insufficient data supplied", x)] <- NA
  x[grepl("NMR", x)] <- NA
  x <- sort(unique(x))
  x
}, how = "list") 

replacedTechnologies <- lapply(replacedTechnologies, function(x) paste0(x, collapse = ";"))


replacedTechnologies.df <- data.frame(technology = unlist(replacedTechnologies))
replacedTechnologies.df$technology[is.na(replacedTechnologies.df$technology) | replacedTechnologies.df$technology == ""] <- "Undefined"
replacedTechnologies.df$technology[grepl(";", replacedTechnologies.df$technology)] <- 'Multiple'  
replacedTechnologies.df$technology[
  replacedTechnologies.df$technology %in% 
    names(which(table(replacedTechnologies.df$technology) < 6))] <- 'Other'  

#metabolomicsDatasets$technology3<- replacedTechnologies.df$technology



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
  technology2 = unlist(replaced2),
  technology3 = replacedTechnologies.df$technology)

"Fill information from met data if not present in file name"
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


metabolomicsTranslation.df[
  metabolomicsTranslation.df$datasetID %in% 
    c("MTBLS103_m_ibanez_02_metabolite_profiling_mass_spectrometry_v2_maf.tsv", 
      "MTBLS794_m_mitochondrial_deficiency_metabolite_profiling_mass_spectrometry_v2_maf.tsv"),]$dataTypeSubgroupsNew <- "metabolomics_MS_LC"

metabolomicsTranslation.df$dataTypeSubgroupsNew2 <- paste0(metabolomicsTranslation.df$dataTypeSubgroupsNew, "_", metabolomicsTranslation.df$technology3)
metabolomicsTranslation.df$dataTypeSubgroupsNew2[grepl("metabolomics_NMR", metabolomicsTranslation.df$dataTypeSubgroupsNew2)] <-
  "metabolomics_NMR"
metabolomicsTranslation.df$dataTypeSubgroupsNew2[
  !grepl("NMR|GC|LC", metabolomicsTranslation.df$dataTypeSubgroupsNew2) |
    grepl("Undefined|Multiple|Other", metabolomicsTranslation.df$dataTypeSubgroupsNew2) |
    (metabolomicsTranslation.df$dataTypeSubgroupsNew2 %in% 
       names(which(table(metabolomicsTranslation.df$dataTypeSubgroupsNew2) < 6)))
  ] <-
  "metabolomics_MS_Other"


dataset <- dataset %>%
  dplyr::mutate(
    dataTypeSubgroups = case_when(
      grepl("metabolomics_", dataType) ~ 
        metabolomicsTranslation.df$dataTypeSubgroupsNew2[match(
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
# "MTBLS555_m_test_131017_metabolite_profiling_mass_spectrometry_v2_maf.tsv"
# "PXD019903_proteinGroups.txt_^LFQ_"                                       
# "PXD019903_proteinGroups.txt_^iBAQ_"                                      
# "PXD019903_proteinGroups.txt_^Intensity_"   

data <- data[!(
  data$datasetID %in% 
    c("PXD006847_CulturedCells_proteinGroups.txt_^LFQ_",
      "PXD006847_CulturedCells_proteinGroups.txt_^iBAQ_",
      "PXD006847_CulturedCells_proteinGroups.txt_^Intensity_",
      "PXD021882_proteinGroups.txt_^LFQ_",
      "PXD021882_proteinGroups.txt_^Intensity_",
      # "MTBLS103_m_ibanez_02_metabolite_profiling_mass_spectrometry_v2_maf.tsv",
      "MTBLS555_m_test_131017_metabolite_profiling_mass_spectrometry_v2_maf.tsv",
      "PXD019903_proteinGroups.txt_^LFQ_",                                      
      "PXD019903_proteinGroups.txt_^iBAQ_",                                      
      "PXD019903_proteinGroups.txt_^Intensity_")),]


# occurrences <- data %>%
#   group_by(across(-c(datasetID, dataTypeSubgroups))) %>%
#   summarise(Occurrences = n(), .groups = "drop")
# 
# result <- data %>%
#   left_join(occurrences, by = setdiff(names(data), c("datasetID", "dataTypeSubgroups")))
# result <- result %>%
#   select(datasetID, dataTypeSubgroups, Occurrences)


##  Duplicate rows removed
data <- data[-which(duplicated(data %>% dplyr::select(-datasetID))), ]

# Add lipidomics MS information
lipidomicsDatasets <- read.csv("lipidomicsMSDatasets.csv")$lipidomicsDatasets

data <- data %>%
  mutate(
    dataType = ifelse(datasetID %in% lipidomicsDatasets, 
                      gsub("metabolomics", "lipidomics", dataType), 
                      dataType),
    dataTypeSubgroups = ifelse(datasetID %in% lipidomicsDatasets, 
                               gsub("metabolomics", "lipidomics", dataTypeSubgroups), 
                               dataTypeSubgroups)
  )

# Add lipidomics NMR information
data[data$datasetID == "MTBLS1894_m_MTBLS1894_organic_NMR_metabolite_profiling_v2_maf.tsv",]$dataType <- "lipidomics_NMR"
data[data$datasetID == "MTBLS1894_m_MTBLS1894_organic_NMR_metabolite_profiling_v2_maf.tsv",]$dataTypeSubgroups <- "lipidomics_NMR"


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
    OldDataTypeNames <- c("lipidomics_MS", "lipidomics_NMR",
                          "metabolomics_MS", "metabolomics_NMR", 
                          "microarray", "microbiome", 
                          "proteomics_expressionatlas_iBAQ", 
                          "proteomics_expressionatlas_Intensity", 
                          "proteomics_pride_LFQ", "proteomics_pride_Intensity",
                          "proteomics_pride_iBAQ", 
                          "RNAseq_fpkms_median", "RNAseq_raw", 
                          "RNAseq_tpms_median", "sc_normalized", 
                          "sc_unnormalized", "scProteomics")
    NewDataTypeNames <- c("Lipidomics (MS)", "Lipidomics (NMR)",
                          "Metabolomics (MS)", "Metabolomics (NMR)", 
                          "Microarray", "Microbiome", 
                          "Proteomics (iBAQ, Expression Atlas)", 
                          "Proteomics (Intensity, Expression Atlas)", 
                          "Proteomics (LFQ, PRIDE)", 
                          "Proteomics (Intensity, PRIDE)", 
                          "Proteomics (iBAQ, PRIDE)", 
                          "RNA-seq (FPKM)", "RNA-seq (raw)", "RNA-seq (TPM)", 
                          "scRNA-seq (normalized)", 
                          "scRNA-seq (unnormalized)", "scProteomics")
    levels <- c('Metabolomics (NMR)', 
                'Metabolomics (MS)',
                
                "Lipidomics (NMR)",
                "Lipidomics (MS)",
                
                'Proteomics (LFQ, PRIDE)',
                'Proteomics (Intensity, PRIDE)',
                'Proteomics (iBAQ, PRIDE)',
                'Proteomics (Intensity, Expression Atlas)',
                'Proteomics (iBAQ, Expression Atlas)', 
                
                'RNA-seq (raw)',
                'RNA-seq (FPKM)',
                'RNA-seq (TPM)',
                'Microarray',
                
                'scProteomics',
                'scRNA-seq (unnormalized)',
                'scRNA-seq (normalized)',
                'Microbiome'
    )
  } else if (dataTypeLevel == "dataTypeSubgroups") {
    OldDataTypeNames <- c(
      "lipidomics_NMR", 
      "lipidomics_MS_LC_TOF", "lipidomics_MS_LC_LTQ", "lipidomics_MS_LC_TQ", 
      "lipidomics_MS_GC_TOF", "lipidomics_MS_GC_Q",
      "lipidomics_MS_Other",
      
      "metabolomics_NMR", 
      "metabolomics_MS_LC_TOF", "metabolomics_MS_LC_LTQ", "metabolomics_MS_LC_TQ",
      "metabolomics_MS_GC_TOF", "metabolomics_MS_GC_TQ",  "metabolomics_MS_GC_Q",                            
      "metabolomics_MS_Other", 
      
      "microarray_Affymetrix", 
      "microarray_Illumina", "microarray_Agilent", "microbiome_16S", 
      "microbiome_WGS", "proteomics_expressionatlas_iBAQ", 
      "proteomics_expressionatlas_Intensity", 
      "proteomics_expressionatlas_LFQ", 
      
      "proteomics_pride_iBAQ_Multiple", "proteomics_pride_iBAQ_Orbitrap", 
      "proteomics_pride_iBAQ_Other", "proteomics_pride_iBAQ_QExactive", 
      "proteomics_pride_iBAQ_timsTOF", "proteomics_pride_iBAQ_TripleTOF", 
      "proteomics_pride_Intensity_maXis", "proteomics_pride_Intensity_Multiple", 
      "proteomics_pride_Intensity_Orbitrap", "proteomics_pride_Intensity_Other", 
      "proteomics_pride_Intensity_QExactive", "proteomics_pride_Intensity_timsTOF", 
      "proteomics_pride_Intensity_TripleTOF", "proteomics_pride_LFQ_maXis", 
      "proteomics_pride_LFQ_Multiple", "proteomics_pride_LFQ_Orbitrap", 
      "proteomics_pride_LFQ_Other", "proteomics_pride_LFQ_QExactive", 
      "proteomics_pride_LFQ_timsTOF", "proteomics_pride_LFQ_TripleTOF", 
      
      "RNAseq_fpkms_median", 
      "RNAseq_raw", 
      "RNAseq_tpms_median", 
      "sc_normalized_SMART-like",
      "sc_normalized_Droplet-based", "sc_unnormalized_SMART-like", 
      "sc_unnormalized_Droplet-based", "scProteomics")
    NewDataTypeNames <- c(
      
      "Lipidomics (NMR)", 
      "Lipidomics (LC-MS, TOF)", "Lipidomics (LC-MS, LTQ)", "Lipidomics (LC-MS, TQ)", 
      "Lipidomics (GC-MS, TOF)", "Lipidomics (GC-MS, Q)", 
      "Lipidomics (MS, Other)",
      
      "Metabolomics (NMR)", 
      "Metabolomics (LC-MS, TOF)", "Metabolomics (LC-MS, LTQ)", "Metabolomics (LC-MS, TQ)",
      "Metabolomics (GC-MS, TOF)", "Metabolomics (GC-MS, TQ)", "Metabolomics (GC-MS, Q)",                             
      "Metabolomics (MS, Other)",
      
      "Microarray (Affymetrix)", 
      "Microarray (Illumina)", "Microarray (Agilent)", "Microbiome (16S)", 
      "Microbiome (WGS)", 
      "Proteomics (iBAQ, Expression Atlas)", 
      "Proteomics (Intensity, Expression Atlas)", 
      "Proteomics (LFQ, Expression Atlas)", 
      
      "Proteomics (iBAQ, PRIDE, Multiple)", "Proteomics (iBAQ, PRIDE, Orbitrap)", 
      "Proteomics (iBAQ, PRIDE, Other)", "Proteomics (iBAQ, PRIDE, Q Exactive)", 
      "Proteomics (iBAQ, PRIDE, timsTOF)", "Proteomics (iBAQ, PRIDE, TripleTOF)", 
      "Proteomics (Intensity, PRIDE, maXis)", "Proteomics (Intensity, PRIDE, Multiple)",
      "Proteomics (Intensity, PRIDE, Orbitrap)", "Proteomics (Intensity, PRIDE, Other)",
      "Proteomics (Intensity, PRIDE, Q Exactive)", "Proteomics (Intensity, PRIDE, timsTOF)",
      "Proteomics (Intensity, PRIDE, TripleTOF)", "Proteomics (LFQ, PRIDE, maXis)",
      "Proteomics (LFQ, PRIDE, Multiple)", "Proteomics (LFQ, PRIDE, Orbitrap)",
      "Proteomics (LFQ, PRIDE, Other)", "Proteomics (LFQ, PRIDE, Q Exactive)",
      "Proteomics (LFQ, PRIDE, timsTOF)", "Proteomics (LFQ, PRIDE, TripleTOF)",
      
      "RNA-seq (FPKM)", 
      "RNA-seq (raw)", 
      "RNA-seq (TPM)", 
      "scRNA-seq (SMART-like, normalized)",
      "scRNA-seq (Droplet-based, normalized)", 
      "scRNA-seq (SMART-like, unnormalized)",
      "scRNA-seq (Droplet-based, unnormalized)", "scProteomics")
    
    levels <- c(
      "Metabolomics (NMR)", 
      "Metabolomics (GC-MS, TOF)", 
      "Metabolomics (GC-MS, TQ)", 
      "Metabolomics (GC-MS, Q)",   
      "Metabolomics (LC-MS, TOF)", 
      "Metabolomics (LC-MS, LTQ)", 
      "Metabolomics (LC-MS, TQ)",
      "Metabolomics (MS, Other)",
      
      "Lipidomics (NMR)", 
      "Lipidomics (GC-MS, TOF)",
      "Lipidomics (GC-MS, Q)", 
      "Lipidomics (LC-MS, TOF)",  
      "Lipidomics (LC-MS, LTQ)", "Lipidomics (LC-MS, TQ)", 
      "Lipidomics (MS, Other)",
      
      "Proteomics (LFQ, PRIDE, Orbitrap)", "Proteomics (Intensity, PRIDE, Orbitrap)", 
      "Proteomics (iBAQ, PRIDE, Orbitrap)", 
      
      "Proteomics (LFQ, PRIDE, Q Exactive)", "Proteomics (Intensity, PRIDE, Q Exactive)", 
      "Proteomics (iBAQ, PRIDE, Q Exactive)",
      
      "Proteomics (LFQ, PRIDE, maXis)", "Proteomics (Intensity, PRIDE, maXis)", 
      
      "Proteomics (LFQ, PRIDE, TripleTOF)", "Proteomics (Intensity, PRIDE, TripleTOF)", 
      "Proteomics (iBAQ, PRIDE, TripleTOF)",  
      
      "Proteomics (LFQ, PRIDE, timsTOF)", "Proteomics (Intensity, PRIDE, timsTOF)",
      "Proteomics (iBAQ, PRIDE, timsTOF)", 
      
      "Proteomics (LFQ, PRIDE, Multiple)", "Proteomics (Intensity, PRIDE, Multiple)",
      "Proteomics (iBAQ, PRIDE, Multiple)", 
      
      "Proteomics (LFQ, PRIDE, Other)", 
      "Proteomics (Intensity, PRIDE, Other)",
      "Proteomics (iBAQ, PRIDE, Other)",  

      "Proteomics (LFQ, Expression Atlas)", 
      "Proteomics (Intensity, Expression Atlas)", 
      "Proteomics (iBAQ, Expression Atlas)",
      
      "RNA-seq (raw)", 
      "RNA-seq (FPKM)", 
      "RNA-seq (TPM)", 
      
      "Microarray (Affymetrix)", "Microarray (Illumina)", 
      "Microarray (Agilent)", 
      
      "scProteomics",
      
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
    !(data$`Data type` %in% c("Metabolomics (MS, Other)",
                              "Lipidomics (MS, Other)",
                              "Proteomics (LFQ, PRIDE, Multiple)", 
                              "Proteomics (LFQ, PRIDE, Other)", 
                              "Proteomics (Intensity, PRIDE, Multiple)",
                              "Proteomics (Intensity, PRIDE, Other)",
                              "Proteomics (iBAQ, PRIDE, Multiple)", 
                              "Proteomics (iBAQ, PRIDE, Other)"
                              )),]
  
  data <- data %>% dplyr::group_by(`Data type`) %>% filter(n() > 5) %>% ungroup
  data$`Data type` <- droplevels(data$`Data type`)
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

  for (var in toBeLog2Transformed) {
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
  # The following datasets are removed in data2.complete compared to data2:
  # [1] "PXD029255_200303_HeLa_IAA-Desthio_10_proteinGroups.txt_^LFQ_"      
  # [2] "PXD029255_200303_HeLa_IAA-Desthio_10_proteinGroups.txt_^Intensity_"
  
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
    height <- 16# 12
    width <- 15# 18
  } else {
    height <- 28
    width <- 19
  }
  
  customColors <- c()
  if (selectedDataTypeLevel == "Data type") {
    customColors <- c("Metabolomics (NMR)" = "#4BADF1", 
                      "Metabolomics (MS)" = "#0033CC",
                      "Lipidomics (MS)" = "#000000",
                      
                      "Proteomics (LFQ, PRIDE)" = "#800000",
                      "Proteomics (Intensity, PRIDE)" =  "#DC143C",
                      "Proteomics (iBAQ, PRIDE)" = "#DE4B7E",
                      "Proteomics (Intensity, Expression Atlas)" = "#FF6600",
                      "Proteomics (iBAQ, Expression Atlas)" = "#FFC0CB",
                      
                      "RNA-seq (raw)" = "#2E5F72",
                      "RNA-seq (FPKM)" = "#467741",
                      "RNA-seq (TPM)" = "#32CD32",
                      "Microarray" = "#5BB3B1",
                      
                      "scProteomics" = "#806FC4",
                      "scRNA-seq (unnormalized)" = "#DDA0DD",
                      "scRNA-seq (normalized)" = "#BA30B5",
                      "Microbiome" = "#6911D3"
    )
  }
  
  plotBoxplots(data2.long, 
               fileNameAddition = paste0(
                 "_allDataTypes_", 
                 gsub(" ", "_", selectedDataTypeLevel)), 
               height = height, width = width,
               customColors = customColors)
  
  plotBoxplots(data2.long, 
               fileNameAddition = paste0("_allDataTypes_ranks_", 
                                         gsub(" ", "_", selectedDataTypeLevel)), 
               height = height, width = width,
               varsToRank = unique(data2.long$variable),
               order = neworder,
               customColors = customColors)
  
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
               order = neworder,
               customColors = customColors)
  
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
  
  if (selectedDataTypeLevel == "Data type") {
    
    # data2.nonVarColsRemoved <- data2 %>% select(-'Data type')
    # data2.nonVarColsRemoved <- 
    #   data2.nonVarColsRemoved[, apply(data2.nonVarColsRemoved, 
    #                                   2, function(x) {length(unique(x)) > 1}) ]
    # 
    # corm <- cor(data2.nonVarColsRemoved)
    # out <- as.data.frame(apply(
    #   corm, 2, function(x) ifelse(abs(x) > 0.9, round(x, 3), "-")))
    # 
    # # Out: min(%NA in samples), max(%NA in samples), %Analytes with NAs, %Samples with NAs
    # remove <- c("max(% NA in analytes)", "min(% NA in samples)", "max(% NA in samples)", 
    #             "% Analytes with NAs", "Median", "Max", "log2(median(Variance of samples))", 
    #             "log2(median(Variance of analytes))")
    
    remove <- c("min(% NA in analytes)", "max(% NA in analytes)", 
                "min(% NA in samples)", "max(% NA in samples)",
                "% Analytes with NAs", 
                "% Samples with NAs",
                "Median", "Max", "log2(median(Variance of samples))", 
                            "log2(median(Variance of analytes))")
    
    plotPCABiplots(df = data2.complete, groupColName = "Data type", 
                   addStr = gsub(" ", "_", selectedDataTypeLevel), 
                   pcaMethod = "svd",
                   ylimLower = -5, ylimUpper = 10.3
                   )
    plotPCABiplots(df = data2, groupColName = "Data type", 
                   addStr = gsub(" ", "_", selectedDataTypeLevel), 
                   pcaMethod = "nipals",
                   ylimLower = -5, ylimUpper = 10.3
                   ) 
    
    umap12 <- plotUMAPplots(df = data2.complete, groupColName = "Data type", 
                  addStr = gsub(" ", "_", selectedDataTypeLevel),
                  pointSize = 0.5,
                  pointAlpha = 0.5)
    
    umap_data <- umap12[["data"]]

    data2.NAsamplesRemoved <- data2
    data2.NAsamplesRemoved <- data2.NAsamplesRemoved[
      !(row.names(data2.NAsamplesRemoved) %in% 
          setdiff(row.names(data2), row.names(data2.complete))),]
    
      
    umapColorPlotsLst <- lapply(
      setdiff(colnames(data2.NAsamplesRemoved), "Data type"),
      function(x, umap_data, data2.NAsamplesRemoved) {
        ggplot(umap_data, aes(x = UMAP1, y = UMAP2, 
                              color = data2.NAsamplesRemoved[, x])) +
          geom_point(aes(stroke = 0.1), 
                     # alpha = 0.25, 
                     size = 0.5, shape = 21) +
          labs(title = x) +
          scale_colour_continuous(type = "viridis") +
          theme_bw() +
          theme(legend.title = element_blank(),
                legend.position = "bottom",
                legend.key.size = unit(0.8, "cm"),
                legend.key.height = unit(0.3, "cm"),
                legend.text = element_text(size = 11),
                plot.title = element_text(size = 15),
                axis.title = element_text(size = 11))
      }, umap_data = umap_data, data2.NAsamplesRemoved = data2.NAsamplesRemoved)
    
    ggsave(filename = "umap_Dim1vsDim2_Data_type_DataCharColors.pdf", 
           patchwork::wrap_plots(umapColorPlotsLst, ncol = 5), 
           height = 28, width = 27)
    
    
    #"% Analytes with NAs", "% Samples with NAs", "Corr(Mean vs. % NA) (Samples)",
    #"% Distinct values", "log2(# Samples)", "log2(# Analytes)" 
    ggsave(filename = "umap_Dim1vsDim2_Data_type_DataCharColors_selected.pdf", 
           patchwork::wrap_plots(umapColorPlotsLst[c(6, 7, 20, 25, 24, 19)], ncol = 3), 
           height = 9, width = 15)
    
  #   supersets <- list(metabolomics = c("Metabolomics (NMR)", 
  #                                      "Metabolomics (MS)"),
  #                     proteomics = c("Proteomics (LFQ, PRIDE)",
  #                                    "Proteomics (Intensity, PRIDE)",
  #                                    "Proteomics (iBAQ, PRIDE)",
  #                                    "Proteomics (Intensity, Expression Atlas)", 
  #                                    "Proteomics (iBAQ, Expression Atlas)"),
  #                     transcriptomics = c("RNA-seq (raw)",
  #                                         "RNA-seq (FPKM)",
  #                                         "RNA-seq (TPM)",
  #                                         "Microarray"),
  #                     sc = c("scProteomics",
  #                            "scRNA-seq (unnormalized)",
  #                            "scRNA-seq (normalized)",
  #                            "Microbiome"))
  #   
 

  }
}

session <- sessionInfo()
sink(paste0("visualizations_sessionInfo.txt"))
print(session)
sink()
