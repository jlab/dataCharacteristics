setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(plyr)
library(dplyr)
library(ggplot2)
library(GGally)
library(ggforce)
theme_set(theme_bw())

dataset <- ldply(list.files("results", pattern = ".csv", full.names = TRUE), read.csv, header=TRUE)
#dataset <- plyr::ldply(list.files(pattern = ".csv", full.names = TRUE), read.csv, header=TRUE)
dataset <- dataset[dataset$nSamples>0,]
dataset <- dataset[dataset$variance>0,]
names(dataset) <- gsub(x = names(dataset), pattern = "\\.log2|\\.wNAs", replacement = "")  
names(dataset)[names(dataset) == 'corColR'] <- 'corSampleMeanNA'
names(dataset)[names(dataset) == 'corRowR'] <- 'corAnalyteMeanNA'

dataset <- dataset %>% 
  mutate(dataType2 = case_when(grepl("\\^iBAQ", datasetID) ~ paste0(dataType, "_iBAQ"),
                               grepl("\\^LFQ", datasetID) ~ paste0(dataType, "_LFQ"),
                               grepl("\\^Intensity", datasetID) ~ paste0(dataType, "_Intensity"),
                               grepl("LC", datasetID, ignore.case = TRUE) & dataType == "metabolomics_MS" ~ paste0(dataType, "_LC"),
                               grepl("GC", datasetID, ignore.case = TRUE) & dataType == "metabolomics_MS" ~ paste0(dataType, "_GC"),
                               TRUE ~ dataType),
         # dataType3 = case_when(grepl("LC", datasetID, ignore.case = TRUE) & 
         #                         grepl("pos", datasetID, ignore.case = TRUE) & 
         #                         dataType == "metabolomics_MS" ~ paste0(dataType, "_LC_pos"),
         #                       grepl("LC", datasetID, ignore.case = TRUE) & 
         #                         grepl("neg", datasetID, ignore.case = TRUE) & 
         #                         dataType == "metabolomics_MS" ~ paste0(dataType, "_LC_neg"),
         #                       grepl("GC", datasetID, ignore.case = TRUE) & 
         #                         grepl("pos", datasetID, ignore.case = TRUE) & 
         #                         dataType == "metabolomics_MS" ~ paste0(dataType, "_GC_pos"),
         #                       grepl("GC", datasetID, ignore.case = TRUE) & 
         #                         grepl("neg", datasetID, ignore.case = TRUE) & 
         #                         dataType == "metabolomics_MS" ~ paste0(dataType, "_GC_neg"),
         #                       TRUE ~ dataType4),
         dataType4 = case_when(grepl("\\^iBAQ", datasetID) ~ paste0(dataType, "_iBAQ"),
                               grepl("\\^LFQ", datasetID) ~ paste0(dataType, "_LFQ"),
                               grepl("\\^Intensity", datasetID) ~ paste0(dataType, "_Intensity"),
                               TRUE ~ dataType),
         dataType5 = case_when(grepl("\\^iBAQ", datasetID) ~ paste0(dataType, "_iBAQ"),
                               grepl("\\^LFQ", datasetID) ~ paste0(dataType, "_LFQ"),
                               grepl("\\^Intensity", datasetID) ~ paste0(dataType, "_Intensity"),
                               dataType == "RNAseq_raw_undecorated" ~ "RNAseq_raw",
                               TRUE ~ dataType) # ,
         # dataTypeMetabolomicsCombined = case_when(grepl("\\^iBAQ", datasetID) ~ paste0(dataType, "_iBAQ"),
         #                                          grepl("\\^LFQ", datasetID) ~ paste0(dataType, "_LFQ"),
         #                                          grepl("\\^Intensity", datasetID) ~ paste0(dataType, "_Intensity"),
         #                                          grepl("metabolomics", dataType, ignore.case = TRUE) ~ "metabolomics",
         #                                          TRUE ~ dataType)
         )

selectedDataTypeCol <- "dataType5"
data <- dataset %>% dplyr::group_by(!!sym(selectedDataTypeCol)) %>% filter(n() > 5) %>% ungroup

write.csv(data, "datasets_results.csv", row.names = FALSE)
write.csv(data.frame(table(data[, selectedDataTypeCol])), "numberOfDatasets.csv", row.names = FALSE)
#############################################

plotNipalsBiplot <- function(df, groups= c(), alpha = 0.5) {
  # See https://stackoverflow.com/a/49788251
  #   # install_github("vqv/ggbiplot")

  library(ggbiplot)
  
  iris_dummy<-df
  iris_dummy[is.na(iris_dummy)]<-7777 #swap out your NAs with a dummy number so prcomp will run
  pca.obj <- prcomp(iris_dummy, center=TRUE, scale.=TRUE)
  
  # scale: One of "UV" (unit variance a=a/\sigma_{a}), 
  # "vector" (vector normalisation b=b/|b|), 
  # "pareto" (sqrt UV) or "none" 
  # to indicate which scaling should be used to scale the matrix with aa variables and b samples. 
  pca.obj2 <- pcaMethods::pca(df, method="nipals", nPcs=4, center=TRUE
                              , scale = "uv"
  )
  
  pca.obj$x<-pca.obj2@scores 
  pca.obj$rotation<-pca.obj2@loadings 
  pca.obj$sdev<-pca.obj2@sDev
  pca.obj$center<-pca.obj2@center
  pca.obj$scale<-pca.obj2@scale
  
  P2 <- ggbiplot::ggbiplot(pca.obj,
                           obs.scale = 1, 
                           var.scale=1,
                           ellipse=T,
                           circle=F,
                           varname.size=3,
                           var.axes=T,
                           groups=groups, 
                           alpha=0)  +
    scale_color_discrete(name = '') +  
    coord_fixed(ratio = 0.4) +
    #theme_bw() +
    ggforce::facet_zoom(xlim = c(-3, 3), ylim = c(-3, 3))

    P2 <- P2 + theme(legend.direction ='horizontal', 
                     legend.position = 'bottom')
    
    P2$layers <- c(geom_point(aes(colour=groups), cex=1, alpha = alpha), P2$layers)
  
  P2
}

colorForCorValues <- function(data, mapping, method="spearman", use="pairwise", ...){
  
  # grab data
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # calculate correlation
  corr <- cor(x, y, method=method, use=use)
  
  # calculate colour based on correlation value
  # Here I have set a correlation of minus one to blue, 
  # zero to white, and one to red 
  # Change this to suit: possibly extend to add as an argument of `my_fn`
  colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
  
  ggally_cor(data = data, mapping = mapping, 
             method = method,
             colour = "black", ...) + 
    theme_void() +
    theme(panel.background = element_rect(fill=fill))
}

plotData <- function(df, groupColName = "", addStr = "") {

  for (group in unlist(unique(df[, groupColName]))){
    print(group)
    df.group <- df[df[, groupColName] == group,]
    pdf(paste0("ggpairs_", group, addStr, ".pdf"), width = 12, height = 12)
    print(ggpairs(df.group, columns = 2:ncol(df.group),
                  lower = list(continuous = wrap("smooth", alpha = 0.3, size = 0.2)),
                  upper = list(continuous = wrap(colorForCorValues, method = "spearman"),
                               continuous = "smooth") #,
                  # mapping=aes(# color = get(groupColName),
                  #   # fill= get(groupColName),
                  #   alpha=0.3)
    ) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8)) +
      ggtitle(group) +
      # theme_bw() +
      theme(strip.placement = "outside", text = element_text(size = 7))
    )
    
    # print(GGally::ggpairs(df.group, columns = 2:ncol(df.group), # ggplot2::aes(colour=get(groupColName)),
    #               lower = list(continuous = wrap("points", alpha = 0.3, size = 0.2),
    #                            continuous = "smooth",
    #                            combo = wrap("dot_no_facet", alpha = 0.4)),
    #               upper = list(continuous = wrap("cor", method = "spearman")),
    #               mapping=aes(# color = get(groupColName),
    #                           # fill= get(groupColName),
    #                           alpha=0.3)) +
    #         ggtitle(group) +
    #         theme_bw() +
    #         theme(strip.placement = "outside", text = element_text(size = 7))
    #       )
    dev.off()
  }
  
  pdf(file = paste0("biplotNipals", addStr, ".pdf"), width = 12, height = 10)
  print(plotNipalsBiplot(df = df %>% dplyr::select(-!!groupColName), 
                         groups= df[[groupColName]],
                         alpha = 0.3))
  dev.off()
  
  pca <- pcaMethods::pca(df %>% dplyr::select(-!!groupColName), method="nipals", nPcs=4, center=TRUE
                         , scale = "uv")
  write.csv(pca@loadings, paste0("loadings", addStr,".csv"))
  dat <- merge(pcaMethods::scores(pca), df, by=0)
  
  library(GGally)
  pdf(paste0("ggpairs_Nipals", addStr, ".pdf"), width = 12, height = 10)
  print(GGally::ggpairs(dat, columns = 2:5, ggplot2::aes(colour=get(groupColName)),
                lower = list(continuous = wrap("smooth", alpha = 0.3, size = 1), 
                             combo = wrap("dot_no_facet", alpha = 0.4)),
                upper=list(continuous = wrap("cor", method = "spearman", size = 3)),
                mapping=aes(color = get(groupColName),
                            fill= get(groupColName), 
                            alpha=0.5)) +
          ggtitle(groupColName) +
          theme_bw())
  dev.off()
  
  library(plotly)
  library(htmlwidgets)
  fig <- plot_ly(dat, x = ~PC1, y = ~PC2, z = ~PC3, 
                 color = ~as.factor(dat[[groupColName]]), 
                 type="scatter3d", mode="markers",
                 #colors = c('#636EFA','#EF553B') , 
                 marker = list(size = 2)
                 # , alpha = 0.75
  ) #%>%
  # add_markers(size = 5, marker=list(sizeref=8, sizemode="area"))
  fig <- fig %>%
    layout(
      title = "3D PCA",
      scene = list(bgcolor = "#e5ecf6"
      )
    )
  
  htmlwidgets::saveWidget(fig, paste0("plotly", addStr,".html"), selfcontained = F, libdir = "lib")
}


margPlot <- ggplot(data, aes(x = corSampleMeanNA, y = corAnalyteMeanNA, colour = get(selectedDataTypeCol))) +
  geom_point(aes(fill = get(selectedDataTypeCol)),  size = 0.8, alpha = 0.5) +
  theme_minimal() +
  theme(legend.title=element_blank())

pdf("marginPlot.pdf", width = 12, height = 10)
ggExtra::ggMarginal(margPlot, groupFill = TRUE, groupColour = TRUE)
dev.off()


# data2 <- data[, c("dataType4", 
#                   "nSamples", 
#                   "nAnalytes", 
#                   "percNATotal.log2", "corColR.log2", "corRowR.log2",  "median.wNAs.log2", 
#                    "medianSampleVariance.wNAs.log2", "medianAnalyteVariance.wNAs.log2", "skewness.wNAs.log2", "prctPC1.wNAs.log2", "prctPC2.wNAs.log2")]
# data2.woMet <- data2[!grepl("metabolomics", data2$dataType4), ]

data2 <- data[, c(selectedDataTypeCol, 
                  "nSamples", 
                  "nAnalytes", 
                  "percNATotal", "corSampleMeanNA", "corAnalyteMeanNA",  "median", 
                  "medianSampleVariance", "medianAnalyteVariance", "skewness", "prctPC1", "prctPC2")]

plotData(df = data2, groupColName = selectedDataTypeCol, addStr = "")

data2$`log2(nSamples)` <- log2(data2$nSamples)
data2$`log2(nAnalytes)` <- log2(data2$nAnalytes)
data2$`log2(medianSampleVariance)` <- log2(data2$medianSampleVariance)
data2$`log2(medianAnalyteVariance)` <- log2(data2$medianAnalyteVariance)
data2[sapply(data2, is.infinite)] <- NA

data2 <- data2[ , !names(data2) %in% 
      c("nSamples", "nAnalytes", "medianSampleVariance", "medianAnalyteVariance")]

# dput(names(data2))

data2.long <- reshape2::melt(data2)
neworder <- c("dataType5", 
              "log2(nSamples)", "log2(nAnalytes)", 
              "log2(medianSampleVariance)", "log2(medianAnalyteVariance)",
              "median", "skewness", "prctPC1", "prctPC2",
              "corSampleMeanNA", "corAnalyteMeanNA", "percNATotal")
data2.long <- dplyr::arrange(dplyr::mutate(data2.long,
                                           variable=factor(variable,levels=neworder)), variable)

# For each data characteristic: Median of the median of each data type
medianValues <- data2.long %>%
  group_by(get(selectedDataTypeCol), variable) %>%  
  summarise(medianValue = median(value, na.rm=T)) %>%
  group_by(variable) %>%  
  summarise(medianValue = median(medianValue, na.rm=T))

write.csv(medianValues, "medianValues.csv", row.names = FALSE)

ggplot.charact <- ggplot(data2.long, aes(forcats::fct_rev(get(selectedDataTypeCol)), value)) +
  geom_boxplot(aes(fill=get(selectedDataTypeCol)), alpha=0.5, outlier.size=0.5) +
  coord_flip() +
  xlab("") +
  ylab("") +
  geom_hline(aes(yintercept = medianValue), medianValues, colour = 'red') +
  facet_wrap( ~ variable, scales = "free_x", ncol=2, strip.position = "bottom") +
  ggplot2::theme_bw() +
  #theme_minimal() + 
  theme(panel.spacing.y=unit(1.5, "lines"),
        legend.title = element_blank(), axis.text.y = element_text(hjust=0, face = "bold"), 
        # legend.position="bottom",
        legend.position = "none",
        legend.justification = "left", legend.direction = "horizontal",
        strip.text = element_text(face="bold"),
        strip.placement = "outside",                      # Place facet labels outside x axis labels.
        strip.background = element_blank(),  # Make facet label background white.
        axis.title = element_blank()) +     
  guides(fill = guide_legend(reverse = TRUE))
ggsave(file=paste0("boxplots.pdf"), ggplot.charact, height=15, width=12)

########################################################################
# getNodeRule <- function(rules, terminal.id, roundDecim=NULL) {
#   rule <- rules[[toString(terminal.id)]]
#   rule.original <- rule
#   #blub <- paste0("rule=", gsub("&", ", rule=", rule))
#   rule <- gsub("&", ", ", rule)
#   
#   rule.split <- strsplit(rule, ",")[[1]]
#   rule.split <- gsub(" $|^ ", "", rule.split)
#   
#   extract <- stringr::str_match(rule.split, "([0-9]{1,2}$)")
#   rule.split <- rule.split[order(rank(strtoi(extract[,2])))]
#   
#   rule.df <- data.frame(rule=rule.split, name=rep(1, length(rule.split)), description=rep(1, length(rule.split)))
#   
#   if (!validatetools::is_infeasible(validate::validator(.data= rule.df))){
#     simplifiedRules <- validatetools::remove_redundancy(validate::validator(.data= rule.df))
#     
#     str <- ""
#     
#     simplifiedRules <- lapply(simplifiedRules$rules, function(rule){rule@expr})
#     for (rule1 in simplifiedRules){
#       str.parts <-  gsub(" ", "", strsplit(toString(rule1), ",")[[1]])
#       
#       newstr <- paste0(str.parts[[2]], str.parts[[1]], str.parts[[3]])
#       newstr <- gsub(" ", "", newstr, fixed = TRUE)
#       
#       if (grepl(paste0(">", toString(strtoi(str.parts[[3]])-1)), str)) {
#         str <- paste0(gsub(paste0(str.parts[[2]], ">", toString(strtoi(str.parts[[3]])-1)), "", str))
#         newstr <- gsub("<=", "==", newstr, fixed = TRUE)
#       } else if (grepl(paste0("<=", toString(strtoi(str.parts[[3]])+1)), str)) {
#         str <- paste0(gsub(paste0("<=", toString(strtoi(str.parts[[3]])+1)), paste0("==", toString(strtoi(str.parts[[3]])+1)), str))
#         newstr <- ""
#       } 
#       
#       if (str == ""){
#         str <- paste0(str, newstr)
#       } else {
#         str <- paste0(str, " & ", newstr) 
#       }
#     }
#     
#     str <- gsub(" &$", "", str)
#     str <- gsub("&  &", "&", str)
#   } else {
#     str <- rule.original
#   }
#   
#   if (!is.null(roundDecim))
#     str <- stringr::str_replace_all(str, "\\b\\d+\\.\\d+\\b", function(x) as.character(round(as.numeric(x), roundDecim)))
#   
#   return(str)
# }
# 
# 
# plotTerminalIdsBarplot <- function(terminal.id, tree, rules, dataTypeFreq.df){
#   
#   node.df <- data_party(tree, id = terminal.id)
#   dataTypeFreq.df.node <- data.frame(table(node.df$dataType4))
#   dataTypeFreq.df.node$totalFreq <- dataTypeFreq.df$Freq
#   dataTypeFreq.df.node$percentage <- dataTypeFreq.df.node$Freq/dataTypeFreq.df$Freq * 100
#   
#   rule <- getNodeRule(rules, terminal.id, roundDecim = 2)
#   
#   gg <- ggplot(data = dataTypeFreq.df.node, aes(x = Var1, y = percentage)) +
#     ggtitle(paste0("Node ", terminal.id, ",\n Rule: ", gsub("&", "&\n", rule))) +
#     geom_bar(stat = "identity", aes(fill = Var1), show.legend = FALSE) +
#     geom_text(aes(label = paste0(round(percentage, 2), " (", Freq, "/", totalFreq,")")), position=position_dodge(width=0.9), hjust=-0.15, angle=90) +
#     xlab("") +
#     ylab("% Total category size") +
#     ylim(0, 100) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1),
#           plot.title = element_text(size = 10))
#   return(gg)
# }


# library(partykit)
# data <- dataset %>% group_by(dataType4) %>% filter(n() > 5) %>% ungroup
# data$dataType4 <- as.factor(data$dataType4)
# 
# dataTypeFreq <- table(data$dataType4)
# weights <- data.frame(1/(dataTypeFreq/sum(dataTypeFreq)))
# colnames(weights) <- c("dataType4", "weights")
# #weights$weights <- weights$weights/sum(weights$weights) 
# data <- merge(x=data, y=weights, 
#              by="dataType4", all.x=TRUE)
# 
# tree <- ctree(dataType4 ~ corColR.log2, data = data, weights = weights, logmincriterion = -1e-200)
# # tree <- ctree(dataType4 ~ corColR.noLog2, data = data, mincriterion = 0.99)
# 
# pdf("ctree_percentages.pdf", width = 24, height = 10)
# plot(tree, margins = c(8, 0, 0, 0), tp_args = list(rot = 45, just = c("right", "top"), text = "vertical"))
# dev.off()
# 
# pvals <- unlist(nodeapply(tree, ids = nodeids(tree), function(n) info_node(n)$p.value))
# pvals <- pvals[pvals <.05]
# 
# library(ggparty)
# dataTypeFreq.df <- data.frame(dataTypeFreq)
# gg <- ggparty(tree) +
#   geom_edge() +
#   geom_edge_label() +
#   geom_node_label(line_list = list(aes(label = splitvar),
#                                    aes(label = paste0("N=", nodesize, ", p", 
#                                                       ifelse(pvals < .001, "<.001", paste0("=", round(pvals, 3)))), 
#                                        size = 10)),
#                   line_gpar = list(list(size = 13), 
#                                    list(size = 10)), 
#                   ids = "inner") +
#   geom_node_label(aes(label = paste0("Node ", id, ", N = ", nodesize)),
#                   ids = "terminal", nudge_y = 0.01, nudge_x = 0.01) # +
#   # geom_node_plot(gglist = list(
#   #   geom_bar(aes(x = "", fill = dataType4),
#   #            position = position_dodge(), color = "black"),
#   #   theme_minimal(),
#   #   theme(panel.grid.major = element_blank(),
#   #         panel.grid.minor = element_blank()),
#   #   # scale_fill_manual(values = c("gray50", "gray80"), guide = FALSE),
#   #   #scale_y_continuous(breaks = seq(0, 100, 20),
#   #   #                   limits = c(0, 100)),
#   #   xlab(""),
#   #   ylab("Frequency"),
#   #   geom_text(aes(x = "", group = dataType4,
#   #                 label = stat(count)),
#   #             stat = "count",
#   #             position = position_dodge(0.9), vjust = -0.7)),
#   #   shared_axis_labels = TRUE)
# 
# ggsave("ggparty_plot.pdf", plot=gg, width=20, height=10)
# 
# 
# terminal.ids <- nodeids(tree, from = NULL, terminal = T)
# rules <-  partykit:::.list.rules.party(tree)
#   
# plots <- lapply(terminal.ids, plotTerminalIdsBarplot, tree=tree, rules=rules, dataTypeFreq.df = dataTypeFreq.df)
# 
# library(gridExtra)
# plots.comb <- do.call("grid.arrange", c(plots, ncol=length(plots)))
# 
# 
# pdf("dataCharacteristics_ctree_ggparty.pdf", height= 15, width = 30)
# 
# #pdf(paste0("lmtree_", additionalString, "_", sett, "_scaled", scaled, "_", eval.measure, "_allDIA", allDIA, 
# #           "_selectedDia", selectedDia, "_selectedSR", selectedSR, "_alpha", alpha, "_minsize", minsize, "_ggparty_plusBoxplot.pdf"), width=plotWidth, height = plotHeight)
# #plot(g)
# gridExtra::grid.arrange(gg,                                    # bar plot spaning two columns
#                         plots.comb,                               # box plot and scatter plot
#                         ncol = 1, nrow = 2, 
#                         layout_matrix = rbind(c(1,1), 2))
# 
# dev.off()


session <- sessionInfo()
sink(paste0("visualizations_sessionInfo.txt"))
print(session)
sink()
