
# Inspired by Reisetter, "Mixture model normalization for non-targeted gas chromatography/mass
# spectrometry metabolomics data", 2017

library(shiny)
library(ggplot2)
library(gridExtra)
library(ggplotify)
library(reshape2)
library(shinydashboardPlus)

library(shinycssloaders)
options(spinner.color="#0275D8", spinner.color.background="#ffffff", spinner.size=2)
##################################################################################
# FUNCTIONS

rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
    matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

getQuantile <- function(vec = c(), prob=0.9){
  vec[is.na(vec)] <- 0
  q <- quantile(vec, probs = c(prob))
  
  propNotNA <- sum(vec>0)/length(vec)
  if(propNotNA < 2*(1-prob))  
    q <- NA
  return(q)
  }

plotBoxplots <- function(mtx, plotTitle, plotSubtitle, yValueMin, yValueMax) {

  quantiles <- reshape2::melt(data.frame(as.list(apply(mtx, 2, getQuantile)), check.names = FALSE))
  
  gg <- ggplot(reshape2::melt(mtx) , aes(x=as.factor(Var2), y=value)) + 
    geom_boxplot() + geom_jitter(color="black", size=0.1, alpha=0.2) + 
    labs(title = plotTitle, subtitle = plotSubtitle, y="Intensity", x = "") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    # stat_summary(fun = getQuantile,
    #              geom = "line",
    #              aes(group = 1),
    #              col = "red") 
    geom_line(data = quantiles,
              aes(x = variable, y = value, group = 1), size = 0.8, color = "red") +
    ylim(c(yValueMin, yValueMax))
  
  print(gg)
}


library(ggpubr)
plotCorrelation <- function(mtx.corr, plotTitle = "", plotSubtitle = "", 
                            ySampleMeanMin = NULL, ySampleMeanMax = NULL) {
    colMeans <- colMeans(mtx.corr, na.rm = TRUE)
    colNaPercentage <- colMeans(is.na(mtx.corr))*100
    # corrPlot <- ggscatter(data.frame(colNaPercentage, colMeans), x = "colNaPercentage", y = "colMeans",
    #     title = plotTitle,
    #     subtitle =  plotSubtitle,
    #     xlab = "NA Percentage",
    #     ylab= "Sample Mean",
    #     add = "reg.line",                                 # Add regression line
    #     conf.int = TRUE,                                  # Add confidence interval
    #     add.params = list(color = "blue",
    #         fill = "lightgray")) +
    #     stat_cor(method = "spearman", label.x = 3) + # Add correlation coefficient
    #     xlim(0, 100) +
    #     theme(plot.title = element_text(margin = margin(10, 0, 10, 0)))
    
    corRes <- cor.test(colMeans, colNaPercentage, method = "spearman")
    corrPlot <- ggplot(data.frame(colNaPercentage, colMeans), aes(colNaPercentage, colMeans)) +
      geom_point() +
      geom_smooth(method = "lm") +
      # ggpubr::theme_pubr() +
      theme_bw() +
      labs(x = "NA Percentage", y = "Sample Mean") +
      ggtitle(plotTitle, subtitle = plotSubtitle) +
      xlim(0, 100) 
    
    if (!is.null(ySampleMeanMin) & !is.null(ySampleMeanMax))
      corrPlot <- corrPlot + ylim(c(ySampleMeanMin, ySampleMeanMax))
    
    if (!is.na(corRes$p.value) & corRes$p.value < 0.05) corrPlot <- corrPlot + ggpubr::stat_cor(method = "spearman", label.x = 3)
    
    print(corrPlot)
}


plotProbability <- function(mtx, plotTitle = "", plotSubtitle = "", withLegend = FALSE) {
  mtx <- mtx[rowSums(is.na(mtx)) != ncol(mtx), ]
  mtx <- mtx[, colSums(is.na(mtx)) != nrow(mtx)]
  
  data.long <- reshape2::melt(mtx)
  colnames(data.long) <- c("Feature", "Sample", "Value")
  data.long$isNA <- as.integer(is.na(data.long$Value))
  
  library(speedglm)
  lmModel <- speedglm::speedlm(Value ~ Feature + Sample, data = data.long,fitted=T)
  predicted <- predict(lmModel,data=data.long)
  
  imputed <- data.long$Value
  imputed[is.na(imputed)] <- lmModel$fitted.values[is.na(imputed)]
  data.long$imputed <- imputed
  
  glmModel <- speedglm::speedglm(isNA ~ imputed + Sample, 
                                 family=binomial(link='logit'), 
                                 data = data.long, fitted = TRUE)

  data.long$prob <- c(predict(glmModel, type="response"))
  
  plot(data.long$imputed, data.long$prob)
  
  library(ggplot2)
  gg <- ggplot(data.long, aes(x = imputed, y = prob, color = Sample)) +
    geom_point(size = 1, alpha = 0.5) +
    theme_bw() +
    xlab("Intensity") +
    ylab("NA Probability") +
    ggtitle(plotTitle, subtitle = plotSubtitle) +
    ylim(c(0, 1))
  
  if (withLegend)  {
    gg <- gg + theme(legend.position="bottom", legend.title=element_blank())
  } else {
    gg <- gg + theme(legend.position = "none")
  }
  
  print(gg)
}


#' @name msb.simulateOmicsNormallyDistributed 
#' @title Simulation of a data matrix with normally distributed noise.
#' 
#' Sample effects, feature effects and batch effects can also be included.
#' 
#' @param nFeatures is the number of rows
#' @param nSamples is the number of columns
#' @param mean overall mean
#' @param sd standard deviation of noise of each data point
#' @param sdSamples standard deviation of noise which affects each sample in the same way
#' @param sdFeatures deviation of noise which affects each feature in the same way
#' @param nBatch is the number of rows
#' @param sdBatch standard deviation of noise which affects each sample in the same batch in the same way
#'
#' @examples
#' dat <- msb.simulateOmicsNormallyDistributed()
#' image(dat)
#' 
#' image(msb.simulateOmicsNormallyDistributed(sdSamples=1))
#' 
#' image(msb.simulateOmicsNormallyDistributed(sdSamples=1, sdBatch=5,nBatch=3))
#' 
#' image(msb.simulateOmicsNormallyDistributed(nFeatures=20, sdSamples=10, sdFeatures=10))

msb.simulateOmicsNormallyDistributed <- function(nFeatures=1000,nSamples=20,mean=0,sd=1,sdSamples=0,sdFeatures=0,nBatch=1,sdBatch=0) {
  r <- rnorm(nFeatures*nSamples,mean=mean,sd=sd) # observational noise
  r <- matrix(r,nrow=nFeatures) # make a matrix
  
  # are samples different?
  if(sdSamples>0){ # each column has a different offset
    sampleOffsets <- rnorm(nSamples,sd=sdSamples)
    r <- r + matrix(rep(sampleOffsets,each=nFeatures), nrow=nFeatures)
  } 
  
  # are features different?
  if(sdFeatures>0){ # each column has a different offset
    featureOffsets <- rnorm(nFeatures,sd=sdFeatures)
    r <- r + matrix(rep(featureOffsets,times=nSamples), ncol=nSamples)
  } 
  
  # Batch Effects?
  if(nBatch>1 && sdBatch>0){ # each column has a different offset
    # assignment of samples to batch (same size if nSamples is a multiple of nBatch):
    indBatch <- round(seq(0.500001,nBatch+0.49999,length.out=nSamples)) 
    batchOffsets <- rnorm(nBatch,sd=sdBatch) 
    for(i in 1:nBatch){
      r[,indBatch==i] <- r[,indBatch==i]+batchOffsets[i] 
    }
  }
  
  colnames(r) <- paste0("Sample ",1:nSamples)
  rownames(r) <- paste0("Feature ",1:nFeatures)
  
  r
}


#' @name msb.applyDetectionLimit
#' @title Add NA for small values, in analogy to a detection limit
#' 
#' The inverse logit is used for the probabilities. Location and width are chosen 
#' to match q50 and q90, i.e. the intensity level where 50% NA and 90% are there.
#' 
#' @param dat data where NA should be added
#' @param q50 the quantile in the data where prob(NA) should be 50%
#' @param q90 the quantile in the data where prob(NA) should be 90%
#' @param qSD the SD to enable sample dependent q50 and q90 in units of the data
#'            q50 and q90 are shifted by a Gaussian random number with this SD
#' @param NAval the number used as NA, default: NA, zero also makes sense
#' @param doPlot should the probablity of NA be plotted? Default: FALSE
#'
#' @examples
#' 
#' # Constant threshold:
#' image(dat <- msb.simulateOmicsNormallyDistributed(mean=6,sdSamples = 1,sdFeatures = 2))
#' dat2 <- msb.applyDetectionLimit(dat,q50=0.2,q90=0.1,doPlot=T)
#' boxplot(cbind(dat2$dat,dat2$datNA))
#' plot(colMeans(dat2$datNA,na.rm=T),colSums(is.na(dat2$datNA),na.rm=T))
#' 
#' # Sample-dependent threshold:
#' image(dat <- msb.simulateOmicsNormallyDistributed(mean=6,sdSamples = 0,sdFeatures = 2))
#' dat2 <- msb.applyDetectionLimit(dat,q50=0.2,q90=0.1,qSD=2,doPlot=T)
#' boxplot(cbind(dat2$dat,dat2$datNA))
#' print(dat2$sampleDependency)
#' plot(colMeans(dat2$datNA,na.rm=T),colSums(is.na(dat2$datNA),na.rm=T))
msb.applyDetectionLimit <- function(dat, q90 = .1, q50 = .2, qSD = 0,
                                    NAval = NA, doPlot = FALSE) {
  
  if (q50 <= q90)
    stop("q50 should be larger than q90")
  
  quants <- quantile(dat, probs = c(q90, q50)) # the quantiles in the data
  logit <- qlogis # ADDED
  scaleX <- (quants[2] - quants[1]) / (logit(0.9) - logit(0.5)) # how the x-axis of the inv logit has to be scaled
  
  # do it column-wise
  sampleDependency <- rnorm(ncol(dat), mean = 0, sd = qSD)
  isNA <- 0 * dat # zeros, same size as dat
  probNoNA <- 0 * dat # zeros, same size as dat
  for (i in 1:ncol(dat)){
    tmpQuants <- quants + sampleDependency[i]
    probNoNA[,i] <- matrix(exp((dat[,i] - tmpQuants[2])/scaleX)/(1+exp((dat[,i]-tmpQuants[2])/scaleX)),nrow=nrow(dat))
    isNA[,i] <- matrix(rbinom(nrow(dat), 1, 1 - probNoNA[,i]), nrow = nrow(dat))
  }
  
  if (doPlot)
    plot(dat, 1 - probNoNA)
  
  datNA <- dat
  datNA[isNA == 1] = NAval
  
  return(list(dat = dat, datNA = datNA, 
              sampleDependency = sampleDependency, 
              quants = quants,
              q50 = q50,
              q90 = q90))
  
}



generateMatrices <- function(input, seed, saveParameters = TRUE) {
  
  n_sample <- input$n_sample #30
  n_metab <- input$n_metab  #3000
  
  q50 <- input$q50 # 0.2
  q90 <- input$q90 # 0.1
  qSD <- input$qSD # 2
  mean <- input$mean # 12
  sd <- input$sd # 1
  sdSamples <- input$sdSamples # 1
  sdFeatures <- input$sdFeatures # 2
  
  if (saveParameters) writeLines(
    paste0("q50: ", q50, "\nq90: ", q90, "\nqSD: ", qSD, "\nmean: ", mean, 
           "\nsd: ", sd, "\nsdSamples: ", sdSamples, "\nsdFeatures: ", sdFeatures), 
    paste0("parameters_Seed", seed, ".txt"))
 
  dat <- msb.simulateOmicsNormallyDistributed(nFeatures=n_metab, 
                                              nSamples=n_sample, 
                                              mean=mean, 
                                              sd = sd,
                                              sdSamples = sdSamples, 
                                              sdFeatures = sdFeatures)
  
  print("Generating matrix")
  
  dat2 <- msb.applyDetectionLimit(dat, q50=q50, q90=q90, qSD=qSD, doPlot=FALSE)
  
  # boxplot(cbind(dat2$dat,dat2$datNA))
  
  idxNA <- sample(seq(length(dat2$dat)), sum(is.na(dat2$datNA)))
  datRandomNA <- dat2$dat
  datRandomNA[idxNA] <- NA
  
  
  list(dat = dat2$dat, datThrNA = dat2$datNA, datRandomNA = datRandomNA)
}

saveData <- function(mtxs, seed){
  saveRDS(mtxs,file=paste0("dataMatrices_Seed", seed, ".RDS"))
}

savePlots <- function(ptlist, seed){
  for (i in seq(length(ptlist))){
    ggsave(filename = paste0(names(ptlist)[i], "_Seed", seed, ".pdf"), 
           plot = ptlist[[i]],
           width = 5,
           height = 4)
  }
}

generatePlots <- function(mtxs, output, seed, savePlots = TRUE) {
  yValueMin <- min(c(c(mtxs$dat), c(mtxs$datThrNA), c(mtxs$datRandomNA)), na.rm = TRUE)
  yValueMax <- max(c(c(mtxs$dat), c(mtxs$datThrNA), c(mtxs$datRandomNA)), na.rm = TRUE)
  
  # ySampleMeanMin <- min(c(colMeans(mtxs$dat, na.rm = TRUE), 
  #                         colMeans(mtxs$datThrNA,na.rm = TRUE), 
  #                         colMeans(mtxs$datRandomNA, na.rm = TRUE)), na.rm = TRUE)
  # ySampleMeanMax <- max(c(colMeans(mtxs$dat, na.rm = TRUE), 
  #                         colMeans(mtxs$datThrNA,na.rm = TRUE), 
  #                         colMeans(mtxs$datRandomNA, na.rm = TRUE)), na.rm = TRUE)
  
  output$plotgraph = renderPlot({
    pt1 <- plotBoxplots(mtx = mtxs$dat, 
                        # plotTitle = "No NAs", 
                        plotTitle = "", 
                        plotSubtitle = "Red line indicates 90th percentile", 
                        yValueMin = yValueMin, yValueMax = yValueMax)
    
    pt2 <- plotBoxplots(mtx = mtxs$datThrNA, 
                        # plotTitle = "Detection limit-related NAs", 
                        plotTitle = "", 
                        plotSubtitle = "Red line indicates 90th percentile", 
                        yValueMin = yValueMin, yValueMax = yValueMax)
    
    pt3 <- plotBoxplots(mtx = mtxs$datRandomNA, 
                        # plotTitle = "Random NAs", 
                        plotTitle = "", 
                        plotSubtitle = "Red line indicates 90th percentile", 
                        yValueMin = yValueMin, yValueMax = yValueMax)
    
    pt4 <- plotCorrelation(mtx.corr = mtxs$dat, 
                           # plotTitle = "No NAs",
                           plotTitle = "",
                           plotSubtitle = "Points correspond to samples"# ,
                           #ySampleMeanMin = ySampleMeanMin, ySampleMeanMax = ySampleMeanMax
    )
    
    pt5 <- plotCorrelation(mtx.corr = mtxs$datThrNA, 
                           # plotTitle = "Detection limit-related NAs",
                           plotTitle = "",
                           plotSubtitle = "Points correspond to samples" # ,
                           # ySampleMeanMin = ySampleMeanMin, ySampleMeanMax = ySampleMeanMax
    )
    
    pt6 <- plotCorrelation(mtx.corr = mtxs$datRandomNA, 
                           # plotTitle = "Random NAs",
                           plotTitle = "",
                           plotSubtitle = "Points correspond to samples" # ,
                           # ySampleMeanMin = ySampleMeanMin, ySampleMeanMax = ySampleMeanMax
    )
    
    pt7 <- plotProbability(mtx = mtxs$dat, 
                           # plotTitle = "No NAs",
                           plotTitle = "",
                           plotSubtitle = "Points correspond to data points")
    
    pt8 <- plotProbability(mtx = mtxs$datThrNA, 
                           # plotTitle = "Detection limit-related NAs",
                           plotTitle = "",
                           plotSubtitle = "Points correspond to data points")
    
    pt9 <- plotProbability(mtx = mtxs$datRandomNA, 
                           # plotTitle = "Random NAs",
                           plotTitle = "",
                           plotSubtitle = "Points correspond to data points")
    
    ptlist <- list(pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8, pt9)
    names(ptlist) <- c("boxplots_NoNA", "boxplots_ThrNA", "boxplots_RandomNA",
                       "correlationPlot_NoNA", "correlationPlot_ThrNA", "correlationPlot_RandomNA",
                       "NAProbabiltyPlot_NoNA", "NAProbabiltyPlot_ThrNA", "NAProbabiltyPlot_RandomNA")
    to_delete <- !sapply(ptlist,is.null)
    ptlist <- ptlist[to_delete] 
    if (length(ptlist)==0) return(NULL)
    
    if (savePlots) savePlots(ptlist, seed)
    
    # grid.arrange(grobs=ptlist,ncol=length(ptlist)/2, nrow = 2)
    # grid.arrange(grobs=ptlist,ncol=length(ptlist)/3, nrow = 3)
    
    
    grid.arrange(arrangeGrob(pt1, pt4, pt7, top=ggpubr::text_grob("No NAs", size = 15, face = "bold")), 
                 arrangeGrob(pt2, pt5, pt8, top=ggpubr::text_grob("Detection limit-related NAs", size = 15, face = "bold")),
                 arrangeGrob(pt3, pt6, pt9, top=ggpubr::text_grob("Random NAs", size = 15, face = "bold")),
                 ncol = 3)
    # grid.arrange(grobs=ptlist,ncol=length(ptlist)/2)
    
  })
}

#####################################################################################
#####################################################################################

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  tags$head(
    # Note the wrapping of the string in HTML()
    tags$style(HTML("
      .btn {
        background-color:  #337ab7;
        color: #fff;
      }
      
      .btn:hover {
        background-color: #1966AB;
        color: #fff;
      }
      
      .btn:active {
        background-color: #011f4b !important;
        color: #fff !important;
      }"))
  ),
  
  
    # # App title ----
    titlePanel(""),
    
    # # Sidebar layout with input and output definitions ----
    sidebarLayout(
      
         # Sidebar panel for inputs ----
         sidebarPanel(
           tags$style(".well {background-color:white;}"),
           accordion(
             id = "accordion1",
             accordionItem(
               title = "Primary simulation parameters",
               status = "danger",
               collapsed = FALSE,
               sliderInput(inputId ="q50", label = "Quantile where prob(NA) should be 50% (q50):",
                           min = 0, max = 1, step = 0.01,
                           value = 0.4),
               sliderInput(inputId ="q90", label = "Quantile where prob(NA) should be 90% (q90):",
                           min = 0, max = 1, step = 0.01,
                           value = 0.2),
               sliderInput(inputId ="qSD", label = "SD to enable sample dependent q50 and q90:",
                           min = 0, max = 10, step = 0.1,
                           value = 3),
               sliderInput(inputId ="sdSamples", label = "SD of noise which affects each sample in the same way:",
                           min = 0, max = 10, step = 0.1,
                           value = 0.5)
              )
             ),
            accordion(
             id = "accordion2",
             accordionItem(
               title = "Additional simulation parameters",
               status = "warning",
               collapsed = TRUE,
               numericInput("n_sample", 
                            label="# Samples:", 
                            value = 20,  min = 1, max = 100,
                            step = 1),
               numericInput("n_metab", 
                            label="# Analytes:", 
                            value = 1000, min = 10, max = 4000,
                            step = 1),
               sliderInput(inputId ="mean", label = "Overall mean:",
                           min = 0, max = 20,
                           step = 0.5, 
                           value = 12),
               sliderInput(inputId ="sd", label = "SD of noise of each data point:",
                           min = 1, max = 10,
                           step = 0.5, 
                           value = 1),
               sliderInput(inputId ="sdFeatures", label = "SD of noise which affects each feature in the same way:",
                           min = 0, max = 10,
                           step = 0.5, 
                           value = 2)
             )
           ),
           actionButton(
             inputId = "submit",
             label = "Plot!"#,
             #style="color: #fff; background-color: #337ab7; border-color: #2e6da4"
           ),
         width = 3),
         # Main panel for displaying outputs ----
         mainPanel(
           conditionalPanel(
             condition = "input.submit > 0",
             style = "display: none;",
             withSpinner(plotOutput(outputId="plotgraph", width="1100px",height="900px"), type = 5)
           )
             
          )
))


server <- shinyServer(function(input, output) {
  observeEvent(
    eventExpr = input[["submit"]],
    handlerExpr = {
      print("PRESSED") #simulation code can go here
      seed <- sample(1:1000000000, 1)
      print(paste0("Seed: ", seed))
      mtxs <- generateMatrices(input, seed)
      # output$thresholds <- renderUI({
      #   HTML(
      #     paste0(sprintf("%s", round(mtxs[[4]], 2)), collapse = "<br>")
      #   )
      # })

      saveData(mtxs, seed)
      generatePlots(mtxs, output, seed, savePlots = TRUE)
    }
  )
})
   
shinyApp(ui, server)