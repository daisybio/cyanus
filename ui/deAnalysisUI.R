
plotbox_height <- "48em"
methods_height <- "35em"

deBody <- function(){
  
  selectionBox <- shinydashboard::box(
    column(
      radioButtons(
        inputId = "da_ds",
        label = span(
          "What type of analysis method do you want to perform?",
          icon("question-circle"),
          id = "da_dsQ"
        ),
        choices = c("Differential Abundance", "Differential States"),
        inline = T
      ),
      bsPopover(
        id = "da_dsQ",
        title = "Analysis type",
        content = HTML(
          "Before doing this, you should have done clustering, preferrably by type. <br> <b>Differential abundance:</b> Differential analysis of cell population abundance. Compares the proportions of cell types across experimental condition and aims to highlight populations that are present at different ratios. <br> <b>Differential States:</b> Differential analysis of the median expression of the state markers in each cell population (i.e. cluster)."
        )
      ),
      uiOutput("deMethodSelection"),
      uiOutput("modelSelection"),
      uiOutput("contrastSelection"),
      uiOutput("deSubselection"),
      width = 6
    ),
    column(
      uiOutput("clusterSelection"),
      uiOutput("markerToTestSelection"),
      uiOutput("extraFeatures"),
      uiOutput("normalizeSelection"),
      width = 6),
    div(
      bsButton(
        "diffExpButton",
        "Start Analysis",
        icon = icon("tools"),
        style = "success"
      ),
      style = "float: right; bottom:5px"
    ),
    title = "Choose Method and Parameters",
    width = 12,
    height = methods_height
  )
  
  plotBox <- shinydashboard::box(
    uiOutput("deBoxPlots"),
    title = span("Boxplots", icon("question-circle"), id = "boxplotPopover"),
    width = 12,
    height = plotbox_height,
  )
  
  boxplotPopover <- bsPopover(
    id = "boxplotPopover",
    title = "Median expression of markers",
    content = "This plots helps to get a rough image of how strong the differences might be."
  )
  
  debody <- tabItem(
    tabName = "de",
    fluidRow(
      shinydashboard::box(
        div(
          "Here, you can compute differential marker expression between your conditions. The boxplots can help to get a rough 
          overview of the different markers in order to see which markers might be differentially expressed between conditions. 
          There are various options for editing the boxplots (e.g. facet by antigen or cluster ID, color by condition, sample ID, ...)"
        ),
        div(
          "There are two different analysis types: Differential abundance compares the proportions of cell types across experimental 
          conditions per cluster and aims to highlight populations that are present at different ratios. The methods for this are"
        ),
        div(
          HTML("<body><a target=\"_blank\" rel=\"noopener noreferrer\" href = https://bioconductor.org/packages/release/bioc/html/edgeR.html>edgeR</a>, <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://rdrr.io/bioc/diffcyt/man/testDA_voom.html>voom</a> and the generalized linear mixed model (GLMM). <br>
          The second option is differential states (differential analysis on the median marker expression of the markers per condition, 
          overall or cluster-wise). This analysis can be performed using <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://www.bioconductor.org/packages/release/bioc/html/limma.html>limma</a> 
               or a linear mixed effect model. </body>")
        ),
        title = h2("Differential Marker Expression Analysis"),
        width = 12
        )
    ), 
    fluidRow(
      plotBox,
      boxplotPopover,
    ),
    
    fluidRow(
      selectionBox,
    ),
    
    fluidRow(
      uiOutput("DEVisualization"),
    ),
    
    fluidRow(
      uiOutput("deTopTable")
    ),
  
  )
  return(debody)
  
}