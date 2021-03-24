
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
          "Before doing this, you should have done clustering, preferrably by type. <br> <b>Differential abundance:</b> Differential analysis of cell population abundance regarding the clusters. Compares the proportions of cell types across experimental condition and aims to highlight populations that are present at different ratios. <br> <b>Differential States:</b> Differential analysis of the marker expression in each cell population (i.e. cluster or overall)."
        )
      ),
      uiOutput("deMethodSelection"),
      uiOutput("modelSelection"),
      uiOutput("contrastSelection"),
      uiOutput("emdInput"),
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
          HTML("<body>There are two different analysis types: Differential abundance compares the proportions of cell types across experimental 
          conditions per cluster and aims to highlight populations that are present at different ratios. The methods for this are 
          <ul><li><a target=\"_blank\" rel=\"noopener noreferrer\" href = https://bioconductor.org/packages/release/bioc/html/edgeR.html>edgeR</a>: The edgeR package is used to fit models and calculate moderated tests at the cluster level. The statistical methods were originally designed for the analysis of gene expression data such as RNA-seq counts. Here, the methods are applied to cell counts.</li>
          <li><a target=\"_blank\" rel=\"noopener noreferrer\" href = https://rdrr.io/bioc/diffcyt/man/testDA_voom.html>voom</a>: Calculates tests for differential abundance of clusters, using functions from the limma package and voom method.  Since count data are often heteroscedastic, the voom method is used to transform raw cluster cell counts and estimate observation-level weights to stabilize the mean-variance relationship. </li>
          <li>a <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://en.wikipedia.org/wiki/Generalized_linear_mixed_model#:~:text=In%20statistics%2C%20a%20generalized%20linear,models%20to%20non%2Dnormal%20data.>generalized linear mixed model</a> (GLMM): a generalized linear model is a flexible generalization of ordinary linear regression that allows for response variables that have error distribution models other than a normal distribution. A GLMM allows to include random effects like patient ID when more than one sample was measured per patient.</li></ul>
          The second option is differential states (differential analysis of marker expression per condition, 
          overall or cluster-wise). This analysis can be performed using 
          <ul><li><a target=\"_blank\" rel=\"noopener noreferrer\" href = https://www.bioconductor.org/packages/release/bioc/html/limma.html>limma</a>: limma is a tool originally designed for detecting differential gene expression, here adapted for differential median marker expression (overall or per cluster). It fits a linear model. </li>
          <li> a <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://en.wikipedia.org/wiki/Mixed_model>linear mixed effect model</a> (LMM): LMMs allow for modelling correlations between different samples in the dataset, i.e. random effects. An example for a random effect is e.g. the patient ID when more than one sample was measured for each patient. When random effects are not specified, the model is equivalent to a linear model. With this method, median marker expression (overall or cluster-wise) is analyzed as well. </li>
          <li> a model based on <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://en.wikipedia.org/wiki/Earth_mover%27s_distance>Earth Mover's Distance</a> (EMD):  With the EMD method, the overall distributions of marker expressions between conditions (overall or cluster-wise) are compared, not just the median marker expression. This is why this method takes longer but can yield more sensitive results. </li>
          </ul></body>")
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