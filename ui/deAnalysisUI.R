
plotbox_height <- "48em"
methods_height <- "35em"

deBody <- function(){
  
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
  plotViolin <- shinydashboard::box(
    uiOutput("deVioPlots"),
    title = span("Violinplots", icon("question-circle"), id = "vioplotPopover"),
    width = 12,
    height = plotbox_height,
  )
  
  vioplotPopover <- bsPopover(
    id = "vioplotPopover",
    title = "Median expression of markers",
    content = "This plots helps to get a rough image of how strong the differences might be."
  )
  
  debody <- tabItem(
    tabName = "de",
    fluidRow(
      shinydashboard::box(
        fluidRow(
          column(
          HTML("Here, you can compute differential marker expression between your conditions. The boxplots can help to get a rough 
          overview of the different markers in order to see which markers might be differentially expressed between conditions. 
          There are various options for editing the boxplots (e.g. facet by antigen or cluster ID, color by condition, sample ID, ...). 
               There are two different analysis types:<br><br>"),
          width=12
        )),
        fluidRow(
        column(
          div(
          HTML("<body><b>Differential Marker Expression</b> (sometimes called differential states) which is a differential analysis of marker expression per condition, 
          overall or cluster-wise. This analysis can be performed using 
          <ul><li><a target=\"_blank\" rel=\"noopener noreferrer\" href = https://www.bioconductor.org/packages/release/bioc/html/limma.html>limma</a>: limma is a tool originally designed for detecting differential gene expression, here adapted for differential median marker expression (overall or per cluster). It fits a linear model and is very fast. </li>
          <li> a <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://en.wikipedia.org/wiki/Mixed_model>linear mixed effect model</a> (LMM): LMMs allow for modelling correlations between different samples in the dataset, i.e. random effects. An example for a random effect is e.g. the patient ID when more than one sample was measured for each patient. When random effects are not specified, the model is equivalent to a linear model. With this method, median marker expression (overall or cluster-wise) is analyzed as well. Very fast.</li>
          <li> a model based on <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://en.wikipedia.org/wiki/Earth_mover%27s_distance>Earth Mover's Distance</a> (CyEMD):  With CyEMD, the overall distributions of marker expressions between conditions (overall or cluster-wise) are compared, not just the median marker expression. This is why this method takes longer but can yield more sensitive results. </li>
          <li> <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04067-x>CytoGLMM</a>: fits a generalized linear mixed model predicting the conditions from the whole expression vectors. Can only handle grouped data. Very fast.</li>
          <li> <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04067-x>CytoGLM</a>: fits a bootstrapped generalized linear model and can therefore be run without grouping variable. High runtime because of the bootstrap replications.</li>
          <li> <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test>Wilcoxon tests</a>: When no grouping variable is specified, a Wilcoxon rank-sum test (= Mann-Whitney U test) is computed, otherwise, a Wilcoxon signed-rank test is performed. Very fast. </li>
          <li> <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://en.wikipedia.org/wiki/Student%27s_t-test>t-test</a>: Computes a paired or unpaired t-test. Very fast. </li>
          </ul>
          </body>")),
          width = 5
        ),
        column(
          div(
          img(src='overview_results.png', align='right', height="100%", width="100%",)
          ),
          width = 7
        )),
        fluidRow(column(
          HTML("<body>The second option is <b>Differential Cluster Abundance</b> (sometimes called differential abundance) compares the proportions of cell types across experimental 
          conditions per cluster and aims to highlight populations that are present at different ratios. The methods for this are 
          <ul><li><a target=\"_blank\" rel=\"noopener noreferrer\" href = https://bioconductor.org/packages/release/bioc/html/edgeR.html>edgeR</a>: The edgeR package is used to fit models and calculate moderated tests at the cluster level. The statistical methods were originally designed for the analysis of gene expression data such as RNA-seq counts. Here, the methods are applied to cell counts.</li>
          <li><a target=\"_blank\" rel=\"noopener noreferrer\" href = https://rdrr.io/bioc/diffcyt/man/testDA_voom.html>voom</a>: Calculates tests for differential abundance of clusters, using functions from the limma package and voom method.  Since count data are often heteroscedastic, the voom method is used to transform raw cluster cell counts and estimate observation-level weights to stabilize the mean-variance relationship. </li>
          <li>a <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://en.wikipedia.org/wiki/Generalized_linear_mixed_model#:~:text=In%20statistics%2C%20a%20generalized%20linear,models%20to%20non%2Dnormal%20data.>generalized linear mixed model</a> (GLMM): a generalized linear model is a flexible generalization of ordinary linear regression that allows for response variables that have error distribution models other than a normal distribution. A GLMM allows to include random effects like patient ID when more than one sample was measured per patient.</li></ul></body>")
        ,width=12)),
        title = h2("Differential Marker Expression Analysis"),
        width = 12
        )
    ), 
    fluidRow(
      tabBox(
        tabPanel(
          fluidRow(plotBox, boxplotPopover),
          value = "plotCounts",
          title = "Boxplots"
        ),
        tabPanel(
          fluidRow(plotViolin, vioplotPopover),
          value = "plotViolins",
          title = "Violinplots"
        ),
        id = "DEplots",
        title = "Visual Analysation",
        width = 12,
        height = plotbox_height
      )
    ),
      
    
    fluidRow(
        uiOutput("selectionBoxDE"),
    ),
    
    fluidRow(
      uiOutput("DEVisualization"),
    ),
    
    fluidRow(
      uiOutput("deTopTable")
    )
  
  )
  return(debody)
  
}