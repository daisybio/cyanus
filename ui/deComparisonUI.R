methods_height_venn <- "51em"

vennBody <- function(){
  
  selectionBoxVenn <- shinydashboard::box(
    column(
      radioButtons(
        inputId = "da_dsVenn",
        label = span(
          "What type of analysis method do you want to perform?",
          icon("question-circle"),
          id = "da_dsVennQ"
        ),
        choices = c("Differential Marker Expression", "Differential Cluster Abundance"),
        inline = T
      ),
      bsPopover(
        id = "da_dsVennQ",
        title = "Analysis type",
        content = HTML(
          "Before doing this, you should have done clustering, preferrably by type. <br> <b>Differential Cluster Abundance:</b> Differential analysis of cell population abundance regarding the clusters. Compares the proportions of cell types across experimental condition and aims to highlight populations that are present at different ratios. <br> <b>Differential Marker Expression:</b> Differential analysis of the marker expression in each cell population (i.e. cluster or overall)."
        )
      ),
      uiOutput("DSVenn"),
      uiOutput("conditionSelectionComp"),
      uiOutput("conditionSelectionPairComp"),
      uiOutput("groupSelectionComp"),
      uiOutput("additionalTermsSelectionComp"),
      uiOutput("emdInputComp"),
      uiOutput("CytoGLM_num_bootComp"),
      uiOutput("deSubselectionComp"),
      width = 6,
      id="column1_comparison"
    ),
    column(
      uiOutput("downsamplingComp"),
      uiOutput("clusterSelectionComp"),
      uiOutput("markerToTestSelectionComp"),
      uiOutput("extraFeaturesComp"),
      uiOutput("normalizeSelectionComp"),
      uiOutput("weightSelectionComp"),
      uiOutput("fdrComp"),
      div(
        bsButton(
          "diffExpButtonVenn",
          "Start Comparison",
          icon = icon("tools"),
          style = "success"
        ),
        style = "float: right; bottom:5px"
      ),
      width = 6,
      id="column2_comparison"),
    title = "Choose Method and Parameters",
    width = 12,
    height = methods_height_venn
  )
  
  
  
  vennBody <- tabItem(
    tabName = "venn",
    fluidRow(
      shinydashboard::box(
        div(
          "Here, you can compare the results of different methods run on the same subset. Choose between Differential Cluster Abundance and Differential Marker Expression methods: "
        ),
        div(
          HTML("<ul><li>Differential Cluster Abundance methods: edgeR, voom, GLMM </li><li>Differential Marker Expression methods: limma, LMM, CyEMD, CytoGLMM, CytoGLM, Wilcoxon rank-sum test, Wilcoxon signed-rank test (for paired data), paired and unpaired t-test.</li></ul><br>For more information, please refer to the DE analysis tab!")
        ),
        title = h2("DE Method Comparison"),
        width = 12
      )
    ),
    fluidRow(
      selectionBoxVenn
    ),
    fluidRow(
      shinydashboard::box(
        div(
          uiOutput("vennTitle"),
          shinycssloaders::withSpinner(plotOutput("vennDiagrams", width = "100%", height = "1000px")),
          uiOutput("downloadVenn")
          ),
        id = "vennDiagramsBox",
        title= "Result Comparison",
        width = 12
      ),
    ),
    fluidRow(
      uiOutput("vennTable")
    )
  )
  return(vennBody)
  
}