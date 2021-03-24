plotbox_height <- "48em"
methods_height <- "35em"

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
        choices = c("Differential Abundance", "Differential States"),
        inline = T
      ),
      bsPopover(
        id = "da_dsVennQ",
        title = "Analysis type",
        content = HTML(
          "Before doing this, you should have done clustering, preferrably by type. <br> <b>Differential abundance:</b> Differential analysis of cell population abundance regarding the clusters. Compares the proportions of cell types across experimental condition and aims to highlight populations that are present at different ratios. <br> <b>Differential States:</b> Differential analysis of the marker expression in each cell population (i.e. cluster or overall)."
        )
      ),
      uiOutput("modelSelectionVenn"),
      uiOutput("contrastSelectionVenn"),
      uiOutput("deSubselectionVenn"),
      width = 6
    ),
    column(
      uiOutput("clusterSelectionVenn"),
      uiOutput("markerToTestSelectionVenn"),
      uiOutput("extraFeaturesVenn"),
      uiOutput("normalizeSelectionVenn"),
      width = 6),
    div(
      bsButton(
        "diffExpButtonVenn",
        "Start Comparison",
        icon = icon("tools"),
        style = "success"
      ),
      style = "float: right; bottom:5px"
    ),
    title = "Choose Method and Parameters",
    width = 12,
    height = methods_height
  )
  
  
  
  vennBody <- tabItem(
    tabName = "venn",
    fluidRow(
      shinydashboard::box(
        div(
          "Here, you can compare the results of different methods run on the same subset. Choose between Differential Abundance and Differential States methods: "
        ),
        div(
          HTML("<ul><li>Differential Abundance methods: edgeR, voom, GLMM </li><li>Differential States methods: limma, LMM</li></ul><br>For more information, please refer to the DE analysis tab!")
        ),
        title = h2("DE Method Comparison"),
        width = 12
      )
    ),
    fluidRow(
      selectionBoxVenn
    ),
    fluidRow(
      shinycssloaders::withSpinner(plotOutput("vennDiagrams", width = "100%", height = "550px"))
    ),
    fluidRow(
      uiOutput("downloadVenn")
    )
  )
  return(vennBody)
  
}