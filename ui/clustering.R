clusteringBody <- tabItem(
  tabName = "clustering",
  fluidRow(shinydashboard::box(
    HTML(
      "Detect and define cell populations for further downstream analysis via
        <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://bioconductor.org/packages/release/bioc/html/FlowSOM.html>FlowSOM</a> clustering &
        <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://bioconductor.org/packages/release/bioc/html/ConsensusClusterPlus.html>ConsensusClusterPlus</a> metaclustering:
        First group cells into clusters with FlowSOM using Self-Organizing Map clustering and Minimal Spanning Trees, and subsequently perform metaclustering with ConsensusClusterPlus which determines cluster count and membership by stability evidence."
    ),
    title = h2("Clustering"),
    width = 12
  )),
  fluidRow(
    shinydashboard::box(
      column(
        selectInput(
          "useFeaturesIn",
          label = "Features to choose from",
          choices = c("Marker by Class",
                      "Marker by Name")
        ),
        uiOutput("featuresOut"),
        uiOutput("assayTypeOut"),
        shinyjs::hidden(div(
          id = "invalidClusteringFeaturesWarning",
          shinydashboard::infoBox(
            title = "Invalid clustering features",
            subtitle = "No features in the current selection.",
            value = "Please select more features.",
            width = 12,
            color = "orange",
            fill = TRUE,
            icon = shiny::icon("exclamation-triangle")
          ))),
        width = 6
      ),
      column(
        sliderInput(
          inputId = "xdim",
          label = "X-dimension of the grid size of the self-organizing map",
          value = 10,
          min = 1,
          max = 100
        ),
        sliderInput(
          inputId = "ydim",
          label = "Y-dimension of the grid size of the self-organizing map",
          value = 10,
          min = 1,
          max = 100
        ),
        sliderInput(
          inputId = "k",
          label = "Maximum Number of Clusters to Evaluate in the Metaclustering",
          value = 20,
          min = 3,
          max = 80
        ),
        shinyjs::hidden(div(
          id = "invalidClusteringParamsWarning",
          shinydashboard::infoBox(
            title = "Invalid clustering dimensions",
            subtitle = "The following has to be TRUE: floor(0.8*(xdim*ydim)) >= maxK.",
            value = "Need higher X/Y-dimension or lower maximum number of clusters.",
            width = 12,
            color = "orange",
            fill = TRUE,
            icon = shiny::icon("exclamation-triangle")
          ))),
        width = 6
      ),
      div(
        bsButton(
          "startClustering",
          "Start Clustering",
          icon = icon("border-none"),
          style = "success"
        ),
        style = "float: right;"
      ),
      id = "clusteringParametersBox",
      title = "Choose Clustering Parameters",
      width = 12,
      collapsible = TRUE,
      collapsed = FALSE
    ),
  ),
  fluidRow(uiOutput("clusteringVisualizationSelection"))
)
