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
        width = 6
      ),
      column(
        sliderInput(
          inputId = "k",
          label = "Maximum Number of Clusters to Evaluate in the Metaclustering",
          value = 20,
          min = 2,
          max = 80
        ),
        sliderInput(
          inputId = "xdim",
          label = "xdim of the grid size of the self-organizing map",
          value = 10,
          min = 2,
          max = 100
        ),
        sliderInput(
          inputId = "ydim",
          label = "ydim of the grid size of the self-organizing map",
          value = 10,
          min = 2,
          max = 100
        ),
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
      title = "Choose Clustering Parameters",
      width = 12,
      collapsible = TRUE,
      collapsed = FALSE
    ),
  ),
  fluidRow(uiOutput("clusteringVisualizationSelection"))
)
