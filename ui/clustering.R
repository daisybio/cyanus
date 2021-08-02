clusteringBody <- tabItem(
  tabName = "clustering",
  fluidRow(shinydashboard::box(
    HTML(
      "Detect and define cell populations for further downstream analysis via
        <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://bioconductor.org/packages/release/bioc/html/FlowSOM.html>FlowSOM</a> clustering &
        <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://bioconductor.org/packages/release/bioc/html/ConsensusClusterPlus.html>ConsensusClusterPlus</a> meta-clustering:<br>
        First group cells into clusters with FlowSOM using Self-Organizing Map clustering and Minimal Spanning Trees, and subsequently perform meta-clustering with ConsensusClusterPlus which determines cluster count and membership by stability evidence.<br>
        A self-organizing map (SOM) is an artificial neural network that reflects the topological information of the input in lower dimensions.<br>
        Meta clustering is carried out based on a consensus hierarchical method. By re-sampling from the original data and clustering the samples, an agreement value between the perturbations is obtained, which leads to meta clusters in a hierarchical structure. The number of meta clusters is manually chosen based on the relative change in area under the cumulative distribution function curve. (Clustering Output, Section 1. Delta Area)<br>
        This clustering approach is fast and performs better than comparable techniques <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://doi.org/10.1002/cyto.a.23030>Weber et al. (2016)</a>."
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
          label = span("X-dimension of the grid size of the self-organizing map", icon("question-circle"), id = "xdimQ"),
          value = 10,
          min = 1,
          max = 100
        ),
        bsPopover(
          id = "xdimQ",
          title = "",
          content = "This specificies the x-dimension of the self-organizing map´s grid. The default 10x10 grid will yield 100 clusters."
        ),
        sliderInput(
          inputId = "ydim",
          label = span("Y-dimension of the grid size of the self-organizing map", icon("question-circle"), id = "ydimQ"),
          value = 10,
          min = 1,
          max = 100
        ),
        bsPopover(
          id = "ydimQ",
          title = "",
          content = "This specificies the y-dimension of the self-organizing map´s grid. The default 10x10 grid will yield 100 clusters."
        ),
        sliderInput(
          inputId = "k",
          label = span("Maximum Number of Clusters to Evaluate in the Meta-Clustering", icon("question-circle"), id = "kQ"),
          value = 20,
          min = 3,
          max = 80
        ),
        bsPopover(
          id = "kQ",
          title = "",
          content = "This specificies the maximum number of clusters to evaluate in the meta-clustering. The default will yield 2 through 20 meta-clusters and one additional cluster containing all cells."
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
