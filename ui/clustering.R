clusteringMethods <-
  c(
    "FlowSOM clustering & ConsensusClusterPlus metaclustering" = "flowSOM",
   # "ClusterX: Fast clustering by automatic search and find of density peaks" = "clusterX",
    "RphenoGraph clustering" = "rphenoGraph"
  )


clusteringBody <- tabItem(
  tabName = "clustering",
  fluidRow(shinydashboard::box(
    div("Here you can cluster your data."),
    title = h2("Clustering"),
    width = 12
  )),
  fluidRow(div(
    id = "clusteringInputs",
    shinydashboard::box(
      selectizeInput("clusteringMethod",
                     "Clustering Methods",
                     clusteringMethods),
      div(
        bsButton(
          "startClustering",
          "Start Clustering",
          icon = icon("border-none"),
          style = "success"
        ),
        style = "float: right;"
      ),
      uiOutput("parameters"),
      title = "Choose Clustering Method",
      width = 12
    )
  )),
  fluidRow(uiOutput("clusteringVisualizationSelection"))
)
