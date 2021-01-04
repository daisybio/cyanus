clusteringMethods <-
  c(
    "FlowSOM clustering & ConsensusClusterPlus metaclustering" = "flowSOM",
   # "ClusterX: Fast clustering by automatic search and find of density peaks" = "clusterX",
    "RphenoGraph clustering" = "rphenoGraph"
  )


clusteringBody <- tabItem(
  tabName = "clustering",
  fluidRow(shinydashboard::box(
    HTML("Choose from two clustering algorithms to detect and define cell populations for further downstream analysis:
        <ul><li> <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://bioconductor.org/packages/release/bioc/html/FlowSOM.html>FlowSOM</a> clustering &
        <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://bioconductor.org/packages/release/bioc/html/ConsensusClusterPlus.html> ConsensusClusterPlus </a> metaclustering:
        First group cells into clusters with FlowSOM using Self-Organizing Map clustering and Minimal Spanning Trees, <br>and subsequently perform metaclustering with ConsensusClusterPlus which determines cluster count and membership by stability evidence. <it>recommended</it></li>
        <li><a target=\"_blank\" rel=\"noopener noreferrer\" href = http://www.cell.com/cell/abstract/S0092-8674(15)00637-6>RphenoGraph</a>:
        Clustering method designed for high-dimensional single-cell data analysis. It works by creating a graph network representing phenotypic similarities between cells by calculating the Jaccard coefficient between nearest-neighbor sets, and then identifying communities using the Louvain method.</li></ul>"),
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
