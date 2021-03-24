markerSelectionBox <- shinydashboard::box(
  uiOutput("useFeaturesInVis"),
  uiOutput("markersVis"),
  id = "visTabs",
  title = "Run your dimensionality reduction", 
  width = 6
)

visbody <- function(){
  
  plotBox <- shinydashboard::box(
    fluidRow(
      column(
        2,
        div(
          uiOutput("methodsVis"),
          uiOutput("radioButtonsColorVis"),
          uiOutput("selectColorBy"),
          uiOutput("assayVis"),
          uiOutput("radioButtonsScale"),
          bsPopover(
            id = "scaleVisQ",
            title = "For coloring by expression",
            content = "Should the counts / the expression data be scaled between 0 and 1 using lower (1%) and upper (99%) expression quantiles?"
          ),
          uiOutput("facet_by"),
          uiOutput("plotDimensions"),
          div(
            bsButton(
              "startDimRed",
              "Start Visualization",
              icon = icon("palette"),
              style = "success", 
              disabled = TRUE
            ),
            style = "float: right;",  
            id = "divStartDimRed"
          ),
        style = "position: relative; height: 500px;"
        ),
      ),
      column(
        9,
        shinycssloaders::withSpinner(plotOutput("visPlot", width = "90%")), 
        ),
      column(
        1,
        uiOutput("plotInfo")
        ), 
      div(
        downloadButton("downloadPlot", "Download Plot"),
        style = "position: absolute; bottom: 10px;right:10px"
      )
    ),
    id = "visPlotBox",
    title = "Dimensionality Reduction",
    width = 12)
  
    
  visbody <- tabItem(
    tabName = "visualization",
    fluidRow(
      shinydashboard::box(
      div("Here, you can apply a dimensionality reduction method to your data and visualize its results. The following methods are currently supported:"
      ),
      div(HTML("<ul><li> <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://arxiv.org/abs/1802.03426> UMAP </a>: Uniform Manifold Approximation and Projection. Strong mathematical background, results don't vary very much (in comparison to e.g. t-SNE), State-of-the-Art nonlinear visualization method. </li> <li> <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://en.wikipedia.org/wiki/T-distributed_stochastic_neighbor_embedding> T-SNE</a>: T-distributed stochastic neighbour embedding. Results vary more because of random initialization. Still very popular. Nonlinear method.</li><li> <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://en.wikipedia.org/wiki/Principal_component_analysis>PCA</a>: Principal Component Analysis. The oldest dimensionality reduction technique (1901) but deterministic and very fast. Still a standard method though outperformed by other dimensionality reduction methods. Linear method.</li><li> <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://en.wikipedia.org/wiki/Multidimensional_scaling> MDS </a>: Multidimensional scaling. Uses pairwise distances. Linear method. </li><li> <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://en.wikipedia.org/wiki/Diffusion_map>Diffusion map</a>: Robust to noise perturbation, computationally inexpensive. Nonlinear method.</li><li> <a target=\"_blank\" rel=\"noopener noreferrer\" href = https://en.wikipedia.org/wiki/Isomap>Isomap</a>: Extends MDS by incorporating geodesic distances imposed by a weighted graph. Nonlinear method. </li></ul>")),
      div("Run the dimensionality reduction with the markers that capture best how dissimilar your cells are. You can color and facet by any marker expression, condition, sample and clustering you wish afterwards. "),
      title = h2("Data Visualization"),
      width = 12)
      ),
    fluidRow(
        uiOutput("runDRparBox"),
        width = 12
    ),
    fluidRow(
      plotBox
    )
  )
  return(visbody)
}




