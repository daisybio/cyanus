
plotbox_height <- "45em"
methods_height <- "35em"

deBody <- function(){
  
  selectionBox <- shinydashboard::box(
      radioButtons(
        inputId = "da_ds",
        label = span("What type of analysis method do you want to perform?", icon("question-circle"), id = "da_dsQ"),
        choices = c("Differential Abundance", "Differential States"), 
        inline = T
      ),
      bsPopover(
        id = "da_dsQ",
        title = "Analysis type",
        content = HTML("Before doing this, you should have done clustering, preferrably by type. <br> <b>Differential abundance:</b> Differential analysis of cell population abundance. Compares the proportions of cell types across experimental condition and aims to highlight populations that are present at different ratios. <br> <b>Differential States:</b> Differential analysis of the median expression of the state markers in each cell population (i.e. cluster).")
      ),
      uiOutput("deMethodSelection"),
      uiOutput("clusterSelection"),
      uiOutput("modelSelection"),
      uiOutput("contrastSelection"),
      
      div(
        bsButton(
          "diffExpButton",
          "Start Analysis",
          icon = icon("tools"),
          style = "success"
        ),
        style = "float: right;"
      ),
      
    title = "Choose Method and Parameters",
    width = 6,
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
          "Here, you can compute differential marker expression between your conditions. "
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
      uiOutput("visDiffExp")
    ),
    
    fluidRow(
             uiOutput("heatmapBox"),
             uiOutput("deTopTable")),
  
  )
  return(debody)
  
}