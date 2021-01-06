

deBody <- function(){
  selectionBox <- shinydashboard::box(
      radioButtons(
        inputId = "da_ds",
        label = span("What type of analysis method do you want?", icon("question-circle"), id = "da_dsQ"),
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
    title = "Select your parameters",
    width = 4
  )
  
  plotBox <- shinydashboard::box(
    uiOutput("deBoxPlots"),
    title = "Boxplots",
    width = 8
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
      selectionBox, 
      plotBox
    )
  )
  
  return(debody)
  
}