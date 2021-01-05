library(diffcyt)

observeEvent(input$deBoxFacet, {
  if(input$deBoxFacet == "cluster_id"){
    shinyjs::show("deBoxK")
  }else{
    shinyjs::hide("deBoxK")
  }
})

output$deMethodSelection <- renderUI({
   methodsDA <- c("edgeR" = "diffcyt-DA-edgeR", "Voom" = "diffcyt-DA-voom", "GLMM" = "diffcyt-DA-GLMM")
   methodsDS <- c("limma" = "diffcyt-DS-limma","LMM" = "diffcyt-DS-LMM")
   if(input$da_ds == "Differential Abundance"){
     choices <- methodsDA
   }else{
     choices <- methodsDS
   }
   div(
     selectizeInput(
      inputId = "chosenDAMethod",
      label = span("Available Methods", icon("question-circle"), id = "deMethodsQ"),
      choices = choices,
      multiple = F
    ),
    bsPopover(
      id = "deMethodsQ",
      title = "Available Methods",
      content = "Depending on what you want to analyse, there are different methods available. Please see their documentation for further explanation"
    )
   )
})

output$clusterSelection <- renderUI({
  selectizeInput(
    inputId = "deCluster",
    label = "Choose the cluster populations you want to compare",
    choices = names(cluster_codes(reactiveVals$sce)), 
    multiple = F
  )
})

output$deBoxPlots <- renderUI({
  uiOutput("deExprsCluster")
  #shinydashboard::tabBox(
  #  tabPanel(
  #    plotOutput("deExprsBoxPlot"), 
  #    title = "Overall marker expression",
  #    value = "deExprsTab", 
  #    with = 12
  #    ),
  #  tabPanel(
  #    uiOutput("deExprsCluster"),
  #    title = "Cluster marker expression",
  #    value = "deExprsClusterTab", 
  #    with = 12
  #  ),
  #  width = 12
  #)
  
})

#output$deExprsBoxPlot <- renderPlot(
#  plotPbExprs(reactiveVals$sce, features = "state", shape_by = "patient_id")
#)

output$deExprsCluster <- renderUI({
  factors <- names(colData(reactiveVals$sce))[!names(colData(reactiveVals$sce)) %in% c("patient_id", "sample_id")]
  fluidRow(column(4,
    dropdownButton(
      tags$h3("Plot Options"),
      selectizeInput("deBoxFacet",
                     "Facet by:",
                     c("antigen", "cluster_id"), 
                     multiple = F, 
                     selected = "antigen"),
      hidden(selectizeInput(
        "deBoxK",
        "Clusters",
        names(cluster_codes(reactiveVals$sce)),
        multiple = F,
        selected = "meta9"
      )),
      selectizeInput("deBoxFeatures",
                     "Markers:",
                     unique(SummarizedExperiment::rowData(reactiveVals$sce)$marker_class), 
                     multiple = F),
      selectizeInput(
        "deBoxColor",
        "Color by:",
        c(names(colData(reactiveVals$sce)), names(cluster_codes(reactiveVals$sce))),
        selected = factors[1],
        multiple = F
      ),
      selectizeInput(
        "deBoxShape",
        "Shape By: ",
        c(names(colData(reactiveVals$sce))),
        multiple = F
      ),
      circle = TRUE,
      status = "info",
      icon = icon("gear"),
      width = "100%",
      tooltip = tooltipOptions(title = "Click to see plot options")
  )),
  column(8, shinycssloaders::withSpinner(
    plotOutput("clusterDEPlot", width = "100%", height = "500px")
  )))
})

output$clusterDEPlot <- renderPlot(
  plotPbExprs(reactiveVals$sce, 
              k = input$deBoxK, 
              features = input$deBoxFeatures, 
              color_by = input$deBoxColor, 
              facet_by = input$deBoxFacet, 
              shape_by = input$deBoxShape)
)





