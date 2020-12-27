server <- function(input, output, session) {
  
  # make reactiveValues and server-wide variables
  tab_ids <- c("welcome", "start", "preprocessing", "visualization", "clustering", "de")
  reactiveVals <- reactiveValues()
  reactiveVals$current_tab <- 1
  
  # read all server files
  sapply(list.files("server", full.names = TRUE), source, environment())
}


