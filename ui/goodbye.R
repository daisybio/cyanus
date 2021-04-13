goodbyeBody <- tabItem(tabName = "goodbye",
                       fluidRow(box(
                         div ("We would like to thank you for using our App for analysing mass cytometry data. We would also be very pleased to receive your feedback under ..."),
                         div ("If you use our CyTOF pipeline for analysing your mass cytometry data please cite ..."),
                         title = h1("Thanks for using our CyTOF Pipeline"),
                         width = 12
                       )),
                       fluidRow(
                         column(
                           12,
                           style="padding-bottom:10px;",
                           bsButton(
                             inputId = "exit",
                             label = "Quit App",
                             icon = icon("times-circle"),
                             style = "warning",
                             block = TRUE,
                           )
                         )
                       )
)
