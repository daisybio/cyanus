welcomeBody <- tabItem(tabName = "welcome",
                       fluidRow(box(
                         div("We are very happy that you came to use our nice tool!"),
                         div(
                           "You can get started with uploading your Data or selecting example datasets by clicking on the button below."
                         ),
                         title = h1("Welcome to the CyTOF Pipeline"),
                         width = 12
                       )))