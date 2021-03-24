welcomeBody <- tabItem(tabName = "welcome",
                       fluidRow(box(
                         div ("CyTOF is an uprising method for discovering biomarkers of the immune system. 
                              Its advantage over Flow Cytometry is that it labels the antibodies with metal 
                              isotopes instead of fluorophores, allowing it to analyse more markers in a single 
                              run while needing fewer cells."),
                         div ("Consequently, CyTOF experiments are becoming a powerful tool to find 
                              immune marker expression differences between different conditions. In order 
                              to facilitate the analysis of CyTOF data for biologists and physicians, 
                              a clear, understandable and user-friendly pipeline is needed."),
                         title = h1("Welcome to the CyTOF Pipeline"),
                         width = 12
                       )))