welcomeBody <- tabItem(tabName = "welcome",
                       fluidRow(box(
                         box(
                           width = 8,
                           div(HTML(
                             'Recently, high-dimensional time-of-flight mass cytometry (CyTOF)
                             has emerged with the ability to identify more than 40 parameters simultaneously.
                             Traditional flow cytometry would require multiple tubes with different
                             antibody panels to cover the same number of markers. Consequently,
                             CyTOF experiments are becoming a powerful tool to unveil new cell subtypes,
                             functions, and biomarkers in many fields, e.g. the discovery of disease-associated
                             immunologic changes in cancer. <br><br>
                             In order to facilitate the analysis of CyTOF data for biologists and physicians,
                             a clear, understandable and user-friendly pipeline is needed.<br><br>
                             Here, we integrated the methods from the CATALYST package for preprocessing,
                             visualization and clustering. For differential abundance detection, we included the
                             diffcyt methods diffcyt-DA-edgeR, diffcyt-DA-voom and diffcyt-DA-GLMM.<br><br>
                             However, many experiments aim to detect differential states
                             within cell populations between samples in different conditions.
                             For this, we integrated the published methods diffcyt-DS-limma,
                             diffcyt-DS-LMM, CytoGLMM and CytoGLM. Additionally, we performed a comprehensive analysis
                             of these existing methods and novel approaches published in <a href="https://doi.org/10.1093/bib/bbab471" target="_blank">Briefings in Bioinformatics</a>. 
                             Since the Wilcoxon rank-sum test and the t-test on sample medians, as well as our novel method CyEMD performed well, 
                             we made them available in this interface.'),
                             style = "font-size: large;"
                           )
                         ),
                         column(4, div(
                           img(src = "cyanus_logo.png", height = "250px", width="auto"), 
                           style = "vertical-align: middle; text-align:center;"
                         )),
                         title = h1("Welcome to CYANUS: CYtometry ANalysis Using Shiny"),
                         width = 12
                       )),
                       fluidRow(
                         shinydashboard::box(
                           width = 12,
                           title = 'Why is there a new logo?',
                           div(
                             HTML(
                               'Our old logo resembled the flower <i>Centaurea cyanus</i>, commonly known as cornflower. 
                               Its strong blue color is utilized as pigment for food and drinks and its leaves can be consumed as tea. 
                               In the past, it has also been used in herbal medicine. <br><br>
                               Tt has recently come to our attention that the flower also has a strong symbolic value. 
                               It is, e.g., the national flower of Estonia. In France, it is a symbol for veterans, victims of war, widows, and orphans. 
                               In Belgium, Sweden, and Finland, it has been used as symbol for liberal parties.
                               In Germany and Austria, however, it has been used by the Nazis and other nationalist, antisemitic movements. <br><br>
                               We explicitly distance ourselves from any connection to this symbolic use and have changed our logo and design. '
                             ),
                             style = "font-size: large;"
                           ),
                           collapsible = TRUE,
                           collapsed = TRUE
                           )
                        ),
                       fluidRow(
                         shinydashboard::box(
                           withSpinner(tableOutput('licenseTable')),
                           title = "Dependencies",
                           width = 12,
                           collapsible = TRUE,
                           collapsed = TRUE
                         )
                       ))
