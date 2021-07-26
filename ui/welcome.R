welcomeBody <- tabItem(tabName = "welcome",
                       fluidRow(box(
                         div ("Recently, high-dimensional time-of-flight mass cytometry (CyTOF) 
                               has emerged with the ability to identify more than 40 parameters simultaneously.
                               Traditional flow cytometry would require multiple tubes with different 
                               antibody panels to cover the same number of markers. Consequently, 
                               CyTOF experiments are becoming a powerful tool to unveil new cell subtypes, 
                               functions, and biomarkers in many fields, e.g. the discovery of disease-associated 
                               immunologic changes in cancer."),
                         div ("In order to facilitate the analysis of CyTOF data for biologists and physicians, 
                               a clear, understandable and user-friendly pipeline is needed."),
                         div("Here, we integrated the methods from the CATALYST package for preprocessing, 
                             visualization and clustering. For differential abundance detection, we included the 
                             diffcyt methods diffcyt-DA-edgeR, diffcyt-DA-voom and diffcyt-DA-GLMM."),
                         div("However, many experiments aim to detect differential states
                             within cell populations between samples in different conditions. 
                             For this, we integrated the published methods diffcyt-DS-limma, 
                             diffcyt-DS-LMM, CytoGLMM and CytoGLM. Additionally, we performed a comprehensive analysis
                             of these existing methods and novel approaches published in …. 
                             Since the Wilcoxon rank-sum test and the t-test on sample medians, as well as
                             our novel method CyEMD performed well, we made them available in this interface as well."),
                         div(img(src="cyanus_logo.png", height="150px", style="float:right; padding:20px;")),
                         title = h1("Welcome to CYANUS: CYtometry ANalysis Using Shiny"),
                         width = 12
                       )),
                       fluidRow(shinydashboard::box(withSpinner(tableOutput('licenseTable')), 
                                    title = "Dependencies",
                                    width = 12, collapsible =TRUE, collapsed = TRUE)))
