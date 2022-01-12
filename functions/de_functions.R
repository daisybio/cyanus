#the EMD functions are in cyEMD.R

#################################### EXAMPLE #################################################### 
#source: cluster_functions, de_functions, de_comparison_functions, cyEMD
#sce <- clusterSCE(sce)
#DA: 
#results <- 
#  runDA(sce, 
#        da_methods = c("diffcyt-DA-edgeR", "diffcyt-DA-voom", "diffcyt-DA-GLMM"), 
#        clustering_to_use = "meta10", contrast_vars = "activated_baseline", 
#        design_matrix_vars = c("patient_id", "activated_baseline"), 
#        fixed_effects = "activated_baseline", 
#        random_effects = "patient_id")
#createVennDiagram(results, DS = F, 0.05)
#DS:
#results <- runDS(sce, 
#                 ds_methods = c("diffcyt-DS-limma","diffcyt-DS-LMM","CyEMD"), 
#                 clustering_to_use = "all", contrast_vars = "activated_baseline", 
#                 design_matrix_vars = c("patient_id", "activated_baseline"), fixed_effects = "activated_baseline", 
#                 random_effects = "patient_id", markers_to_test = "state", 
#                 cyEMD_condition = "activated_baseline", binSize = 0, nperm = 100)
#createVennDiagram(results, DS = T, 0.05)
###################################################################################################

effectSize <- function(sce, condition, group = NULL, k=NULL, use_assay="exprs", use_mean=FALSE) {
  
  use_assay <- match.arg(use_assay, SummarizedExperiment::assayNames(sce))
  CATALYST:::.check_cd_factor(sce, condition)
  stopifnot(nlevels(sce[[condition]]) == 2)
  
  dt <- data.table::as.data.table(t(SummarizedExperiment::assay(sce, 'exprs')))
  dt$condition <- sce[[condition]]
  
  if (!is.null(group)) {
    CATALYST:::.check_cd_factor(sce, group)
    dt$group <- sce[[group]]
  }
  message("calculating effect-sizes....")  
  if (is.null(k)) {
    return(calculateEffectSize(dt, use_mean))
  } else {
    k <- CATALYST:::.check_k(sce, k)
    c_ids <- CATALYST::cluster_ids(sce, k)
    return(data.table::rbindlist(sapply(levels(c_ids), function(c_id) calculateEffectSize(dt[c_ids==c_id], use_mean), simplify = FALSE), idcol ='cluster_id'))
  }
}

calculateEffectSize <- function(dt, use_mean=FALSE, hedges.correction = TRUE){
  res <- list()
  if ('group' %in% names(dt)) {
    if (use_mean)
      agg_fun <- mean
    else 
      agg_fun <- median
    dt_melt <- data.table::melt(dt, id.vars=c('condition', 'group'))
    dt_melt <- dt_melt[, .(value=agg_fun(value)), by=list(condition, variable, group)]
    # data.table::dcast(dt_melt, variable + group ~ condition)
    dt_melt[, comb := paste(variable, condition, sep='::')]
    comp <- split(dt_melt[, sort(unique(comb))], ceiling(seq_along(dt_melt[, sort(unique(comb))])/2))
    res$grouped <- data.table::as.data.table(rstatix::cohens_d(dt_melt, value ~ comb, comparisons = comp, paired=T, hedges.correction = hedges.correction))[, -'.y.']
    dt <- dt[, -'group']
  }
  dt_melt <- data.table::melt(dt, id.vars='condition')
  dt_melt[, comb := paste(variable, condition, sep='::')]
  comp <- split(dt_melt[, sort(unique(comb))], ceiling(seq_along(dt_melt[, sort(unique(comb))])/2))
  res$overall <- data.table::as.data.table(rstatix::cohens_d(dt_melt, value ~ comb, comparisons = comp))[, -'.y.']
  return(data.table::rbindlist(res, idcol = 'overall_group'))
}

runDA <- function(sce, parameters = NULL,
                  da_methods = c("diffcyt-DA-edgeR", "diffcyt-DA-voom", "diffcyt-DA-GLMM"),
                  clustering_to_use = "all", normalize = T, trend_edgeR = "none", blockID = NULL, 
                  contrast_vars = NULL, design_matrix_vars = NULL, fixed_effects = NULL, random_effects = NULL){
  
  #either: first call prepDiffExp by yourself and give runDA the result as list of parameters
  #or: specify contrast_vars + design_matrix_vars/fixed_effects+random_effects and we make it for you
  
  da_methods <- match.arg(da_methods, several.ok = TRUE)
  if(is.null(parameters)){
    parameters <- list()
  }
  
  results <- list()
  for(method in da_methods){
    message(paste("Calculating results for", method, "..."))
    
    if(is.null(parameters[[method]]) & method == "diffcyt-DA-GLMM"){
      message(paste("Making specific parameters for:", method ,"\nfixed effects:", toString(fixed_effects), ", random effects:", toString(random_effects), ", contrast:", toString(contrast_vars)))
      parameters[[method]] <- prepDiffExp(sce = sce, contrastVars = contrast_vars, colsFixed = fixed_effects, colsRandom = random_effects, method = method
      )
    }else if(is.null(parameters[[method]])){
      message(paste("Making specific parameters for:", method,"\ninclude in design matrix:", toString(design_matrix_vars), ", contrast:", toString(contrast_vars)))
      parameters[[method]] <- prepDiffExp(sce = sce, contrastVars = contrast_vars, colsDesign = design_matrix_vars, method = method)
    }
    
    out <- diffcyt::diffcyt(
      d_input = sce,
      design = parameters[[method]][["design"]],
      formula = parameters[[method]][["formula"]],
      contrast = parameters[[method]][["contrast"]],
      analysis_type = "DA",
      method_DA = method,
      clustering_to_use = clustering_to_use,
      normalize = normalize,
      trend_method = trend_edgeR,
      block_id = blockID
    )
    results[[method]] <- rowData(out$res)
  }
  return(results)
}

runDS <- function(sce,
                  clustering_to_use,
                  contrast_vars,
                  markers_to_test = "state",
                  ds_methods = c("diffcyt-DS-limma","diffcyt-DS-LMM","CyEMD","BEZI", "ZAGA","ZAIG",
                                 "hurdleBeta", "CytoGLMM", "CytoGLM", "logRegression", "wilcoxon_median", "kruskal_median", "t_test"),
                  design_matrix_vars = NULL, fixed_effects = NULL, random_effects = NULL,
                  parallel = FALSE, parameters = NULL, cyEMD_nperm = 500, cyEMD_binsize = NULL,
                  include_weights = TRUE, trend_limma = TRUE, blockID = NULL, 
                  cytoGLMM_num_boot = 500, time_methods = FALSE) {

  
  #for limma and LMM: 
  ##for parameters:
  ###either: first call prepDiffExp by yourself and give runDA the result as list of parameters
  ###or: specify contrast_vars + design_matrix_vars/fixed_effects+random_effects and we make it for you
  if (is.null(trend_limma)) trend_limma <- FALSE
  if (is.null(include_weights)) include_weights <- FALSE
  stopifnot(is.logical(include_weights), is.logical(trend_limma))
  CATALYST:::.check_sce(sce, TRUE)
  k <- CATALYST:::.check_k(sce, clustering_to_use)
  CATALYST:::.check_cd_factor(sce, contrast_vars)
  
  ds_methods <- match.arg(ds_methods, several.ok = TRUE)
  #extra_args <- list(...)
  
  if (is.null(parameters)) {
    parameters <- list()
  }
  
  results <- list()
  timings <- list()
  
  for (method in c("diffcyt-DS-limma", "diffcyt-DS-LMM")) {
    if (method %in% ds_methods) {
      message(paste("Calculating results for", method, "..."))
      if (is.null(parameters[[method]]) &
          method == "diffcyt-DS-LMM") {
        # message(
        #   paste(
        #     "Making specific parameters for:",
        #     method ,
        #     "fixed effects:",
        #     fixed_effects,
        #     ", random effects:",
        #     random_effects,
        #     ", contrast:",
        #     contrast_vars
        #   )
        # )
        parameters[[method]] <-
          prepDiffExp(
            sce = sce,
            contrastVars = contrast_vars,
            colsFixed = fixed_effects,
            colsRandom = random_effects,
            method = method
          )
      } else if (is.null(parameters[[method]])) {
        # message(
        #   paste(
        #     "Making specific parameters for:",
        #     method,
        #     "include in design matrix:",
        #     design_matrix_vars,
        #     ", contrast:",
        #     contrast_vars
        #   )
        # )
        parameters[[method]] <-
          prepDiffExp(
            sce = sce,
            contrastVars = contrast_vars,
            colsDesign = design_matrix_vars,
            method = method
          )
      }
      #extra args: blockID, trend_limma, markersToTest, includeWeights
      #trend <- ifelse(is.null(extra_args$trend_limma), TRUE, extra_args$trend_limma)
      
      t <- system.time(
        out <- suppressMessages(diffcyt::diffcyt(
          d_input = sce,
          design = parameters[[method]][["design"]],
          formula = parameters[[method]][["formula"]],
          contrast = parameters[[method]][["contrast"]],
          analysis_type = "DS",
          method_DS = method,
          clustering_to_use = clustering_to_use,
          block_id = blockID,
          trend = trend_limma,
          markers_to_test = getMarkersToTest(sce, method, markers_to_test),
          weights = include_weights
        ))
      )
      ds_methods <- ds_methods[ds_methods != method]
      results[[method]] <- rowData(out$res)
      timings[[method]] <- t
    }
  }
  if (length(ds_methods) == 0)
    if(time_methods){
      return(list(results = results,
                  times = timings))
    }else{
      return(results)
    }
  
  if (time_methods) {
    for (method in ds_methods) {
      t <- system.time(
        result <-
          timeMethod(
            method,
            sce,
            markers_to_test,
            clustering_to_use,
            contrast_vars,
            random_effects,
            cyEMD_binsize,
            cyEMD_nperm,
            parallel,
            include_weights,
            cytoGLMM_num_boot,
            fixed_effects
          )
      )
      results[[method]] <- result
      timings[[method]] <- t
      gc()
    }
    return(list(results = results,
                times = timings))
  } else {
    for (method in ds_methods) {
      result <-
        timeMethod(
          method,
          sce,
          markers_to_test,
          clustering_to_use,
          contrast_vars,
          random_effects,
          cyEMD_binsize,
          cyEMD_nperm,
          parallel,
          include_weights,
          cytoGLMM_num_boot,
          fixed_effects
        )
      results[[method]] <- result
      gc()
    }
    return(results)
  }
  
}

timeMethod<- function(method, sce, markers_to_test, clustering_to_use, 
                      contrast_vars, random_effects,
                      cyEMD_binsize, cyEMD_nperm, parallel, 
                      include_weights, cytoGLMM_num_boot, cytoGLMM_additional_covariates){
  # get the cluster_ids for the current meta cluster and iterate over the clusters
  if ("type" %in% markers_to_test | "state" %in% markers_to_test) {
    sce <- CATALYST::filterSCE(sce, CATALYST::marker_classes(sce) %in% markers_to_test)
  } else {
    sce <- CATALYST::filterSCE(sce, rownames(sce) %in% markers_to_test)
  }
  cluster_ids <- CATALYST::cluster_ids(sce, clustering_to_use)
  res <- sapply(levels(cluster_ids), function(curr_cluster_id) {
    cluster_results <- list()
    #subset sce to the desired cluster_id
    sce_cluster <- CATALYST::filterSCE(sce, cluster_id == curr_cluster_id, k = clustering_to_use)
    if("CyEMD" == method){
      message(sprintf("calculating CyEMD for cluster %s", curr_cluster_id))
      out <-
        cyEMD(
          sce = sce_cluster,
          condition = contrast_vars,
          binSize = cyEMD_binsize,
          nperm = cyEMD_nperm,
          parallel = parallel
        )
      cluster_results[["CyEMD"]] <- out
      
    } else if ("BEZI" == method){
      message(sprintf("calculating BEZI for cluster %s", curr_cluster_id))
      # call BEZI from gamlss
      out <- sceGAMLSS(
        sce = sce_cluster,
        method = c("BEZI"),
        condition = contrast_vars,
        random_effect = random_effects,
        parallel = parallel
      )
      cluster_results[["BEZI"]] <- out
      
    } else if ("ZAGA" == method){
      message(sprintf("calculating ZAGA for cluster %s", curr_cluster_id))
      # call ZAGA from gamlss
      out <- sceGAMLSS(
        sce = sce_cluster,
        method = c("ZAGA"),
        condition = contrast_vars,
        random_effect = random_effects,
        parallel = parallel
      )
      cluster_results[["ZAGA"]] <- out
    } else if ("ZAIG" == method){
      message(sprintf("calculating ZAIG  for cluster %s", curr_cluster_id))
      # call ZAIG from gamlss
      out <- sceGAMLSS(
        sce = sce_cluster,
        method = c("ZAIG"),
        condition = contrast_vars,
        random_effect = random_effects,
        parallel = parallel
      )
      cluster_results[["ZAIG"]] <- out
      
    } else if ("hurdleBeta" == method){
      message(sprintf("calculating betaHurdle for cluster %s", curr_cluster_id))
      # call hurdleBeta
      out <- hurdleBeta(
        sce = sce_cluster,
        condition = contrast_vars,
        random_effect = random_effects,
        weighted = include_weights,
        parallel = parallel
      )
      cluster_results[["hurdleBeta"]] <- out
      
    }else if ("CytoGLMM" == method){
      message(sprintf("calculating CytoGLMM for cluster %s", curr_cluster_id))
      # call CytoGLMM
      out <- runCytoGLMM(
        sce = sce_cluster,
        condition = contrast_vars,
        method = "cytoglmm",
        random_effect = random_effects,
        additional_covariates = cytoGLMM_additional_covariates,
        parallel = parallel
      )
      cluster_results[["CytoGLMM"]] <- out
      
    } else if ("CytoGLM" == method) {
      message(sprintf("calculating CytoGLM for cluster %s", curr_cluster_id))
      # call CytoGLMM
      out <- runCytoGLMM(
        sce = sce_cluster,
        condition = contrast_vars,
        method = "cytoglm",
        random_effect = random_effects,
        additional_covariates = cytoGLMM_additional_covariates,
        parallel = parallel,
        num_boot = cytoGLMM_num_boot
      )
      cluster_results[["CytoGLM"]] <- out
    }else if("logRegression" == method){
      message(sprintf("calculating logstic regression model for cluster %s", curr_cluster_id))
      # call logistic regression
      out <- logistic_regression(
        sce = sce,
        condition = contrast_vars,
        random_effect = random_effects,
        parallel = parallel
      )
      cluster_results[["logRegression"]] <- out
      
    }else if("wilcoxon_median" == method){
      message(sprintf("calculating wilcoxon median test for cluster %s", curr_cluster_id))
      out <- median_test(test = "Wilcoxon", 
                         sce,
                         condition = contrast_vars,
                         random_effect = random_effects,
                         features = SummarizedExperiment::rowData(sce)$marker_name,
                         assay_to_use = "exprs", 
                         parallel = parallel)
      cluster_results[["wilcoxon_median"]] <- out
      
    }else if("kruskal_median" == method){
      message(sprintf("calculating kruskal median test for cluster %s", curr_cluster_id))
      out <- median_test(test = "Kruskal_Wallis", 
                         sce,
                         condition = contrast_vars,
                         random_effect = random_effects,
                         features = SummarizedExperiment::rowData(sce)$marker_name,
                         assay_to_use = "exprs", 
                         parallel = parallel)
      cluster_results[["kruskal_median"]] <- out
    }else if("t_test" == method){
      message(sprintf("calculating t-test test for cluster %s", curr_cluster_id))
      out <- median_test(test = "T_Test", 
                         sce,
                         condition = contrast_vars,
                         random_effect = random_effects,
                         features = SummarizedExperiment::rowData(sce)$marker_name,
                         assay_to_use = "exprs", 
                         parallel = parallel)
      cluster_results[["t_test"]] <- out
    }
    return(data.table::rbindlist(cluster_results, use.names = TRUE, fill = TRUE)) #[, variable :=c('Min.', '1st Qu.', 'Median', 'Mean', '3rd Qu.', 'Max.')]
  }, simplify=FALSE)
  other_res <- data.table::rbindlist(res, idcol = 'cluster_id')
  other_res[, cluster_id:=as.factor(cluster_id)]
  other_res[, p_adj := p.adjust(p_val, "BH")] #, by="method"
  return(other_res)
}






# get appropriate vector for each method containing the markers that want to be tested
getMarkersToTest <- function(sce, ds_method, markers){
  if("type" %in% markers | "state" %in% markers){
    by <- "Marker by Class"
  }else{
    by <- "Marker by Name"
  }
  
  if (ds_method %in% c("diffcyt-DS-limma", "diffcyt-DS-LMM")){
    # vector for diffcyt methods to test on specific markers
    is_marker <- SummarizedExperiment::rowData(sce)$marker_class %in% c("type","state")
    if(by == "Marker by Class"){
      markers_to_test <- (SummarizedExperiment::rowData(sce)$marker_class %in% markers)[is_marker]
    }else{
      markers_to_test <- row.names(sce)[is_marker] %in% markers
    }
  } else if (ds_method %in% c("SigEMD", "DEsingle")){
    if(by == "Marker by Class"){
      markers_to_test <- SummarizedExperiment::rowData(sce)[SummarizedExperiment::rowData(sce)$marker_class %in% markers,]$marker_name
    }else{
      markers_to_test <- markers
    }
  }
  return (markers_to_test)
}

prepDiffExp <- function(sce, contrastVars, colsDesign, colsFixed, colsRandom=NULL,
                        method = c( "diffcyt-DA-edgeR", "diffcyt-DA-voom","diffcyt-DA-GLMM", "diffcyt-DS-limma", "diffcyt-DS-LMM") ){
  match.arg(method)
  parameters <- list()
  parameters[["ei"]] <- metadata(sce)$experiment_info
  
  if(method %in% c("diffcyt-DA-edgeR", "diffcyt-DA-voom", "diffcyt-DS-limma")){
    
    parameters[["design"]] <- diffcyt::createDesignMatrix(parameters[["ei"]], cols_design = colsDesign)
    parameters[["contrast"]] <- createCustomContrastMatrix(sce, contrastVars, parameters[["design"]], designMatrix = T)
    
  }else if(method == "diffcyt-DS-LMM"){
    parameters[["formula"]] <- diffcyt::createFormula(parameters[["ei"]], cols_fixed = colsFixed, cols_random = colsRandom)
    parameters[["contrast"]] <- createCustomContrastMatrix(sce, contrastVars, diffcyt::createDesignMatrix(parameters[["ei"]], cols_design = colsFixed), designMatrix = T)
  }else{
    parameters[["formula"]] <- diffcyt::createFormula(parameters[["ei"]], cols_fixed = colsFixed, cols_random = colsRandom)
    parameters[["contrast"]] <- createCustomContrastMatrix(sce, contrastVars, colsFixed, designMatrix = F)
  }
  return(parameters) 
}

#Then run: 
#diffcyt::diffcyt(d_input = sce,
#                 design = parameters[["design"]],
#                 formula = parameters[["formula"]],
#                 contrast = parameters[["contrast"]],
#                 analysis_type = c("DA", "DS"),
#                 method_DS = c("diffcyt-DS-limma", "diffcyt-DS-LMM"),
#                 method_DA = c("diffcyt-DA-edgeR", "diffcyt-DA-voom", "diffcyt-DA-GLMM"),
#                 clustering_to_use = "all",
#                 markers_to_test = markersToTest, 
#                 trend = trend, 
#                 normalize = T, ...)

createCustomContrastMatrix <- function(sce, contrastVars, matrix, designMatrix = T){
  if(designMatrix){
    #the entries have to correspond to the columns of the design matrix
    cnames <- colnames(matrix)
    bool <- getBools(cnames, contrastVars)
    contrast <- diffcyt::createContrast(bool)
    print(contrast)
    return(contrast)
  }else{
    #the entries have to correspond to the levels of the fixed effect terms in the model formula
    lvlList <- lapply(matrix, function(x){levels(colData(sce)[[x]])})
    names(lvlList) <- matrix
    bool <- getBools(matrix, contrastVars)
    names(bool) <- matrix
    contrast <- unlist(lapply(names(lvlList), function(x){
      return( c( 0, rep(bool[x], length(lvlList[[x]]) -1 )) ) 
    }))
    print(contrast)
    return(diffcyt::createContrast(unname(contrast)))
  }
}

getBools <- function(names, contrastVars){
  bool <- sapply(names, function(x){
    any(sapply(contrastVars, function(y){
      grepl(y,x, fixed = T )
    }))
  })
  bool <- as.numeric(bool)
  return(bool)
}

doConditionSubselection <- function(sce, subselected_categories, mappedCategories, contrastVars, exclusionList, method){
  catCount <- sapply(colnames(metadata(sce)$experiment_info)[!colnames(metadata(sce)$experiment_info) %in% c("n_cells", "sample_id")], function(x){
    return(0)
  })
  names(catCount) <- colnames(metadata(sce)$experiment_info)[!colnames(metadata(sce)$experiment_info) %in% c("n_cells", "sample_id")]
  exclude <- list()
  for(s in subselected_categories){
    category <- mappedCategories[[s]]
    catCount[[category]] <- catCount[[category]] + 1 
    if(!is.null(exclude[[category]])){
      exclude[[category]] <- c(exclude[[category]], s)
    }else{
      exclude[[category]] <- s
    }
  }
  for(x in names(catCount)){
    if(x %in% contrastVars & catCount[[x]] == 1){
      message("You want to analyse a condition you subsetted. That is not meaningful. Try again.", type = "error")
      return(NULL)
    }
  }
  #for rendering the heatmap
  exclusionList[[method]] <- exclude
  for(cat in names(exclude)){
    sub <- exclude[[cat]]
    print(sprintf("only using %s from the condition %s", sub, cat))
    sce <- CATALYST::filterSCE(sce, get(cat) %in% sub)
  }
  return(list(sce = sce, exclusionList = exclusionList))
}

runDS_old <- function(sce, condition, random_effect = NULL, 
                  ds_methods = c("diffcyt-DS-limma", "diffcyt-DS-LMM", "CyEMD","SigEMD", "DEsingle"), 
                  k = "all", parallel = TRUE, features = c("type","state"), ...) {
  ds_methods <- match.arg(ds_methods, several.ok = TRUE)
  features <- match.arg(features, several.ok = TRUE)
  
  extra_args <- list(...)
  
  result <- list()
  if ("limma" %in% ds_methods) {
    message("Using limma")
    parameters <- prepDiffExp(sce, contrastVars = c(condition), colsDesign = c(condition), method = "diffcyt-DS-limma")
    
    #blockID <- metadata(sce_dual_ab)$experiment_info[["patient_id"]]
    
    markers_to_test <- getMarkersToTest(sce,"limma",features)
    limma_res <- diffcyt::diffcyt(d_input = sce,
                                  design = parameters[["design"]],
                                  contrast = parameters[["contrast"]],
                                  analysis_type = "DS",
                                  method_DS = "diffcyt-DS-limma",
                                  clustering_to_use = k,
                                  markers_to_test = markers_to_test)
    
    result$limma <- limma_res
    
  }
  if ("diffcyt-DS-LMM" %in% ds_methods) {
    message("Using LMM")
    parameters <-  prepDiffExp(sce, contrastVars = c(condition), colsFixed = c(condition), colsRandom=c(random_effect), method = "diffcyt-DS-LMM")
    
    markers_to_test <- getMarkersToTest(sce,"LMM",features)
    LMM_res <- diffcyt::diffcyt(d_input = sce,
                                formula = parameters[["formula"]],
                                contrast = parameters[["contrast"]],
                                analysis_type = "DS",
                                method_DS = "diffcyt-DS-LMM",
                                clustering_to_use = k,
                                markers_to_test = markers_to_test)
    
    result$LMM <- LMM_res
  }
  if ("SigEMD" %in% ds_methods) {
    source("SigEMD/sigEMD.R")
    message("Using SigEMD")
    nperm <- ifelse(is.null(extra_args$nperm), 100, extra_args$nperm)
    permute_samples <- ifelse(is.null(extra_args$permute_samples), FALSE, extra_args$permute_samples)
    
    markers_to_test <- getMarkersToTest(sce,"SigEMD",features)
    
    # make subselection
    sce_SigEMD <- filterSCE(sce, marker_name %in% markers_to_test)
    
    SigEMD_res <- SigEMD(sce_SigEMD, k, condition, Hur_gene=rownames(sce_SigEMD), nperm=nperm, parallel=parallel, permute_samples=permute_samples)
    SigEMD_res$overall <- data.table::rbindlist(lapply(SigEMD_res, function(x) x$emdall))
    
    result$SigEMD <- SigEMD_res
  }
  if ("DEsingle" %in% ds_methods) {
    message("Using DEsingle")
    
    markers_to_test <- getMarkersToTest(sce,"DEsingle",features)
    
    # make subselection
    sce_DEsingle <- filterSCE(sce, marker_name %in% markers_to_test)
    
    DEsingle_res <- DEsingleSCE(sce_DEsingle, condition, k, parallel=parallel)
    
    
    result$DEsingle <- DEsingle_res
  }
  result
}


plotPbExprsMod <- function (x, k = "meta20", features = "state", assay = "exprs", 
                            fun = c("median", "mean", "sum"), facet_by = c("antigen", 
                                                                           "cluster_id"), color_by = "condition", group_by = color_by, 
                            shape_by = NULL, size_by = FALSE, geom = c("both", "points", 
                                                                       "boxes"), jitter = TRUE, ncol = NULL) 
{
  fun <- match.arg(fun)
  geom <- match.arg(geom)
  facet_by <- match.arg(facet_by)
  stopifnot(is.logical(jitter), length(jitter) == 1)
  if (!is.null(ncol)) 
    stopifnot(is.numeric(ncol), length(ncol) == 1, ncol%%1 == 
                0)
  if (facet_by == "cluster_id") {
    CATALYST:::.check_sce(x, TRUE)
    k <- CATALYST:::.check_k(x, k)
  }
  else CATALYST:::.check_sce(x)
  CATALYST:::.check_assay(x, assay)
  CATALYST:::.check_cd_factor(x, color_by)
  CATALYST:::.check_cd_factor(x, group_by)
  CATALYST:::.check_cd_factor(x, shape_by)
  shapes <- CATALYST:::.get_shapes(x, shape_by)
  if (is.null(shapes)) 
    shape_by <- NULL
  x <- x[CATALYST:::.get_features(x, features), ]
  if (any(c(facet_by, group_by) == "cluster_id")) {
    x$cluster_id <- cluster_ids(x, k)
    by <- c("cluster_id", "sample_id")
  }
  else by <- "sample_id"
  ms <- CATALYST:::.agg(x, by, fun, assay)
  df <- reshape2::melt(ms, varnames = c("antigen", by[length(by)]))
  df[[by[length(by)]]] <- factor(df[[by[length(by)]]], levels(colData(x)[[by[length(by)]]]))
  if (length(by) == 2) 
    names(df)[ncol(df)] <- "cluster_id"
  x_var <- ifelse(facet_by == "antigen", group_by, "antigen")
  if (!is.null(df$cluster_id)) 
    df$cluster_id <- factor(df$cluster_id, levels(x$cluster_id))
  i <- match(df$sample_id, x$sample_id)
  j <- setdiff(names(colData(x)), c(names(df), "cluster_id"))
  df <- cbind(df, colData(x)[i, j, drop=FALSE])
  ncs <- table(as.list(colData(x)[by]))
  ncs <- rep(c(t(ncs)), each = nrow(x))
  if (size_by) {
    size_by <- "n_cells"
    df$n_cells <- ncs
  }
  else size_by <- NULL
  df <- df[ncs > 0, , drop = FALSE]
  ggplot(df, aes_string(x_var, "value", col = color_by)) + 
    facet_wrap(facet_by, ncol = ncol, scales = "free_y") + 
    (if (geom != "boxes") 
      geom_point(alpha = 0.8, position = (if (jitter) {
        position_jitterdodge(jitter.width = 0.2, jitter.height = 0)
      }
      else "identity"), aes_string(fill = color_by, size = size_by, 
                                   shape = shape_by))) + (if (geom != "points") 
                                     geom_boxplot(alpha = 0.4, width = 0.8, fill = NA, outlier.color = NA, 
                                                  show.legend = FALSE)) + scale_shape_manual(values = shapes) + 
    scale_size_continuous(range = c(0.5, 3)) + guides(fill = FALSE, 
                                                      size = guide_legend(order = 3), shape = guide_legend(order = 2, 
                                                                                                           override.aes = list(size = 3)), col = guide_legend(order = 1, 
                                                                                                                                                              override.aes = list(alpha = 1, size = 3))) + ylab(paste(fun, 
                                                                                                                                                                                                                      ifelse(assay == "exprs", "expression", assay))) + theme_bw() + 
    theme(legend.key.height = unit(0.8, "lines"), axis.text = element_text(color = "black"), 
          strip.text = element_text(face = "bold"), strip.background = element_rect(fill = NA, 
                                                                                    color = NA), panel.grid.minor = element_blank(), 
          panel.grid.major = element_line(color = "grey", 
                                          size = 0.2)) + if (length(unique(c(x_var, color_by, 
                                                                             group_by))) == 1) {
                                            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
                                          }else {
    theme(axis.text.x = element_text(angle = 45, hjust = 1, 
                                     vjust = 1))
  }
}


