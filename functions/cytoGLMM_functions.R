simulateSCE <- function (n_samples = 16,
                         n_cells = 1000,
                         n_markers = 10,
                         rho_b = 0.1,
                         rho_u = 0.1,
                         sigma_b = 1,
                         sigma_u = 1,
                         beta_treatment = log(24.7),
                         beta_control = log(22.9),
                         n_true = 3)
{
  library(data.table)
  n_donors <- n_samples / 2
  donor <- rep(rep(seq_len(n_donors), each = n_cells), 2)
  samples <- rep(seq_len(n_samples), each = n_cells)
  condition <- c(rep("treatment", length(donor) / 2), rep("control",
                                                          length(donor) /
                                                            2))
  md <- data.frame(
    sample_id = seq_len(n_samples),
    patient_id = rep(sprintf("p%02d", seq_len(n_donors)), 2),
    condition = factor(c(
      rep("treatment", n_donors), rep("control", n_donors)
    ), levels = c("control", "treatment"))
  )
  md$file_name <- sprintf("%s.fcs", md$sample_id)
  
  protein_names <- sprintf("m%02d", seq_len(n_markers))
  rcov <- function(rho, sigma) {
    corr <- rho ^ toeplitz(0:(n_markers - 1))
    sigma_vec <- rep(sigma, n_markers)
    diag(sigma_vec) %*% corr %*% diag(sigma_vec)
  }
  rcov_block <- function(rho, sigma) {
    corr <- diag(1, nrow = n_markers)
    corr_act <- rho ^ toeplitz(0:(n_true - 1))
    corr_notact <- rho ^ toeplitz(0:(n_markers - n_true -
                                       1))
    corr[seq_len(n_true), seq_len(n_true)] <- corr_act
    corr[(n_true + 1):n_markers, (n_true + 1):n_markers] <-
      corr_notact
    sigma_vec <- rep(sigma, n_markers)
    diag(sigma_vec) %*% corr %*% diag(sigma_vec)
  }
  Sigma_b <- rcov(rho_b, sigma_b)
  Sigma_u <- rcov(rho_u, sigma_u)
  b <-
    MASS::mvrnorm(n = n_samples * n_cells, mu = rep(0, n_markers), Sigma_b)
  u <- MASS::mvrnorm(n = n_donors, mu = rep(0, n_markers), Sigma_u)
  u <- u[donor,]
  beta <- matrix(beta_control, nrow = nrow(b), ncol = n_markers)
  beta[, seq_len(n_true)] <- ifelse(condition == "treatment",
                                    beta_treatment, beta_control)
  log_lambda <- beta + b + u
  lambda <- exp(log_lambda)
  y <- rpois(length(lambda), lambda)
  dim(y) <- dim(lambda)
  colnames(y) <- protein_names
  panel <- data.frame(
    fcs_colname = protein_names,
    antigen = protein_names,
    marker_class = "state",
    differential = c(rep(TRUE, n_true), rep(FALSE, n_markers - n_true))
  )
  
  fs <- flowCore::flowSet(mapply(
    function(file_name, s_id)
      flowCore::flowFrame(y[s_id == samples, ], description = list(GUID = file_name)),
    md$file_name,
    md$sample_id
  ))
  CATALYST::prepData(
    x = fs,
    panel = panel,
    md = md,
    panel_cols = list(differential = "differential")
  )
}


runCytoGLMM <-
  function(sce,
           condition,
           group,
           method = "cytoglmm",
           protein_names = SingleCellExperiment::rowData(sce)$marker_name,
           assay_to_use = "exprs",
           num_cores = 1,
           num_boot = 100) {
    match.arg(method, c("cytoglmm", "cytoglm"))
    match.arg(assay_to_use, names(SummarizedExperiment::assays(sce)))
    data <-
      as.data.frame(t(SummarizedExperiment::assay(sce, assay_to_use)))
    marker_names <-
      match.arg(protein_names, rowData(sce)$marker_name, several.ok = TRUE)
    marker_names <- sapply(marker_names, function(marker) {
      gsub("[^[:alnum:]]", "", marker)
    })
    colnames(data) <- marker_names
    data$donor <- SingleCellExperiment::colData(sce)[[group]]
    data$condition <- SingleCellExperiment::colData(sce)[[condition]]
    args <-
      list(
        df_samples_subset = data,
        protein_names = marker_names,
        condition = "condition",
        group = "donor",
        num_cores = num_cores
      )
    if (method == "cytoglmm") {
      fit <- do.call(CytoGLMM::cytoglmm, args = args)
    } else if (method == "cytoglm") {
      args$num_boot <- num_boot
      fit <- do.call(CytoGLMM::cytoglm, args = args)
    } else
      stop("invalid method")
    return(fit)
  }
