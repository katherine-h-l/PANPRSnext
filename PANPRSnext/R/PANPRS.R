#' Run the package on the provided data set
#' @export
test_pkg <- function(
  debug_output = FALSE
) {
  data("summaryZ")
  data("Nvec")
  data("plinkLD")
  data("funcIndex")
  output <- gsfPEN_R(
    summary_z = summaryZ, # nolint: object_usage_linter.
    n_vec = Nvec, # nolint: object_usage_linter.
    plinkLD = plinkLD, # nolint: object_usage_linter.
    func_index = funcIndex, # nolint: object_usage_linter.
    debug_output = debug_output
  )

  return(output)
}

#' Run the package on the provided data set
#' DEBUG VERSION
#' @export
debug_test_pkg <- function(
  summaryZ,
  funcIndex,
  debug_output = FALSE,
  sub_tuning = 1,
  lambda_vec_limit_len = c(1.5, 1)
) {
  data("Nvec")
  data("plinkLD")
  output <- gsfPEN_R(
    summary_z = summaryZ, # nolint: object_usage_linter.
    n_vec = Nvec, # nolint: object_usage_linter.
    plinkLD = plinkLD, # nolint: object_usage_linter.
    func_index = funcIndex, # nolint: object_usage_linter.
    sub_tuning = sub_tuning,
    lambda_vec_limit_len = lambda_vec_limit_len,
    debug_output = debug_output
  )

  return(output)
}

#' Main Function
#' @export
gsfPEN_R <- function(
    summary_z,
    n_vec,
    plinkLD,
    func_index,
    n_iter = 1000,
    upper_val = NULL,
    breaking = 1,
    z_scale = 1,
    tuning_matrix = NULL,
    p_threshold = NULL,
    p_threshold_params = c(0.5, 10^-4, 4),
    tau_factor = c(1 / 25, 1, 3),
    len_lim_lambda = 4,
    sub_tuning = 4,
    lim_lambda = c(0.5, 0.9),
    len_lambda = 4,
    lambda_vec = NULL,
    lambda_vec_limit_len = c(1.5, 3),
    df_max = NULL,
    debug_output = FALSE
) {

  time <- proc.time()

  if (z_scale != 1) {
    error("Tuning values set-up for multiple traits analysis requires z_scale=1.") # nolint: object_usage_linter.
  }

  num_func <- ncol(func_index)

  P <- nrow(summary_z)
  Q <- length(n_vec)

  summary_betas <- matrix(0, nrow = P, ncol = Q)
  SD_vec <- matrix(0, nrow = P, ncol = Q)

  for (ii in 1:Q) {
    summary_betas[, ii] <- summary_z[, ii] / sqrt(n_vec[ii])
    SD_vec[, ii] <- 1 / sqrt(n_vec[ii])
  }

  rownames(summary_betas) <- rownames(summary_z)

  if (is.null(df_max)) {
    df_max <- ceiling(0.7 * P)
  }

  if (is.null(p_threshold)) {
    p_threshold <- seq(p_threshold_params[1], p_threshold_params[2], length.out = p_threshold_params[3])
  }

  if (any(c(is.null(tuning_matrix), is.null(lambda_vec)))) {
    median_val <- median(apply(abs(summary_betas), 1, sum), na.rm = TRUE)
    tau_vec <- sort(median_val * tau_factor)
    lim_lambda <- quantile(abs(summary_z[, 1]), lim_lambda)

    output <- Tuning_setup_group_func(
      lambda_vec,
      lambda_vec_limit_len,
      p_threshold,
      num_func,
      tau_vec,
      sub_tuning,
      lim_lambda,
      len_lambda,
      len_lim_lambda,
      median_val
    )

    func_lambda <- output$func_lambda
    lambda_vec <- output$lambda_vec
    tuning_matrix <- output$tuning_matrix

    rm(output)
  } else {
    func_lambda <- permutations( # nolint: object_usage_linter.
      length(lambda_vec),
      num_func,
      repeats.allowed = TRUE
    ) - 1
  }

  beta_index <- c(seq_len(nrow(summary_betas))) - 1
  SNP_names <- rownames(summary_betas)

  # Takes only the SNPs that are present in both the PLINK data set and the GWAs
  ld_J <- PlinkLD_transform(plinkLD, SNP_names)

  J_id_matrix <- matrix(nrow = nrow(ld_J), ncol = 2)

  mat1 <- match(ld_J[, 1], SNP_names)
  mat2 <- match(ld_J[, 2], SNP_names)

  J_id_matrix[, 1] <- beta_index[mat1]
  J_id_matrix[, 2] <- beta_index[mat2]

  ld_J[, 1] <- J_id_matrix[, 1]
  ld_J[, 2] <- J_id_matrix[, 2]

  od <- order(J_id_matrix[, 1], J_id_matrix[, 2], decreasing = FALSE)
  ld_J <- ld_J[od, ]


  wind <- which(!beta_index %in% ld_J[, 1])
  num_indices <- length(wind)

  index_J <- -1
  if (num_indices > 0) {
    index_J <- beta_index[wind]
  }

  counts <- table(ld_J[, 1])
  num_SNP <- length(counts)

  index_S <- c(0, cumsum(counts)[-num_SNP])
  index_E <- cumsum(counts) - 1

  index_matrix <- matrix(nrow = num_SNP, ncol = 3)
  index_matrix[, 1] <- as.numeric(names(counts))

  index_matrix[, 2] <- index_S
  index_matrix[, 3] <- index_E

  ld_vec <- ld_J[, 3]
  ld_J <- ld_J[, 2]


  if (is.null(upper_val)) {
    upper_val <- ceiling(max(abs(summary_betas), na.rm = TRUE) * 50)
  }

  if (nrow(func_index) != nrow(summary_betas)) {
    stop("nrow of summary_betas and row of func_index do not match.")
  }


  z_matrix <- matrix(1 - unlist(func_index), nrow = nrow(func_index), ncol = ncol(func_index), byrow = TRUE)
  rownames(z_matrix) <- rownames(func_index)

  sum_func_index <- apply(z_matrix, 1, sum)
  Ifunc_SNP <- rep(0, P)
  Ifunc_SNP[which(sum_func_index != 0)] <- 1

  lambda0_vec <- abs(-qnorm(p_threshold / 2))


  nrow_index_matrix <- nrow(index_matrix)
  ncol_index_matrix <- ncol(index_matrix)

  nrow_z_matrix <- nrow(z_matrix)
  ncol_z_matrix <- ncol(z_matrix)

  nrow_func_lambda <- nrow(func_lambda)
  ncol_func_lambda <- ncol(func_lambda)

  nrow_tuning_matrix <- nrow(tuning_matrix)
  ncol_tuning_matrix <- ncol(tuning_matrix)

  nrow_all_tuning_matrix <- nrow_tuning_matrix * length(p_threshold) * nrow_func_lambda
  ncol_all_tuning_matrix <- (num_func + 1) + 2

  nrow_beta_matrix <- nrow_all_tuning_matrix
  ncol_beta_matrix <- P * Q

  dims <- c(
    num_SNP,                # 1
    P,                      # 2
    Q,                      # 3
    nrow_index_matrix,      # 4
    ncol_index_matrix,      # 5
    nrow_z_matrix,          # 6
    ncol_z_matrix,          # 7
    nrow_func_lambda,       # 8
    ncol_func_lambda,       # 9
    nrow_tuning_matrix,     # 10
    ncol_tuning_matrix,     # 11
    nrow_all_tuning_matrix, # 12 == num_tuning
    ncol_all_tuning_matrix, # 13
    nrow_beta_matrix,       # 14
    ncol_beta_matrix        # 15
  )

  params <- c(
    upper_val,              # 1
    n_iter,                 # 2
    breaking,               # 3
    z_scale,                # 4
    df_max,                 # 5
    length(p_threshold),    # 6
    num_indices             # 7
  )

  print(paste0("Number of total tuning combinations = ", nrow_all_tuning_matrix))

  Z <- gsfPEN_cpp(
    summary_betas,
    ld_J,
    index_matrix,
    index_J,
    ld_vec,
    SD_vec,
    tuning_matrix,
    lambda0_vec,
    z_matrix,
    lambda_vec,
    func_lambda,
    Ifunc_SNP,
    dims,
    params
  )

  beta_matrix <- Z$beta_matrix
  colnames(beta_matrix) <- paste0(
    rep(SNP_names, times = Q),
    ".trait",
    rep(c(1:Q), each = P)
  )

  all_tuning_matrix <- Z$all_tuning_matrix
  colnames(all_tuning_matrix) <- c(
    "lambda0",
    paste0("lambdaf", c(1:num_func)),
    "lamda2",
    "tau2"
  )

  num_iter_vec <- Z$num_iter_vec

  # Remove the tuning combinations that did not converge (correspons to -2 in num_iter_vec)
  if (!debug_output) {
    print("Removing tuning combinations that did not converge.")
    converge_index <- which(num_iter_vec > 0)
    num_iter_vec <- num_iter_vec[converge_index]
    beta_matrix <- beta_matrix[converge_index, ]
    all_tuning_matrix <- all_tuning_matrix[converge_index, ]
  }

  output <- Clean_results(
    beta_matrix = beta_matrix,
    num_iter_vec = num_iter_vec,
    all_tuning_matrix = all_tuning_matrix
  )

  time <- proc.time() - time
  print(paste0("Time elapsed: ", time[3], " seconds"))

  return(output)
}