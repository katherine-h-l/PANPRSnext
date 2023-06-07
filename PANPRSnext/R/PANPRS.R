#' Testing function
#' @export
test_pkg <- function() {
  data("summaryZ")
  data("Nvec")
  data("plinkLD")
  data("funcIndex")
  output <- gsfPEN_R(
    summary_z = summaryZ,
    n_vec = Nvec,
    plinkLD = plinkLD,
    func_index = funcIndex
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
    p.Threshold = NULL,
    p.Threshold_params = c(0.5, 10^-4, 4),
    tau_factor = c(1 / 25, 1, 3),
    llim_length = 4,
    sub_tuning = 4,
    lim_lambda = c(0.5, 0.9),
    len_lambda = 4,
    lambda_vec_func = NULL,
    lambda_vec_func_limit_len = c(1.5, 3),
    df_max = NULL,
    debug_output = 0
) {

  time <- proc.time()

  if (z_scale != 1) {
    error("Tuning values set-up for multiple traits analysis requires z_scale=1.")
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


  if (is.null(p.Threshold)) {
    p.Threshold <- seq(p.Threshold_params[1], p.Threshold_params[2], length.out = p.Threshold_params[3])
  }


  if (any(c(is.null(tuning_matrix), is.null(lambda_vec_func)))) {
    median_val <- median(apply(abs(summary_betas), 1, sum), na.rm = T)
    tau_vec <- sort(median_val * tau_factor)
    lim_lambda <- quantile(abs(summary_z[, 1]), lim_lambda)

    output <- Tuning_setup_group_func(
      lambda_vec_func, lambda_vec_func_limit_len, p.Threshold, num_func, tau_vec, sub_tuning, lim_lambda,
      len_lambda, llim_length, median_val
    )

    func_lambda <- output[[1]]
    lambda_vec_func <- output[[2]]
    tuning_matrix <- output[[3]]
  } else {
    func_lambda0 <- permutations(length(lambda_vec_func), num_func, repeats.allowed = T)
    func_lambda <- func_lambda0 - 1
  }



  beta_index <- c(1:nrow(summary_betas)) - 1
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

  od <- order(J_id_matrix[, 1], J_id_matrix[, 2], decreasing = F)
  ld_J <- ld_J[od, ]

  wind <- which(!beta_index %in% ld_J[, 1])

  index_J <- -1

  if (length(wind) > 0) {
    index_J <- beta_index[wind]
  }

  counts <- table(ld_J[, 1])
  num_SNP <- length(counts)


  index_S <- c(0, cumsum(counts)[-num_SNP])
  index_E <- cumsum(counts) - 1

  index_matrix <- matrix(nrow = num_SNP, ncol = 3)

  index_matrix[, 1] <- as.numeric(names(counts))

  nrow_index_matrix <- num_SNP
  ncol_index_matrix <- ncol(index_matrix)


  index_matrix[, 2] <- index_S
  index_matrix[, 3] <- index_E


  ld_vec <- ld_J[, 3]

  ld_J <- ld_J[, 2]

  length_ldJ <- length(ld_J)

  if (is.null(upper_val)) {
    upper_val <- ceiling(max(abs(summary_betas), na.rm = T) * 50)
  }


  if (nrow(func_index) != nrow(summary_betas)) {
    stop("nrow of summary_betas and row of func_index do not match.")
  }

  z_matrix <- matrix(1 - unlist(func_index), nrow = nrow(func_index), ncol = ncol(func_index), byrow = TRUE)
  rownames(z_matrix) <- rownames(func_index)

  sum_func_index <- apply(z_matrix, 1, sum)
  Ifunc_SNP <- rep(0, P)
  Ifunc_SNP[which(sum_func_index != 0)] <- 1


  nrow_func_lambda <- nrow(func_lambda)
  ncol_func_lambda <- ncol(func_lambda)
  leng.p.Threshold <- length(p.Threshold)

  nrow_z_matrix <- nrow(z_matrix)
  ncol_z_matrix <- ncol(z_matrix)

  num_tuning <- nrow(tuning_matrix) * leng.p.Threshold * nrow_func_lambda
  print(paste0("Number of total tuning combinations = ", num_tuning))


  nrow_all_tuning_matrix <- num_tuning
  ncol_all_tuning_matrix <- (num_func + 1) + 2
  all_tuning_matrix <- matrix(0, nrow = nrow_all_tuning_matrix, ncol = ncol_all_tuning_matrix)

  lambda0_vec <- abs(-qnorm(p.Threshold / 2))

  nrow_beta_matrix <- num_tuning
  ncol_beta_matrix <- P * Q

  dims <- c(
    P,                      # 1
    Q,                      # 2
    nrow_index_matrix,      # 3
    ncol_index_matrix,      # 4
    nrow_z_matrix,          # 5
    ncol_z_matrix,          # 6
    nrow_func_lambda,       # 7
    ncol_func_lambda,       # 8
    nrow(tuning_matrix),    # 9
    ncol(tuning_matrix),    # 10
    nrow_all_tuning_matrix, # 11 == num_tuning
    ncol_all_tuning_matrix, # 12
    nrow_beta_matrix,       # 13
    ncol_beta_matrix        # 14
  )

  params <- c(
    n_iter,           # 1
    breaking,         # 2
    z_scale,          # 3
    df_max,           # 4
    leng.p.Threshold, # 5
    length(wind)      # 6
  )

  num_iter_vec <- rep(0, num_tuning)
  beta_matrix <- matrix(0, nrow_beta_matrix, ncol_beta_matrix)

  Z <- gsfPEN_cpp(
    summary_betas,
    ld_J,
    num_iter_vec,
    index_matrix,
    index_J,
    ld_vec,
    upper_val,
    SD_vec,
    tuning_matrix,
    beta_matrix,
    lambda0_vec,
    z_matrix,
    all_tuning_matrix,
    lambda_vec_func,
    func_lambda,
    Ifunc_SNP,
    dims,
    params
  )

  return(dims)

  # Z <- .C("gsfFunc",
  #         as.double(t(summary_betas)),
  #         as.integer(ld_J),
  #         as.integer(dims1),
  #         num_iter_vec = as.integer(num_iter_vec),
  #         as.integer(t(index_matrix)),
  #         as.integer(index_J),
  #         as.double(ld_vec),
  #         as.integer(ChrIndexBeta),
  #         as.double(upper_val),
  #         as.double(Init_summary_betas),
  #         as.double(t(SD_vec)),
  #         as.double(t(tuning_matrix)),
  #         beta_matrix = as.double(t(beta_matrix)),
  #         penalty,
  #         as.double(lambda0_vec),
  #         as.double(t(z_matrix)),
  #         as.integer(dims2),
  #         all_tuning_matrix = as.double(t(all_tuning_matrix)),
  #         as.double(lambda_vec_func),
  #         as.integer(t(func_lambda)),
  #         as.integer(Ifunc_SNP),
  #         PACKAGE = "PANPRS"
  # )
  #
  # beta_matrix <- matrix(Z$beta_matrix, nrow = num_tuning, ncol = ncol_beta_matrix, byrow = TRUE)
  # colnames(beta_matrix) <- paste0(rep(SNP_names, times = Q), ".trait", rep(c(1:Q), each = P))
  #
  # all_tuning_matrix <- matrix(Z$all_tuning_matrix, nrow = nrow_all_tuning_matrix, ncol = ncol_all_tuning_matrix, byrow = TRUE)
  # num_iter_vec <- Z$num_iter_vec
  #
  # colnames(all_tuning_matrix) <- c("lam0", paste0("lam.f", c(1:num_func)), "lam2", "tau")
  #
  # output <- Cleaning(beta_matrix = beta_matrix, num_iter_vec = num_iter_vec, all_tuning_matrix = all_tuning_matrix)
  # num_iter_vec <- output[[1]]
  # beta_matrix <- output[[2]]
  # all_tuning_matrix <- output[[3]]
  # rm(output)
  # if (debug_output == 0) {
  #   convergeIndex <- which(num_iter_vec > 0)
  #   num_iter_vec <- num_iter_vec[convergeIndex]
  #   beta_matrix <- beta_matrix[convergeIndex, ]
  #   all_tuning_matrix <- all_tuning_matrix[convergeIndex, ]
  # }
  #
  # ll <- list(beta_matrix = beta_matrix, num_iter_vec = num_iter_vec, all_tuning_matrix = all_tuning_matrix)
  #
  # time <- proc.time() - time
  # print(paste0("Time elapsed: ", time[3], " seconds"))
  #
  # return(ll)

}
