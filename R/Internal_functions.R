Tuning_setup_group_only <- function(
    tau_vec,
    sub_tuning,
    lim_lambda,
    len_lambda,
    len_lim_lambda,
    median_val) {
  lambda_vec <- seq(min(lim_lambda), max(lim_lambda), len = len_lambda)
  tuning_matrix <- cbind(lambda_vec, 1, 0, 1)

  top_vec <- seq(min(lim_lambda), max(lim_lambda), len = len_lim_lambda)

  for (i in seq_along(tau_vec)) {
    tau <- tau_vec[i]

    for (t in 2:length(top_vec)) {
      top <- top_vec[t]
      lambda_vec <- seq(min(lim_lambda), (top - 0.05), len = sub_tuning)

      temp <- (top - lambda_vec) * (median_val + tau)
      tuning_matrix <- rbind(tuning_matrix, cbind(lambda_vec, tau, temp, tau))
    }
  }
  return(tuning_matrix)
}

Tuning_setup_group_func <- function(
    lambda_vec,
    lambda_vec_limit_len,
    p_threshold,
    num_func,
    tau_vec,
    sub_tuning,
    lim_lambda,
    len_lambda,
    len_lim_lambda,
    median_val) {
  lim_lambda_vec <- seq(min(lim_lambda), max(lim_lambda), len = len_lambda)

  top_vec <- seq(min(lim_lambda), max(lim_lambda), len = len_lim_lambda)

  start <- 1
  for (i in seq_along(tau_vec)) {
    tau <- tau_vec[i]

    for (t in 2:length(top_vec)) {
      top <- top_vec[t]
      lim_lambda_vec <- seq(min(lim_lambda), (top - 0.05), len = sub_tuning)

      temp <- (top - lim_lambda_vec) * (median_val + tau)
      if (start == 1) {
        tuning_matrix <- cbind(lim_lambda_vec, tau, temp, tau)
        start <- 0
      } else {
        tuning_matrix <- rbind(tuning_matrix, cbind(lim_lambda_vec, tau, temp, tau))
      }
    }
  }

  if (is.null(lambda_vec)) {
    lambda_vec <- seq(0, lambda_vec_limit_len[1], length.out = lambda_vec_limit_len[2])
  }

  func_lambda <- permutations( # nolint: object_usage_linter.
    length(lambda_vec),
    num_func,
    repeats.allowed = TRUE
  ) - 1

  output <- list(
    func_lambda = func_lambda,
    lambda_vec = lambda_vec,
    tuning_matrix = tuning_matrix
  )

  return(output)
}

Non_zero <- function(xx) {
  return(length(which(xx != 0)))
}

Clean_results <- function(
    beta_matrix,
    num_iter_vec,
    all_tuning_matrix) {
  # Selects only the unique rows, corresponding to unique tuning parameters
  tuning_vec <- apply(all_tuning_matrix, 1, paste0, collapse = ":")
  unique_vec <- unique(tuning_vec)
  mat <- match(unique_vec, tuning_vec)

  num_iter_vec <- num_iter_vec[mat]
  beta_matrix <- beta_matrix[mat, ]
  all_tuning_matrix <- all_tuning_matrix[mat, ]

  # Reorders the rows by the number of non-zero coefficients
  num_counts <- apply(beta_matrix, 1, Non_zero)
  ord <- order(num_counts)

  num_iter_vec <- num_iter_vec[ord]
  beta_matrix <- beta_matrix[ord, ]
  all_tuning_matrix <- all_tuning_matrix[ord, ]

  output <- list(
    beta_matrix = beta_matrix,
    num_iter_vec = num_iter_vec,
    all_tuning_matrix = all_tuning_matrix
  )
  return(output)
}
