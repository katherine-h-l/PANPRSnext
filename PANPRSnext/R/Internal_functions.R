Tuning_setup_group_only <- function(
    tau_vec,
    sub_tuning,
    lim_lambda,
    len_lambda,
    llim_length,
    median_val
) {
  lambda_vec <- seq((min(lim_lambda)), (max(lim_lambda)), len = len_lambda)
  tuning_matrix <- cbind(lambda_vec, 1, 0, 1)

  topvec <- seq(from = min(lim_lambda), to = max(lim_lambda), len = llim_length)

  for (i1 in 1:length(tau_vec)) {
      tau <- tau_vec[i1]

      for (t1 in 2:length(topvec)) {
          top <- topvec[t1]
          lambda_vec <- seq(min(lim_lambda), (top - 0.05), len = sub_tuning)

          temp <- (top - lambda_vec) * (median_val + tau)
          tuning_matrix <- rbind(tuning_matrix, cbind(lambda_vec, tau, temp, tau))
      }
  }
  return(tuning_matrix)
}




Tuning_setup_group_func <- function(
    lambda_vec_func,
    lambda_vec_func_limit_len,
    p.Threshold,
    num_func,
    tau_vec,
    sub_tuning,
    lim_lambda,
    len_lambda,
    llim_length,
    median_val
) {
  lambda_vec <- seq((min(lim_lambda)), (max(lim_lambda)), len = len_lambda)

  topvec <- seq(from = min(lim_lambda), to = max(lim_lambda), len = llim_length)

  start <- 1
  for (i1 in 1:length(tau_vec)) {
      tau <- tau_vec[i1]

      for (t1 in 2:length(topvec)) {
          top <- topvec[t1]
          lambda_vec <- seq(min(lim_lambda), (top - 0.05), len = sub_tuning)

          temp <- (top - lambda_vec) * (median_val + tau)
          if (start == 1) {
              tuning_matrix <- cbind(lambda_vec, tau, temp, tau)
              start <- 0
          } else {
              tuning_matrix <- rbind(tuning_matrix, cbind(lambda_vec, tau, temp, tau))
          }
      }
  }
  if (is.null(lambda_vec_func)) {
      lambda_vec_func <- seq(0, lambda_vec_func_limit_len[1], length.out = lambda_vec_func_limit_len[2])
  }
  n.threshold <- length(p.Threshold)
  leng.lambda <- length(lambda_vec_func)
  func_lambda0 <- permutations(length(lambda_vec_func), num_func, repeats.allowed = T)
  func_lambda <- func_lambda0 - 1

  output <- list()
  output[[1]] <- func_lambda
  output[[2]] <- lambda_vec_func
  output[[3]] <- tuning_matrix

  return(output)
}

Nonzero <- function(xx) {
  return(length(which(xx != 0)))
}

Cleaning <- function(
    beta_matrix,
    num_iter_vec,
    all_tuning_matrix
) {
  tuningvec <- apply(all_tuning_matrix, 1, paste0, collapse = ":")
  uniqvec <- unique(tuningvec)
  mat <- match(uniqvec, tuningvec)
  num_iter_vec <- num_iter_vec[mat]
  beta_matrix <- beta_matrix[mat, ]
  all_tuning_matrix <- all_tuning_matrix[mat, ]
  num_counts <- apply(beta_matrix, 1, Nonzero)
  od <- order(num_counts)
  num_iter_vec <- num_iter_vec[od]
  beta_matrix <- beta_matrix[od, ]
  all_tuning_matrix <- all_tuning_matrix[od, ]
  output <- list()
  output[[1]] <- num_iter_vec
  output[[2]] <- beta_matrix
  output[[3]] <- all_tuning_matrix
  return(output)
}
