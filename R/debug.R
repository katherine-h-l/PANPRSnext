#' Run the package on the provided data set
#' @export
test_pkg <- function(
    debug_output = FALSE) {
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
    lambda_vec_limit_len = c(1.5, 1)) {
  data("summaryZ")
  data("Nvec")
  data("plinkLD")
  data("funcIndex")
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
