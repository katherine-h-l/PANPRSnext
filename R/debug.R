#' Run gsPEN on the provided data set
#' @export
test_gsPEN <- function(...) {
  data("summaryZ")
  data("Nvec")
  data("plinkLD")
  data("funcIndex")

  output <- gsPEN_R(
    summary_z = summaryZ, # nolint: object_usage_linter.
    n_vec = Nvec, # nolint: object_usage_linter.
    plinkLD = plinkLD, # nolint: object_usage_linter.
    ...
  )

  return(output)
}

#' Run gsfPEN on the provided data set
#' @export
test_gsfPEN <- function(...) {
  data("summaryZ")
  data("Nvec")
  data("plinkLD")
  data("funcIndex")

  output <- gsfPEN_R(
    summary_z = summaryZ, # nolint: object_usage_linter.
    n_vec = Nvec, # nolint: object_usage_linter.
    plinkLD = plinkLD, # nolint: object_usage_linter.
    func_index = funcIndex, # nolint: object_usage_linter.
    ...
  )

  return(output)
}
