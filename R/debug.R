#' Run gsPEN on the provided data set
#' @param ... Additional arguments to pass to gsPEN_R
#' @importFrom utils data
#' @export
test_gsPEN <- function(...) {
  summaryZ <- Nvec <- plinkLD <- funcIndex <- NULL
  data("summaryZ", envir = environment())
  data("Nvec", envir = environment())
  data("plinkLD", envir = environment())
  data("funcIndex", envir = environment())

  output <- gsPEN_R(
    summary_z = summaryZ, # nolint: object_usage_linter.
    n_vec = Nvec, # nolint: object_usage_linter.
    plinkLD = plinkLD, # nolint: object_usage_linter.
    ...
  )

  return(output)
}

#' Run gsfPEN on the provided data set
#' @param ... Additional arguments to pass to gsfPEN_R
#' @importFrom utils data
#' @export
test_gsfPEN <- function(...) {
  summaryZ <- Nvec <- plinkLD <- funcIndex <- NULL
  data("summaryZ", envir = environment())
  data("Nvec", envir = environment())
  data("plinkLD", envir = environment())
  data("funcIndex", envir = environment())

  output <- gsfPEN_R(
    summary_z = summaryZ, # nolint: object_usage_linter.
    n_vec = Nvec, # nolint: object_usage_linter.
    plinkLD = plinkLD, # nolint: object_usage_linter.
    func_index = funcIndex, # nolint: object_usage_linter.
    ...
  )

  return(output)
}
