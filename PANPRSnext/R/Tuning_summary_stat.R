Tuning_summary_stat <- function(
    beta_vec,
    family = "gaussian",
    penalty = "LOG",
    n.tau = 6,
    n.lambdas = 100,
    lambda.min = NULL,
    pfactor = 0.1,
    min_to_max = FALSE
) {
    output <- list()
    output[["lambda"]] <- NULL
    output[["tau"]] <- NULL

    if (family == "gaussian") {
        if (penalty == "LOG") {
            l1.max <- max(abs(beta_vec))
            if (is.null(lambda.min)) {
                min_lambda <- l1.max * pfactor
            } else {
                min_lambda <- lambda.min
            }
            maxTau <- max(abs(beta_vec))
            if (min_to_max) {
                Thresholds <- c(exp(seq(log(min_lambda), log(l1.max), len = n.lambdas)))
            } else {
                Thresholds <- c(exp(seq(log(l1.max), log(min_lambda), len = n.lambdas)))
            }
            tauset <- c(exp(seq(log(1e-6), log(maxTau), len = n.tau)))

            tau <- lambda <- c()

            for (t in 1:length(Thresholds)) {
                thres <- Thresholds[t]
                slambda <- tauset * thres

                tau <- c(tau, tauset)
                lambda <- c(lambda, slambda)
            }

            if (length(tau) != length(lambda)) {
                stop("error")
            }

            output[["lambda"]] <- lambda
            output[["tau"]] <- tau
        }
    }
    return(output)
}
