#' @importFrom stats rnorm
#'
rNoParents <- function (N,
                       b0,
                       s) {

  return (rnorm(n = N,
                mean = b0,
                sd = s))

}

#' @importFrom stats as.formula rnorm
#'
rMParents <- function (N,
                       mParents,
                       parentData,
                       b0,
                       b1,
                       s) {

  lmFormula <- c()

  lmFormula[[1]] <- '~ b0'

  for (e in 1:mParents) {

    lmFormula[[e + 1]] <- paste0('b1[[',
                                 e,
                                 ']] * parentData[[',
                                 e,
                                 ']]')

  }

  # Create a formula object to pass to the mean argument of the rnorm function.
  # We subset using [[2]] to get only the right hand side of the formula.
  f <- as.formula(paste(lmFormula,
                        collapse = ' + '))[[2]]

  means <- eval(f)

  return (rnorm(n = N,
                mean = means,
                sd = s))

}
