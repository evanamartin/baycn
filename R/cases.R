#' rNoParents
#'
#' Creates data from a normal distribution. This function is used for generating
#' data for nodes with no parents.
#'
#' @param N The number of observations to simulate.
#'
#' @param b0 The mean of the normal distribution.
#'
#' @param s The standard deviation of the normal distribution.
#'
#' @return A vector containing the simulated data.
#'
#' @importFrom stats rnorm
#'
rNoParents <- function (N,
                       b0,
                       s) {

  return (rnorm(n = N,
                mean = b0,
                sd = s))

}

#' rMParents
#'
#' @param N The number of observations to simulate.
#'
#' @param mParents The number of parents of the current gene.
#'
#' @param parentData A list containing the data of the parents of the current
#' gene.
#'
#' @param b0 The slope of the linear model
#' b0 + b1[[1]] * parentData[[1]] + b1[[2]] * parentData[[2]] + ...
#'
#' @param b1 A vector containing the regression coefficients of the parents of
#' the gene for which data is being simulated.
#'
#' @param s The standard deviation of the gene for which data is being
#' simulated.
#'
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
