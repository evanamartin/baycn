#' makeFormula
#'
#' Takes the estimates calculated in the linear model from the estimates
#' function and creates a linear model to pass as the mean to the dnorm
#' function.
#'
#' @param estimates The estimates of the linear model calculated in the
#' estimates function.
#'
#' @param nParents The number of parents of the current node
#'
#' @param parentData The data matrix containing only the columns of the
#' current node and the parents of the current node.
#'
#' @return A vector of means calculated by using the estimates from the
#' estimates function.
#'
#' @importFrom stats as.formula
#'
#' @export
#'
makeFormula <- function (estimates,
                         nParents,
                         parentData) {

  lmFormula <- c()

  lmFormula[[1]] <- '~ estimates[1]'

  for (e in 1:nParents) {

    lmFormula[[e + 1]] <- paste0('estimates[[',
                                 e + 1,
                                 ']] * parentData[, ',
                                 e + 1,
                                 ']')

  }

  # Create a formula object to pass to the mean argument of the dnorm function.
  # We subset using [[2]] to get only the right hand side of the formula.
  f <- as.formula(paste(lmFormula,
                        collapse = ' + '))[[2]]

  return (eval(f))

}
