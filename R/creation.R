# sNorm
#
# Generates normal data for source nodes (nodes with no parents).
#
# @param N The number of observations to simulate.
#
# @param b0 The mean of the normal distribution.
#
# @param s The standard deviation of the normal distribution.
#
# @return A vector containing the simulated data.
#
#' @importFrom stats rnorm
#'
sNorm <- function (N,
                   b0,
                   s) {

  return (rnorm(n = N,
                mean = b0,
                sd = s))

}

# sBinom
#
# Generates binomial data for source nodes (nodes with no parents).
#
# @param N The number of observations to simulate.
#
# @param p The probability of success.
#
# @return A vector containing the simulated data.
#
sBinom <- function (N,
                    p) {

  return (sample(x = 0:1,
                 size = N,
                 replace = TRUE,
                 prob = c(p, 1 - p)))

}

# sMulti
#
# Generates multinomial data for source nodes (nodes with no parents).
#
# @param N The number of observations to simulate.
#
# @param q The frequency of the reference allele.
#
# @return A vector containing the simulated data.
#
#' @importFrom stats rnorm
#'
sMulti <- function (N,
                    q) {

  return (sample(x = 0:2,
                 size = N,
                 replace = TRUE,
                 prob = c((1 - q)^2,
                          2 * q * (1 - q),
                          q^2)))

}

# cNorm
#
# @param N The number of observations to simulate.
#
# @param mParents The number of parent nodes.
#
# @param parentData A list containing the data of the parents nodes.
#
# @param b0 The slope of the linear model
# b0 + b1[[1]] * parentData[[1]] + b1[[2]] * parentData[[2]] + ...
#
# @param b1 A vector containing the regression coefficients of the parent nodes.
#
#' @importFrom stats as.formula rnorm
#'
cNorm <- function (N,
                   mParents,
                   parentData,
                   b0,
                   b1,
                   s) {

  # Create a vector that will hold a linear model as a character string.
  lmFormula <- vector(length = mParents + 1)

  # The first element is the slope.
  lmFormula[[1]] <- '~ b0'

  # Loop through the parents for the current node.
  for (e in 1:mParents) {

    # Fill in the remaining elements with the regression coefficient and the
    # data for the corresponding parent node.
    lmFormula[[e + 1]] <- paste0('b1[[',
                                 e,
                                 ']] * parentData[[',
                                 e,
                                 ']]')

  }

  # Create a formula object and evaluate it. The subset ([[2]]) extracts just
  # the right hand side of the formula.
  means <- eval(as.formula(paste(lmFormula,
                                 collapse = ' + '))[[2]])

  # Simulate data from a normal distribution and return it.
  return (rnorm(n = N,
                mean = means,
                sd = s))

}

# cBinom
#
# @param N The number of observations to simulate.
#
# @param mParents The number of parent nodes.
#
# @param parentData A list containing the data of the parents nodes.
#
# @param b0 The slope of the linear model
# b0 + b1[[1]] * parentData[[1]] + b1[[2]] * parentData[[2]] + ...
#
# @param b1 A vector containing the regression coefficients of the parent nodes.
#
#' @importFrom stats as.formula
#'
cBinom <- function (N,
                    mParents,
                    parentData,
                    b0,
                    b1) {

  # Create a vector that will hold a linear model as a character string.
  lmFormula <- vector(length = mParents + 1)

  # The first element is the slope.
  lmFormula[[1]] <- '~ b0'

  # Loop through the parents for the current node.
  for (e in 1:mParents) {

    # Fill in the remaining elements with the regression coefficient and the
    # data for the corresponding parent node.
    lmFormula[[e + 1]] <- paste0('b1[[',
                                 e,
                                 ']] * parentData[[',
                                 e,
                                 ']]')

  }

  # Create a formula object and evaluate it. The subset ([[2]]) extracts just
  # the right hand side of the formula.
  f <- eval(as.formula(paste(lmFormula,
                             collapse = ' + '))[[2]])

  # Inverse logit using the parent nodes.
  prob <- exp(f) / (1 + exp(f))

  # Generate a uniform random variable to determine the state (0 or 1).
  rU <- runif(N, 0, 1)

  # Return the binomial rv.
  return (ifelse(rU < prob, 0, 1))

}
