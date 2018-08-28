#' estimates
#'
#' Calculates the estimates of the mean and standard deviation for the current
#' node. If the current node doesn't have any parent then the mean and standard
#' deviation are simply calculated from the data. If the current node has at
#' least one parent the mean is estimated by a linear model with the current
#' node as the independent variable and parent nodes as independent variables.
#' The standard deviation is estimated by the standard error of the residuals.
#' These estimates are later used to calculate the log likelihood of the current
#' node.
#'
#' @param adjMatrix The adjacency matrix created from the DNA of the current
#' individual.
#'
#' @param data A matrix with the nodes across the columns and the observations
#' along the rows.
#'
#' @param node The column number of the current node, the node for which we will
#' calculate the estimates for the mean and standard deviation.
#'
#' @return A list containing the number of parents of the current node, a vector
#' with the estimates of the beta coefficients and standard error of the
#' residuals, and the data matrix containing only the columns of the current
#' node and its parents.
#'
#' @importFrom stats sd lm
#'
#' @export
#'
estimates <- function (adjMatrix,
                       node,
                       data) {

  nParents <- sum(adjMatrix[, node])

  if (nParents == 0) {

    # If there are no parents the mean and standard deviation of the data
    # are used to calculate the log likelihood
    estimates <- c()

    estimates[[1]] <- mean(data[, node])

    estimates[[2]] <- sd(data[, node])

    parentData <- data.frame(y = data[, node])

  } else {

    # If a node has one or more parents a linear model is used for the mean of
    # the data and the standard error of the residuals is used to estimate the
    # standard deviation.

    # Take only the columns of the data that pertain to the child node and its
    # parents.
    parentData <- data.frame(y = data[, node],
                             data[, which(adjMatrix[, node] != 0)])

    lmSummary <- summary(lm(y ~ .,
                            data = parentData))

    estimates <- as.vector(lmSummary$coefficients[, 1])

    estimates <- append(estimates, lmSummary$sigma)

  }

  return (list(nParents = nParents,
               estimates = estimates,
               parentData = parentData))

}
