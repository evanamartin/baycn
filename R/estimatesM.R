#' estimatesM
#'
#' Calculates the estimates of the probability of being in each group of the
#' current genetic variant.
#'
#' @param adjMatrix The adjacency matrix created from the edge directions of the
#' current graph.
#'
#' @param data A matrix with the nodes across the columns and the observations
#' down the rows.
#'
#' @param node The column number of the current node, the node for which we will
#' calculate the estimates for the mean and standard deviation.
#'
#' @return A list containing the number of parents of the current node, a vector
#' with the estimates of the beta coefficients, and the data matrix containing
#' only the columns of the current node and its parents.
#'
#' @importFrom MASS polr
#'
#' @export
#'
estimatesM <- function (adjMatrix,
                        data,
                        node) {

  # Get the number of parents for the current node.
  nParents <- sum(adjMatrix[, node])

  if (nParents == 0) {

    # Calculate the counts for each level of the current genetic variant.
    counts <- as.vector(table(data[, node]))

    # If there are no parents the probability of each level of the genetic
    # variant vector will be calculated directly from the data.
    estimates <- list(counts = counts,
                      probs = counts / sum(counts))

    # parentData is part of the returned list and must be assigned a value even
    # if it won't be used in any calculations.
    parentData <- NULL

  } else {

    # If a node has one or more parents multinomial logistic regression is used
    # to calculate the probability of each level occuring.

    # Take only the columns of the data matrix that are the parents of the child
    # node.
    parentData <- data.frame(y = data[, node],
                             data[, which(adjMatrix[, node] != 0)])

    # Run ordinal logistic regression with the function polr from the package
    # MASS the response variable must be a factor. There is one beta coefficient
    # and a slope for the k - 1 levels of the response.
    olrSummary <- summary(polr(as.factor(y) ~ .,
                               data = parentData))

    # Extract the coefficients from the summary. The beta coefficient is first
    # followed by the intercepts.
    estimates <- as.vector(olrSummary$coefficients[, 1])


  }

  return (list(nParents = nParents,
               estimates = estimates,
               parentData = parentData))

}
