#' cllNormal
#'
#' Calculates the log likelihood for a node with normal data.
#'
#' @param adjMatrix The adjacency matrix according to the edge directions of the
#' current graph.
#'
#' @param data The data matrix.
#'
#' @param logLikelihood A vector containing the log likelihood for each genetic
#' variant node.
#'
#' @param nGV The number of genetic variants in the graph.
#'
#' @param nNodes The total number of nodes in the graph.
#'
#' @return A vector with the log likelihood for each node.
#'
#' @importFrom stats lm logLik sd
#'
cllNormal <- function (adjMatrix,
                       data,
                       logLikelihood,
                       nGV,
                       nNodes){

  # The log likelihood for the gene expression nodes will start from nGV+1 and
  # go to the value of nNodes.
  for (e in (nGV + 1):nNodes) {

    # Get the number of parents for the current node.
    nParents <- sum(adjMatrix[, e])

    # Loop through each of the gene expression nodes and calculate the log
    # likelihood given the edge directions of the current graph.
    if (nParents == 0) {

      # Store the log likelihood for the current node.
      logLikelihood[[e]] <- sum(dnorm(x = data[, e],
                                      mean = mean(data[, e]),
                                      sd = sd(data[, e]),
                                      log = TRUE))

    } else {

      # Run a linear model on the current node given its parents and calculate
      # the log likelihood.
      logLikelihood[[e]] <- logLik(lm(data[, e] ~ data[, which(adjMatrix[, e] != 0)]))[1]

    }

  }

  return (logLikelihood)

}
