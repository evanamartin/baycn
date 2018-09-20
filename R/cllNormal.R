#' cllNormal
#'
#' Calculates the log likelihood for a node with normal data.
#'
#' @param adjMatrix The adjacency matrix according to the edge directions of the
#' current graph.
#'
#' @param data The data matrix.
#'
#' @param k A vector containing the number of parents for each genetic variant
#' node. This vector will be used if the scoreFun argument is set to 'BIC'.
#'
#' @param logLikelihood A vector containing the log likelihood for each genetic
#' variant node.
#'
#' @param nGV The number of genetic variants in the graph.
#'
#' @param nNodes The total number of nodes in the graph.
#'
#' @return A list containing two elements. The first element is a vector with
#' the log likelihood for each node. The second element is a vector with the
#' number of parents for each node.
#'
#' @importFrom stats lm logLik sd
#'
cllNormal <- function (adjMatrix,
                       data,
                       k,
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

      ll <- sum(dnorm(x = data[, e],
                      mean = mean(data[, e]),
                      sd = sd(data[, e]),
                      log = TRUE))

    } else {

      # Take only the columns of the data that pertain to the child node and its
      # parents.
      parentData <- data.frame(y = data[, e],
                               data[, which(adjMatrix[, e] != 0)])

      model <- lm(y ~ .,
                  data = parentData)

      ll <- logLik(model)[1]

    }

    # Store the number of parents for each node. This will be used if the
    # scoreFun argument is set to 'BIC'.
    k[[e]] <- nParents

    # Store the log likelihood for the current node.
    logLikelihood[[e]] <- ll

  }

  return (list(logLikelihood = logLikelihood,
               nParents = k))

}
