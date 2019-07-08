#' cllMultinom
#'
#' Calculates the log likelihood for a node with multinomial data.
#'
#' @param adjMatrix The adjacency matrix according to the edge directions of the
#' current graph.
#'
#' @param data The data matrix.
#'
#' @param nGV The number of genetic variants in the graph.
#'
#' @param nNodes The total number of nodes in the graph.
#'
#' @return A vector with the log likelihood for each node.
#'
#' @importFrom MASS polr
#'
#' @importFrom stats dmultinom logLik
#'
cllMultinom <- function (adjMatrix,
                         data,
                         nGV,
                         nNodes){

  # A vector to hold the log likelihood for each node in the graph.
  logLikelihood <- vector(mode = 'numeric',
                          length = nNodes)

  # If there are no genetic variants in the data return an empty list for the
  # log likelihood and the number of parents for each node. This will be used by
  # the cllNormal function.
  if (nGV == 0) {

    return(logLikelihood)

  }

  # The log likelihood for the genetic variant nodes will start from 1 and go to
  # the value of nGV.
  for (e in 1:nGV) {

    # Get the number of parents for the current node.
    nParents <- sum(adjMatrix[, e])

    # Loop through each of the genetic variant nodes and calculate the log
    # likelihood given the edge directions of the current graph.
    if (nParents == 0) {

      # Calculate the counts for each level of the current genetic variant.
      counts <- as.vector(table(data[, e]))

      # Calculate the probability of each level of counts.
      lprob <- log(counts / sum(counts))

      # Store the log likelihood for the current node.
      logLikelihood[[e]] <- sum(lprob * counts)

    } else {

      # Take only the columns of the data that pertain to the child node and its
      # parents.
      parentData <- data.frame(y = data[, e],
                               data[, which(adjMatrix[, e] != 0)])

      # Run ordered logistic regression on the current genetic variant given its
      # parents.
      model <- polr(as.factor(y) ~ .,
                    data = parentData)

      # Store the log likelihood for the current node.
      logLikelihood[[e]] <- logLik(model)[1]

    }

  }

  return (logLikelihood)

}
