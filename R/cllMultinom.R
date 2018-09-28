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
#' @return A list containing two elements. The first element is a vector with
#' the log likelihood for each node. The second element is a vector with the
#' number of parents for each node.
#'
#' @importFrom MASS polr
#'
#' @importFrom stats dmultinom logLik
#'
#' @export
#'
cllMultinom <- function (adjMatrix,
                         data,
                         nGV,
                         nNodes){

  # A vector to hold the number of parents for each node in the graph.
  k <- c()

  # A vector to hold the log likelihood for each node in the graph.
  logLikelihood <- c()

  # If there are no genetic variants in the data return an empty list for the
  # log likelihood and the number of parents for each node. This will be used by
  # the cllNormal function.
  if (nGV == 0) {

    return(list(logLikelihood = logLikelihood,
                nParents = k))

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

      # Create a vector to hold the log of the probability of seeing each value
      # the observed number of times.
      lp <- c()

      # Loop through each of the levels of the genetic variant to calculate the
      # log likelihood.
      for (v in 1:length(counts)) {

        lp[[v]] <- lprob[[v]] * counts[[v]]

      }

      ll <- sum(lp)

    } else {

      # Take only the columns of the data that pertain to the child node and its
      # parents.
      parentData <- data.frame(y = data[, e],
                               data[, which(adjMatrix[, e] != 0)])

      # Run ordered logistic regression on the current genetic variant given its
      # parents.
      model <- polr(as.factor(y) ~ .,
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
