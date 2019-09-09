#' @importFrom stats dnorm
#'
cllEdge <- function (coordinates,
                     data,
                     graph,
                     m,
                     nGV,
                     nEdges,
                     nNodes) {

  # Convert the DNA of the current graph to an adjacency matrix
  adjMatrix <- toAdjMatrix(coordinates = coordinates,
                           graph = graph,
                           nEdges = nEdges,
                           nNodes = nNodes)

  # Calculate the log likelihood for the genetic variant nodes.
  cllM <- cllMultinom(adjMatrix = adjMatrix,
                      data = data,
                      nGV = nGV,
                      nNodes = nNodes)

  # Calculate the log likelihood for the gene expression nodes.
  cll <- cllNormal(adjMatrix = adjMatrix,
                   data = data,
                   logLikelihood = cllM,
                   nGV = nGV,
                   nNodes = nNodes)

  # Add the log likelihood to the end of the graph vector.
  graph[[m]] <- sum(cll)

  return (graph)

}

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

      # Run ordered logistic regression on the current genetic variant given its
      # parents and calculate the log likelihood.
      logLikelihood[[e]] <- logLik(polr(as.factor(data[, e]) ~ data[, which(adjMatrix[, e] != 0)]))[1]

    }

  }

  return (logLikelihood)

}

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
