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
#' @importFrom stats dmultinom
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

    # Get the number of parents, estimates for the mean and standard
    # deviation of the current node, and the data of the parent nodes
    ep <- estimatesM(adjMatrix = adjMatrix,
                     data = data,
                     node = e)

    # Loop through each of the genetic variant nodes and calculate the log
    # likelihood given the edge directions of the current graph.
    if (ep$nParents == 0) {

      ll <- dmultinom(x = ep$estimates$counts,
                      prob = ep$estimates$probs,
                      log = TRUE)

    } else {

      ll <- 0

    }

    # Store the number of parents for each node. This will be used if the
    # scoreFun argument is set to 'BIC'.
    k[[e]] <- ep$nParents

    # Store the log likelihood for the current node.
    logLikelihood[[e]] <- ll

  }

  return (list(logLikelihood = logLikelihood,
               nParents = k))

}
