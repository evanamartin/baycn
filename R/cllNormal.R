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
#' node.
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
cllNormal <- function (adjMatrix,
                       data,
                       k,
                       logLikelihood,
                       nGV,
                       nNodes){

  # The log likelihood for the gene expression nodes will start from nGV+1 and
  # go to the value of nNodes.
  for (e in (nGV + 1):nNodes) {

    # Get the number of parents, estimates for the mean and standard
    # deviation of the current node, and the data of the parent nodes
    ep <- estimatesN(adjMatrix = adjMatrix,
                     data = data,
                     node = e)

    # Loop through each of the gene expression nodes and calculate the log
    # likelihood given the edge directions of the current graph.
    if (ep$nParents == 0) {

      ll <- sum(dnorm(x = data[, e],
                      mean = ep$estimates[[1]],
                      sd = ep$estimates[[2]],
                      log = TRUE))

    } else {

      ll <- sum(dnorm(x = data[, e],
                      mean = makeFormula(estimates = ep$estimates,
                                         nParents = ep$nParents,
                                         parentData = ep$parentData),
                      sd = ep$estimates[[length(ep$estimates)]],
                      log = TRUE))

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
