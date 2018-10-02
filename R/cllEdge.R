#' cllEdge
#'
#' Calculates the log likelihood of the current node when running a genetic
#' algorithm on the direction of the edges in the graph.
#'
#' @param coordinates A matrix of the row and column numbers of the nonzero
#' elements in the adjacency matrix. The row numbers make up the first row of
#' the coordinates matrix and the column numbers make up the second row.
#'
#' @param data A matrix with the nodes across the columns and the observations
#' along the rows.
#'
#' @param graph The edge directions of the current graph.
#'
#' @param m The length of the graph vector. This is the number of edges in the
#' graph plus the fitness of the graph.
#'
#' @param nGV The number of genetic variants in the graph.
#'
#' @param nNodes The number of nodes in the graph.
#'
#' @param pmr Logical. If true the Metropolis-Hastings algorithm will use the
#' Principle of Mendelian Randomization, PMR. This prevents the direction of an
#' edge pointing from a gene expression node to a genetic variant node.
#'
#' @param scoreFun A character string indicating what method to use for
#' calculating the fitness.
#'
#' @return The log likelihood of the graph according to the orientation of the
#' edges denoted by the DNA of the current graph.
#'
#' @importFrom stats dnorm
#'
#' @export
#'
cllEdge <- function (coordinates,
                     data,
                     graph,
                     m,
                     nGV,
                     nNodes,
                     pmr,
                     scoreFun) {

  # Convert the DNA of the current graph to an adjacency matrix
  adjMatrix <- toAdjMatrix(coordinates = coordinates,
                           graph = graph,
                           nNodes = nNodes)

  # If there are genetic variants in the data and if pmr is TRUE then the
  # adjacency matrix needs to be changed so there are no directed edges pointing
  # from a gene expression node to a genetic variant node. If there was an edge
  # from a ge node to a gv node it will be removed. This essentially increases
  # the prior probability of not seeing an edge between the two nodes because
  # when the edge directions are changed in the mutate function the direction
  # from a ge node to a gv node will be changed to the edge not being present.
  if (nGV != 0 & pmr == TRUE & any(adjMatrix[(nGV + 1):nNodes, 1:nGV] == 1)) {

    # Remove any edges where a ge node is the parent of a gv node.
    graph <- toED(adjMatrix = adjMatrix,
                  coordinates = coordinates,
                  graph = graph,
                  nGV = nGV,
                  nNodes = nNodes)

    # change the adjacency matrix to match the edge directions.
    adjMatrix[(nGV + 1):nNodes, 1:nGV] <- 0

  }

  # Calculate the log likelihood for the genetic variant nodes.
  cllM <- cllMultinom(adjMatrix = adjMatrix,
                      data = data,
                      nGV = nGV,
                      nNodes = nNodes)

  # Calculate the log likelihood for the gene expression nodes.
  cll <- cllNormal(adjMatrix = adjMatrix,
                   data = data,
                   k = cllM$nParents,
                   logLikelihood = cllM$logLikelihood,
                   nGV = nGV,
                   nNodes = nNodes)

  switch(scoreFun,

         'logLikelihood' = {

           # Add the log likelihood the the end of the graph vector.
           graph[[m]] <- sum(cll$logLikelihood)

           return (graph)

         },

         'BIC' = {

           # Calculate the sample size. This will be used in calculating BIC.
           n <- dim(data)[2]

           # Vector to hold the BIC for each node in the graph.
           bic <- c()

           # Calculate BIC for each node in the graph.
           for (e in 1:nNodes) {

             bic[[e]] <- (log(n) *
                            (cll$nParents[[e]] + 2) -
                            2 * cll$logLikelihood[[e]])

           }

           graph[[m]] <- sum(bic)

           return (graph)

         })

}
