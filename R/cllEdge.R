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
#' @param individual The edge directions of the current graph.
#'
#' @param nGV The number of genetic variants in the graph.
#'
#' @param nNodes The number of nodes in the graph.
#'
#' @param scoreFun A character string indicating what method to use for
#' calculating the fitness.
#'
#' @return The log likelihood of the graph according to the orientation of the
#' edges denoted by the DNA of the current individual.
#'
#' @importFrom stats dnorm
#'
#' @export
#'
cllEdge <- function (coordinates,
                     data,
                     individual,
                     nGV,
                     nNodes,
                     scoreFun) {

  # Convert the DNA of the current individual to an adjacency matrix
  adjMatrix <- toAdjMatrix(coordinates = coordinates,
                           individual = individual,
                           nNodes = ncol(data))

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

           return (sum(cll$logLikelihood))

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

           return (sum(bic))

         })

}
