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
#' @param nEdges The number of edges in the graph.
#'
#' @param nGV The number of genetic variants in the graph.
#'
#' @param nNodes The number of nodes in the graph.
#'
#' @return The log likelihood of the graph according to the orientation of the
#' edges denoted by the DNA of the current graph.
#'
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
