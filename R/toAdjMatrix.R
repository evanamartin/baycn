#' toAdjMatrix
#'
#' Creates an adjacency matrix from the DNA of the current individual. The DNA
#' is converted to an adjacency matrix to determine the parents of each node
#' in order to calculate its log likelihood.
#'
#' @param individual A vector containing the directions of the edges in the
#' graph and the log likelihood of the individual.
#'
#' @param coordinates A matrix with the row and column coordinates of the
#' edges from the adjacency matrix from the MRPC function. The row numbers of
#' each nonzero element make up the first row in the matrix and the column
#' numbers of each nonzero element make up the second row of the matrix.
#'
#' @param nNodes The number of nodes in the graph.
#'
#' @return An adjacency matrix.
#'
#' @export
#'
toAdjMatrix <- function (individual,
                         coordinates,
                         nNodes) {

  adjMatrix <- matrix(rep(0, nNodes^2),
                      nrow = nNodes,
                      ncol = nNodes)

  # create the adjacency matrix with the edges oriented according to the DNA.
  # This will be used in determining which function to use when calculating the
  # log likelihood for each node.
  for(e in 1:(length(individual) - 1)) {

    if (individual[[e]] == 0) {

      adjMatrix[coordinates[1, e], coordinates[2, e]] <- 1

    } else if (individual[[e]] == 1) {

      adjMatrix[coordinates[2, e], coordinates[1, e]] <- 1

    }

  }

  return(adjMatrix)

}
