toAdjMatrix <- function (coordinates,
                         graph,
                         nEdges,
                         nNodes) {

  # Create a pxp matrix of 0s to be converted to an adjacency matrix.
  adjMatrix <- matrix(0,
                      nrow = nNodes,
                      ncol = nNodes)

  # create the adjacency matrix with the edges oriented according to the DNA.
  # This will be used in determining which function to use when calculating the
  # log likelihood for each node.
  for(e in 1:nEdges) {

    if (graph[[e]] == 0) {

      adjMatrix[coordinates[1, e], coordinates[2, e]] <- 1

    } else if (graph[[e]] == 1) {

      adjMatrix[coordinates[2, e], coordinates[1, e]] <- 1

    }

  }

  return(adjMatrix)

}
