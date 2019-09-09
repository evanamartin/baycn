detType <- function (coord,
                     nEdges,
                     nGV,
                     pmr) {

  if (pmr == TRUE) {

    # Vector to hold the edge type.
    forbidden <- vector(mode = 'integer',
                        length = nEdges)

    # Loop through the coordinates matrix and determine the edge type based on
    # where the coordinates fall in the adjacency matrix.
    for (e in 1:nEdges) {

      # Look for edges in the lower left matrix. Where the adjacency matrix is
      # divided into four submatrices by the number of genetic variants present.
      if (coord[1, e] <= nGV & coord[2, e] > nGV) {

        # The edge is between a gv and ge
        forbidden[[e]] <- 1

      } else {

        # Any other type of edge (gv-gv or ge-ge)
        forbidden[[e]] <- 0

      }

    }

  } else {

    # If pmr is set to false all edges will be edge type 0 because we don't care
    # to distinguish between the different edge types.
    forbidden <- c(rep(0, nEdges))

  }

  return (forbidden)

}
