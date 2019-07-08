#' coordinates
#'
#' Takes the adjacency matrix and extracts the row and column numbers for each
#' nonzero element. It also takes the nonzero elements and collapses them into a
#' vector. This vector is the DNA of the individual. The length of the DNA is
#' the number of edges in the graph.
#'
#' @param adjMatrix The adjacency matrix returned by the MRPC algorithm. The
#' adjacency matrix is a matrix of zeros and ones. The ones represent an edge
#' and also indicates the direction of that edge.
#'
#' @return A matrix. The first row in the matrix holds the row locations of the
#' nonzero elements in the adjacency matrix and the second row holds the column
#' locations of the nonzero elements in the adjacency matrix.
#'
coordinates <- function (adjMatrix) {

  # add the adjacency matrix to its transpose. I do this to make sure all of the
  # edges in the graph are in the upper triangle of the adjacency matrix.
  AtA <- adjMatrix + t(adjMatrix)

  # I need the location of the edges in the adjacency matrix. I will use this
  # information in calculating the log likelihood.
  # Use only the upper matrix to extract the positions of the edges in the
  # graph and to create the DNA. The diagonal and the lower triangle will be
  # zeros.
  AtA[lower.tri(AtA)] <- 0

  # Get the rows and columns of the nonzero entries of the upper AtA matrix
  Rows <- c()
  Columns <- c()

  i <- 1

  for (e in 1:nrow(AtA)) {

    for (v in 1:ncol(AtA)) {

      if (AtA[e, v] != 0) {

        Rows[[i]] <- e

        Columns[[i]] <- v

        i <- i + 1

      }

    }

  }

  # Save the coordinates of the nonzero cells of the upper AtA matrix to use in
  # determining which equation to use when calculating the log likelihood.
  coordinates <- rbind(Rows, Columns)

  return (coordinates)

}
