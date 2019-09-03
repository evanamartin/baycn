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

  Rows <- c()
  Columns <- c()

  i <- 1

  for (e in 1:nrow(adjMatrix)) {

    for (v in 1:ncol(adjMatrix)) {

      if (e < v && (adjMatrix[e, v] != 0 || adjMatrix[v, e] != 0)) {

        # Get the rows and columns of the nonzero entries of the upper
        # triangular adjacency matrix.
        Rows[[i]] <- e
        Columns[[i]] <- v

        i <- i + 1

      }

    }

  }

  return (rbind(Rows, Columns))

}
