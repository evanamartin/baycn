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
