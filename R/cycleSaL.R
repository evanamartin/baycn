#' cycleSaL
#'
#' The first step in identifying potential directed cycles in the graph is to
#' find the nodes that could form a cycle. Since a node must be connected to at
#' least two other nodes to be a part of a cycle we start by summing all of the
#' rows of the undirected adjacency matrix and creating a new matrix that only
#' has the rows of the undirected adjacency matrix with a sum of two or higher.
#' The function cycleSaL creates the matrix of nodes that could form a directed
#' cycle.
#'
#' Creates a stem and leaf matrix from the adjacency matrix with the nodes that
#' could potentially form a cycle. The row names of the matrix created by the
#' cycleSaL function are the stems and the 1s in the matrix are the leaves.
#'
#' @param adjMatrix The adjacency matrix of the graph.
#'
#' @return  A binary matrix where the 1s represent an edge between the row node
#' and the column node. If there are no cycles in the graph the function returns
#' NULL.
#'
cycleSaL <- function (adjMatrix) {

  # Add the adjacency matrix to its transpose to make all edges undirected.
  undirected <- adjMatrix + t(adjMatrix)

  # If a graph with undirected edges is passed to the removeCycles function then
  # the row sums could be 2 or greater even if that node is not in a potential
  # cycle because the cell in the undirected matrix for those nodes will be a 2.
  # To fix this problem we change any 2 in the undirected matrix to a 1.
  undirected[undirected == 2]  <- 1

  # Sum each row across all columns. The nodes that can potentially form a cycle
  # will have a sum greater than 2.
  sums <- apply(X = undirected,
                MARGIN = 1,
                FUN = sum)

  # Create a matrix whose rows are the nodes that could potentially form a
  # cycle.
  SaL <- undirected[which(sums >= 2), ]

  # If there is only one row in the SaL matrix then it is turned into a vector
  # and we cannot give row names to a vector. The following if statement checks
  # if SaL is a vector or if it is a matrix with a dimension of two. If the
  # number of rows is 2 or fewer there cannot be a cycle in the graph.
  if (is.null(dim(SaL))) {

    return (NULL)

    # I had to break up the if statement because if SaL is a vector then the
    # statement in else if would return logical(0)
  } else if (dim(SaL)[1] <= 2) {

    return(NULL)

  }

  # Name the rows according to their placement in the undirected matrix. The
  # name of each row is the row number of the undirected adjacency matrix where
  # there was a sum of 2 or greater.
  rownames(SaL) <- which(sums >= 2)

  # The column names are also the column numbers of the undirected adjacency
  # matrix. At this point in the function we keep all of the columns in the
  # undirected adjacency matrix.
  colnames(SaL) <- c(1:dim(undirected)[2])

  # Get the number of rows and columns in the matrix of potential cycles. These
  # will be used later to determine which entries in the cycles matrix need to
  # remain a one and which need to be changed to zero. The cells that will be
  # changed to a zero are the column numbers that don't also appear in the rows.
  nRows <- dim(SaL)[1]
  nCols <- dim(SaL)[2]

  # Go through each element in the cycles matrix and turn the ones to zeros if
  # the column number the one occurs in does not appear in the rownames vector.
  for (e in 1:nRows) {

    for (v in 1:nCols) {

      if (SaL[e, v] == 1) {

        # If the column number doesn't show up in the row names change the one
        # to a zero in the cycles matrix.
        if (!any(as.numeric(rownames(SaL)) == v)) {

          SaL[e, v] <- 0

        }

      }

    }

  }

  # get the new row sums after the nodes that can't make up a cycle are removed
  # from the cycles matrix.
  newSums <- apply(X = SaL,
                   MARGIN = 1,
                   FUN = sum)

  # The following while lookp will continue to trim down the SaL matrix when
  # there is a string of nodes that are connected to at least two nodes but none
  # of the nodes form a cycle.
  while (any(newSums <= 1)) {

    # Only keep the rows with a sum >= 2
    SaL <- SaL[which(newSums >= 2), ]

    # Get the new row dimensions
    nRows <- dim(SaL)[1]

    # If there is only one row left then the dim function will return NULL
    if (is.null(nRows)) {

      return (NULL)

      # If there are 2 or fewer rows with a row sum of 2 or greater there are no
      # cycles in the graph and cycleSaL returns a zero.
    } else if (nRows <= 2) {

      return (NULL)

    } else {

      # Go through each element in the cycles matrix and turn the ones to zeros
      # if the column number the one occurs in does not appear in the rownames
      # vector.
      for (e in 1:nRows) {

        for (v in 1:nCols) {

          if (SaL[e, v] == 1) {

            # If the column number doesn't show up in the row names change the
            # one to a zero in the cycles matrix.
            if (!any(as.numeric(rownames(SaL)) == v)) {

              SaL[e, v] <- 0

            }

          }

        }

      }

      # This will be used to check if there are any rows that have a sum that is
      # less than 2. If there are then the while loop will continue. If not then
      # the while loop finishes.
      newSums <- apply(X = SaL,
                       MARGIN = 1,
                       FUN = sum)

    }

  }

  # Remove any columns that do not contain any 1s
  SaL <- SaL[, apply(X = SaL,
                     MARGIN = 2,
                     FUN = sum) != 0]

  return (SaL)

}
