#' cycleCED
#'
#' Extracts the coordinates of the adjacency matrix for each edge that is part
#' of a directed cycle and takes the names of the nodes in each cycle and
#' creates a vector of edge directions that form a directed cycle.
#'
#' @param nBranches An integer. The number of branches in the tree.
#'
#' @param tBranches A list of cycles. Each element in the list is a vector of
#' nodes that could form a directed cycle.
#'
#' @return A list with the coordinates of the adjacency matrix for each edge in
#' each cycle as the first element in the list, a vector of columns visited by
#' the current cycle as the second element, and the direction of the edges that
#' form a directed cycle as the third element.
#'
cycleCED <- function (nBranches,
                      tBranches) {

  # Create a list that will have a list for each edge in each directed cycle.
  # The first level of lists are for each cycle and the sublists are for each
  # edge in the cycle.
  cCoord <- vector(mode = 'list',
                   length = nBranches)

  # A list that will hold the column indices of each column visited in the SaL
  # matrix by the current cycle.
  Columns <- vector(mode = 'list',
                    length = nBranches)

  # The edgeDir list will hold the vectors of the edge directions that will form
  # a cycle in the graph.
  edgeDir <- vector(mode = 'list',
                    length = nBranches)

  # Loop through each trimmed branch or cycle present in the
  # graph.
  for (e in 1:nBranches) {

    # Get the nubmer of edges in the current cycle
    nEdges <- length(tBranches[[e]]) - 1

    # A sublist that holds the columns indices for each column visited in the
    # current cycle.
    Columns[[e]] <- vector(mode = 'integer',
                           length = nEdges)

    # A vector to hold the directions of each edge for the directed cycle.
    edgeDir[[e]] <- numeric(nEdges)

    # This for loop gets the row and column indices for each edge in the cycle.
    for (v in 1:nEdges) {

      # The row index is for the parent node and the column index is for the
      # child node. These indices will be used to determine if a 0 or 1 should
      # be present in the adjacency matrix in order for the nodes to form a
      # directed cycle.
      #                        Row index            Column index
      cCoord[[e]][[v]] <- c(tBranches[[e]][[v]], tBranches[[e]][[v + 1]])

      # Keep the columns visited by each edge in the current cycle. This will be
      # used to determine if there are disjoint cycles in the graph.
      Columns[[e]][[v]] <- as.numeric(tBranches[[e]][[v + 1]])

      # If the value of the first node is less than the value of the second node
      # the direction of the edge points from the node with a smaller index to
      # the node with a larger index.
      if (as.numeric(tBranches[[e]][[v]]) <
          as.numeric(tBranches[[e]][[v + 1]])) {

        edgeDir[[e]][[v]] <- 0

        # If the value of the first node is greater than the value of the second
        # node the direction of the edge points from the node with a larger
        # index to the node with a smaller index.
      } else {

        edgeDir[[e]][[v]] <- 1

      }

    }

  }

  # Get the unique column numbers for each column visited
  uColumns <- unique(unlist(Columns))

  return (list(cCoord = cCoord,
               columns = uColumns,
               edgeDir = edgeDir))

}
