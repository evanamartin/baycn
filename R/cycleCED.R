#' cycleCED
#'
#' Extracts the coordinates of the adjacency matrix for each edge that is part
#' of a directed cycle and takes the names of the nodes in each cycle and
#' creates a vector of edge directions that form a directed cycle.
#'
#' @param tBranches A list of cycles. Each element in the list is a vector of
#' nodes that could form a directed cycle.
#'
#' @return A list with the coordinates of the adjacency matrix for each edge in
#' each cycle as the first element and the direction of the edges that form a
#' directed cycle as the second element.
#'
#' @export
#'
cycleCED <- function (tBranches) {

  # Create a list that will have a list for each edge in each directed cycle.
  # The first level of lists are for each cycle and the sublists are for each
  # edge in the cycle.
  cCoord <- list()

  # The edgeDir list will hold the vectors of the edge directions that will form
  # a cycle in the graph.
  edgeDir <- list()

  # This for loop loops through each trimmed branch or cycle present in the
  # graph.
  for (e in 1:length(tBranches)) {

    # A sublist that holds the coordinates of the directed edge for all the
    # edges in the eth cycle.
    cCoord[[e]] <- list()

    # A vector to hold the directions of each edge for the directed cycle.
    edgeDir[[e]] <- list()

    # This for loop gets the row and column indices for each edge in the cycle.
    for (v in 1:(length(tBranches[[e]]) - 1)) {

      # The row index is for the parent node and the column index is for the
      # child node. These indices will be used to determine if a 0 or 1 should
      # be present in the adjacency matrix in order for the nodes to form a
      # directed cycle.
      #                        Row index            Column index
      cCoord[[e]][[v]] <- c(tBranches[[e]][[v]], tBranches[[e]][[v + 1]])

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

    edgeDir[[e]] <- unlist(edgeDir[[e]])

  }

  return (list(cCoord = cCoord,
               edgeDir = edgeDir))

}
