#' cycleCoord
#'
#' Extracts the coordinates from the adjacency matrix of each edge that is part
#' of a directed cycle.
#'
#' @param tBranches A list cycles. Each element in the list is a vector of nodes
#' that could form a directed cycle.
#'
#' @return A list of row and column coordinates for each edge in a directed
#' cycle. The coordinates are to the adjacency matrix.
#'
#' @export
#'
cycleCoord <- function (tBranches) {

  # Create a list that will have a list for each edge in each directed cycle.
  # The first level of lists are for each cycle and the sublists are for each
  # edge in the cycle.
  cCoord <- list()

  # This for loop loops through each trimmed branch, or cycle, present in the
  # graph.
  for (e in 1:length(tBranches)) {

    # A sublist that holds the coordinates of the directed edge for all the
    # edges in the eth cycle.
    cCoord[[e]] <- list()

    # This for loop gets the row and column indices for each edge in the cycle.
    for (v in 1:(length(tBranches[[e]]) - 1)) {

      # The row index is for the parent node and the column index is for the
      # child node. These indices will be used to determine if a 0 or 1 should
      # be present in the adjacency matrix in order for the nodes to form a
      # directed cycle.
      #                        Row index            Column index
      cCoord[[e]][[v]] <- c(tBranches[[e]][[v]], tBranches[[e]][[v + 1]])

    }

  }

  return (cCoord)

}
