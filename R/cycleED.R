#' cycleED
#'
#' Takes the names of the nodes that could form a cycle and creates a vector of
#' edge directions that form a directed cycle.
#'
#' @param tBranches A list of trimmed branches. The branches in this list
#' contain only the nodes that could form a cycle.
#'
#' @return A list of edge directions that forms a directed cycle for each
#' trimmed branch passed to the cycleED function.
#'
#' @export
#'
cycleED <- function (tBranches) {

  # Get the Number of cycles created by the cycleTB function
  nCycle <- length(tBranches)

  # The edgeDir list will hold the vectors of the edge directions that will form
  # a cycle in the graph.
  edgeDir <- list()

  # Loop through all the cycles
  for (e in 1:nCycle) {

    edgeDir[[e]] <- list()

    # Loop through each element in the cycle to determine which direction the
    # edge should be facing to create a cycle.
    for (v in 1:(length(tBranches[[e]]) - 1)) {

      # If the value of the first node is less than the value of the second node
      # the direction of the edge points from the node with a smaller index to
      # the node with a larger index.
      if (as.numeric(tBranches[[e]][[v]]) < as.numeric(tBranches[[e]][[v+1]])) {

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

  return (edgeDir)

}
