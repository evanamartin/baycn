#' cycleED
#'
#' Takes the names of the nodes that could form a cycle and creates a vector of
#' DNA that forms a cycle between the nodes.
#'
#' @param tBranches A list of trimmed branches. The branches in this list
#' contain only the nodes that could form a cycle.
#'
#' @return A list of DNA that creates a cycle.
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
      # the direction of the edge points from the first to the second.
      if (as.numeric(tBranches[[e]][[v]]) < as.numeric(tBranches[[e]][[v+1]])) {

        edgeDir[[e]][[v]] <- 0

        # If the value of the first node is greater than the value of the second
        # node the direction of the edge points from the first to the second.
      } else {

        edgeDir[[e]][[v]] <- 1

      }

    }

    edgeDir[[e]] <- unlist(edgeDir[[e]])

  }

  return (edgeDir)

}
