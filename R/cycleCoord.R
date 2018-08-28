#' cycleCoord
#'
#' Gets the coordinates in the adjacency matrix of the nodes that could form a
#' cycle.
#'
#' @param tBranches A list of the nodes that could form a cycle.
#'
#' @return A list of row column coordinates for each node that could form a
#' cycle. The coordinates are to the adjacency matrix.
#'
#' @export
#'
cycleCoord <- function (tBranches) {

  # From the tBranches list create the pairs of nodes between each edge for each
  # list in the tBranches list
  cCoord <- list()

  for (e in 1:length(tBranches)) {

    cCoord[[e]] <- list()

    for (v in 1:(length(tBranches[[e]]) - 1)) {

      cCoord[[e]][[v]] <- c(tBranches[[e]][[v]], tBranches[[e]][[v + 1]])

    }

  }

  return (cCoord)

}
