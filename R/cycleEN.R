#' cycleEN
#'
#' Extracts the edge number of each edge in the graph.
#'
#' @param adjMatrix The adjacency matrix of the graph.
#'
#' @param cCoord A list containing the coordinates in the adjacency matrix of
#' the nodes that could form a cycle.
#'
#' @return A list containing the edge numbers for the edges that could form a
#' directed cycle for each cycle in the graph.
#'
#' @export
#'
cycleEN <- function (adjMatrix,
                     cCoord) {

  # This line returns the coordinates of the adjacency matrix in a matrix where
  # the column numbers correspond to the position in the DNA of each edge. This
  # will be used to subset the individual to remove the possibility of a cycle
  # forming.
  position <- coordinates(adjMatrix)

  # Number of directed cycles in the graph.
  nCycles <- length(cCoord)

  # Stores the location of each edge for each cycle in the tBranches list
  edgeNum <- vector(mode = 'list',
                    length = nCycles)

  # The first for loop loops through each of the cycles present in the graph.
  for (e in 1:nCycles) {

    nEdges <- length(cCoord[[e]])

    # Create a list to store the edge number for each edge in the cycle.
    edgeNum[[e]] <- vector(mode = 'numeric',
                           length = nEdges)

    # This for loop is the length of the current cycle, the number of edges in
    # the cycle.
    for (v in 1:nEdges) {

      # This loops through all of the coordinates for each edge in the graph and
      # compares the coordinates of the current edge from the cycle to the
      # coordinates of each edge in the graph until it finds a match. When it
      # finds a match it stores the column number in the coordinates matrix
      # which is the edge number of the current edge in the graph.
      for (a in 1:dim(position)[2]) {

        if (setequal(position[, a], cCoord[[e]][[v]])) {

          edgeNum[[e]][[v]] <- a

          # Stop checking for the edge number once it is found.
          break

        }

      }

    }

  }

return (edgeNum)

}
