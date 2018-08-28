#' cycleEN
#'
#' Gets the number of each edge in the DNA that could form a cycle.
#'
#' @param adjMatrix The adjacency matrix of the graph.
#'
#' @param cCoord A list containing the coordinates in the adjacency matrix of
#' the nodes that could form a cycle
#'
#' @return A list containing the edge numbers for the edges that could form a
#' cycle for each cycle in the graph.
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

  # Stores the location of each edge for each cycle in the tBranches list
  edgeNum <- list()

  # Create a for loop that loops through the position matrix for each element
  # in the cycleCoor list. This loop will be within another loop that will
  # create a new matrix of DNA position for each edge in each cycle.
  # the for loop starts across the columns then moves down the rows.
  for (e in 1:length(cCoord)) {

    edgeNum[[e]] <- list()

    for (v in 1:length(cCoord[[e]])) {

      for (a in 1:dim(position)[2]) {

        if (setequal(position[, a], cCoord[[e]][[v]])) {

          edgeNum[[e]][[v]] <- a

        }

      }

    }

    edgeNum[[e]] <- unlist(edgeNum[[e]])

  }

  return (edgeNum)

}
