#' decimal
#'
#' Turns the directions and locations of the cycles in the graph in to a decimal
#' number.
#'
#' @param edgeDir A list containing the directions of the edges for each cycle
#' in the graph.
#'
#' @param edgeNum A list containing the edge index of each edge in each cycle
#' of the graph.
#'
#' @return A vector containing the decimal number of each cycle.
#'
decimal <- function(edgeDir,
                    edgeNum) {

  # Number of cycles in the graph
  nCycles <- length(edgeDir)

  # Store the edge directions and indices that form a cycle as a decimal number.
  dn <- vector(mode = 'integer',
               length = nCycles)

  # The length of edgeDir and edgeNum will always be the same. This length will
  # be used to loop through each list and turn the edge directions and numbers
  # of each cycle into a decimal number.
  for (e in 1:nCycles) {

    dn[[e]] <- sum(edgeDir[[e]] * 3^edgeNum[[e]]) + sum(edgeNum[[e]])

  }

  return (dn)

}
