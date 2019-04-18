#' decimalCI
#'
#' Calculates the decimal number of the cycles of an individual according to its
#' current edge directions.
#'
#' @param edgeNum A list containing the index of each edge for each cycle in
#' the graph.
#'
#' @param individual A vector containing the DNA and the log likelihood of the
#' current individual.
#'
#' @return The decimal numbers for each potential cycle given the edge
#' directions of the current individual.
#'
#' @export
#'
decimalCI <- function (edgeNum,
                       individual) {

  # Get the directions of the edges between the nodes that could create a cycle.
  # These directions along with the edge numbers will be used to calculate a
  # decimal number for each potential cycle given the directions of the current
  # individual.
  curDir <- vector(mode = 'list',
                   length =  length(edgeNum))

  for (e in 1:length(edgeNum)) {

    # Directions of the edges that could form a cycle from the current
    # individual.
    curDir[[e]] <- individual[edgeNum[[e]]]

  }

  # Calculate the decimal number for each potential cycle for the current
  # individual. This will be used to determine if there is a cycle and which
  # edges are involved in the cycle.
  dnCI <- decimal(edgeDir = curDir,
                  edgeNum = edgeNum)

  return (dnCI)

}
