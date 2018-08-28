#' decimalCI
#'
#' Calculates the decimal number of the cycles of an individual according to its
#' current edge directions.
#'
#' @param edgeNum A list containing the number of each edge for each cycle in
#' the graph.
#'
#' @param wCycle A vector of the locations in the cycleDN vector where each
#' unique number occurs. The unique numbers refer to each of the unique cycles
#' in the graph. This argument is NULL if there are no cycles in the graph.
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
                       wCycle,
                       individual) {

  # Get the directions of the edges between the nodes that could create a cycle.
  # These directions along with the edge numbers will be used to calculate a
  # decimal number for each potential cycle given the directions of the current
  # individual.
  curDir <- list()
  for (e in 1:length(wCycle)) {

    # Directions of the edges that could form a cycle from the current
    # individual.
    curDir[[e]] <- individual[edgeNum[[wCycle[[e]]]]]

  }

  # Calculate the decimal number for each potential cycle for the current
  # individual. This will be used to determine if there is a cycle and which
  # edges are involved in the cycle.
  dnCI <- decimal(edgeDir = curDir,
                  edgeNum = edgeNum[wCycle])

  return (dnCI)

}
