#' rmCycle
#'
#' Changes the edge direction of the current graph until there are no directed
#' cylces.
#'
#' @param dnUnique A vector of unique decimal numbers for each potential cycle
#' in the graph.
#'
#' @param edgeNum A list containing the numbers for each of the edges in each of
#' the potential cycles in the graph.
#'
#' @param edgeType A 0 or 1 indicating whether the edge is a gv-ge edge (1) or
#' a gv-gv or ge-ge edge (1).
#'
#' @param individual A vector of the edge directions, and the log likelihood
#' when the function is called within the populate function, of an individual.
#'
#' @param pmr Logical. If true the Metropolis-Hastings algorithm will use the
#' Principle of Mendelian Randomization, PMR. This prevents the direction of an
#' edge pointing from a gene expression node to a genetic variant node.
#'
#' @param prior A vector containing the prior probability of seeing each edge
#' direction.
#'
#' @param wCycle A vector of the locations in the cycleDN vector where each
#' unique number occurs. The unique numbers refer to each of the unique cycles
#' in the graph. This argument is NULL if there are no cycles in the graph.
#'
#' @return A vector of the DNA of the individual with the cycles removed.
#'
#' @export
#'
rmCycle <- function (dnUnique,
                     edgeNum,
                     edgeType,
                     individual,
                     pmr,
                     prior,
                     wCycle) {

  # Get the decimal numbers for each potential cycle in the current individual.
  dnCI <- decimalCI(edgeNum,
                    wCycle,
                    individual)

  # If there is a cycle change the direction of an edge invovled in the cycle.
  while (any(dnCI == dnUnique)) {

    # Get the edges that create a cycle in the current individual.
    whichCycle <- which(dnCI == dnUnique)

    for (e in 1:length(whichCycle)) {

      # Get the edges involved in the eth cycle
      curEdges <- edgeNum[[whichCycle[[e]]]]

      # Choose one of these edges and change it.
      whichEdge <- sample(x = curEdges,
                          size = 1)

      if(individual[[whichEdge]] == 0) {

        # Calculate the probability of moving to edge state 1 or 2 (location 2
        # or 3 in the prior vector)
        probability <- cPrior(edges = c(2, 3),
                              edgeType = edgeType[[whichEdge]],
                              pmr = pmr,
                              prior = prior)

        individual[[whichEdge]] <- sample(x = c(1, 2),
                                          size = 1,
                                          prob = probability)

      } else {

        # Calculate the probability of moving to edge state 0 or 2 (location 1
        # or 3 in the prior vector)
        probability <- cPrior(edges = c(1, 3),
                              edgeType = edgeType[[whichEdge]],
                              pmr = pmr,
                              prior = prior)

        individual[[whichEdge]] <- sample(x = c(0, 2),
                                          size = 1,
                                          prob = probability)

      }

    }

    # Get the decimal numbers for each potential cycle in the current individual
    # after an edge in each cycle has been changed.
    dnCI <- decimalCI(edgeNum,
                      wCycle,
                      individual)

  }

  return (individual)

}
