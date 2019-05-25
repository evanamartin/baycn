#' mutate
#'
#' Takes the current graph and changes at least one of the edges from its
#' current edge state to a different edge state.
#'
#' @param edgeType A 0 or 1 indicating whether the edge is a gv-ge edge (1) or
#' a gv-gv or ge-ge edge (0).
#'
#' @param graph A vector containing the current edge states for each edge.
#'
#' @param nEdges The number of edges in the graph.
#'
#' @param pmr Logical. If true the Metropolis-Hastings algorithm will use the
#' Principle of Mendelian Randomization, PMR. This prevents the direction of an
#' edge pointing from a gene expression node to a genetic variant node.
#'
#' @param prior A vector containing the prior probability of seeing each edge
#' direction.
#'
#' @param ztbProb A vector of probabilities for 1:nEdges from a zero truncated
#' binomial distribution.
#'
#' @return A vector of the edge directions and log likelihood.
#'
#' @importFrom stats rbinom
#'
#' @export
#'
mutate <- function (edgeType,
                    graph,
                    nEdges,
                    pmr,
                    prior,
                    ztbProb) {

  # Create a new graph. This graph will be a copy of the original and
  # it will be subject to mutation.
  newGraph <- graph

  # Generate a zero truncated binomial random variable to determine the number
  # of edges to change states.
  nChange <- sample(x = 1:nEdges,
                    size = 1,
                    prob = ztbProb)

  # Select which edges will be changed from their current state to a different
  # state.
  wEdges <- sample(x = 1:nEdges,
                   size = nChange,
                   replace = FALSE)

  # Loop through the edges to be mutated
  for (e in 1:nChange) {

    if (newGraph[[wEdges[[e]]]] == 0) {

      # Calculate the probability of moving to edge state 1 or 2 (location 2 or
      # 3 in the prior vector)
      probability <- cPrior(edges = c(2, 3),
                            edgeType = edgeType[[wEdges[[e]]]],
                            pmr = pmr,
                            prior = prior)

      newGraph[[wEdges[[e]]]] <- sample(x = c(1, 2),
                                        size = 1,
                                        prob = probability)

    } else if (newGraph[[wEdges[[e]]]] == 1) {

      # Calculate the probability of moving to edge state 0 or 2 (location 1 or
      # 3 in the prior vector)
      probability <- cPrior(edges = c(1, 3),
                            edgeType = edgeType[[wEdges[[e]]]],
                            pmr = pmr,
                            prior = prior)

      newGraph[[wEdges[[e]]]] <- sample(x = c(0, 2),
                                        size = 1,
                                        prob = probability)

    } else {

      # Calculate the probability of moving to edge state 0 or 1 (location 1 or
      # 2 in the prior vector)
      probability <- cPrior(edges = c(1, 2),
                            edgeType = edgeType[[wEdges[[e]]]],
                            pmr = pmr,
                            prior = prior)

      newGraph[[wEdges[[e]]]] <- sample(x = c(0, 1),
                                        size = 1,
                                        prob = probability)

    }

  }

  return (newGraph)

}
