#' mutate
#'
#' Takes the edge directions of the current graph mutates it according to the
#' mutation rate passed to the function.
#'
#' @param edgeType A 0 or 1 indicating whether the edge is a gv-ge edge (1) or
#' a gv-gv or ge-ge edge (1).
#'
#' @param graph A vector containing the DNA of the graph who was
#' selected by the selectBest function.
#'
#' @param m The length of the graph vector. This is the number of edges in the
#' graph plus the fitness of the graph.
#'
#' @param mutationRate A number between 0 and 1. As a rule of thumb the mutation
#' rate should be 1 / length of the DNA.
#'
#' @param pmr Logical. If true the Metropolis-Hastings algorithm will use the
#' Principle of Mendelian Randomization, PMR. This prevents the direction of an
#' edge pointing from a gene expression node to a genetic variant node.
#'
#' @param prior A vector containing the prior probability of seeing each edge
#' direction.
#'
#' @return A vector of the edge directions and log likelihood.
#'
#' @importFrom stats rbinom
#'
#' @export
#'
mutate <- function (edgeType,
                    graph,
                    m,
                    mutationRate,
                    pmr,
                    prior) {

  # Create a new graph. This graph will be a copy of the original and
  # it will be subject to mutation.
  newGraph <- graph

  # Generate a binomial random variable to determine the number of edges to
  # mutate.
  nEdges <- rbinom(n = 1,
                   size = (m - 1),
                   prob = mutationRate)

  # If there was no mutation return the original graph.
  if (nEdges == 0) {

    return (graph)

  }

  # Select which edges will be changed from their current state to a different
  # state.
  wEdges <- sample(x = 1:(m - 1),
                   size = nEdges,
                   replace = FALSE)

  # Loop through the edges to be mutated
  for (e in 1:nEdges) {

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
