#' mutate
#'
#' Takes the edge directions of the current graph mutates it according to the
#' mutation rate passed to the function.
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
#' @param prior A vector containing the prior probability of seeing each edge
#' direction.
#'
#' @return A vector of the edge directions and log likelihood.
#'
#' @importFrom stats rbinom
#'
#' @export
#'
mutate <- function (graph,
                    m,
                    mutationRate,
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
  for (e in 1:length(wEdges)) {

    if (newGraph[[wEdges[[e]]]] == 0) {

      newGraph[[wEdges[[e]]]] <- sample(x = c(1, 2),
                                             size = 1,
                                             prob = c(prior[[2]],
                                                      prior[[3]]))

    } else if (newGraph[[wEdges[[e]]]] == 1) {

      newGraph[[wEdges[[e]]]] <- sample(x = c(0, 2),
                                             size = 1,
                                             prob = c(prior[[1]],
                                                      prior[[3]]))

    } else {

      newGraph[[wEdges[[e]]]] <- sample(x = c(0, 1),
                                             size = 1,
                                             prob = c(prior[[1]],
                                                      prior[[2]]))

    }

  }

  return (newGraph)

}
