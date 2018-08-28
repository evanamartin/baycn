#' mutate
#'
#' Takes the DNA of the cloned winner and mutates it according to the mutation
#' rate passed to the function.
#'
#' @param coordinates A matrix of the row and column numbers of the nonzero
#' elements in the adjacency matrix. The row numbers make up the first row of
#' the coordinates matrix and the column numbers make up the second row.
#'
#' @param data A matrix with the nodes across the columns and the observations
#' along the rows.
#'
#' @param dnUnique A vector of unique decimal numbers. These numbers represent
#' each cycle in the graph. This argument is NULL if there are no cycles in the
#' graph.
#'
#' @param edgeNum A list containing the number of each edge in each of the
#' cycles in the graph. This argument is NULL if there are no cycles in the
#' graph.
#'
#' @param individual A vector containing the DNA of the individual who was
#' selected by the selectBest function.
#'
#' @param m The length of the individual vector. This is the number of edges in
#' the graph and the fitness of the graph.
#'
#' @param mutationRate A number between 0 and 1. As a rule of thumb the mutation
#' rate should be 1 / length of the DNA.
#'
#' @param prior A vector containing the prior probability of seeing each edge
#' direction.
#'
#' @param scoreFun A character string indicating what method to use for
#' calculating the fitness.
#'
#' @param wCycle A vector of the locations in the cycleDN vector where each
#' unique number occurs. The unique numbers refer to each of the unique cycles
#' in the graph. This argument is NULL if there are no cycles in the graph.
#'
#' @return A vector of the DNA and the log likelihood.
#'
#' @importFrom stats rbinom
#'
#' @export
#'
mutate <- function (coordinates,
                    data,
                    dnUnique,
                    edgeNum,
                    individual,
                    m,
                    mutationRate,
                    prior,
                    scoreFun,
                    wCycle) {

  # Create a new individual. This individual will be a copy of the original and
  # it will be subject to mutation.
  newIndividual <- individual

  # Generate a binomial random variable to determine the number of edges to
  # mutate.
  nEdges <- rbinom(n = 1,
                   size = (m - 1),
                   prob = mutationRate)

  # If there was no mutation return the original individual.
  if (nEdges == 0) {

    return (individual)

  }

  # Select which edges, if any, will be changed from their current state to a
  # different state.
  wEdges <- sample(x = 1:(m - 1),
                   size = nEdges,
                   replace = FALSE)

  # Loop through the edges to be mutated
  for (e in 1:length(wEdges)) {

    if (newIndividual[[wEdges[[e]]]] == 0) {

      newIndividual[[wEdges[[e]]]] <- sample(x = c(1, 2),
                                             size = 1,
                                             prob = c(prior[[2]],
                                                      prior[[3]]))

    } else if (newIndividual[[wEdges[[e]]]] == 1) {

      newIndividual[[wEdges[[e]]]] <- sample(x = c(0, 2),
                                             size = 1,
                                             prob = c(prior[[1]],
                                                      prior[[3]]))

    } else {

      newIndividual[[wEdges[[e]]]] <- sample(x = c(0, 1),
                                             size = 1,
                                             prob = c(prior[[1]],
                                                      prior[[2]]))

    }

  }

  # If there are potential cycles in the graph check if the mutated individual
  # has a cycle in it.
  if(!is.null(dnUnique)) {

    # Check for cycles and remove them if they appear in the individual.
    newIndividual <- rmCycle(dnUnique = dnUnique,
                             edgeNum = edgeNum,
                             prior = prior,
                             individual = newIndividual,
                             wCycle = wCycle)

    # It is possible that after removing any directed cycles in the graph the
    # new individual is the same as the original individual. If that is the case
    # return the original individual.
    if (identical(newIndividual[1:(m - 1)], individual[1:(m - 1)])) {

      return (individual)

    }

    # calculate the fitness then return the new individual if the likelihood
    # is higher. If the likelihood is lower return the new individual with
    # probability of new likelihood / old likelihood.
    newIndividual[[m]] <- cllEdge(individual = newIndividual,
                                  coordinates = coordinates,
                                  data = data,
                                  scoreFun = scoreFun)

    # Determine whether to keep the original individual, individual, or the
    # new individual, newIndividual. The original individual is the individual
    # that was passed to the mutate function and the new individual is the
    # individual that was subject to mutation.
    oldOrNew <- caRatio(m = m,
                        newIndividual = newIndividual,
                        oldIndividual = individual,
                        prior = prior,
                        scoreFun = scoreFun)

    return (oldOrNew)

  } else {

    # calculate the fitness then return the new individual if the likelihood is
    # higher. If the likelihood is lower return the new individual with
    # probability of new likelihood / old likelihood.
    newIndividual[[m]] <- cllEdge(individual = newIndividual,
                                  coordinates = coordinates,
                                  data = data,
                                  scoreFun = scoreFun)

    # Determine whether to keep the original individual, individual, or the new
    # individual, newIndividual. The original individual is the individual that
    # was passed to the mutate function and the new individual is the individual
    # that was subject to mutation.
    oldOrNew <- caRatio(m = m,
                        newIndividual = newIndividual,
                        oldIndividual = individual,
                        prior = prior,
                        scoreFun = scoreFun)

    return (oldOrNew)

  }

}
