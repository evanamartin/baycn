#' evolveEdges
#'
#' @param N An integer for the number of individuals to generate.
#'
#' @param adjMatrix The adjacency matrix returned by the MRPC algorithm. The
#' adjacency matrix is a matrix of zeros and ones. The ones represent an edge
#' and also indicates the direction of that edge.
#'
#' @param data A matrix with the nodes across the columns and the observations
#' along the rows.
#'
#' @param tournamentSize An integer. The number of individuals that will be
#' passed to the selectBest and selectWorst functions.
#'
#' @param mutationRate A number between 0 and 1. As a rule of thumb the mutation
#' rate should be 1 / length of the DNA.
#'
#' @param generations The number of times the individual with the highest log
#' likelihood is selected, cloned, and mutated.
#'
#' @param scoreFun A character string indicating what method to use for
#' calculating the fitness.
#'
#' @return A matrix where the first m columns are the DNA and the m + 1 column
#' is the log likelihood. Each row represents one individual.
#'
#' @export
#'
evolveEdges <- function (N,
                         adjMatrix,
                         data,
                         tournamentSize,
                         mutationRate,
                         generations,
                         scoreFun) {

  # Check for potential cycles in the graph and return the edge directions, edge
  # numbers, and the decimal numbers for each cycle if any exist.
  cp <- cyclePrep(adjMatrix)

  # Check the output of the cyclePrep function if there aren't any cycles run
  # the populate function as normal.
  if (is.null(cp)) {

    # Randomly generate the starting population.
    pop <- populate(N = N,
                    adjMatrix = adjMatrix,
                    data = data,
                    dnUnique = NULL,
                    edgeNum = NULL,
                    wCycle = NULL,
                    scoreFun = scoreFun)

    for (e in 1:generations) {

      # select the location of the individual with the highest fitness among the
      # individuals selected for the tournament to reproduce
      bestI <- selectBest(population = pop$population,
                          tournamentSize = tournamentSize,
                          scoreFun = scoreFun)

      # from the individuals selected for the losing tournament select the
      # location of the individual with the lowest fitness
      worstI <- selectWorst(population = pop$population,
                            tournamentSize = tournamentSize,
                            scoreFun = scoreFun)

      # repalce the loser DNA with the winner DNA
      pop$population[worstI, ] <- pop$population[bestI, ]

      # mutate the winner DNA
      pop$population[worstI, ] <- mutate(individual = pop$population[worstI, ],
                                         mutationRate = mutationRate,
                                         coordinates = pop$coordinates,
                                         data = data,
                                         dnUnique = NULL,
                                         edgeNum = NULL,
                                         wCycle = NULL,
                                         scoreFun = scoreFun)
    }

  } else {

    # Randomly generate the starting population and remove any cycles that show
    # up in any of the individuals.
    pop <- populate(N = N,
                    adjMatrix = adjMatrix,
                    data = data,
                    dnUnique = cp$dnUnique,
                    edgeNum = cp$edgeNum,
                    wCycle = cp$wCycle,
                    scoreFun = scoreFun)

    for (e in 1:generations) {

      # select the location of the individual with the highest fitness among the
      # individuals selected for the tournament to reproduce
      bestI <- selectBest(population = pop$population,
                          tournamentSize = tournamentSize,
                          scoreFun = scoreFun)

      # from the individuals selected for the losing tournament select the
      # location of the individual with the lowest fitness
      worstI <- selectWorst(population = pop$population,
                            tournamentSize = tournamentSize,
                            scoreFun = scoreFun)

      # repalce the loser DNA with the winner DNA
      pop$population[worstI, ] <- pop$population[bestI, ]

      # mutate the winner DNA
      pop$population[worstI, ] <- mutate(individual = pop$population[worstI, ],
                                         mutationRate = mutationRate,
                                         coordinates = pop$coordinates,
                                         data = data,
                                         dnUnique = cp$dnUnique,
                                         edgeNum = cp$edgeNum,
                                         wCycle = cp$wCycle,
                                         scoreFun = scoreFun)
    }

  }

  return (pop$population)

}
