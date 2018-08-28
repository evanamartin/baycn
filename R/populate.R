#' populate
#'
#' Creates a matrix where the first m columns are the DNA and the m + 1 column
#' is the log likelihood. Each row represents one individual.
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
#' @param dnUnique A vector of unique decimal numbers. These numbers represent
#' each cycle in the graph. This argument is NULL if there are no cycles in the
#' graph.
#'
#' @param edgeNum A list containing the number of each edge in each of the
#' cycles in the graph. This argument is NULL if there are no cycles in the
#' graph.
#'
#' @param wCycle A vector of the locations in the cycleDN vector where each
#' unique number occurs. The unique numbers refer to each of the unique cycles
#' in the graph. This argument is NULL if there are no cycles in the graph.
#'
#' @param scoreFun A character string indicating what method to use for
#' calculating the fitness.
#'
#' @return A list where the first element of the list is the population matrix
#' and the second element is the coordinates matrix created by the
#' coordinates function.
#'
#' @export
#'
populate <- function (N,
                      adjMatrix,
                      data,
                      dnUnique,
                      edgeNum,
                      wCycle,
                      scoreFun) {

  # Call the coordinates function to get the coordinates of all the nonzero
  # elements in the adjacency matrix and to create the DNA of the individuals.
  coord <- coordinates(adjMatrix = adjMatrix)

  # get the length of the DNA
  n <- dim(coord)[2]

  # create an empty matrix to hold the DNA of each individual and its log
  # likelihood with respect to its DNA.
  population <- matrix(nrow = N,
                       ncol = n + 1)

  # If there aren't any potential cycles in the graph generate the population
  # with out any constraints.
  if (is.null(dnUnique)) {

    # create a population whose DNA are the edges found by MRPC.
    for (e in 1:N) {

      population[e, 1:n] <- sample(x = 0:2,
                                   size = n,
                                   replace = TRUE)

      population[e, n + 1] <- cllEdge(individual = population[e, ],
                                      coordinates = coord,
                                      data = data,
                                      scoreFun = scoreFun)

    }

  } else {

    # Create a population whose DNA are the edges found by MRPC and with all the
    # cycles removed.
    for (e in 1:N) {

      population[e, 1:n] <- sample(x = 0:2,
                                   size = n,
                                   replace = TRUE)

      population[e, 1:n] <- rmCycle(dnUnique,
                                    edgeNum,
                                    wCycle,
                                    individual = population[e, 1:n])

      population[e, n + 1] <- cllEdge(individual = population[e, ],
                                      coordinates = coord,
                                      data = data,
                                      scoreFun = scoreFun)

    }

  }

  return (list(population = population,
               coordinates = coord))

}
