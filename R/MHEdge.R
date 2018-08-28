#' MHEdge
#'
#' @param adjMatrix The adjacency matrix returned by the MRPC algorithm. The
#' adjacency matrix is a matrix of zeros and ones. The ones represent an edge
#' and also indicates the direction of that edge.
#'
#' @param data A matrix with the nodes across the columns and the observations
#' along the rows.
#'
#' @param iterations An integer for the number of iterations to run the MH
#' algorithm.
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
#' @return A matrix where the first m columns are the DNA and the m + 1 column
#' is the log likelihood. Each row represents one individual.
#'
#' @export
#'
MHEdge <- function (adjMatrix,
                    data,
                    iterations,
                    mutationRate,
                    prior = c(0.33,
                              0.33,
                              0.33),
                    scoreFun) {

  # Check for potential cycles in the graph and return the edge directions, edge
  # numbers, and the decimal numbers for each cycle if any exist.
  cp <- cyclePrep(adjMatrix)

  # Call the coordinates function to get the coordinates of all the nonzero
  # elements in the adjacency matrix and to create the DNA of the individuals.
  coord <- coordinates(adjMatrix = adjMatrix)

  # Get the number of edges in the graph.
  nEdges <- dim(coord)[2]

  # Create the length of the individual which is the number of edges plus the
  # fitness of the individual.
  m <- nEdges + 1

  # Create a matrix to hold the accepted graph at each iteration of the MH
  # algorithm.
  mcGraph <- matrix(nrow = iterations,
                    ncol = m)

  # Check the output of the cyclePrep function if there aren't any cycles run
  # the populate function as normal.
  if (is.null(cp)) {

    # Randomly generate a starting graph
    graph <- sample(x = c(0, 1, 2),
                    size = nEdges,
                    replace = TRUE,
                    prob = prior)

    graph[[nEdges + 1]] <- cllEdge(individual = graph,
                                   coordinates = coord,
                                   data = data,
                                   scoreFun = scoreFun)

    for (e in 1:iterations) {

      # Mutate the graph at each iteration. A new graph is proposed and accepted
      # or rejected within the mutate function.
      graph <- mutate(individual = graph,
                      m = m,
                      mutationRate = mutationRate,
                      coordinates = coord,
                      data = data,
                      prior = prior,
                      dnUnique = NULL,
                      edgeNum = NULL,
                      wCycle = NULL,
                      scoreFun = scoreFun)

      # A matrix to hold the accepted graph at each iteration
      mcGraph[e, ] <- graph

    }

  } else {

    # Randomly generate a starting graph
    graph <- sample(x = c(0, 1, 2),
                    size = nEdges,
                    replace = TRUE,
                    prob = prior)

    # Remove any directed cycles in the graph.
    graph <- rmCycle(dnUnique = cp$dnUnique,
                     edgeNum = cp$edgeNum,
                     wCycle = cp$wCycle,
                     prior = prior,
                     individual = graph)

    # Calculate the log likelihood of the graph.
    graph[[nEdges + 1]] <- cllEdge(individual = graph,
                                   coordinates = coord,
                                   data = data,
                                   scoreFun = scoreFun)

    for (e in 1:iterations) {

      # Mutate the graph at each iteration. A new graph is proposed and accepted
      # or rejected within the mutate function.
      graph <- mutate(individual = graph,
                      m = m,
                      mutationRate = mutationRate,
                      coordinates = coord,
                      data = data,
                      prior = prior,
                      dnUnique = cp$dnUnique,
                      edgeNum = cp$edgeNum,
                      wCycle = cp$wCycle,
                      scoreFun = scoreFun)

      # A matrix to hold the accepted graph at each iteration
      mcGraph[e, ] <- graph

    }

  }

  return (mcGraph)

}
