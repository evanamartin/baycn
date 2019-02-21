#' mhEdge
#'
#' The main function for the Metropolis-Hastings algorithm. It returns the
#' posterior distribution of the edge directions.
#'
#' @param adjMatrix The adjacency matrix returned by the MRPC algorithm. The
#' adjacency matrix is a matrix of zeros and ones. The ones represent an edge
#' and also indicates the direction of that edge.
#'
#' @param data A matrix with the nodes across the columns and the observations
#' down the rows.
#'
#' @param iterations An integer for the number of iterations to run the MH
#' algorithm.
#'
#' @param mutationRate A number between 0 and 1. As a rule of thumb the mutation
#' rate should be 1 / (number of edges in the graph).
#'
#' @param nGV The number of genetic variants in the graph.
#'
#' @param pmr Logical. If true the Metropolis-Hastings algorithm will use the
#' Principle of Mendelian Randomization, PMR. This prevents the direction of an
#' edge pointing from a gene expression node to a genetic variant node.
#'
#' @param prior A vector containing the prior probability of seeing each edge
#' direction.
#'
#' @return A matrix where the first m columns are the edge directions and the
#' m + 1 column is the log likelihood. Each row represents one graph.
#'
#' @examples
#' # Generate data under topology M1.
#' data_m1 <- m1gv(N = 200, p = 0.27, ss = 1)
#'
#' # Use the true edges as the inferred edges.
#' adjmatrix_m1 <- matrix(c(0, 1, 0,
#'                          0, 0, 1,
#'                          0, 0, 0),
#'                        byrow = TRUE,
#'                        nrow = 3)
#'
#' # Run the Metropolis-Hastings algorithm with the Principle of Mendelian
#' # Randomization (PMR).
#' mh_m1_pmr <- mhEdge(adjMatrix = adjmatrix_m1,
#'                     data = data_m1,
#'                     iterations = 500,
#'                     mutationRate = 1/2,
#'                     nGV = 1,
#'                     pmr = TRUE,
#'                     prior = c(0.05,
#'                               0.05,
#'                               0.9))
#'
#' # Run the Metropolis-Hastings algorithm on the same data as the previous
#' # example but with pmr set to FALSE.
#' mh_m1 <- mhEdge(adjMatrix = adjmatrix_m1,
#'                 data = data_m1,
#'                 iterations = 500,
#'                 mutationRate = 1/2,
#'                 nGV = 1,
#'                 pmr = FALSE,
#'                 prior = c(0.05,
#'                           0.05,
#'                           0.9))
#'
#' @export
#'
mhEdge <- function (adjMatrix,
                    data,
                    iterations,
                    mutationRate,
                    nGV,
                    pmr = FALSE,
                    prior = c(0.05,
                              0.05,
                              0.9)) {

  # Check for potential cycles in the graph and return the edge directions, edge
  # numbers, and the decimal numbers for each cycle if any exist.
  cp <- cyclePrep(adjMatrix)

  # Call the coordinates function to get the coordinates of all the nonzero
  # elements in the adjacency matrix and to create the DNA of the individuals.
  coord <- coordinates(adjMatrix = adjMatrix)

  # Get the edge types for each edge in the graph. The two edge types are gv-ge
  # or gv-gv, ge-ge.
  edgeType <- detType(coord = coord,
                      nGV = nGV,
                      pmr = pmr)

  # Get the number of edges in the graph.
  nEdges <- dim(coord)[2]

  # Get the number of nodes in the graph. This will be used in for loops
  # throughout the function.
  nNodes <- ncol(adjMatrix)

  # Create the length of the individual which is the number of edges plus the
  # fitness of the individual.
  m <- nEdges + 1

  # Create a matrix to hold the accepted graph at each iteration of the MH
  # algorithm.
  mcGraph <- matrix(nrow = iterations,
                    ncol = m)

  # Randomly generate a starting graph.
  graph <- sample(x = c(0, 1, 2),
                  size = nEdges,
                  replace = TRUE,
                  prob = prior)

  # Check the output of the cyclePrep function. If there are cycles in the graph
  # run the rmCycle function to remove any directed cycles.
  if (!is.null(cp$dnUnique)) {

    # Remove any directed cycles in the graph.
    graph <- rmCycle(dnUnique = cp$dnUnique,
                     edgeNum = cp$edgeNum,
                     edgeType = edgeType,
                     individual = graph,
                     pmr = pmr,
                     prior = prior,
                     wCycle = cp$wCycle)

  }

  # Calculate the log likelihood for the starting graph given the data.
  graph <- cllEdge(coordinates = coord,
                   data = data,
                   graph = graph,
                   m = m,
                   nGV = nGV,
                   nNodes = nNodes,
                   pmr = pmr)

  # The initial graph is the first graph in the Markov chain.
  mcGraph[1, ] <- graph

  # Run through the iterations of the Metropolis Hastings algorithm.
  for (e in 2:iterations) {

    # Mutate the graph at each iteration.
    mGraph <- mutate(edgeType = edgeType,
                     graph = mcGraph[e - 1, ],
                     m = m,
                     mutationRate = mutationRate,
                     pmr = pmr,
                     prior = prior)

    # If there are potential directed cycles in the graph check if they are
    # present and remove them if they are.
    if (!is.null(cp$dnUnique)) {

      # Remove any directed cycles in the graph.
      mGraph <- rmCycle(dnUnique = cp$dnUnique,
                        edgeNum = cp$edgeNum,
                        edgeType = edgeType,
                        individual = mGraph,
                        pmr = pmr,
                        prior = prior,
                        wCycle = cp$wCycle)

    }

    # If the graph after mutation and removing directed cycles is the same as
    # the graph at the previous iteration then return that graph without
    # calculating the log likelihood.
    if (identical(mGraph[1:nEdges], mcGraph[e - 1, 1:nEdges])) {

      # Add the e - 1 graph to the eth row of the mcGraph matrix.
      mcGraph[e, ] <- mcGraph[e - 1, ]

      # Skip the rest of the for loop and go to the next iteration
      next

    }

    # Calculate the log likelihood of the mutated graph
    mGraph <- cllEdge(coordinates = coord,
                      data = data,
                      graph = mGraph,
                      m = m,
                      nGV = nGV,
                      nNodes = nNodes,
                      pmr = pmr)

    # When pmr is set to true sometimes the graph that is passed to cllEdge
    # is different from the graph that is returned by cllEdge. This happens when
    # a gv node is the child of a ge node . It is possible that the changes made
    # to the graph in the cllEdge function revert the graph back to the graph at
    # iteration e - 1. We need to check if the two graphs are equal and return
    # the graph at iteration e - 1 if they are.
    if (identical(mGraph[1:nEdges], mcGraph[e - 1, 1:nEdges])) {

      # Add the e - 1 graph to the eth row of the mcGraph matrix.
      mcGraph[e, ] <- mcGraph[e - 1, ]

      # Skip the rest of the for loop and go to the next iteration.
      next

    }

    # Determine whether to accept the new graph or the original graph.
    acceptedGraph <- caRatio(current = mcGraph[e - 1, ],
                             edgeType = edgeType,
                             m = m,
                             pmr = pmr,
                             proposed = mGraph,
                             prior = prior)

    # Store the accepted graph at each iteration in the mcGraph matrix.
    mcGraph[e, ] <- acceptedGraph

  }

  return (mcGraph)

}
