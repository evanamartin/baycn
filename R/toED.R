toED <- function (adjMatrix,
                  coordinates,
                  graph,
                  nGV,
                  nNodes) {

  # Create a matix of NAs.
  unholy <- matrix(nrow = nNodes,
                   ncol = nNodes)

  # Fill the values only in the lower left martix. This section of the adjacency
  # matrix are the gene expression nodes that could be the parent of the genetic
  # variant nodes.
  unholy[(nGV + 1):nNodes, 1:nGV] <- adjMatrix[(nGV + 1):nNodes, 1:nGV]

  # Get the row and column indicies for the nonzero elements. The elements that
  # are a one indicate a gene expression node being a parent of a genetic
  # variant node. This edge will be removed.
  wCoord <- which(unholy == 1,
                  arr.ind = TRUE)

  # Create an empty vector to store the locations of the edges that will be
  # removed from the graph.
  wEdge <- c()

  # Change the edges where the gene expression nodes are the parents of the
  # genetic variant nodes to a 2, indicating there is no edge between the two
  # nodes.
  for (e in 1:dim(wCoord)[1]) {

    # Get the edge location.
    wEdge[[e]] <- which(apply(X = coordinates,
                              MARGIN = 2,
                              FUN = function(x, want) setequal(x, want),
                              wCoord[e, ])) # wCoord[e, ] is passed to want in
    # the function defined in FUN. It is the row and column indicies of the
    # nonzero elements in the unholy matrix.

    # Remove the edge from the graph.
    graph[[wEdge[[e]]]] <- 2

  }

  return (graph)

}
