#' efMH
#'
#' Calculates the number of times each edge is oriented in each of the two
#' possible directions.
#'
#' @param dEdge A matrix from the output of the MHEdges function.
#'
#' @return A matrix. The number of rows is equal to the number of edges in the
#' graph. It contains three columns. The first column is the frequency of the
#' edge whose directions goes from the node with a smaller index to the node
#' with a larger index. The second column contains the frequencies of the edges
#' whose direction is opposite from the first column. The third column is the
#' frequency of not having that edge in the graph.
#'
#' @export
#'
efMH <- function (dEdge) {

  # Get the number of edges to loop through when calculating the frequency of
  # each edge direction for each node.
  nEdges <- ncol(dEdge) - 1

  # Get the number of rows in the dEdge matrix
  nRow <- nrow(dEdge)

  # Create a matrix to hold the frequencies of each direction for each edge.
  allFrequencies <- matrix(nrow = nEdges,
                           ncol = 3)

  # 'one' represents an edge going from the node with a smaller index to the
  # node with a larger index. 'two' represents an edge going from the node with
  # a larger index to the node with a smaller index.
  colnames(allFrequencies) <- c('zero', 'one', 'two')

  # Tally each 0, 1, or 2 for each edge in all evolved populations in the list
  # passed to the edgeFrequency function. These counts will then be divided by
  # the total number of 0s, 1s, and 2s.
  for (e in 1:nEdges) {

    edge_0 <- sum(dEdge[, e] == 0)
    edge_1 <- sum(dEdge[, e] == 1)
    edge_2 <- sum(dEdge[, e] == 2)

    allFrequencies[e, ] <- c(round(edge_0 / nRow, 4),
                             round(edge_1 / nRow, 4),
                             round(edge_2 / nRow, 4))

  }

  return(allFrequencies)

}
