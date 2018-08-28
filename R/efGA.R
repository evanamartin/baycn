#' efGA
#'
#' Calculates the number of times each edge is oriented in each of the two
#' possible directions.
#'
#' @param evolvedPop A list of evolved populations.
#'
#' @return A matrix. The number of rows is equal to the number of edges in the
#' graph. It contains two columns. The first column is the frequency of the edge
#' whose directions goes from the node with a smaller index to the node with a
#' larger index. The second column contains the frequencies of the edges whose
#' direction is opposite from the first column.
#'
#' @export
#'
efGA <- function (evolvedPop) {

  # Get the number of edges to loop through when calculating the frequency of
  # each edge direction for each node.
  nColumns <- ncol(evolvedPop[[1]]) - 1

  # Create a matrix to hold the frequencies of each direction for each edge.
  allFrequencies <- matrix(nrow = nColumns,
                           ncol = 3)

  # 'one' represents an edge going from the node with a smaller index to the
  # node with a larger index. 'two' represents an edge going from the node with
  # a larger index to the node with a smaller index.
  colnames(allFrequencies) <- c('zero', 'one', 'two')

  # Tally each 0, 1, or 2 for each edge in all evolved populations in the list
  # passed to the edgeFrequency function. These counts will then be divided by
  # the total number of 0s, 1s, and 2s.
  for (e in 1:nColumns) {

    edge_0 <- 0
    edge_1 <- 0
    edge_2 <- 0

    for(v in 1:length(evolvedPop)) {

      edge_0 <- edge_0 + sum(evolvedPop[[v]][, e] == 0)
      edge_1 <- edge_1 + sum(evolvedPop[[v]][, e] == 1)
      edge_2 <- edge_2 + sum(evolvedPop[[v]][, e] == 2)

    }

    total <- edge_0 + edge_1 + edge_2

    allFrequencies[e, ] <- c(round(edge_0 / total, 4),
                             round(edge_1 / total, 4),
                             round(edge_2 / total, 4))

  }

  return(allFrequencies)

}
