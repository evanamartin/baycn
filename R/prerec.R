#' prerec
#'
#' Calculates the precision and recall of the inferred graph.
#'
#' @param iMatrix The undirected adjacency matrix of the inferred graph.
#'
#' @param tMatrix The undirected adjacency matrix of the true graph.
#'
#' @return A list. The first element is the precision and the second element is
#' the recall of the inferred graph.
#'
#' @export
#'
prerec <- function (iMatrix,
                    tMatrix) {

  # Extract the upper triangular matrix of the true graph's adjacency matrix.
  tEdges <- tMatrix[upper.tri(tMatrix)]

  # Extract the upper triangular matrix of the inferred graph's adjacency
  # matrix.
  iEdges <- iMatrix[upper.tri(iMatrix)]

  # Determine which edges were correctly identified.
  nCorrect <- sum(iEdges == 1 & tEdges == 1)

  # Calculate the total number of edges inferred.
  iTotal <- sum(iEdges)

  # Calculate the total number of edges in the true graph.
  tTotal <- sum(tEdges)

  # Calculate the precision of the graph.
  precision <- ifelse(test = iTotal == 0,
                      yes = 0,
                      no = nCorrect / iTotal)

  # Calculate the recall of the graph.
  recall <- nCorrect / tTotal

  return (list(precision = precision,
               recall = recall))

}
