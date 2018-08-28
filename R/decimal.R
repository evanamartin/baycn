#' decimal
#'
#' Turns the directions and locations of the cycles in the graph in to a decimal
#' number.
#'
#' @param edgeDir A list containing the directions of the edges for each cycle
#' in the graph.
#'
#' @param edgeNum A list containing the edge number of each edge in each cycle
#' of the graph.
#'
#' @return A vector containing the decimal number of each cycle.
#'
#' @export
#'
decimal <- function(edgeDir,
                    edgeNum) {

  # Store the edge directions and numbers that form a cycle as a decimal number.
  dn <- c()

  # The length of edgeDir and edgeNum will always be the same. This length will
  # be used to loop through each list and turn the edge directions and numbers
  # of each cycle into a decimal number.
  for (e in 1:length(edgeDir)) {

    # Hold a decimal number for the ith edge direction and number.
    iDecimal <- c()

    for (v in 1:length(edgeNum[[e]])) {

      # Take each edge direction for each edge number and turn it into a decimal
      # number.
      iDecimal[[v]] <- eval(quote(edgeDir[[e]][[v]] * 3^edgeNum[[e]][[v]]))


    }

    dn[[e]] <- sum(iDecimal, edgeNum[[e]])

  }

  return (dn)

}
