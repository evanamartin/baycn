#' toDecimal
#'
#' Turns the edge directions of the graph into a decimal number.
#'
#' @param edgeDir A matrix of the posterior distribution sample.
#'
#' @return A vector containing the decimal number of each graph.
#'
toDecimal <- function(edgeDir) {

  # Store the edge directions and numbers that form a cycle as a decimal number.
  dn <- c()

  # The length of edgeDir and edgeNum will always be the same. This length will
  # be used to loop through each list and turn the edge directions and numbers
  # of each cycle into a decimal number.
  for (e in 1:dim(edgeDir)[1]) {

    # Hold a decimal number for the ith edge direction and number.
    iDecimal <- c()

    # Create a counter that will start at zero to use as the exponent in
    # calculating the decimal for each graph.
    counter <- 0

    for (v in dim(edgeDir)[2]:1) {

      # Take each edge direction for each edge number and turn it into a decimal
      # number.
      iDecimal[[v]] <- eval(quote(edgeDir[e, v] * 3^counter))

      counter <- counter + 1

    }

    dn[[e]] <- sum(iDecimal)

  }

  return (dn)

}
