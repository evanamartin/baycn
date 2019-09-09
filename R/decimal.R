decimal <- function(edgeDir,
                    edgeNum) {

  # Number of cycles in the graph
  nCycles <- length(edgeDir)

  # Store the edge directions and indices that form a cycle as a decimal number.
  dn <- vector(mode = 'integer',
               length = nCycles)

  # The length of edgeDir and edgeNum will always be the same. This length will
  # be used to loop through each list and turn the edge directions and numbers
  # of each cycle into a decimal number.
  for (e in 1:nCycles) {

    dn[[e]] <- sum(edgeDir[[e]] * 3^edgeNum[[e]]) + sum(edgeNum[[e]])

  }

  return (dn)

}

decimalCI <- function (edgeNum,
                       individual) {

  # Get the directions of the edges between the nodes that could create a cycle.
  # These directions along with the edge numbers will be used to calculate a
  # decimal number for each potential cycle given the directions of the current
  # individual.
  curDir <- vector(mode = 'list',
                   length =  length(edgeNum))

  for (e in 1:length(edgeNum)) {

    # Directions of the edges that could form a cycle from the current
    # individual.
    curDir[[e]] <- individual[edgeNum[[e]]]

  }

  # Calculate the decimal number for each potential cycle for the current
  # individual. This will be used to determine if there is a cycle and which
  # edges are involved in the cycle.
  dnCI <- decimal(edgeDir = curDir,
                  edgeNum = edgeNum)

  return (dnCI)

}
