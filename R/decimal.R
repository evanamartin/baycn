# decimal
#
# Turns the directions and locations of the cycles in the graph in to a decimal
# number.
#
# @param edgeDir A list containing the directions of the edges for each cycle
# in the graph.
#
# @param edgeID A list containing the edge index of each edge in each cycle
# of the graph.
#
# @return A list containing the decimal number of each cycle in each cycle size.
#
decimal <- function(edgeDir,
                    edgeID) {

  # Get the number of cycle sizes.
  ncs <- length(edgeDir)

  # Store the edge directions and indices that form a cycle as a decimal number.
  dn <- vector(mode = 'list',
               length = ncs)

  # Loop through each cycle size and calculate the decimal for each cycle.
  for (e in 1:ncs) {

    # Number of cycles in the graph
    nCycles <- length(edgeDir[[e]])

    # Create a sub vector to store the decimal for each cycle
    dn[[e]] <- vector(mode = 'integer',
                      length = nCycles)

    # Loop through each cycle in the current cycle size and calculate the
    # decimal number.
    for (v in 1:nCycles) {

      # Calculate the decimal for the current edge directions and indices.
      dn[[e]][[v]] <- (sum(edgeDir[[e]][[v]] * 3^edgeID[[e]][[v]]) +
                         sum(edgeID[[e]][[v]]))

    }

  }

  return (dn)

}

# decimalCI
#
# Calculates the decimal number of the cycles of an individual according to its
# current edge directions.
#
# @param edgeID A list containing the index of each edge for each cycle in
# the graph.
#
# @param individual A vector containing the DNA and the log likelihood of the
# current individual.
#
# @return The decimal numbers for each potential cycle given the edge
# directions of the current individual.
#
decimalCI <- function (edgeID,
                       individual) {

  # Get the number of cycle sizes
  ncs <- length(edgeID)

  # Get the directions of the edges between the nodes that could create a cycle.
  # These directions along with the edge numbers will be used to calculate a
  # decimal number for each potential cycle given the directions of the current
  # individual.
  curDir <- vector(mode = 'list',
                   length =  ncs)

  # Loop through each cycle size.
  for (e in 1:ncs) {

    # Get the number of cycles in the eth cycle size.
    nCycles <- length(edgeID[[e]])

    # Loop through each cycle in the eth cycle size.
    for (v in 1:nCycles) {

      # Extract the directions of the edges that could form a cycle from the
      # current individual.
      curDir[[e]][[v]] <- individual[edgeID[[e]][[v]]]

    }

  }

  # Calculate the decimal number for each potential cycle for the current
  # individual. This will be used to determine if there is a cycle and which
  # edges are involved in the cycle.
  dnCI <- decimal(edgeDir = curDir,
                  edgeID = edgeID)

  return (dnCI)

}
