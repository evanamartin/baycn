# cll
#
# Calculates the log likelihood of the current node.
#
# @param data A matrix with the nodes across the columns and the observations
# along the rows.
#
# @param nodeIndex The index of the for loop for the nodes
#
# @param v The index of the for loop for the parent combinations
#
# @param nGV The number of genetic variants in the graph.
#
# @param nNodes The number of nodes in the graph.
#
# @param pVector The vector of parent nodes for the current node.
#
# @return The log likelihood of the graph according to the orientation of the
# edges denoted by the DNA of the current graph.
#
#' @importFrom stats dnorm
#'
cll <- function (data,
                 nodeIndex,
                 nGV,
                 nNodes,
                 pVector,
                 v) {

  # Check if the current node is a genetic variant. If e is less than the number
  # of genetic variants then the current node is a genetic variant.
  if (nodeIndex <= nGV) {

    # Calculate the log likelihood for the multinomial variables.
    logLikelihood <- cllMultinom(data = data,
                                 nodeIndex = nodeIndex,
                                 pVector = pVector,
                                 v = v)

  } else {

    # Calculate the log likelihood for the normal variables.
    logLikelihood <- cllNormal(data = data,
                               nodeIndex = nodeIndex,
                               pVector = pVector,
                               v = v)

  }

  return (logLikelihood)

}

# cllMultinom
#
# Calculates the log likelihood for a node with multinomial data.
#
# @param data The data matrix.
#
# @param nodeIndex The index of the for loop for the nodes
#
# @param pVector The vector of parent nodes for the current node.
#
# @param v The index of the for loop for the parent combinations
#
# @return The log likelihood of the current node.
#
#' @importFrom MASS polr
#'
#' @importFrom stats dmultinom logLik
#'
cllMultinom <- function (data,
                         nodeIndex,
                         v,
                         pVector){

  # Loop through each of the genetic variant nodes and calculate the log
  # likelihood given the edge directions of the current graph.
  if (v == 1) {

    # Calculate the counts for each level of the current genetic variant.
    counts <- as.vector(table(data[, nodeIndex]))

    # Calculate the probability of each level of counts.
    lprob <- log(counts / sum(counts))

    # Store the log likelihood for the current node.
    ll <- sum(lprob * counts)

  } else {

    # Store the log likelihood for the current node.
    ll <- logLik(polr(as.factor(data[,nodeIndex])~data[,which(pVector!=0)]))[1]

  }

  return (ll)

}

# cllNormal
#
# Calculates the log likelihood for a node with normal data.
#
# @param data The data matrix.
#
# @param nodeIndex The index of the for loop for the nodes
#
# @param pVector The vector of parent nodes for the current node.
#
# @param v The index of the for loop for the parent combinations
#
# @return The log likelihood for the current node.
#
#' @importFrom stats lm logLik sd
#'
cllNormal <- function (data,
                       nodeIndex,
                       v,
                       pVector){

  if (v == 1) {

    # Store the log likelihood for the current node.
    ll <- sum(dnorm(x = data[, nodeIndex],
                    mean = mean(data[, nodeIndex]),
                    sd = sd(data[, nodeIndex]),
                    log = TRUE))

  } else {

    # Store the log likelihood for the current node.
    ll <- logLik(lm(data[, nodeIndex] ~ data[, which(pVector != 0)]))[1]

  }

  return (ll)

}
