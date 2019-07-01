#' show
#'
#' @param object An object of class mcmc.
#'
#' @importFrom methods show
#'
#' @export
setMethod('show',
          signature = 'mcmc',
          definition = function (object) {

            cat('Object of class: ',
                class(object),
                '\n',
                sep = '')

          })

#' summary
#'
#' @param object An object of class mcmc.
#'
#' @param ... Other Arguments passed to methods.
#'
#' @export
setMethod('summary',
          signature(object = 'mcmc'),
          definition = function (object, ...) {

            returned <- list()

            # Get the number of edges to loop through when calculating the
            # probability of each edge state.
            nEdges <- ncol(object@chain)

            # Get the number of samples kept.
            nRow <- nrow(object@chain)

            # Create a matrix to hold the posterior probability of each edge.
            posteriorProb <- matrix(nrow = nEdges,
                                    ncol = 3)

            # zero - the edge points from the node with a smaller index to the
            # node with a larger index
            # one - the edge points from the node with a larger index to the
            # node with a smaller index
            # two - the edge is absent between the two nodes
            colnames(posteriorProb) <- c('zero', 'one', 'two')

            rownames(posteriorProb) <- paste('edge', 1:nEdges, sep = '')

            # Calculate the proportion of 0, 1, and 2 for each edge in the
            # network.
            for (e in 1:nEdges) {

              edge_0 <- sum(object@chain[, e] == 0)
              edge_1 <- sum(object@chain[, e] == 1)
              edge_2 <- sum(object@chain[, e] == 2)

              posteriorProb[e, ] <- c(round(edge_0 / nRow, 4),
                                      round(edge_1 / nRow, 4),
                                      round(edge_2 / nRow, 4))

            }

            returned$posterior <- posteriorProb

            # Create a matrix for the summary of the likelihood.
            logLikSummary <- matrix(nrow = 1,
                                    ncol = 5)

            colnames(logLikSummary) <- c('Min', '1Q', 'Median', '3Q', 'Max')
            rownames(logLikSummary) <- ''

            logLikSummary[1, ] <- round(summary(object@likelihood)[c(1:3,
                                                                     5:6)],
                                        2)

            returned$likelihood <- logLikSummary

            # Calculate the number of unique graphs
            nGraphs <- length(unique(object@decimal))

            returned$graphs <- nGraphs

            returned$runtime <- object@time

            class(returned) <- 'summary.mcmc'

            returned

          })

#' @method print summary.mcmc
#'
#' @export
print.summary.mcmc <- function (x, ...) {

  # Display the posterior probability for each edge.
  cat('Posterior probability: \n')
  print(x$posterior)

  # Display the min, 1q, median, 3q, and max.
  cat('\n',
      'log likelihood: \n')
  print(x$likelihood)

  # Display the number of unique graphs
  cat('\n',
      'number of unique graphs: \n',
      x$graphs,
      '\n')

  # Display the amount of time it took to complete.
  cat('\n',
      'Run time in seconds: \n',
      x$runtime)

}

#' plot
#'
#' @param x An object of classs mcmc.
#'
#' @param y Optional if x is the appropriate structure.
#'
#' @param ... Other Arguments passed to methods.
#'
#' @import ggplot2
#'
#' @export
setMethod('plot',
          signature(x = 'mcmc'),
          definition = function (x, y, ...) {

            # Number of samples kept
            nSamples <- length(x@likelihood)

            # Convert the likelihood to a data frame to pass to ggplot.
            logLikelihood <- data.frame(likelihood = x@likelihood)

            # Convert the decimal numbers to a data frame to pass to ggplot.
            decimalNum <- data.frame(decimal = x@decimal)

            # Create the trace plot for the likelihood.
            p <- ggplot(logLikelihood, aes(x = 1:nSamples,
                                           y = likelihood)) +
              geom_line(color = '#5500cc',
                        size = 1) +
              labs(x = '') +
              labs(y = 'log likelihood') +
              ylim(min(logLikelihood), max(logLikelihood))

            # Create the trace plot for the decimal number.
            g <- ggplot(decimalNum, aes(x = 1:nSamples,
                                        y = decimal)) +
              geom_line(color = '#cc5500',
                        size = 1) +
              labs(x = '') +
              labs(y = 'decimal number') +
              ylim(min(decimalNum), max(decimalNum))

            # print the two plots with one column and two rows.
            multiplot(p,
                      g,
                      cols = 1)

          })
