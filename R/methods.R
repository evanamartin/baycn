#' show
#'
#' @param object An object of class baycn.
#'
#' @importFrom methods show
#'
#' @export
setMethod('show',
          signature = 'baycn',
          definition = function (object) {

            cat('Object of class: ',
                class(object),
                '\n',
                sep = '')

          })

#' summary
#'
#' @param object An object of class baycn.
#'
#' @param ... Other Arguments passed to methods.
#'
#' @export
setMethod('summary',
          signature(object = 'baycn'),
          definition = function (object, ...) {

            returned <- list()

            # Get the number of edges to loop through when calculating the
            # probability of each edge state.
            nEdges <- ncol(object@chain)

            # Get the number of samples kept.
            nRow <- nrow(object@chain)

            returned$posterior <- object@posteriorES

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

            # Add the runtime in seconds to the returned list.
            returned$runtime <- object@time

            # Add the iterations, burn in, and step size to the returned list.
            returned$iterations <- object@iterations
            returned$burnIn <- object@burnIn
            returned$stepSize <- object@stepSize

            class(returned) <- 'summary.baycn'

            returned

          })

#' @method print summary.baycn
#'
#' @export
print.summary.baycn <- function (x, ...) {

  # Display the posterior probability for each edge.
  cat('Posterior probability: \n')
  print(x$posterior)

  # Display the min, 1q, median, 3q, and max.
  cat('\n',
      'Log likelihood: \n',
      sep = '')
  print(x$likelihood)

  # Display the number of unique graphs
  cat('\n',
      'Number of unique graphs: ',
      x$graphs,
      sep = '')

  # Display the amount of time it took to complete.
  cat('\n',
      'Run time in seconds: ',
      x$runtime,
      sep = '')

  # Display the number of iterations, burn in, and step size
  cat('\n',
      'Iterations: ',
      x$iterations,
      sep = '')

  cat('\n',
      'Burn in: ',
      x$burnIn,
      '%',
      sep = '')

  cat('\n',
      'Step size: ',
      x$stepSize,
      sep = '')

}

#' plot
#'
#' @param x An object of classs baycn.
#'
#' @param y Optional if x is the appropriate structure.
#'
#' @param ... Other Arguments passed to methods.
#'
#' @aliases plot,baycn-method
#'
#' @import ggplot2
#'
#' @importFrom egg ggarrange
#'
#' @export
setMethod('plot',
          signature(x = 'baycn'),
          definition = function (x, y, ...) {

            likelihood <- NULL

            # Number of samples kept
            nSamples <- length(x@likelihood)

            # Convert the likelihood to a data frame to pass to ggplot.
            logLikelihood <- data.frame(likelihood = x@likelihood)

            # Convert the decimal numbers to a data frame to pass to ggplot.
            decimalNum <- data.frame(decimal = x@decimal)

            # Create the trace plot for the likelihood.
            p <- ggplot(logLikelihood, aes(x = 1:nSamples,
                                           y = likelihood)) +
              ggtitle('log likelihood') +
              theme(panel.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    plot.title = element_text(hjust = 0.5),
                    axis.line = element_line(color = 'black')) +
              geom_line(color = '#5500cc',
                        size = 1) +
              labs(x = '') +
              labs(y = '') +
              ylim(min(logLikelihood), max(logLikelihood))

            # Create the trace plot for the decimal number.
            g <- ggplot(decimalNum, aes(x = 1:nSamples,
                                        y = decimal)) +
              ggtitle('decimal number') +
              theme(panel.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    plot.title = element_text(hjust = 0.5),
                    axis.line = element_line(color = 'black')) +
              geom_line(color = '#cc5500',
                        size = 1) +
              labs(x = '') +
              labs(y = '') +
              ylim(min(decimalNum), max(decimalNum))

            # print the two plots with one column and two rows.
            egg::ggarrange(p, g, ncol = 1)

          })
