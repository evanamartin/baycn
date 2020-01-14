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

            # Print the class name.
            cat('Object of class: ',
                class(object),
                '\n',
                sep = '')

            # Print the number of iterations.
            cat('Number of iterations: ',
                object@iterations,
                '\n',
                sep = '')

            # Print the stepsize.
            cat('Step size: ',
                object@stepSize,
                '\n',
                sep = '')

            # Print the number of nodes.
            cat('Number of nodes in the network: ',
                dim(object@posteriorPM)[1],
                '\n',
                sep = '')

            # Print the number of edges.
            cat('Number of edges considered: ',
                dim(object@posteriorES)[1],
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

            # Get the number of edges to loop through when calculating the
            # probability of each edge state.
            nEdges <- ncol(object@chain)

            # Get the number of samples kept.
            nRow <- nrow(object@chain)

            # Create a matrix for the summary of the likelihood.
            logLikSummary <- matrix(nrow = 1,
                                    ncol = 5)

            # Name the rows and columns of the log likelihood summary matrix.
            colnames(logLikSummary) <- c('Min', '1Q', 'Median', '3Q', 'Max')
            rownames(logLikSummary) <- ''

            # Calculate the summary of the log likelihood.
            logLikSummary[1, ] <- round(summary(object@likelihood)[c(1:3,
                                                                     5:6)],
                                        2)

            # Display the posterior probability for each edge.
            cat('Posterior probability: \n',
                sep = '')
            print(object@posteriorES)

            # Display the min, 1q, median, 3q, and max.
            cat('\n',
                'Log likelihood: \n',
                sep = '')
            print(logLikSummary)

            # Display the number of unique graphs
            cat('\n',
                'Number of unique graphs: ',
                length(unique(object@decimal)),
                sep = '')

            # Display the amount of time it took to complete.
            cat('\n',
                'Run time in seconds: ',
                object@time,
                sep = '')

            # Display the number of iterations, burn in, and step size
            cat('\n',
                'Iterations: ',
                object@iterations,
                '\n',
                'Burn in: ',
                object@burnIn,
                '%',
                '\n',
                'Step size: ',
                object@stepSize,
                sep = '')

          })

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

            # The following lines are to avoid the note 'Undefined global
            # functions or variables: likelihood, decimal'.
            likelihood <- NULL
            decimal <- NULL

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
