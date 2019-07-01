#' mcmc class
#'
#' @slot chain A matrix where the rows contain the vector of edge states for the
#' graph that is accepted at each iteration.
#'
#' @slot decimal A vector where each element is the decimal number of the graph
#' accepted at each iteration of the MH algorithm.
#'
#' @slot likelihood A vector where each element is the likelihood of the graph
#' accepted at each iteration of the MH algorithm.
#'
#' @slot time The time in seconds it took the MH algorithm to run.
#'
setClass(Class = 'mcmc',
         slots = c(chain = 'matrix',
                   decimal = 'numeric',
                   likelihood = 'numeric',
                   time = 'numeric'))
