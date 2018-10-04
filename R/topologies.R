#' g2
#'
#' Simulates data under topology G2.
#'
#' @param N The number of observations to generate.
#'
#' @param ss The value of the beta_1 coefficient in the linear model
#' beta_0 + beta_1 * x_1. This term is referred to as the signal strength.
#'
#' @export
g2 <- function (N,
                ss) {

  T1 <- rNoParents(N = N,
                   b0 = 0,
                   s = 1)

  T2 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T1),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T4 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T2),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T3 <- rMParents(N = N,
                  mParents = 2,
                  parentData = list(T1, T4),
                  b0 = 0,
                  b1 = list(ss, ss),
                  s = 1)

  return (cbind(T1, T2, T3, T4))

}

#' h2
#'
#' Simulates data under topology H2.
#'
#' @param N The number of observations to generate.
#'
#' @param ss The value of the beta_1 coefficient in the linear model
#' beta_0 + beta_1 * x_1. This term is referred to as the signal strength.
#'
#' @export
h2 <- function (N,
                ss) {

  T1 <- rNoParents(N = N,
                   b0 = 0,
                   s = 1)

  T2 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T1),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T3 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T1),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T4 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T2),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T5 <- rMParents(N = N,
                  mParents = 2,
                  parentData = list(T3, T4),
                  b0 = 0,
                  b1 = list(ss, ss),
                  s = 1)

  return (cbind(T1, T2, T3, T4, T5))

}

#' layer
#'
#' Simulates data under the layer topology.
#'
#' @param N The number of observations to generate.
#'
#' @param p The frequency of the reference allele.
#'
#' @param ss The value of the beta_1 coefficient in the linear model
#' beta_0 + beta_1 * x_1. This term is referred to as the signal strength.
#'
#' @export
layer <- function (N,
                   p,
                   ss) {

  V <- sample(x = 0:2,
              size = N,
              replace = TRUE,
              prob = c((1 - p)^2, 2 * p * (1 - p), p^2))

  T1 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(V),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T2 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(V),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T3 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(V),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T4 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T1),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T5 <- rMParents(N = N,
                  mParents = 2,
                  parentData = list(T1, T2),
                  b0 = 0,
                  b1 = list(ss, ss),
                  s = 1)

  T6 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T2),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T7 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T3),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  return (cbind(V, T1, T2, T3, T4, T5, T6, T7))

}

#' m1ge
#'
#' Simulates data under topology M1 with only gene expression nodes.
#'
#' @param N The number of observations to generate.
#'
#' @param ss The value of the beta_1 coefficient in the linear model
#' beta_0 + beta_1 * x_1. This term is referred to as the signal strength.
#'
#' @export
m1ge <- function (N,
                  ss) {

  T1 <- rNoParents(N = N,
                   b0 = 0,
                   s = 1)

  T2 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T1),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T3 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T2),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  return (cbind(T1, T2, T3))

}

#' m1gv
#'
#' Simulates data under topology M1 with one genetic variant node.
#'
#' @param N The number of observations to generate.
#'
#' @param p The frequency of the reference allele.
#'
#' @param ss The value of the beta_1 coefficient in the linear model
#' beta_0 + beta_1 * x_1. This term is referred to as the signal strength.
#'
#' @export
m1gv <- function (N,
                  p,
                  ss) {

  V <- sample(x = 0:2,
              size = N,
              replace = TRUE,
              prob = c((1 - p)^2, 2 * p * (1 - p), p^2))

  T1 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(V),
                  b0 = 0,
                  b1 = list(1),
                  s = 1)

  T2 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T1),
                  b0 = 0,
                  b1 = list(1),
                  s = 1)

  return (cbind(V, T1, T2))

}

#' m2ge
#'
#' Simulates data under topology M2 with only gene expression nodes.
#'
#' @param N The number of observations to generate.
#'
#' @param ss The value of the beta_1 coefficient in the linear model
#' beta_0 + beta_1 * x_1. This term is referred to as the signal strength.
#'
#' @export
m2ge <- function (N,
                  ss) {

  T1 <- rNoParents(N = N,
                   b0 = 0,
                   s = 1)

  T3 <- rNoParents(N = N,
                   b0 = 0,
                   s = 1)

  T2 <- rMParents(N = N,
                  mParents = 2,
                  parentData = list(T1, T3),
                  b0 = 0,
                  b1 = list(ss, ss),
                  s = 1)

  return (cbind(T1, T2, T3))

}

#' m2gv
#'
#' Simulates data under topology M2 with one genetic variant node.
#'
#' @param N The number of observations to generate.
#'
#' @param p The frequency of the reference allele.
#'
#' @param ss The value of the beta_1 coefficient in the linear model
#' beta_0 + beta_1 * x_1. This term is referred to as the signal strength.
#'
#' @export
m2gv <- function (N,
                  p,
                  ss) {

  V <- sample(x = 0:2,
              size = N,
              replace = TRUE,
              prob = c((1 - p)^2, 2 * p * (1 - p), p^2))

  T2 <- rNoParents(N = N,
                   b0 = 0,
                   s = 1)

  T1 <- rMParents(N = N,
                  mParents = 2,
                  parentData = list(V, T2),
                  b0 = 0,
                  b1 = list(1, 1),
                  s = 1)

  return (cbind(V, T1, T2))

}

#' m3gv
#'
#' Simulates data under topology M3 with one genetic variant node.
#'
#' @param N The number of observations to generate.
#'
#' @param p The frequency of the reference allele.
#'
#' @param ss The value of the beta_1 coefficient in the linear model
#' beta_0 + beta_1 * x_1. This term is referred to as the signal strength.
#'
#' @export
m3gv <- function (N,
                  p,
                  ss) {

  V <- sample(x = 0:2,
              size = N,
              replace = TRUE,
              prob = c((1 - p)^2, 2 * p * (1 - p), p^2))

  T1 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(V),
                  b0 = 0,
                  b1 = list(1),
                  s = 1)

  T2 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(V),
                  b0 = 0,
                  b1 = list(1),
                  s = 1)

  return (cbind(V, T1, T2))

}

#' mpge
#'
#' Simulates data under the multi-parent topology.
#'
#' @param N The number of observations to generate.
#'
#' @param ss The value of the beta_1 coefficient in the linear model
#' beta_0 + beta_1 * x_1. This term is referred to as the signal strength.
#'
#' @export
mpge <- function (N,
                  ss) {

  T1 <- rNoParents(N = N,
                   b0 = 0,
                   s = 1)

  T2 <- rNoParents(N = N,
                   b0 = 0,
                   s = 1)

  T3 <- rNoParents(N = N,
                   b0 = 0,
                   s = 1)

  T4 <- rMParents(N = N,
                  mParents = 3,
                  parentData = list(T1, T2, T3),
                  b0 = 0,
                  b1 = list(ss, ss, ss),
                  s = 1)

  return (cbind(T1, T2, T3, T4))

}

#' mpgv
#'
#' Simulates data under the multi-parent topology.
#'
#' @param N The number of observations to generate.
#'
#' @param p The frequency of the reference allele.
#'
#' @param ss The value of the beta_1 coefficient in the linear model
#' beta_0 + beta_1 * x_1. This term is referred to as the signal strength.
#'
#' @export
mpgv <- function (N,
                  p,
                  ss) {

  V <- sample(x = 0:2,
              size = N,
              replace = TRUE,
              prob = c((1 - p)^2, 2 * p * (1 - p), p^2))

  T1 <- rNoParents(N = N,
                   b0 = 0,
                   s = 1)

  T2 <- rNoParents(N = N,
                   b0 = 0,
                   s = 1)

  T3 <- rMParents(N = N,
                  mParents = 3,
                  parentData = list(V, T1, T2),
                  b0 = 0,
                  b1 = list(ss, ss, ss),
                  s = 1)

  return (cbind(V, T1, T2, T3))

}

#' me
#'
#' Simulates data under topology ME.
#'
#' @param N The number of observations to generate.
#'
#' @param ss The value of the beta_1 coefficient in the linear model
#' beta_0 + beta_1 * x_1. This term is referred to as the signal strength.
#'
#' @export
me <- function (N,
                ss) {

  T1 <- rNoParents(N = N,
                   b0 = 0,
                   s = 1)

  T2 <- rNoParents(N = N,
                   b0 = 0,
                   s = 1)

  T3 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T1),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T4 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T2),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T5 <- rMParents(N = N,
                  mParents = 2,
                  parentData = list(T3, T4),
                  b0 = 0,
                  b1 = list(ss, ss),
                  s = 1)

  return (cbind(T1, T2, T3, T4, T5))

}

#' nc11
#'
#' Simulates data under topology NC11.
#'
#' @param N The number of observations to generate.
#'
#' @param ss The value of the beta_1 coefficient in the linear model
#' beta_0 + beta_1 * x_1. This term is referred to as the signal strength.
#'
#' @export
nc11 <- function (N,
                  ss) {

  T1 <- rNoParents(N = N,
                   b0 = 0,
                   s = 1)

  T7 <- rNoParents(N = N,
                   b0 = 0,
                   s = 1)

  T2 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T1),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T3 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T2),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T4 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T3),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T5 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T4),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T8 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T7),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T9 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T8),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T10 <- rMParents(N = N,
                   mParents = 1,
                   parentData = list(T9),
                   b0 = 0,
                   b1 = list(ss),
                   s = 1)

  T11 <- rMParents(N = N,
                   mParents = 1,
                   parentData = list(T10),
                   b0 = 0,
                   b1 = list(ss),
                   s = 1)

  T6 <- rMParents(N = N,
                  mParents = 2,
                  parentData = list(T5, T7),
                  b0 = 0,
                  b1 = list(ss, ss),
                  s = 1)

  return (cbind(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11))

}

#' nc4
#'
#' Simulates data under topology NC4.
#'
#' @param N The number of observations to generate.
#'
#' @param ss The value of the beta_1 coefficient in the linear model
#' beta_0 + beta_1 * x_1. This term is referred to as the signal strength.
#'
#' @export
nc4 <- function (N,
                 ss) {

  T1 <- rNoParents(N = N,
                   b0 = 0,
                   s = 1)

  T2 <- rNoParents(N = N,
                   b0 = 0,
                   s = 1)

  T3 <- rMParents(N = N,
                  mParents = 2,
                  parentData = list(T1, T2),
                  b0 = 0,
                  b1 = list(ss, ss),
                  s = 1)

  T4 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T3),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  return (cbind(T1, T2, T3, T4))

}

#' pc
#'
#' Simulates data under topology PC.
#'
#' @param N The number of observations to generate.
#'
#' @param ss The value of the beta_1 coefficient in the linear model
#' beta_0 + beta_1 * x_1. This term is referred to as the signal strength.
#'
#' @export
pc <- function (N,
                ss) {

  T1 <- rNoParents(N = N,
                   b0 = 0,
                   s = 1)

  T2 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T1),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T3 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T2),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T4 <- rNoParents(N = N,
                   b0 = 0,
                   s = 1)

  T5 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T2),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T6 <- rMParents(N = N,
                  mParents = 2,
                  parentData = list(T1, T5),
                  b0 = 0,
                  b1 = list(ss, ss),
                  s = 1)

  T7 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T6),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T8 <- rMParents(N = N,
                  mParents = 2,
                  parentData = list(T1, T5),
                  b0 = 0,
                  b1 = list(ss, ss),
                  s = 1)

  return (cbind(T1, T2, T3, T4, T5, T6, T7, T8))

}

#' star
#'
#' Simulates data under the star topology.
#'
#' @param N The number of observations to generate.
#'
#' @param p The frequency of the reference allele.
#'
#' @param ss The value of the beta_1 coefficient in the linear model
#' beta_0 + beta_1 * x_1. This term is referred to as the signal strength.
#'
#' @export
star <- function (N,
                  p,
                  ss) {

  V <- sample(x = 0:2,
              size = N,
              replace = TRUE,
              prob = c((1 - p)^2, 2 * p * (1 - p), p^2))

  T1 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(V),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T2 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T1),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T3 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T1),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T4 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T1),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  T5 <- rMParents(N = N,
                  mParents = 1,
                  parentData = list(T1),
                  b0 = 0,
                  b1 = list(ss),
                  s = 1)

  return (cbind(V, T1, T2, T3, T4, T5))

}
