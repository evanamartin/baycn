################################################################################
# Topology M1 with one genetic variant
# V -> T1 -> T2
################################################################################

m1gv <- function (N, p, ss) {

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

################################################################################
# Topology M3 with one genetic variant
# T1 <- V -> T2
################################################################################

m3gv <- function (N, p, ss) {

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
