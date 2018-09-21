################################################################################
# Generate data under topology M1 with only gene expression nodes:
# T1 -> T2 -> T3
################################################################################

m1ge <- function (N, ss) {

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

################################################################################
# Generate data under topology H2:
# T1 -> T2 -> T4 -> T5 <- T3 <- T1
################################################################################

h2 <- function (N, ss) {

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
