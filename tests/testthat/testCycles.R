context('Identify cycles')

test_that('cyclePrep returns the correct cycles',{

  # The following graph has three cycles two of which are nested in the largest
  # cycle.
  # T1 - T2 - T3 - T1
  # T1 - T3 - T4 - T5 - T6 - T1
  # T1 - T2 - T3 - T4 - T5 - T6 - T1
  adjMatrix <- matrix(c(0, 1, 1, 0, 0, 1,
                        1, 0, 1, 0, 0, 0,
                        1, 1, 0, 1, 0, 0,
                        0, 0, 1, 0, 1, 0,
                        0, 0, 0, 1, 0, 1,
                        1, 0, 0, 0, 1, 0),
                      byrow = TRUE,
                      nrow = 6)

  cp <- cyclePrep(adjMatrix)

  expect_identical(cp$dnUnique, c(91, 3269, 16, 3191, 50, 53))
  expect_equal(cp$edgeNum, list(c(2, 4, 1),
                                c(3, 7, 6, 5, 4, 1),
                                c(1, 4, 2),
                                c(3, 7, 6, 5, 2),
                                c(2, 5, 6, 7, 3),
                                c(1, 4, 5, 6, 7, 3)))
  expect_equal(length(cp$wCycle), 6)

  # The following adjacency matrix is from topology H2
  adjMatrix <- matrix(c(0, 1, 1, 0, 0,
                        0, 0, 0, 1, 0,
                        0, 0, 0, 0, 1,
                        0, 0, 0, 0, 1,
                        0, 0, 0, 0, 0),
                      byrow = TRUE,
                      nrow = 5)

  cp <- cyclePrep(adjMatrix)

  expect_identical(cp$dnUnique, c(288, 105))
  expect_equal(cp$edgeNum, list(c(2, 4, 5, 3, 1),
                                c(1, 3, 5, 4, 2)))
  expect_equal(length(cp$wCycle), 2)
})

test_that('cyclePrep returns NULL if there are no cycles in the graph', {

  # adjacency matrix for topology NC4
  #           T4
  #           ^
  #     T1 -> T3 <- T2
  adjMatrix <- matrix(c(0, 0, 1, 0,
                        0, 0, 1, 0,
                        0, 0, 0, 1,
                        0, 0, 0, 0),
                      byrow = TRUE,
                      nrow = 4)

  cp <- cyclePrep(adjMatrix)

  expect_true(is.null(cp$dnUnique))
  expect_true(is.null(cp$edgeNum))
  expect_true(is.null(cp$wCycle))

  # The adjacency matrix for topology NC8
  # T1 - T3 - T5 - T7
  # |    |    |    |
  # T2   T4   T6   T8
  adjMatrix <- matrix(c(0, 1, 1, 0, 0, 0, 0, 0,
                        1, 0, 0, 0, 0, 0, 0, 0,
                        1, 0, 0, 1, 1, 0, 0, 0,
                        0, 0, 1, 0, 0, 0, 0, 0,
                        0, 0, 1, 0, 0, 1, 1, 0,
                        0, 0, 0, 0, 1, 0, 0, 0,
                        0, 0, 0, 0, 1, 0, 0, 1,
                        0, 0, 0, 0, 0, 0, 1, 0),
                      byrow = TRUE,
                      nrow = 8)

  cp <- cyclePrep(adjMatrix)

  expect_true(is.null(cp$dnUnique))
  expect_true(is.null(cp$edgeNum))
  expect_true(is.null(cp$wCycle))

})
