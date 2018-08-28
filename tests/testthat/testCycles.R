context('Check cycles')

test_that('cyclePrep returns the correct number of cycles',{

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

  expect_equal(length(cp$dnUnique), 6)
  expect_equal(length(cp$edgeNum), 6)
  expect_equal(length(cp$wCycle), 6)

  # The following adjacency matrix comes from the data from the pcalg package
  adjMatrix <- matrix(c(0, 1, 0, 0, 0, 1, 0, 1,
                        0, 0, 1, 0, 1, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 1, 0, 1,
                        0, 0, 0, 0, 0, 0, 1, 0,
                        0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0),
                      byrow = TRUE,
                      nrow = 8)

  cp <- cyclePrep(adjMatrix)

  expect_equal(length(cp$dnUnique), 6)
  expect_equal(length(cp$edgeNum), 6)
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

  expect_equal(length(cp$dnUnique), 2)
  expect_equal(length(cp$edgeNum), 2)
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

  expect_true(is.null(cp))

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

  expect_true(is.null(cp))

})

test_that('cyclePrep returns the correct decimal for each cycle', {

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

})
