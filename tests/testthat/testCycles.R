context('Identify cycles')

test_that('cyclePrep returns the correct cycles',{

  # Nested cycles --------------------------------------------------------------

  # Fully connected four node graph.
  adjMatrix<- matrix(1,
                     nrow = 4,
                     ncol = 4)

  diag(adjMatrix) <- 0

  cp <- cyclePrep(adjMatrix = adjMatrix,
                  nGV = 0,
                  pmr = FALSE)

  # Check that cycleDN and edgeID have the same length.
  expect_true(length(cp$cycleDN) == length(cp$edgeID))
  expect_true(length(cp$cycleDN[[1]]) == length(cp$edgeID[[1]]))
  expect_true(length(cp$cycleDN[[2]]) == length(cp$edgeID[[2]]))
  # Check that the correct number of cycle sizes is returned.
  expect_true(length(cp$cycleDN) == 2)
  # Check that the correct number of cycles for each cycle size is correct.
  expect_true(length(cp$cycleDN[[1]]) == 8)
  expect_true(length(cp$cycleDN[[2]]) == 6)
  # Check that the correct cycles are found.
  expect_setequal(cp$cycleDN, list(c(91, 255, 825, 258, 16, 749, 36, 38),
                                   c(827, 260, 266, 752, 122, 41)))

  # Disjoint cycles ------------------------------------------------------------

  # Two disjoint cycles (one three node cycle and one four node cycle).
  adjMatrix <- matrix(c(0, 1, 1, 0, 0, 0, 0,
                        1, 0, 1, 0, 0, 0, 0,
                        1, 1, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 1, 1, 0,
                        0, 0, 0, 1, 0, 1, 1,
                        0, 0, 0, 1, 1, 0, 1,
                        0, 0, 0, 0, 1, 1, 0),
                      nrow = 7,
                      byrow = TRUE)

  cp <- cyclePrep(adjMatrix = adjMatrix,
                  nGV = 0,
                  pmr = FALSE)

  # Check that cycleDN and edgeID have the same length.
  expect_true(length(cp$cycleDN) == length(cp$edgeID))
  expect_true(length(cp$cycleDN[[1]]) == length(cp$edgeID[[1]]))
  expect_true(length(cp$cycleDN[[2]]) == length(cp$edgeID[[2]]))
  # Check that the correct number of cycle sizes is returned.
  expect_true(length(cp$cycleDN) == 2)
  # Check that the correct number of cycles for each cycle size is correct.
  expect_true(length(cp$cycleDN[[1]]) == 6)
  expect_true(length(cp$cycleDN[[2]]) == 2)
  # Check that all the cycles are found.
  expect_setequal(cp$cycleDN, list(c(36, 15, 825, 7311, 2208, 258),
                                   c(2292, 6828)))

  # Single cycle ---------------------------------------------------------------

  # The following adjacency matrix is from topology GN5
  adjMatrix <- matrix(c(0, 1, 1, 0, 0,
                        0, 0, 0, 1, 0,
                        0, 0, 0, 0, 1,
                        0, 0, 0, 0, 1,
                        0, 0, 0, 0, 0),
                      byrow = TRUE,
                      nrow = 5)

  cp <- cyclePrep(adjMatrix = adjMatrix,
                  nGV = 0,
                  pmr = FALSE)

  # Check that cycleDN and edgeID have the same length.
  expect_true(length(cp$cycleDN) == length(cp$edgeID))
  expect_true(length(cp$cycleDN[[1]]) == length(cp$edgeID[[1]]))
  # Check that the correct number of cycle sizes is returned.
  expect_true(length(cp$cycleDN) == 1)
  # Check that the correct number of cycles for each cycle size is correct.
  expect_true(length(cp$cycleDN[[1]]) == 2)
  # Check that all the cycles are found.
  expect_setequal(cp$cycleDN, list(c(288, 105)))

})

test_that('cyclePrep returns NULL if there are no cycles in the graph', {

  #           T4
  #           ^
  #     T1 -> T3 <- T2
  adjMatrix <- matrix(c(0, 0, 1, 0,
                        0, 0, 1, 0,
                        0, 0, 0, 1,
                        0, 0, 0, 0),
                      byrow = TRUE,
                      nrow = 4)

  cp <- cyclePrep(adjMatrix = adjMatrix,
                  nGV = 0,
                  pmr = FALSE)

  # Check that NULL is returned when there are no cycles in the graph.
  expect_null(cp$cycleDN)
  expect_null(cp$edgeNum)

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

  cp <- cyclePrep(adjMatrix = adjMatrix,
                  nGV = 0,
                  pmr = FALSE)

  # Check that NULL is returned when there are no cycles in the graph.
  expect_null(cp$cycleDN)
  expect_null(cp$edgeNum)

})

test_that('cyclePrep correctly uses the PMR', {

  # Fully connected three node graph.
  adjMatrix<- matrix(1,
                     nrow = 3,
                     ncol = 3)

  diag(adjMatrix) <- 0

  cp <- cyclePrep(adjMatrix = adjMatrix,
                  nGV = 1,
                  pmr = TRUE)

  # Check that NULL is returned for a three node graph when using the PMR.
  expect_null(cp$cycleDN)
  expect_null(cp$edgeNum)

  # Fully connected five node graph.
  adjMatrix<- matrix(1,
                     nrow = 5,
                     ncol = 5)

  diag(adjMatrix) <- 0

  cp <- cyclePrep(adjMatrix = adjMatrix,
                  nGV = 1,
                  pmr = TRUE)

  # Check that cycleDN and edgeID have the same length.
  expect_true(length(cp$cycleDN) == length(cp$edgeID))
  expect_true(length(cp$cycleDN[[1]]) == length(cp$edgeID[[1]]))
  expect_true(length(cp$cycleDN[[2]]) == length(cp$edgeID[[2]]))
  # Check that the correct number of cycle sizes is returned.
  expect_true(length(cp$cycleDN) == 2)
  # Check that the correct number of cycles for each cycle size is correct.
  expect_true(length(cp$cycleDN[[1]]) == 8)
  expect_true(length(cp$cycleDN[[2]]) == 6)
  # Check that the correct cycles are found.
  expect_setequal(cp$cycleDN, list(c(6823, 19947, 65637, 19710, 748, 59801,
                                     2208, 2210),
                                   c(65883, 19956, 20442, 59808, 8778, 2217)))

})
