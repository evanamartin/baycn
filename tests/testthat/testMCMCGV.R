context('baycn with genetic variants')

test_that('the mhEdge function infers the correct graph using PMR',{

  # Adjacency matrix for M1 - M4 -----------------------------------------------

  # Fully connected adjacency matrix for M1 - M4
  am_m <- matrix(c(0, 1, 1,
                   0, 0, 1,
                   0, 0, 0),
                 byrow = TRUE,
                 nrow = 3)

  # Run baycn on data simulated for M1 with GV ---------------------------------

  set.seed(1)

  data_m1 <- simdata(b0 = 0,
                     N = 500,
                     s = 1,
                     ss = 1,
                     graph = 'm1_gv')

  baycn_m1 <- mhEdge(adjMatrix = am_m,
                     burnIn = 0.2,
                     data = data_m1,
                     iterations = 1000,
                     nGV = 1,
                     pmr = TRUE,
                     prior = c(0.05,
                               0.05,
                               0.9),
                     progress = FALSE,
                     thinTo = 200)

  # Run baycn on data simulated for M3 with GV ---------------------------------

  data_m3 <- simdata(b0 = 0,
                     N = 500,
                     s = 1,
                     ss = 1,
                     graph = 'm3_gv')

  baycn_m3 <- mhEdge(adjMatrix = am_m,
                     burnIn = 0.2,
                     data = data_m3,
                     iterations = 1000,
                     nGV = 1,
                     pmr = TRUE,
                     prior = c(0.05,
                               0.05,
                               0.9),
                     progress = TRUE,
                     thinTo = 200)

  # Calculate the MSE for M1 and M3 --------------------------------------------

  ep_m1 <- matrix(c(1, 0, 0,
                    0, 0, 1,
                    1, 0, 0),
                  byrow = TRUE,
                  nrow = 3)

  mse_m1 <- sum((baycn_m1@posteriorES[, 2:4] - ep_m1)^2)

  ep_m3 <- matrix(c(1, 0, 0,
                    1, 0, 0,
                    0, 0, 1),
                  byrow = TRUE,
                  nrow = 3)

  mse_m3 <- sum((baycn_m3@posteriorES[, 2:4] - ep_m3)^2)

  expect_true(mse_m1 < 0.1)
  expect_true(mse_m3 < 0.1)

})
