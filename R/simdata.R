#' simdata
#'
#' A function for simulating data under various topologies for continuous and
#' mixed data.
#'
#' @param b0 The mean of the variable if it does not have any parents. If the
#' variable has one or more parents it is the slope in the linear model that is
#' the mean of the normally distributed variables.
#'
#' @param N The number of observations to simulate.
#'
#' @param s The standard deviation of the normal distribution.
#'
#' @param graph A character string of the graph for which data will be
#' simulated. The graphs that can be chosen are m1_ge, m1_gv, m1_cp, m1_cc,
#' m1_iv, m2_ge, m2_gv, m2_cp, m2_cc, m2_iv, mp_ge, mp_gv, gn4, gn5, gn8, gn11,
#' layer, and star.
#'
#' The following figures show the graph for each of the topologies listed above.
#' The nodes with a circle around the name are normally distributed and the
#' nodes with a diamond around the name are distributed multinomial. The nodes
#' labeled with a C represent confounding variables. The Principle of Mendelian
#' Randomization (PMR) can be used on graphs with discrete and continuous
#' random variables. This introduces the constraint that the continuous random
#' variables cannot be parents of discrete random variables and edges between
#' these types of variables only have two states: absent and directed with the
#' discrete random variable as the parent.
#'
#' m1_ge - Topology M1 with three continuous random variables. In this case M1
#' has two other Markov equivalent graphs.
#'
#' m1_gv - Topology M1 with one discrete random variable U and two continuous
#' random variables. When using the PMR this graph does not have any other
#' Markov equivalent graphs.
#'
#' We consider three types of confounding variables (Yang et al., 2017):
#'
#' m1_cp - Topology M1 with n common parent confounding variables.
#'
#' m1_cc - Topology M1 with n common child confounding variables.
#'
#' m1_iv - Topology M1 with n intermediate confounding variables.
#'
#' \if{html}{\figure{m1.png}{options: width="75\%" alt="Figure: m1.png"}}
#' \if{latex}{\figure{m1.pdf}{options: width=12cm}}
#'
#' m2_ge - Topology M2 with three continuous random variables. This graph is a v
#' structure and does not have any other Markov equivalent graphs.
#'
#' m2_gv - Topolog M2 with one discrete random variable U and two continuous
#' random variables. This graph is a v structure and does not have any other
#' Markov equivalent graphs.
#'
#' m2_cp - Topology M2 with n common parent confounding variables.
#'
#' m2_cc - Topology M2 with n common child confounding variables.
#'
#' m2_iv - Topology M2 with n intermediate confounding variables.
#'
#' \if{html}{\figure{m2.png}{options: width="75\%" alt="Figure: m2.png"}}
#' \if{latex}{\figure{m2.pdf}{options: width=12cm}}
#'
#' mp_ge - The multi-parent topology with continuous random variables. This
#' graph is made up of multiple v structures and has no other Markov equivalent
#' graphs.
#'
#' mp_gv - The multi-parent topology with one discrete random variable. This
#' graph is made up of multiple v structures and has no other Markov equivalent
#' graphs.
#'
#' \if{html}{\figure{mp.png}{options: width="55\%" alt="Figure: mp.png"}}
#' \if{latex}{\figure{mp.pdf}{options: width=7cm}}
#'
#' gn4 - Topology GN4 is formed by combining topologies M1 and M2. The Markov
#' equivalence class for this topology is made up of three different graphs.
#'
#' gn5 - Topology GN5 has three other Markov equivalent graphs.
#'
#' gn8 - Topology GN8 has three overlapping cycles, two v structures, and two
#' other Markov equivalent graphs.
#'
#' gn11 - Topology GN11 has two sub-graphs separated by a v structure at node
#' T6.
#'
#' \if{html}{\figure{gn.png}{options: width="55\%" alt="Figure: gn.png"}}
#' \if{latex}{\figure{gn.pdf}{options: width=7cm}}
#'
#' layer - The layer topology has no other Markov equivalent graphs when using
#' the PMR and is made up of multiple M1 topologies.
#'
#' \if{html}{\figure{layer.png}{options: width="25\%" alt="Figure: layer.png"}}
#' \if{latex}{\figure{layer.pdf}{options: width=2.5cm}}
#'
#' star - The star topology has no other Markov equivalent graphs when using the
#' PMR and is made up of multiple M1 topologies.
#'
#' \if{html}{\figure{star.png}{options: width="25\%" alt="Figure: star.png"}}
#' \if{latex}{\figure{star.pdf}{options: width=2.5cm}}
#'
#' @param ss The coefficient of the parent nodes (if there are any) in the
#' linear model that is the mean of the normally distributed variables. This
#' coefficient is referred to as the signal strength.
#'
#' @param p The frequency of the reference allele.
#'
#' @param ssc The signal strength of the confounding variables.
#'
#' @param nConfounding The number of confounding variables to simulate.
#'
#' @return A matrix with the variables across the columns and the observations
#' down the rows.
#'
#' @references Yang, F., Wang, J., The GTEx Consortium, Pierce, B. L., and
#' Chen, L. S. (2017).
#' Identifying cis-mediators for trans-eQTLs across many human tissues using
#' genomic mediation analysis. \emph{Genome Res.} 27, 1859-1871.
#'
#' @examples
#' # Generate data under topology GN4.
#' data_gn4 <- simdata(b0 = 1,
#'                     N = 500,
#'                     s = 1,
#'                     graph = 'gn4',
#'                     ss = 1)
#'
#' # Display the first few rows of the data.
#' data_gn4[1:5, ]
#'
#' # Generate data under topology M1 with 3 intermediate confounding variables.
#' data_m1_iv <- simdata(b0 = 0,
#'                       N = 500,
#'                       s = 1,
#'                       graph = 'm1_iv',
#'                       ss = 1,
#'                       p = 0.1,
#'                       ssc = 0.2,
#'                       nConfounding = 3)
#'
#' # Show the first few rows of the data.
#' data_m1_iv[1:5, ]
#'
#' @export
#'
simdata <- function (b0 = 0,
                     N = 500,
                     s = 1,
                     ss = 1,
                     graph = 'gn4',
                     p = 0.1,
                     ssc = 0.2,
                     nConfounding = 2) {

  switch(graph,

         # m1_ge ---------------------------------------------------------------

         'm1_ge' = {

           T1 <- rNoParents(N = N,
                            b0 = b0,
                            s = s)

           T2 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T1),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T3 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T2),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           return (cbind(T1, T2, T3))

         },

         # m1_gv ---------------------------------------------------------------

         'm1_gv' = {

           U <- sample(x = 0:2,
                       size = N,
                       replace = TRUE,
                       prob = c((1 - p)^2,
                                2 * p * (1 - p),
                                p^2))

           T1 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(U),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T2 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T1),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           return (cbind(U, T1, T2))

         },

         # m1_cp ---------------------------------------------------------------

         'm1_cp' = {

           U <- sample(0:2,
                       size = N,
                       replace = TRUE,
                       prob = c(p^2,
                                2 * p * (1 - p),
                                (1 - p)^2))

           # create a vector with the signal strength of the parents for T1
           ssT1 <- vector(mode = 'numeric',
                          length = nConfounding + 1)

           ssT1[[1]] <- ss

           # Create a list with the data of the parents for T1 as the elements
           # of the list.
           parT1 <- vector(mode = 'list',
                           length = nConfounding + 1)

           parT1[[1]] <- U

           # create a vector with the signal strength of the parents for T2
           ssT2 <- vector(mode = 'numeric',
                          length = nConfounding + 1)

           ssT2[[1]] <- ss

           # Create a list with the data of the parents for T2 as the elements
           # of the list.
           parT2 <- vector(mode = 'list',
                           length = nConfounding + 1)

           # create a list to hold the data for the confounding variables
           con <- matrix(nrow = N,
                         ncol = nConfounding)

           # add column names to the matrix of confounding variables
           colnames(con) <- paste0('C', 1:nConfounding)

           # Simulate the data for the hidden variables
           for (a in 1:nConfounding) {

             con[, a] <- rNoParents(N = N,
                                    b0 = b0,
                                    s = s)

             # Add the signal strength of the confounding variable to the ssT1
             # vector.
             ssT1[[a + 1]] <- ssc

             # Add the data of the current confounding variable to the parT2
             # list.
             parT1[[a + 1]] <- con[, a]

             # Add the signal strength of the confounding variable to the ssT2
             # vector.
             ssT2[[a + 1]] <- ssc

             # Add the data of the current confounding variable to the parT2
             # list.
             parT2[[a + 1]] <- con[, a]

           }

           T1 <- rMParents(N = N,
                           mParents = nConfounding + 1,
                           parentData = parT1,
                           b0 = b0,
                           b1 = ssT1,
                           s = s)

           # Fill in the first element of parT2 with data generated from T1.
           parT2[[1]] <- T1

           T2 <- rMParents(N = N,
                           mParents = nConfounding + 1,
                           parentData = parT2,
                           b0 = b0,
                           b1 = ssT2,
                           s = s)

           return (cbind(U, T1, T2, con))

         },

         # m1_cc ---------------------------------------------------------------

         'm1_cc' = {

           U <- sample(x = 0:2,
                       size = N,
                       replace = TRUE,
                       prob = c((1 - p)^2,
                                2 * p * (1 - p),
                                p^2))

           T1 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(U),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T2 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T1),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           # create a list to hold the data for the confounding variables
           con <- matrix(nrow = N,
                         ncol = nConfounding)

           # add column names to the matrix of confounding variables
           colnames(con) <- paste0('C', 1:nConfounding)

           # Simulate the data for the hidden variables
           for (a in 1:nConfounding) {

             con[, a] <- rMParents(N = N,
                                   mParents = 2,
                                   parentData = list(T1, T2),
                                   b0 = b0,
                                   b1 = c(ssc, ssc),
                                   s = s)

           }

           return (cbind(U, T1, T2, con))

         },

         # m1_iv ---------------------------------------------------------------

         'm1_iv' = {

           U <- sample(0:2,
                       size = N,
                       replace = TRUE,
                       prob = c(p^2,
                                2 * p * (1 - p),
                                (1 - p)^2))

           T1 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(U),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           # create a list to hold the data for the confounding variables
           con <- matrix(nrow = N,
                         ncol = nConfounding)

           # add column names to the matrix of confounding variables
           colnames(con) <- paste0('C', 1:nConfounding)

           # create a vector with the signal strength of the parents for T2
           ssT2 <- vector(mode = 'numeric',
                          length = nConfounding + 1)

           ssT2[[1]] <- ss

           # Create a list with the data of the parents for T2 as the elements
           # of the list.
           parT2 <- vector(mode = 'list',
                           length = nConfounding + 1)

           parT2[[1]] <- T1

           # Simulate the data for the hidden variables
           for (a in 1:nConfounding) {

             con[, a] <- rMParents(N = N,
                                   mParents = 1,
                                   parentData = list(T1),
                                   b0 = b0,
                                   b1 = c(ssc),
                                   s = s)

             # Add the signal strength of the confounding variable to the ssT2
             # vector.
             ssT2[[a + 1]] <- ssc

             # Add the data of the current confounding variable to the parT2
             # list.
             parT2[[a + 1]] <- con[, a]

           }

           T2 <- rMParents(N = N,
                           mParents = nConfounding + 1,
                           parentData = parT2,
                           b0 = b0,
                           b1 = ssT2,
                           s = s)

           return (cbind(U, T1, T2, con))

         },

         # m2_ge ---------------------------------------------------------------

         'm2_ge' = {

           T1 <- rNoParents(N = N,
                            b0 = b0,
                            s = s)

           T3 <- rNoParents(N = N,
                            b0 = b0,
                            s = s)

           T2 <- rMParents(N = N,
                           mParents = 2,
                           parentData = list(T1, T3),
                           b0 = b0,
                           b1 = c(ss, ss),
                           s = s)

           return (cbind(T1, T2, T3))

         },

         # m2_gv ---------------------------------------------------------------

         'm2_gv' = {

           U <- sample(x = 0:2,
                       size = N,
                       replace = TRUE,
                       prob = c((1 - p)^2,
                                2 * p * (1 - p),
                                p^2))

           T2 <- rNoParents(N = N,
                            b0 = b0,
                            s = s)

           T1 <- rMParents(N = N,
                           mParents = 2,
                           parentData = list(U, T2),
                           b0 = b0,
                           b1 = c(ss, ss),
                           s = s)

           return (cbind(U, T1, T2))

         },

         # m2_cp ---------------------------------------------------------------

         'm2_cp' = {

           U <- sample(0:2,
                       size = N,
                       replace = TRUE,
                       prob = c(p^2,
                                2 * p * (1 - p),
                                (1 - p)^2))

           # create a vector with the signal strength of the parents for T1
           ssT1 <- vector(mode = 'numeric',
                          length = nConfounding + 2)

           # Fill in the first and second elements of the ssT1 vector. Both
           # elements will be ss becasue the first two parents of T1 are U and
           # T2.
           ssT1[[1]] <- ss
           ssT1[[2]] <- ss

           # Create a list with the data of the parents for T1 as the elements
           # of the list.
           parT1 <- vector(mode = 'list',
                           length = nConfounding + 2)

           # The first parent of T1 is U.
           parT1[[1]] <- U

           # create a vector with the signal strength of the parents for T2
           ssT2 <- vector(mode = 'numeric',
                          length = nConfounding)

           # Create a list with the data of the parents for T2 as the elements
           # of the list.
           parT2 <- vector(mode = 'list',
                           length = nConfounding)

           # create a list to hold the data for the confounding variables
           con <- matrix(nrow = N,
                         ncol = nConfounding)

           # add column names to the matrix of confounding variables
           colnames(con) <- paste0('C', 1:nConfounding)

           # Simulate the data for the hidden variables
           for (a in 1:nConfounding) {

             con[, a] <- rNoParents(N = N,
                                    b0 = b0,
                                    s = s)

             # Add the signal strength of the confounding variable to the ssT1
             # vector.
             ssT1[[a + 2]] <- ssc

             # Add the data of the current confounding variable to the parT2
             # list.
             parT1[[a + 2]] <- con[, a]

             # Add the signal strength of the confounding variable to the ssT2
             # vector.
             ssT2[[a]] <- ssc

             # Add the data of the current confounding variable to the parT2
             # list.
             parT2[[a]] <- con[, a]

           }

           T2 <- rMParents(N = N,
                           mParents = nConfounding,
                           parentData = parT2,
                           b0 = b0,
                           b1 = ssT2,
                           s = s)

           # Fill in the second element of the parT1 list with the data from T2
           parT1[[2]] <- T2

           T1 <- rMParents(N = N,
                           mParents = nConfounding + 2,
                           parentData = parT1,
                           b0 = b0,
                           b1 = ssT1,
                           s = s)

           return (cbind(U, T1, T2, con))

         },

         # m2_cc ---------------------------------------------------------------

         'm2_cc' = {

           U <- sample(x = 0:2,
                       size = N,
                       replace = TRUE,
                       prob = c((1 - p)^2,
                                2 * p * (1 - p),
                                p^2))

           T2 <- rNoParents(N = N,
                            b0 = b0,
                            s = s)

           T1 <- rMParents(N = N,
                           mParents = 2,
                           parentData = list(U, T2),
                           b0 = b0,
                           b1 = c(ss, ss),
                           s = s)

           # create a list to hold the data for the confounding variables
           con <- matrix(nrow = N,
                         ncol = nConfounding)

           # add column names to the matrix of confounding variables
           colnames(con) <- paste0('C', 1:nConfounding)

           # Simulate the data for the hidden variables
           for (a in 1:nConfounding) {

             con[, a] <- rMParents(N = N,
                                   mParents = 2,
                                   parentData = list(T1, T2),
                                   b0 = b0,
                                   b1 = c(ssc, ssc),
                                   s = s)

           }

           return (cbind(U, T1, T2, con))

         },

         # m2_iv ---------------------------------------------------------------

         'm2_iv' = {

           U <- sample(0:2,
                       size = N,
                       replace = TRUE,
                       prob = c(p^2,
                                2 * p * (1 - p),
                                (1 - p)^2))

           T2 <- rNoParents(N = N,
                            b0 = b0,
                            s = s)

           # create a list to hold the data for the confounding variables
           con <- matrix(nrow = N,
                         ncol = nConfounding)

           # add column names to the matrix of confounding variables
           colnames(con) <- paste0('C', 1:nConfounding)

           # create a vector with the signal strength of the parents for T1
           ssT1 <- vector(mode = 'numeric',
                          length = nConfounding + 2)

           # The first two elements will be the signal strength of U and T2.
           ssT1[[1]] <- ss
           ssT1[[2]] <- ss

           # Create a list with the data of the parents for T1 as the elements
           # of the list.
           parT1 <- vector(mode = 'list',
                           length = nConfounding + 2)

           # The first to elements will be the data from U and T2.
           parT1[[1]] <- U
           parT1[[2]] <- T2

           # Simulate the data for the hidden variables
           for (a in 1:nConfounding) {

             con[, a] <- rMParents(N = N,
                                   mParents = 1,
                                   parentData = list(T2),
                                   b0 = b0,
                                   b1 = c(ssc),
                                   s = s)

             # Add the signal strength of the confounding variable to the ssT1
             # vector.
             ssT1[[a + 2]] <- ssc

             # Add the data of the current confounding variable to the parT1
             # list.
             parT1[[a + 2]] <- con[, a]

           }

           T1 <- rMParents(N = N,
                           mParents = nConfounding + 2,
                           parentData = parT1,
                           b0 = b0,
                           b1 = ssT1,
                           s = s)

           return (cbind(U, T1, T2, con))

         },

         # mp_ge ---------------------------------------------------------------

         'mp_ge' = {

           T1 <- rNoParents(N = N,
                            b0 = b0,
                            s = s)

           T2 <- rNoParents(N = N,
                            b0 = b0,
                            s = s)

           T3 <- rNoParents(N = N,
                            b0 = b0,
                            s = s)

           T4 <- rMParents(N = N,
                           mParents = 3,
                           parentData = list(T1, T2, T3),
                           b0 = b0,
                           b1 = c(ss, ss, ss),
                           s = s)

           return (cbind(T1, T2, T3, T4))

         },

         # mp_gv ---------------------------------------------------------------

         'mp_gv' = {

           U <- sample(x = 0:2,
                       size = N,
                       replace = TRUE,
                       prob = c((1 - p)^2,
                                2 * p * (1 - p),
                                p^2))

           T1 <- rNoParents(N = N,
                            b0 = b0,
                            s = s)

           T2 <- rNoParents(N = N,
                            b0 = b0,
                            s = s)

           T3 <- rMParents(N = N,
                           mParents = 3,
                           parentData = list(U, T1, T2),
                           b0 = b0,
                           b1 = c(ss, ss, ss),
                           s = s)

           return (cbind(U, T1, T2, T3))

         },

         # gn4 -----------------------------------------------------------------

         'gn4' = {

           T1 <- rNoParents(N = N,
                            b0 = b0,
                            s = s)

           T2 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T1),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T4 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T2),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T3 <- rMParents(N = N,
                           mParents = 2,
                           parentData = list(T1, T4),
                           b0 = b0,
                           b1 = c(ss, ss),
                           s = s)

           return (cbind(T1, T2, T3, T4))

         },

         # gn5 -----------------------------------------------------------------

         'gn5' = {

           T1 <- rNoParents(N = N,
                            b0 = b0,
                            s = s)

           T2 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T1),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T3 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T1),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T4 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T2),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T5 <- rMParents(N = N,
                           mParents = 2,
                           parentData = list(T3, T4),
                           b0 = b0,
                           b1 = c(ss, ss),
                           s = s)

           return (cbind(T1, T2, T3, T4, T5))

         },

         # gn8 -----------------------------------------------------------------

         'gn8' = {

           T1 <- rNoParents(N = N,
                            b0 = b0,
                            s = s)

           T2 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T1),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T3 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T2),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T4 <- rNoParents(N = N,
                            b0 = b0,
                            s = s)

           T5 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T2),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T6 <- rMParents(N = N,
                           mParents = 2,
                           parentData = list(T1, T5),
                           b0 = b0,
                           b1 = c(ss, ss),
                           s = s)

           T7 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T6),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T8 <- rMParents(N = N,
                           mParents = 2,
                           parentData = list(T1, T5),
                           b0 = b0,
                           b1 = c(ss, ss),
                           s = s)

           return (cbind(T1, T2, T3, T4, T5, T6, T7, T8))

         },

         # gn11 ----------------------------------------------------------------

         'gn11' = {

           T1 <- rNoParents(N = N,
                            b0 = b0,
                            s = s)

           T7 <- rNoParents(N = N,
                            b0 = b0,
                            s = s)

           T2 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T1),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T3 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T2),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T4 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T3),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T5 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T4),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T8 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T7),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T9 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T8),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T10 <- rMParents(N = N,
                            mParents = 1,
                            parentData = list(T9),
                            b0 = b0,
                            b1 = c(ss),
                            s = s)

           T11 <- rMParents(N = N,
                            mParents = 1,
                            parentData = list(T10),
                            b0 = b0,
                            b1 = c(ss),
                            s = s)

           T6 <- rMParents(N = N,
                           mParents = 2,
                           parentData = list(T5, T7),
                           b0 = b0,
                           b1 = c(ss, ss),
                           s = s)

           return (cbind(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11))

         },

         # layer ---------------------------------------------------------------

         'layer' = {

           U <- sample(x = 0:2,
                       size = N,
                       replace = TRUE,
                       prob = c((1 - p)^2,
                                2 * p * (1 - p),
                                p^2))

           T1 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(U),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T2 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(U),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T3 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(U),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T4 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T1),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T5 <- rMParents(N = N,
                           mParents = 2,
                           parentData = list(T1, T2),
                           b0 = b0,
                           b1 = c(ss, ss),
                           s = s)

           T6 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T2),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T7 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T3),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           return (cbind(U, T1, T2, T3, T4, T5, T6, T7))

         },

         # layer_cp ------------------------------------------------------------

         'layer_cp' = {

           U <- sample(x = 0:2,
                       size = N,
                       replace = TRUE,
                       prob = c((1 - p)^2,
                                2 * p * (1 - p),
                                p^2))

           C1 <- rNoParents(N = N,
                            b0 = b0,
                            s = s)

           C2 <- rNoParents(N = N,
                            b0 = b0,
                            s = s)

           T1 <- rMParents(N = N,
                           mParents = 2,
                           parentData = list(U, C1),
                           b0 = b0,
                           b1 = c(ss, ssc),
                           s = s)

           T2 <- rMParents(N = N,
                           mParents = 2,
                           parentData = list(U, C2),
                           b0 = b0,
                           b1 = c(ss, ssc),
                           s = s)

           T3 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(U),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T4 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T1),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T5 <- rMParents(N = N,
                           mParents = 3,
                           parentData = list(T1, T2, C1),
                           b0 = b0,
                           b1 = c(ss, ss, ssc),
                           s = s)

           T6 <- rMParents(N = N,
                           mParents = 2,
                           parentData = list(T2, C2),
                           b0 = b0,
                           b1 = c(ss, ssc),
                           s = s)

           T7 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T3),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           return (cbind(U, T1, T2, T3, T4, T5, T6, T7, C1, C2))

         },

         # layer_iv ------------------------------------------------------------

         'layer_iv' = {

           U <- sample(x = 0:2,
                       size = N,
                       replace = TRUE,
                       prob = c((1 - p)^2,
                                2 * p * (1 - p),
                                p^2))

           T1 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(U),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           C1 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T1),
                           b0 = b0,
                           b1 = c(ssc),
                           s = s)

           T2 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(U),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           C2 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T2),
                           b0 = b0,
                           b1 = c(ssc),
                           s = s)

           T3 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(U),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T4 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T1),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T5 <- rMParents(N = N,
                           mParents = 3,
                           parentData = list(T1, T2, C1),
                           b0 = b0,
                           b1 = c(ss, ss, ssc),
                           s = s)

           T6 <- rMParents(N = N,
                           mParents = 2,
                           parentData = list(T2, C2),
                           b0 = b0,
                           b1 = c(ss, ssc),
                           s = s)

           T7 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T3),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           return (cbind(U, T1, T2, T3, T4, T5, T6, T7, C1, C2))

         },

         # star ----------------------------------------------------------------

         'star' = {

           U <- sample(x = 0:2,
                       size = N,
                       replace = TRUE,
                       prob = c((1 - p)^2,
                                2 * p * (1 - p),
                                p^2))

           T1 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(U),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T2 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T1),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T3 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T1),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T4 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T1),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           T5 <- rMParents(N = N,
                           mParents = 1,
                           parentData = list(T1),
                           b0 = b0,
                           b1 = c(ss),
                           s = s)

           return (cbind(U, T1, T2, T3, T4, T5))

         })

}
