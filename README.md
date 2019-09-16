[![Build Status](https://travis-ci.org/evanamartin/baycn.svg?branch=master)](https://travis-ci.org/evanamartin/baycn)

# baycn <img src='man/figures/logo.png' align="right" height="139" /></a>

baycn is a Bayesian hybrid approach for inferring Directed Acyclic Graphs (DAGs) for continuous, discrete, and mixed data. The algorithm can use the graph inferred by another more efficient graph inference method as input; the input graph may contain false edges or undirected edges, but can help reduce the search space to a more manageable size. A Bayesian Markov chain Monte Carlo algorithm is then used to infer the probability of direction and absence for the edges in the network.

## Installation

To install baycn from github use the following function from the devtools package.

``` r
devtools::install_github('evanamartin/baycn')
```

## Examples - continuous data

The following R code simulates and analyzes data from topology GN4 which has four continuous random variables. The following figure shows the true graph for topology GN4.

<img src='man/figures/gn4.jpg' align="center" width="100" />

```r
# Generate data under topology gn4.
data_gn4 <- simdata(b0 = 0,
                    N = 200,
                    s = 1,
                    graph = 'gn4',
                    ss = 1)

# Use the true edges for the input.
am_gn4 <- matrix(c(0, 1, 1, 0,
                   0, 0, 0, 1,
                   0, 0, 0, 0,
                   0, 0, 1, 0),
                 byrow = TRUE,
                 nrow = 4)

mh_gn4 <- mhEdge(data = data_gn4,
                 adjMatrix = am_gn4,
                 burnIn = 0.2,
                 iterations = 1000,
                 nGV = 0,
                 pmr = FALSE,
                 prior = c(0.05,
                           0.05,
                           0.9),
                 thinTo = 200)

summary(mh_gn4)
```

## Examples - mixed data

The following R code simulates and analyzes data from topology M1 which has one discrete random variable U and two continuous random variables. The following figure shows the true graph for topology M1.

<img src='man/figures/m1_gv.jpg' align="center" width="100" />

```r
# Generate data under topology m1_gv.
data_m1 <- simdata(b0 = 0,
                   N = 200,
                   s = 1,
                   graph = 'm1_gv',
                   ss = 1,
                   p = 0.27)

# Use the true edges for the input.
am_m1 <- matrix(c(0, 1, 0,
                  0, 0, 1,
                  0, 0, 0),
                byrow = TRUE,
                nrow = 3)

# Run the Metropolis-Hastings algorithm with the Principle of Mendelian
# Randomization (PMR).
mh_m1_pmr <- mhEdge(data = data_m1,
                    adjMatrix = am_m1,
                    burnIn = 0.2,
                    iterations = 1000,
                    nGV = 1,
                    pmr = TRUE,
                    prior = c(0.05,
                              0.05,
                              0.9),
                    thinTo = 200)

summary(mh_m1_pmr)
```
