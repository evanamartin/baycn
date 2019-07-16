# baycn <img src='man/figures/logo.png' align="right" height="139" /></a>

baycn is a Bayesian hybrid approach for inferring a Directed Acyclic Graph (DAG) for continuous, discrete, and mixed data. To speed up the inference of the DAG the algorithm uses a constraint-based method (e.g., the PC algorithm) to reduce the search space; a score-based search method is then used to infer the probability of direction for the edges in the network.

# Installation

To install baycn from github use the following function from the devtools package

``` r
devtools::install_github('Evatar/baycn')
```

# Examples - continuous data

Simulate data under topology GN4.
![gn4](gn4.jpg)

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
