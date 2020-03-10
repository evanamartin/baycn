This release fixes several bugs in the cycle finding and removing functions,
updates unit tests, adds the geuvadis data set, and includes the ability to
analyze binomial data.

-------

## Test environments

* local OS X install, R 3.6.1
* local Windows install, R 3.6.1
* Ubuntu 16.04.6 (on travis-ci), R 3.6.2

## R CMD check results

0 errors | 0 warnings | 0 notes

## Note to CRAN maintainers

In response to an email from Martin Maechler I have changed all instances of c()
to numeric() with the numeric vector initialized to full length.
