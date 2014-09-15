## bare-bones/skeleton testing
library(betararef)
library(testthat)

set.seed(101)
b1 <- betasim()
expect_equal(dim(b1),c(3,30))
d1 <- calcbeta(b1)
expect_equal(var(d1),0)
expect_equal(mean(d1),sqrt(3)/3)
d2 <- calcbeta(b1,method="manhattan")
expect_equal(mean(d2),3.2974,tol=1e-5)
d3 <- calcbeta(b1,method="manhattan",trans="identity")
expect_equal(mean(d3),115.47,tol=1e-5)

rr <- betararef:::rtrpois(100000,lambda=4.5)
expect_equal(mean(rr),4.5,tol=5e-3)
