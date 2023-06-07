Rcpp::compileAttributes("~/Desktop/RTestDev/TestPKG")

library(devtools)

build("~/Desktop/RTestDev/TestPKG")
install("~/Desktop/RTestDev/TestPKG")
load_all("~/Desktop/RTestDev/TestPKG")

abc <- rep(0, 10)

test_ublas(as.integer(abc))

