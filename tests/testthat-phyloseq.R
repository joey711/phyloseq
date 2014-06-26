library("testthat")
packageVersion("phyloseq")
# As suggested for opt-out option on testing by users, recommended by CRAN
# http://adv-r.had.co.nz/Testing.html
# Previously, best practice was to put all test files in inst/tests and ensure that R CMD check ran them by putting the following code in tests/test-all.R:  
# library(testthat)
# library(yourpackage)
# test_package("yourpackage")
# Now, recommended practice is to put your tests in tests/testthat, and ensure R CMD check runs them by putting the following code in tests/test-all.R:
# library(testthat)
# test_check("yourpackage")
# The advantage of this new structure is that the user has control over whether or not tests are installed using the –install-tests parameter to R CMD install, or INSTALL_opts = c(“–install-tests”) argument to install.packages(). I’m not sure why you wouldn’t want to install the tests, but now you have the flexibility as requested by CRAN maintainers.
test_check("phyloseq")
