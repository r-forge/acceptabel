\name{AcceptABEL-package}
\alias{AcceptABEL-package}
\alias{AcceptABEL}
\docType{package}
\title{
Estimate false positive rate (FPR) in your GWAS study
}
\description{
Package allows you to estimate FPR in your GWAS study using an approach alternative to randomization tests. Estimation is based upon simulating several artificial phenotypes and running association tests to assess the false positive rate at different p-value thresholds. 
}
\details{
\tabular{ll}{
Package: \tab AcceptABEL\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-02-26\cr
License: \tab GPL (>= 2)\cr
}
Package contains several functions that form a pipeline for estimating FPR. There are three main steps in the FPR analysis:
1/ Preparing simulated phenotypes
2/ Running association tests
3/ Determining FPR for different thresholds
In the example provided below, FPR is determined for default p-values of: 0.001, 0.01 and 0.05. First, an example dataset stored in 
test_data structure is loaded and quality control is run using GenABEL check.marker function. Following these initial preparations, 
artificial phenotypes are created for 0.001 percent of all the markers, association tests are run and, finally, FPR are determined. 
For more details on each of the functions included in the package, see appropriate man page.
}
\author{
Marcin Kierczak, Xia Shen

Maintainer: Marcin Kierczak <marcin.kierczak@slu.se>, Xia Shen <xia.shen@slu.se>
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
# Load GenABEL functions
#require(GenABEL)
# Load example data
#data("test_data")
# Run basic quality control
#qc1 <- check.marker(data=data.orig, callrate=0.95, perid.call=0.95, ibs.exclude="lower", maf=0.01, p.level=1e-10)
#data.qc1 <- data.orig[qc1$idok, qc1$snpok]
# Simulate phenotypes for 0.001% of the markers
#data.simulated <- simulatePhenotypes(data=data.qc1, frac=0.001, H2 = .8, means = c(-1, 0, 1), variances = c(.8, 0, .3), verbose = TRUE)
# Perform association tests and record their outcomes
#data.models <- fitModels(simulatedData = data.simulated, verbose = TRUE)
# Use results from the previous step to determine FPR for different 
# p-value thresholds (by default: 0.001, 0.01 and 0.05)
#fpr <- computeFPR(results = data.models, plot = TRUE, verbose = TRUE)
}
