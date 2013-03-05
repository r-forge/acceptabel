\name{getSimulatedPhenotype}
\alias{getSimulatedPhenotype}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getSimulatedPhenotype(data, snp.name, trait.name, covars = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{snp.name}{
%%     ~~Describe \code{snp.name} here~~
}
  \item{trait.name}{
%%     ~~Describe \code{trait.name} here~~
}
  \item{covars}{
%%     ~~Describe \code{covars} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, snp.name, trait.name, covars = "") 
{
    variances <- c(1, 1, 1)
    snp.num <- which(snpnames(data) == snp.name)
    gtypes <- as.numeric(data@gtdata[, snp.num])
    means <- getMeans(data, snp.name, trait.name)
    if (sum(is.na(gtypes)) > 0) {
        warning("Missing genotypes found. Filling with the most frequent allele.")
        gtypes[is.na(gtypes)] <- 0
    }
    freqs <- as.vector(table(gtypes)/length(gtypes))
    phtypes <- rep(0, length(gtypes))
    phtypes[gtypes == 0] <- rnorm(length(gtypes[gtypes == 0]), 
        mean = means[1], sd = sqrt(variances[1]))
    phtypes[gtypes == 1] <- rnorm(length(gtypes[gtypes == 1]), 
        mean = means[2], sd = sqrt(variances[2]))
    phtypes[gtypes == 2] <- rnorm(length(gtypes[gtypes == 2]), 
        mean = means[3], sd = sqrt(variances[3]))
    H2 <- getH2(data, snp.name, trait.name, covars)
    phtypes <- phtypes + getError(H2, freqs, means, length(phtypes))
    phtypes
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
