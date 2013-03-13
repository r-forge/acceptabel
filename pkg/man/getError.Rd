\name{getError}
\alias{getError}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get residuals coming from a given distribution
}
\description{
Get residuals coming from a distribution with a given set of parameters
}
\usage{
getError(H2, freqs, means, ninds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{H2}{
Heritability
}
  \item{freqs}{
Allele frequencies
}
  \item{means}{
Phenotype means per genotype
}
  \item{ninds}{
Number of individuals
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
function (H2, freqs, means, ninds) 
{
    if (length(freqs) < 3) {
        freqs[3] <- 0
    }
    VG <- freqs[1] * (means[1])^2 + freqs[2] * (means[2])^2 + 
        freqs[3] * (means[3])^2
    sigmaSq <- (VG - H2 * VG)/H2
    error <- rnorm(ninds, 0, sqrt(sigmaSq))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
