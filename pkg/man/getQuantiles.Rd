\name{getQuantiles}
\alias{getQuantiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getQuantiles(data, snp.name, trait.name, cov.names, n = 100, progress.bar = T)
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
  \item{cov.names}{
%%     ~~Describe \code{cov.names} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{progress.bar}{
%%     ~~Describe \code{progress.bar} here~~
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
function (data, snp.name, trait.name, cov.names, n = 100, progress.bar = T) 
{
    quantiles <- c()
    data.tmp <- data
    phdata(data.tmp) <- cbind(phdata(data.tmp), simpheno = rep(0, 
        nids(data)))
    snp.num <- which(snpnames(data) == snp.name)
    if (progress.bar) {
        pb <- txtProgressBar(style = 3)
    }
    for (i in 1:n) {
        simulated.phenotype <- getSimulatedPhenotype(data, snp.name, 
            trait.name, covars = cov.names)
        phdata(data.tmp)$simpheno <- simulated.phenotype
        form <- paste("simpheno ~ ", paste(cov.names, collapse = "+"), 
            sep = "")
        an0 <- try(qtscore(as.formula(form), data.tmp))
        if (typeof(an0) != "try-error") {
            pvals <- an0@results$P1df
            quantile <- rank(pvals)[snp.num]/nsnps(data)
            quantiles <- c(quantiles, quantile)
        }
        if (progress.bar) {
            setTxtProgressBar(pb, i/n)
        }
    }
    quantiles
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
