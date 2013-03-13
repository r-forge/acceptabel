\name{getMeans}
\alias{getMeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get means
}
\description{
get means
}
\usage{
getMeans(data, snp.name, trait.name)
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
function (data, snp.name, trait.name) 
{
    snp.num <- which(snpnames(data) == snp.name)
    gtypes <- as.numeric(data@gtdata[, snp.num])
    m1 <- mean(na.omit(data@phdata[which(gtypes == 0), trait.name]))
    m2 <- mean(na.omit(data@phdata[which(gtypes == 1), trait.name]))
    m3 <- mean(na.omit(data@phdata[which(gtypes == 2), trait.name]))
    c(m1, m2, m3)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
