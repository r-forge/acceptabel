\name{getH2}
\alias{getH2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getH2}
\description{
Function for estimating broad-sense heritability.
}
\usage{
getH2(data, snp.name, trait.name, covars = "")
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
\author{Marcin Kierczak, Xia Shen}
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
    snp.num <- which(snpnames(data) == snp.name)
    gtypes <- as.factor(as.character(data@gtdata[, snp.num]))
    phtypes <- data@phdata[, trait.name]
    if (covars[1] == "") {
        covars <- c("1")
        dat <- data.frame(pheno = phtypes, geno = gtypes)
    }
    else {
        dat <- data.frame(pheno = phtypes, geno = gtypes, data@phdata[, 
            covars])
        names(dat) <- c("pheno", "geno", covars)
    }
    names.kept <- covars
    lm.form <- paste("pheno ~ ", paste(names.kept, collapse = "+"), 
        sep = "")
    lm.obj <- lm(lm.form, data = dat)
    r2.1 <- summary(lm.obj)$r.squared
    names.kept <- c("geno", covars)
    lm.form <- paste("pheno ~ ", paste(names.kept, collapse = "+"), 
        sep = "")
    lm.obj <- lm(lm.form, data = dat)
    r2.2 <- summary(lm.obj)$r.squared
    H2 <- abs(r2.2 - r2.1)
    H2
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
