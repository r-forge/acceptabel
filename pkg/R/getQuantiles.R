getQuantiles <-
function(data, snp.name, trait.name, cov.names, n=100, progress.bar=T) {
quantiles <- c()
#### Get Chi2 for your SNP
data.tmp <- data
phdata(data.tmp) <- cbind(phdata(data.tmp), simpheno=rep(0, nids(data)))
snp.num <- which(snpnames(data) == snp.name)
# Get Chi2 values from simulations
if (progress.bar) { pb <- txtProgressBar(style=3) }
for (i in 1:n) {
simulated.phenotype <- getSimulatedPhenotype(data, snp.name, trait.name, covars=cov.names)
phdata(data.tmp)$simpheno <- simulated.phenotype
form <- paste("simpheno ~ ", paste(cov.names, collapse="+"), sep = "")
an0 <- try(qtscore(as.formula(form), data.tmp))
if (typeof(an0) != "try-error") {
pvals <- an0@results$P1df
quantile <- rank(pvals)[snp.num]/nsnps(data)
quantiles <- c(quantiles, quantile)
}
if (progress.bar) { setTxtProgressBar(pb, i/n) }
}
quantiles
}
