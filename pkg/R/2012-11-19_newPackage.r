setwd("~/Desktop/AcceptABEL/")
load("data/test_data.rda")

data <- data.orig
snp.name <- "BICF2P176848"
trait.name <- "igaValue"
cov.names <- c("sex", "cad")
n <- 5
q <- getAccScore(data, snp.name, trait.name, cov.names, n=1, progress.bar=T)

# Arabidopsis data
data <- load.gwaa.data(phenofile="~/Research/ArabiStrat/call_method_32/arabidopsis_phenotype_published.csv", genofile="~/Research/ArabiStrat/call_method_32/arabidopsis_genotype_published.raw")
an0 <- qtscore(X43_FLC,data=data)
plot(an0, pch=4, cex=.5)
q <- getAccScore(data=data, snp.name="3_2227823", trait.name="X43_FLC", n=10, cov.names=F)

getAccScore <- function(data, snp.name, trait.name, n=100, progress.bar=T, cov.names=F) {
#  snp.name="3_2227823"; trait.name="X33_avrRpm1"; n=1; cov.names=F; progress.bar = T
  quantiles <- rep(NA, n)
  lambdas <- rep(NA, n)
  lambdas.se <- rep(NA, n)
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
    if (!inherits(an0, 'try-error')) {
      pvals <- an0@results$P1df
      quantile <- rank(pvals)[snp.num]/nsnps(data)
      quantiles[i] <- quantile
      lambdas[i] <- lambda(an0)$estimate
      lambdas.se[i] <- lambda(an0)$se
    }
    if (progress.bar) { setTxtProgressBar(pb, i/n) }
  }
  df <- data.frame(quantile=quantiles, lambda=lambdas, se=lambdas.se)
  score <- exp(1 - mean(lambdas))
  mean.lambda <- mean(lambdas)
  mean.quantile <- mean(quantiles)
  SE <- sqrt(var(lambdas))/sqrt(length(lambdas))
  setClass("result", representation(score="numeric", mean.lambda="numeric", lambda.SE="numeric", mean.quantile="numeric", details="data.frame"))
  result <- new("result", score = score, mean.lambda = mean.lambda, lambda.SE = SE, mean.quantile = mean.quantile, details = df)
  result
}
  
getH2 <- function(data, snp.name, trait.name, covars=F) {
  # Get H2
  snp.num <- which(snpnames(data) == snp.name)
  gtypes <- as.factor(as.character(data@gtdata[,snp.num]))
  phtypes <- data@phdata[,trait.name]
  if (is.logical(covars)) {
    covars <- "1"
    dat <- data.frame(pheno=phtypes, geno=gtypes)
  }
  else {
    dat <- data.frame(pheno=phtypes, geno=gtypes, data@phdata[,covars])
    names(dat) <- c("pheno", "geno", covars)
  }
  # Without SNP
  names.kept <- covars
  lm.form <- paste("pheno ~ ", paste(names.kept, collapse="+"), sep = "")
  lm.obj <- lm(lm.form, data=dat)
  r2.1 <- summary(lm.obj)$r.squared
  # With SNP
  names.kept <- c("geno", covars)
  lm.form <- paste("pheno ~ ", paste(names.kept, collapse="+"), sep = "")
  lm.obj <- lm(lm.form, data=dat)
  r2.2 <- summary(lm.obj)$r.squared
  # Compute h2
  H2 <- abs(r2.2 - r2.1)
  H2
}
 
getMeans <- function(data, snp.name, trait.name) {
  snp.num <- which(snpnames(data) == snp.name)
  gtypes <- as.numeric(data@gtdata[,snp.num])
  m1 <- mean(data@phdata[which(gtypes==0) ,trait.name], na.rm=T)
  m2 <- mean(data@phdata[which(gtypes==1) ,trait.name], na.rm=T)
  m3 <- mean(data@phdata[which(gtypes==2) ,trait.name], na.rm=T)
  result <- c(m1,m2,m3)
  result[is.na(result)] <- 0
  result
}

getVariances <- function(data, snp.name, trait.name) {
  snp.num <- which(snpnames(data) == snp.name)
  gtypes <- as.numeric(data@gtdata[,snp.num])
  v1 <- var(data@phdata[which(gtypes==0) ,trait.name], na.rm=T)
  v2 <- var(data@phdata[which(gtypes==1) ,trait.name], na.rm=T)
  v3 <- var(data@phdata[which(gtypes==2) ,trait.name], na.rm=T)
  result <- c(v1,v2,v3)
  result[is.na(result)] <- 0
  result
}

getError <- function(H2, freqs, means, ninds) {
  VG <- freqs[1] * (means[1])^2 + freqs[2] * (means[2])^2 + freqs[3] * (means[3])^2 
  sigmaSq <- (VG - H2 * VG) / H2
  error <- rnorm(ninds, 0, sqrt(sigmaSq))
  error
}

getSimulatedPhenotype <- function(data, snp.name, trait.name, covars=F, means="actual", variances="actual") {
  if (means == "actual") { means <- getMeans(data, snp.name, trait.name) }
  if (variances == "actual") { variances <- getVariances(data, snp.name, trait.name) }
  snp.num <- which(snpnames(data) == snp.name)
  gtypes <- as.numeric(data@gtdata[ ,snp.num])
  # IMPORTANT!!! Assumes that missing have the most frequent genotype!
  if (sum(is.na(gtypes)) > 0) {
    warning("Missing genotypes found. Filling with the most frequent allele.")
    gtypes[is.na(gtypes)] <- 0
  }
  freqs <- c(length(gtypes[gtypes==0]), length(gtypes[gtypes==1]), length(gtypes[gtypes==2]))/length(gtypes)
  phtypes <- rep(0, length(gtypes))
  phtypes[gtypes==0] <- rnorm(length(gtypes[gtypes==0]), mean=means[1], sd=sqrt(variances[1]))
  phtypes[gtypes==1] <- rnorm(length(gtypes[gtypes==1]), mean=means[2], sd=sqrt(variances[2]))
  phtypes[gtypes==2] <- rnorm(length(gtypes[gtypes==2]), mean=means[3], sd=sqrt(variances[3]))
  H2 <- getH2(data, snp.name, trait.name, covars)
  phtypes <- phtypes + getError(H2, freqs, means, length(phtypes))
  phtypes
}