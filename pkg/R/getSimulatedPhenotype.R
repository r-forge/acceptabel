getSimulatedPhenotype <-
function(data, snp.name, trait.name, covars="") {
  variances <- c(1,1,1)
  snp.num <- which(snpnames(data) == snp.name)
  gtypes <- as.numeric(data@gtdata[ ,snp.num])
  means <- getMeans(data, snp.name, trait.name)
  # IMPORTANT!!! Assumes that missing have the most frequent genotype!
  if (sum(is.na(gtypes)) > 0) {
    warning("Missing genotypes found. Filling with the most frequent allele.")
    gtypes[is.na(gtypes)] <- 0
  }
  freqs <- as.vector(table(gtypes)/length(gtypes))
  phtypes <- rep(0, length(gtypes))
  phtypes[gtypes==0] <- rnorm(length(gtypes[gtypes==0]), mean=means[1], sd=sqrt(variances[1]))
  phtypes[gtypes==1] <- rnorm(length(gtypes[gtypes==1]), mean=means[2], sd=sqrt(variances[2]))
  phtypes[gtypes==2] <- rnorm(length(gtypes[gtypes==2]), mean=means[3], sd=sqrt(variances[3]))
  H2 <- getH2(data, snp.name, trait.name, covars)
  phtypes <- phtypes + getError(H2, freqs, means, length(phtypes))
  phtypes
}
