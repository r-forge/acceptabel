getH2 <-
function(data, snp.name, trait.name, covars="") {
  # Get H2
  snp.num <- which(snpnames(data) == snp.name)
  gtypes <- as.factor(as.character(data@gtdata[,snp.num]))
  phtypes <- data@phdata[,trait.name]
  if (covars[1] == "") {
    covars <- c("1")
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
