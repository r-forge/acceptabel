getMeans <-
function(data, snp.name, trait.name) {
  snp.num <- which(snpnames(data) == snp.name)
  gtypes <- as.numeric(data@gtdata[,snp.num])
  m1 <- mean(na.omit(data@phdata[which(gtypes==0) ,trait.name]))
  m2 <- mean(na.omit(data@phdata[which(gtypes==1) ,trait.name]))
  m3 <- mean(na.omit(data@phdata[which(gtypes==2) ,trait.name]))
  c(m1,m2,m3)
}
