getError <-
function(H2, freqs, means, ninds) {
  # If only two variants, e.g. inbred lines, make sure the third frequency is 0
  if (length(freqs) < 3) {
    freqs[3] <- 0
  }
  VG <- freqs[1]*(means[1])^2 + freqs[2]*(means[2])^2 + freqs[3]*(means[3])^2 
  sigmaSq <- (VG - H2 * VG) / H2
  error <- rnorm(ninds, 0, sqrt(sigmaSq))
}
