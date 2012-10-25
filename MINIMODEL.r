ministep = function(Vcurrent, t, A) {
  # Initialize new vector
  Vnext = Vcurrent
  # Clear out deaths
  Vnext[12] = 0
  # Add new entrants based on birth rate
  Vnext[1] = Vnext[1]	+ NewEntt[t] * 1e6
  # Initialize the rate matrix
  RateMat = matrix(0, 12, 12)
  dimnames(RateMat) = list(StatNam, StatNam)
  # Include deaths
  RateMat[,12] = mubt[t] + alpha * c(0, VmuHIV[1:3], 0, VmuHIV, 0)
  # Calculate guideline readiness
  RateMat[cbind(1:4, 5:8)] = A
  # Calculate HIV status changes
  RateMat[cbind(1:7, 2:8)] = c(rHIVt[t], H1toH2, H2toH3, 0, rHIVt[t], H1toH2, H2toH3)
  # Calculate treatment availability
  RateMat[cbind(6:8, 9:11)] = c(D1t[t], D2t[t], D3t[t])
  # Scale and correct the diagonal entries to make the matrix stochastic
  RateMat	= RateMat/12
  diag(RateMat) = 1 - rowSums(RateMat)
  # Calculate the new vector
  Vnext = Vnext %*% RateMat
  # Calculate the number of treated people
  Treated = sum(Vnext[9:11])
  # Return value
  list(Vector = Vnext, Treated = Treated)
}
