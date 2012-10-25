miniSim = function(Setting, Avector, C = 1)  {
    source("MINIPARAM.r")
    source("MINIMODEL.r")
    # Initializing the population - this will need to change!!!
    V1950 = dget("V1950.Rdata")
    Nums = sapply(split(V1950[-505], rep(0:6, each = 72)), sum)
    perm = c(1, 2, 4, 6, 3, 5, 7)
    V1980 = Nums[perm]
    V1980 = c(V1980[1:4], rep(0, 4), V1980[5:7], 0)
  	Vcurrent = (V1980/sum(V1980)) * InitPop[Setting] * 1e6
    # Run the model until the end of 2012
    OutMat = matrix(NA, 32*12, 12)
    Diffs = rep(NA, 32*12)
  	for (t in (1:(32*12))) {
      realT = t + 30 * 12
      Output = ministep(Vcurrent, realT, Avector[t])
      Vnext = Output$Vector
      OutMat[t,] = Vnext
      Diffs[t] = Output$Treated - ARTVolt[t]
      Vcurrent = Vnext
  	}
  	Penalty = sum(Diffs^2, na.rm = TRUE) + C * sum(diff(Avector)^2)
    list(objval = Penalty, diffs = Diffs, output = OutMat)
}