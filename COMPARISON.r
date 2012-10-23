# This function is used for testing the time taken by different timestep functions

rm(list = ls())
Vparam <- read.csv(file="MeanParam.csv",row.names=1); Setting <- "SouthAfrica"
V1950=dget("V1950.Rdata")
source("PARAMv5.r")

myProc = function() {
  DIAG <- 2
  OutMat1 <- OutMatInit <- OutMat
	Vcurrent <- V1950/sum(V1950)*InitPop[Setting]*10^6
	for(t in 1:(62*12)) {
    Out <- timestep(Vcurrent,t,ArtNdCov11,DIAG,OutMat1)
		Vcurrent <- Out$Vnext
    OutMatInit[t,] <- Out$Vout
		if (t == 732) {
      ArtNdCov11 <- Out$Vout["ArtNdCov"]
      break
    }
  }
	V2012 <- Vcurrent
	V2012
}

myTest = function(scriptName = "TIMESTEPv7.r") {
  source(scriptName)
  ptm = proc.time()
  V = as.vector(myProc())
  Time = proc.time()-ptm
  Vgold = as.vector(dget("TimestepV6Output.RData"))
  stopifnot(all.equal(Vgold, V, check.attributes = FALSE))
  Time
}

