solveForA = function(numInit, Targets, t) {
  multiplier = D1t[t + 1] * numInit[2] + D2t[t + 1] * numInit[3] + D3t[t + 1] * numInit[4]
  constant1  = rHIVt[t] * D1t[t + 1] * numInit[5] + (12 - (H1toH2 + D1t[t] + mubt[t] + alpha * VmuHIV[1])) * D1t[t  + 1] * numInit[6] + D1t[t] * (12 - (mubt[t + 1] + alpha * VmuHIV[4])) * numInit[6] + (12 - (mubt[t] + alpha * VmuHIV[4])) * (12 - (mubt[t + 1] + alpha * VmuHIV[4])) * numInit[9]
  constant2  = H1toH2   * D2t[t + 1] * numInit[6] + (12 - (H2toH3 + D2t[t] + mubt[t] + alpha * VmuHIV[2])) * D2t[t  + 1] * numInit[7] + D2t[t] * (12 - (mubt[t + 1] + alpha * VmuHIV[5])) * numInit[7] + (12 - (mubt[t] + alpha * VmuHIV[5])) * (12 - (mubt[t + 1] + alpha * VmuHIV[5])) * numInit[10]
  constant3  = H2toH3   * D3t[t + 1] * numInit[7] + (12 - (         D3t[t] + mubt[t] + alpha * VmuHIV[3])) * D3t[t  + 1] * numInit[8] + D3t[t] * (12 - (mubt[t + 1] + alpha * VmuHIV[6])) * numInit[8] + (12 - (mubt[t] + alpha * VmuHIV[6])) * (12 - (mubt[t + 1] + alpha * VmuHIV[6])) * numInit[11]
  constant   = constant1 + constant2 + constant3
  A = (Targets[t + 2] - constant)/multiplier
  A
}

fullSolve = function(Setting, Targets) {
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
    	curA = solveForA(Vcurrent, Targets, realT)
      Output = ministep(Vcurrent, realT, curA)
    }
    
}