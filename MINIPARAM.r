  # Initialize the parameter vector
  Vparam = read.csv(file="MeanParam.csv", row.names = 1)
  tempx		<- if(class(Vparam)=="data.frame") {rownames(Vparam)} else {names(Vparam)}
	Vparam 	<- as.vector(as.matrix(Vparam))
	names(Vparam) <- tempx
  
  # External inputs
  NewEnt	 = read.csv("NewEnt.csv")
	BgMort	 = read.csv("BgMort.csv")
	HIVIncid = read.csv("HIVIncid.csv")
  
  # Initial population
	InitPop = c(0.244, 0.436, 0.297, 8.404, 0.155)
  names(InitPop) = c("Botswana","Lesotho","Namibia","SouthAfrica","Swaziland")
  
  # New entrants
	NewEntt = rep(NA, 92 * 12)
	for (i in 0:91) {
    NewEntt[i*12+(1:12)] = seq(NewEnt[i+1,Setting], NewEnt[i+2,Setting], length.out = 13)[-13]
  }
	NewEntt = NewEntt/12		# indexed by t, absolute number of new adult entrants over time
	
	# Mortality rates
	muH1		= exp(Vparam["muH1"])
	muH2		= exp(Vparam["muH2"])
	muH3		= exp(Vparam["muH3"])
	muT1		= exp(Vparam["muT1"])
	muT2		= exp(Vparam["muT2"])
	muT3		= exp(Vparam["muT3"])
	muTBH		= exp(Vparam["muTBH"])
	VmuHIV	= c(muH1, muH2, muH3, muT1, muT2, muT3)  # Note the order difference wrt PARAMv5.r!
	
	Tunmub = 1.0 		# Multiplier to conduct SA on background mortality rates
	mubt = rep(NA, 92 * 12)
	for (i in 0:91) { 
    mubt[i*12+(1:12)] = seq(BgMort[i+1,Setting]*Tunmub, BgMort[i+2,Setting]*Tunmub, length.out = 13)[-13]  
  }
  
  # HIV incidence
	TunrHIV = HIVIncid[63,Setting]/HIVIncid[62,Setting] + exp(Vparam["TunrHIV"]) -1 # Tuning parameter for adjusting HIV incidence post 2011
	HIVIncid2 = HIVIncid[,Setting]
  HIVIncid2[63:93] = HIVIncid[62,Setting]*TunrHIV^(1:31)
	rHIVt	= rep(NA, 92 * 12)
  for (i in 0:91) { 
    rHIVt[(i*12+1):(i*12+12)] = seq(HIVIncid2[i+1], HIVIncid2[i+2], length.out = 13)[-13]  
  }
  
  # HAART uptake rate
	ArtHistt	<- rep(NA,61*12);  for (i in 0:60){ ArtHistt[(i*12+1):(i*12+12)] <- seq(ArtHist[i+1,Setting],ArtHist[i+2,Setting],length.out=13)[-13]  }
	rArtSU	<- exp(Vparam["rArtSU"])   # Annual factor increase in ART volume 
	ARTVolt	<- rep(NA, 31*12)
	for (i in 0:30) { 
    ARTVolt[(i*12+1):(i*12+12)] = seq(ArtHist[62,Setting]*rArtSU^i, ArtHist[62,Setting]*rArtSU^(i+1), length.out = 13)[-13]  
  }
	
	# Transition probabilities
	H1toH2	<- exp(Vparam["H1toH2"]) # Rate of trans from CD4 >350 to CD4 350-200
	H2toH3	<- exp(Vparam["H2toH3"]) # Rate of trans from CD4 350-200 to CD4 200-0
	  
  # Treatment enrolment rates
  D1t = rep(0, 92 * 12)
  D2t = rep(0, 92 * 12)
  D2t[60*12 + (0:24)] = seq(0, 2, length.out = 25)
  D3t = rep(2, 92 * 12)
	
	# Compartment names
	StatNam = c("H-G-", "H1G-", "H2G-", "H3G-", "H-G+", "H1G+", "H2G+", "H3G+", "H1Tx", "H2Tx", "H3Tx", "Dead")
	
	# Adjustment factor to net out HIV mortality
	alpha = 1.25