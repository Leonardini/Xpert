##########   TB XPERT DIAGNOSTIC MODEL 2012  ##########   
##########            PARAM FILE             ##########
## NOTES
# a.	This gets run once per simulation, and sets up all the data structures for the simulation.
# b.	Takes as its arguments (i) model parameter inputs, either point estimates a random parameter set, (ii) indicators for the option above, and if a random parameter set , the index number of the parameter set
# c.	Defines indices for pulling out different health states
# d.	Loads in external data, incl. starting vector
# e.	Retransform all variables
# f.	Creates static model matrix
# g.	Creates other structures used in TIMESTEP

## A function for retransforming logits
	logitT <- function(x) { exp(x)/(1+exp(x)) }

## SOME VECTORS USED TO PULL OUT PARTS OF MATRIX:
	Vtemp1	<- rep(0:4*7,14)+rep(0:13*36,each=5)				 	# Pulls out basic states in state vector
	Vtemp2	<- rep(2:8,14)+rep(0:13*36,each=7)						# Pulls out one Strain from RateMat
	Vtemp4 	<- Vtemp1[rep(1:5,7)+rep(0:6*10,each=5)]					# Pulls out tx naive basic states from RateMat
	Vtemp6	<- rep(rep(5:8,5)+rep(0:4*7,each=4),14)+rep(0:13*36,each=20)	# Pulls out tx states from RateMat
	Vtemp7	<- rep(rep(3:8,5)+rep(0:4*7,each=6),14)+rep(0:13*36,each=30)	# Pulls out Active disease states from RateMat
	Vtemp8 	<- rep(c(1,2+0:4*7),14)+rep(0:13*36,each=6)				# Pulls out Susceptible and Latent states from RateMat
	Vtemp9	<- rep(rep(3:4,5)+rep(0:4*7,each=2),14)+rep(0:13*36,each=10) 	# Pulls out untreated active disease from RateMat

######### PARAMETER DEFINITIONS ######### 
	tempx		<- if(class(Vparam)=="data.frame") {rownames(Vparam)} else {names(Vparam)}
	Vparam 	<- as.vector(as.matrix(Vparam))
	names(Vparam) <- tempx
	DIAG <- 2  

## EXT INPUTs
	NewEnt	<- read.csv("NewEnt.csv")
	BgMort	<- read.csv("BgMort.csv")
	HIVIncid	<- read.csv("HIVIncid.csv")
	DotsAR	<- read.csv("DotsAR.csv")
	ArtHist	<- read.csv("ArtHist.csv")	# Only goes through 2010
	TxSuc		<- read.csv("TxSuc.csv")	
	TxDef		<- read.csv("TxDef.csv")

### SETTING-SPECIFIC PARAMETERS

## NEW ENTRANTS & INIT POP
	InitPop <- c(0.244,0.436,0.297,8.404,0.155); names(InitPop) <- c("Botswana","Lesotho","Namibia","SouthAfrica","Swaziland")
	NewEntt 	<- rep(NA,92*12)
				for (i in 0:91){ NewEntt[(i*12+1):(i*12+12)] <- seq(NewEnt[i+1,Setting],NewEnt[i+2,Setting],length.out=13)[-13]  }
	NewEntt 	<- NewEntt/12		# indexed by t, abs number of new adult entrants over time 

## TRANSMISSION PARAMETER (per year, by basic state category)
	CR1950init  <- c(11,13,15,15,13); names(CR1950init) <- c("Botswana","Lesotho","Namibia","SouthAfrica","Swaziland")
	CR1950 	<- exp(Vparam["CR1950"])
	CRr		<- exp(Vparam["CRr"])/12 	
	CRt		<- CR1950*exp(-CRr*1:(92*12))

## CURRENT RATE OF ATTENDING HEALTH CENTER, fulfilling criteria for TB suspect, and getting tested 
	NDTst1990 	<- exp(Vparam["NDTunTst1990"])  
	NDTst2010 	<- exp(Vparam["NDTunTst2010"])  
	DTst1990 	<- DotsAR[43,Setting]*exp(Vparam["TunTst1990"])  
	DTst2010 	<- DotsAR[62,Setting]*exp(Vparam["TunTst2010"])

	Dtest 	<-  c(rep(0,37),seq(0,DTst1990,length.out=6),seq(DTst1990,DTst2010,length.out=19)[-1],rep(DTst2010,32))
	NDtest 	<-  c(rep(0,20),seq(0,NDTst1990,length.out=18)[-18],seq(NDTst1990,NDTst2010,length.out=6),rep(NDTst2010,50))

	NDTestt <- DTestt <- rep(NA,92*12)
		for (i in 0:91){ 	 DTestt[(i*12+1):(i*12+12)] <- seq(Dtest[i+1],  Dtest[i+2],length.out=13)[-13] 
 					NDTestt[(i*12+1):(i*12+12)] <- seq(NDtest[i+1],NDtest[i+2],length.out=13)[-13]   }

## CURRENT RATE OF ATTENDING HEALTH CENTER, fulfilling criteria for TB suspect, and getting tested 
	rTstIn	<- 1		# Rate of HS attendence for In, relative to Ip
	rTstSL	<- exp(Vparam[paste("rTstSL",Setting,sep="")])	# CHANGE? # Rate of HS attendence for Su and Ls, relative to Ip
	Vtestfreq 	<- rep(c(rTstSL,rep(c(rTstSL,rTstIn,1,rep(0,4)),5)),14)  # vector of testing rates, relative to Ip rate

## Probability of Culture and DST in base-case DOTS algorithms
	pDstU		<- 0									# P(DST) in positive TB diagnoses, treatment naive
	pDstR		<- logitT(Vparam[paste("pDstR",Setting,sep="")])	# P(DST) in positive TB diagnoses, treatment experienced
	pCulU		<- logitT(Vparam[paste("pCulU",Setting,sep="")])	# P(DST) in positive TB diagnoses, treatment naive
	pCulR		<- logitT(Vparam[paste("pCulR",Setting,sep="")])	# P(DST) in positive TB diagnoses, treatment experienced

## LOSS TO FOLLOW-UP BETWEEN TEST AND TX, by short/long delay
	pLtfuS	<- logitT(Vparam["pLtfuS"])	# Short delay
	pLtfuL	<- logitT(Vparam["pLtfuL"])	# Long delay

## TREATMENT PARAMETERS
	TunTxSuc	<- exp(Vparam["TunTxSuc"]) 	
	TunTxDef	<- exp(Vparam["TunTxDef"]) 	
	pReTx		<- logitT(Vparam["pReTx"])	# Probability Tx failure identified, patient initiated on tx experienced reg (may be same)
	TunTxMort	<- exp(Vparam["TunTxMort"])	# Multiplier to adjust mortality rates while on treatment into reasonable range (based on observed data)
	TxEft	<- rep(NA,92*12)
		for (i in 0:91){ TxEft[(i*12+1):(i*12+12)] <- seq(TxSuc[i+1,Setting],TxSuc[i+2,Setting],length.out=13)[-13]  }
					# indexed by t, treatment success (cures/(cure+default)) for 1st line regimen with pansensitive disease
	TxEft <- logitT(log(TxEft/(1-TxEft))*TunTxSuc)
	pDefND	<- logitT(Vparam["pDefND"])/(1-logitT(Vparam["pDefND"]))
	pDeft		<- rep(NA,92*12)
			for (i in 0:91){ pDeft[(i*12+1):(i*12+12)] <- seq(TxDef[i+1,Setting],TxDef[i+2,Setting],length.out=13)[-13]  }
					# indexed by t, treatment success (cures/(cure+default)) for 1st line regimen with pansensitive disease
	pDeft 	<- logitT(log(pDeft/(1-pDeft))*TunTxDef)
	pDeft		<- pDeft/(1-pDeft)	# To properly calculate exit rate for default

## HIV INCIDENCE, HAART UPTAKE RATE
	TunrHIV 	<- HIVIncid[63,Setting]/HIVIncid[62,Setting]+exp(Vparam["TunrHIV"])-1 # Tuning parameter for adjusting HIV incidence post 2011
	HIVIncid2 	<- HIVIncid[,Setting]; HIVIncid2[63:93] <- HIVIncid[62,Setting]*TunrHIV^(1:31)
	rHIVt		<- rep(NA,92*12);  for (i in 0:91){ rHIVt[(i*12+1):(i*12+12)] <- seq(HIVIncid2[i+1],HIVIncid2 [i+2],length.out=13)[-13]  }
	ArtHistt	<- rep(NA,61*12);  for (i in 0:60){ ArtHistt[(i*12+1):(i*12+12)] <- seq(ArtHist[i+1,Setting],ArtHist[i+2,Setting],length.out=13)[-13]  }

	ARTConstr	<- 0  # Indicator = 1 if future ART demand constrained (by below), =0 if 80% (universal) coverage, = 2 if Xpert constrained at Status-quo
	ArtNdCov11	<- 0 # Just initializing this object
	ArtFutCov	<- logitT(Vparam["ArtFutCov"])
	rArtSU	<- exp(Vparam["rArtSU"])   # Annual factor increase in ART volume 
	ARTVolt	<- rep(NA,31*12)
		for (i in 0:30){ ARTVolt[(i*12+1):(i*12+12)] <- seq(ArtHist[62,Setting]*rArtSU^i,ArtHist[62,Setting]*rArtSU^(i+1),length.out=13)[-13]  }

	PriCD4200t	<- c(rep(1,60*12),seq(1,0,length.out=5*12),rep(0,27*12))   
			# Factor for transitioning between CD4<200 proritization and equal priority in ART enrollment
 
	H1toH2	<- exp(Vparam["H1toH2"]) # Rate of trans from CD4 >350 to CD4 350-200
	H2toH3	<- exp(Vparam["H2toH3"]) # Rate of trans from CD4 350-200 to CD4 200-0

## UNIT COSTS
	CXpt		<- 0									# Per test Xpert cost
	CSmr		<- exp(Vparam[paste("CSmr",Setting,sep="")])		# Per test Smear cost 
	CXry		<- exp(Vparam[paste("CXry",Setting,sep="")])		# Per test Xray cost
	CCltr		<- exp(Vparam[paste("CCltr",Setting,sep="")])		# Per test Culture cost
	CDst		<- exp(Vparam[paste("CDst",Setting,sep="")])		# Per test DST cost
	CVisS		<- exp(Vparam[paste("CVisS",Setting,sep="")])		# Per visit Clinic visit cost (short, tx check-in)
	CVisL		<- exp(Vparam[paste("CVisL",Setting,sep="")])		# Per visit Clinic visit cost (long, diagnosis)
	CArt		<- exp(Vparam[paste("CArt",Setting,sep="")])		# Per month HAART cost
	CIpd		<- exp(Vparam[paste("CIpd",Setting,sep="")])		# Per day inpatient cost	 #### NEW ###
	TunRegCost	<- exp(Vparam["TunRegCost"])					# Multiplier on regimen costs

### GLOBAL PARAMETERS

## MORTALITY RATES
	Tunmub	<- 1.0 		# Multiplier to conduct SA on background mortality rates
	mubt		<-rep(NA,92*12)
				for (i in 0:91){ mubt[(i*12+1):(i*12+12)] <- seq(BgMort[i+1,Setting]*Tunmub,BgMort[i+2,Setting]*Tunmub,length.out=13)[-13]  }
						# indexed by t, background mortality rate over time 
	muIn		<- exp(Vparam["muIn"])		# excess mort for smear-neg active disease  
	muIp		<- exp(Vparam["muIp"])		# excess mort for smear-pos active disease    
	muH1		<- exp(Vparam["muH1"])
	muH2		<- exp(Vparam["muH2"])
	muH3		<- exp(Vparam["muH3"])
	muT1		<- exp(Vparam["muT1"])
	muT2		<- exp(Vparam["muT2"])
	muT3		<- exp(Vparam["muT3"])
	muTBH		<- exp(Vparam["muTBH"])
	VmuHIV	<- c(muH1,muT1,muH2,muT2,muH3,muT3)

## STRAIN FITNESS (with pansensitive=1.0) N.B. applied to contact rates
	TunRelFit	<- logitT(Vparam["TunRelFit"])*3.7	# Tuning parameter for RelFit => acts as a multiplier on the absolute difference from 1. Valid range: 0 - 3.7 (=1/(1-0.73))
	RelFit	<- 1-(1-c(1,0.95,0.85,0.73,0.73))*TunRelFit  

## ACQUIRED RESISTANCE TUNING
	TunAR		<- exp(Vparam["TunAR"])	# Multiplier applied to rates of acquired resistance

## TRANSMISSION PARAMETER (per year, by basic state category)
	TrIn		<- exp(Vparam["TrIn"])	# Contact rate for In as a fraction of Ip

## TUNING PARAMETER for HAART Parameter Values (0=same as HIV-pos no HAART, 1 = same as HIV-neg)
	TunHAART	<- logitT(Vparam["TunHAART"])

## PROBABILITY OF FAST BREAKDOWN, by HIV status
	pfastN	<- logitT(Vparam["pfastN"])
	pfastH1	<- logitT(Vparam["pfastH1"])
	pfastH2	<- logitT(Vparam["pfastH2"])
	pfastH3	<- logitT(Vparam["pfastH3"])
	pfastT1	<- (1-TunHAART)*pfastH1 + TunHAART*pfastN
	pfastT2	<- (1-TunHAART)*pfastH2 + TunHAART*pfastN
	pfastT3	<- (1-TunHAART)*pfastH3 + TunHAART*pfastN
	Vpfast	<- c(pfastN,pfastH1,pfastT1,pfastH2,pfastT2,pfastH3,pfastT3)

## PROBABILITY OF GOING INTO In FOLLOWING BREAKDOWN, by HIV status
	pToIpN	<- logitT(Vparam["pToIpN"])
	pToIpH1	<- logitT(Vparam["pToIpH1"])
	pToIpH2	<- logitT(Vparam["pToIpH2"])
	pToIpH3	<- logitT(Vparam["pToIpH3"])
	pToIpT1	<- (1-TunHAART)*pToIpH1 + TunHAART*pToIpN
	pToIpT2	<- (1-TunHAART)*pToIpH2 + TunHAART*pToIpN
	pToIpT3	<- (1-TunHAART)*pToIpH3 + TunHAART*pToIpN
	VpToIp	<- c(pToIpN,pToIpH1,pToIpT1,pToIpH2,pToIpT2,pToIpH3,pToIpT3)

## Rate OF BREAKDOWN FROM Ls (Ls -> In , Ip), by HIV status
	rBreakDN	<- exp(Vparam["rBreakDN"])
	rBreakDH1	<- exp(Vparam["rBreakDH1"])
	rBreakDH2	<- exp(Vparam["rBreakDH2"])
	rBreakDH3	<- exp(Vparam["rBreakDH3"])
	rBreakDT1	<- (1-TunHAART)*rBreakDH1 + TunHAART*rBreakDN
	rBreakDT2	<- (1-TunHAART)*rBreakDH2 + TunHAART*rBreakDN
	rBreakDT3	<- (1-TunHAART)*rBreakDH3 + TunHAART*rBreakDN
	VrBreakD	<- c(rBreakDN,rBreakDH1,rBreakDT1,rBreakDH2,rBreakDT2,rBreakDH3,rBreakDT3)

## RATE OF CONVERSION FROM In TO Ip
	rNtoP		<- exp(Vparam["rNtoP"])	

## RATE OF SPONTANEOUS CURE, by smear status and HIV status
	rIToLsN	<- exp(Vparam["rIToLsN"])	
	rIToLsH1	<- exp(Vparam["rIToLsH1"])
	rIToLsH2	<- exp(Vparam["rIToLsH2"])
	rIToLsH3	<- exp(Vparam["rIToLsH3"])
	rIToLsT1	<- (1-TunHAART)*rIToLsH1 + TunHAART*rIToLsN
	rIToLsT2	<- (1-TunHAART)*rIToLsH2 + TunHAART*rIToLsN
	rIToLsT3	<- (1-TunHAART)*rIToLsH3 + TunHAART*rIToLsN
	VrIToLs	<- c(rIToLsN,rIToLsH1,rIToLsT1,rIToLsH2,rIToLsT2,rIToLsH3,rIToLsT3)

## PARTIAL IMMUNITY afforded by current/prior infection (a value of 1 = full immunity)
	PartImN	<- logitT(Vparam["PartImN"])	
	PartImH1	<- logitT(Vparam["PartImH1"])	
	PartImH2	<- logitT(Vparam["PartImH2"])	
	PartImH3	<- logitT(Vparam["PartImH3"])	
	PartImT1	<- (1-TunHAART)*PartImH1 + TunHAART*PartImN
	PartImT2	<- (1-TunHAART)*PartImH2 + TunHAART*PartImN
	PartImT3	<- (1-TunHAART)*PartImH3 + TunHAART*PartImN
	VPartIm	<- c(PartImN,PartImH1,PartImT1,PartImH2,PartImT2,PartImH3,PartImT3)

####################################################################

## TEST CHARACTERISTICS
	SensSpIn	<- 0.0  	# Sensitivity of Sputum smear in Sputum negative TB suspects
	SensSpIp	<- 1.0  	# Sensitivity of Sputum smear in Sputum positive TB suspects
	SpecSpSL	<- logitT(Vparam["SpecSpSL"])		# Specificity of Sputum smear in Su & Ls

	SensXrIn	<- logitT(Vparam["SensXrIn"])		# Sensitivity of chest xray in Sputum negative TB suspects
	SensXrIp	<- logitT(Vparam["SensXrIp"])		# Sensitivity of chest xray in Sputum positive TB suspects
	SpecXrSL	<- logitT(Vparam["SpecXrSL"])		# Specificity of chest xray in Su & Ls

	SensCuIn	<- logitT(Vparam["SensCuIn"])  	# Sensitivity of culture in TB suspects
	SensCuIp	<- logitT(Vparam["SensCuIp"])		# Sensitivity of culture in Sputum positive TB suspects
	SpecCuSL	<- logitT(Vparam["SpecCuSL"])		# Specificity of culture in Su & Ls
	
	SensXpIn	<- logitT(Vparam["SensXpIn"])		# Sensitivity of Xpert in Sputum negative TB suspects
	SensXpIp	<- logitT(Vparam["SensXpIp"])		# Sensitivity of Xpert in Sputum positive TB suspects
	SpecXpSL	<- logitT(Vparam["SpecXpSL"])		# Specificity of Xpert in Su & Ls

	SensXpRIF	<- logitT(Vparam["SensXpRIF"])	# Sensitivity of Xpert for RIF resistance
	SpecXpRIF	<- logitT(Vparam["SpecXpRIF"])	# Specificity of Xpert for RIF resistance

## REGIMEN PARAMETERS (similar yet different!)
 # Regimen Duration (no uncertainty?)
	Dur1st <- 6;  DurInh <- 9;  DurRif <- 18;  DurMdr <- 21; DurXdr <- 21; DurND <- 18
	VDur	<- c(Dur1st,DurInh,DurRif,DurMdr,DurXdr,DurND)

 # Tx Efficacy
	TxEf1		<- logitT(Vparam["TxEf1"]) 
	TxEf2		<- logitT(Vparam["TxEf2"]) 
	TxEf3		<- logitT(Vparam["TxEf3"]) 
	TxEf4		<- logitT(Vparam["TxEf4"]) 
	TxEf5		<- logitT(Vparam["TxEf5"]) 
	TxEf6		<- logitT(Vparam["TxEf6"]) 
	TxEfMat 	<- rbind(c(1    ,TxEf1,TxEf1,TxEf2,TxEf2),
				   c(TxEf3,TxEf3,TxEf1,TxEf4,TxEf4), 
				   c(TxEf3,TxEf1,TxEf3,TxEf4,TxEf4), 
				   c(TxEf3,TxEf3,TxEf3,TxEf3,TxEf4), 
				   c(TxEf3,TxEf3,TxEf3,TxEf3,TxEf3), 
				   c(TxEf5,TxEf5,TxEf5,TxEf6,TxEf6)) 

 # Rate of acquired resistance
	ARr1		<- exp(Vparam["ARr1"]) 
	ARr2		<- exp(Vparam["ARr2"]) 
	ARr3		<- exp(Vparam["ARr3"]) 
	ARr4		<- exp(Vparam["ARr4"]) 
	ARr5		<- exp(Vparam["ARr5"]) 
	ARr6		<- exp(Vparam["ARr6"]) 
	ARr7		<- exp(Vparam["ARr7"]) 
	ARrND		<- exp(Vparam["ARrND"]) 

	ARMat 	<- rbind(c(ARr1,ARr2,ARr3,ARr5,ARr5,0   ),
				   c(0   ,ARr2,0   ,ARr4,0   ,ARr7), 
				   c(ARr1,0   ,0   ,0   ,ARr4,ARr7), 
				   c(0   ,0   ,0   ,0   ,0   ,ARr6), 
				   c(0   ,0   ,0   ,0   ,0   ,ARr6), 
				   c(ARr1*ARrND,ARr2*ARrND,ARr3*ARrND,ARr5,ARr5,ARr7))*TunAR

# Regimen Costs and other features
	Creg1		<- exp(Vparam["Creg1"]) 
	Creg2		<- exp(Vparam["Creg2"]) 
	Creg3		<- exp(Vparam["Creg3"]) 
	Creg4		<- exp(Vparam["Creg4"]) 
	Creg5		<- exp(Vparam["Creg5"]) 
	Nvis1L	<- exp(Vparam["Nvis1L"]) 	
	Nvis2L	<- exp(Vparam["Nvis2L"]) 
	Nsmr1L	<- exp(Vparam["Nsmr1L"]) 
	Nsmr2L	<- exp(Vparam["Nsmr2L"]) 
	Ncul1L	<- exp(Vparam["Ncul1L"]) 
	Ncul2L	<- exp(Vparam["Ncul2L"]) 
	Nxry1L	<- exp(Vparam["Nxry1L"]) 
	Nxry2L	<- exp(Vparam["Nxry2L"]) 
	NipdMDR	<- exp(Vparam["NipdMDR"]) 
	OthMat 	<- rbind(c(Creg1,Nvis1L,Nsmr1L,Ncul1L,Nxry1L,0),
				   c(Creg2,Nvis2L,Nsmr2L,Ncul2L,Nxry2L,0), 
				   c(Creg3,Nvis2L,Nsmr2L,Ncul2L,Nxry2L,0), 
				   c(Creg4,Nvis2L,Nsmr2L,Ncul2L,Nxry2L,NipdMDR), 
				   c(Creg5,Nvis2L,Nsmr2L,Ncul2L,Nxry2L,NipdMDR), 
				   c(Creg1,Nvis1L,Nsmr1L,Ncul1L,Nxry1L,0)) 
	OthMat[,1]	<- OthMat[,1]*TunRegCost	

## DALY WEIGHTS
	DwtTB		<- logitT(Vparam["DwtTB"])		# DALY Weight for active TB
	DwtH1		<- logitT(Vparam["DwtH1"])
	DwtH2		<- logitT(Vparam["DwtH2"])
	DwtH3		<- logitT(Vparam["DwtH3"])
	DwtT1		<- logitT(Vparam["DwtT1"])
	DwtT2		<- logitT(Vparam["DwtT2"])
	DwtT3		<- logitT(Vparam["DwtT3"])
	VDwt 		<- rep(0,504);  VDwt[Vtemp7] <- DwtTB
	VDwt 		<- 1-(1-VDwt)*(1-rep(c(0,DwtH1,DwtH2,DwtH3,DwtT1,DwtT2,DwtT3),each=72))

## OTHER
	dr		<- 0.03 # exp(Vparam["dr"])		  	# Discount rate
	Drt <- rep(1,30*12); for(i in 1:length(Drt)) { Drt[i] <- (1+dr)^(-(i-1)/12) }

	PhaseIn1 	<- c(rep(0,41*12-1),seq(0,1,length.out=20*12),rep(1,31*12+1))	
					# Vector for phase-in from primitive diagnostics to current diagnostics, linear phase in over past 20 years (1991-2011)
	PIPrd		<- 36		# No months phase in period for new diagnostic technology
	PhaseIn2 	<- c(rep(0,62*12-1),seq(0,1,length.out=PIPrd),rep(1,30*12-PIPrd+1))	
					# Vector for implementing phase-in of new algorithm

#######################################################
####### DIAGNOSTIC ALGORITHMS ####################################
#######################################################
	Nalg <- 3 # no. algorithms

# Non-DOTS Algorithm (constant across all strategies)

	RegMat0 <- matrix(NA,6,504); rownames(RegMat0) <- c("1stL","2ndLMonoINH","2ndLMonoRIF","2ndLMDR","2ndLXDR","NonDots")
	RegMat0[,c(Vtemp1+6,Vtemp1+8)] <- c(0,0,0,0,0,1)

	TxMatAlg0 <- matrix(0,14,504);	rownames(TxMatAlg0) <- c("Duration","TxEff","AR1-2","AR1-3","AR1-4","AR2-4",
		"AR3-4","AR4-5","RegCost","ClinVis","Smears","Cultures","Xrays","Inpt Days")

	# Vector of all Active disease states
	TruPosND	<- rep(c(rep(c(SensSpIn*(1-pLtfuS),SensSpIp*(1-pLtfuS)),5),rep(c(SensSpIn*(1-pLtfuS),SensSpIp*(1-pLtfuS)),5)),7)
											# True positive diagnoses starting treatment
	TruPosNDB	<- rep(c(rep(c(SensSpIn,SensSpIp),5),rep(c(SensSpIn,SensSpIp),5)),7)
											# True positive diagnoses identified (whether started treatment or not)
	# Vector of all susceptible and latent states
	FalsPosND	<- rep(c(rep((1-SpecSpSL)*(1-pLtfuS),6),rep((1-SpecSpSL)*(1-pLtfuS),6)),7) 
											# False positive diagnoses starting treatment
	FalsPosNDB	<- rep(c(rep((1-SpecSpSL),6),rep((1-SpecSpSL),6)),7) 
											# False positive diagnoses identified (whether started treatment or not)

	CTstSnU		<- CSmr+CVisL+CVisS
	CTstSnR		<- CSmr+CVisL+CVisS
	CTstSpU		<- CSmr+CVisL+CVisS	
	CTstSpR		<- CSmr+CVisL+CVisS		
	CTstSnUH		<- CSmr+CVisL+CVisS		
	CTstSnRH		<- CSmr+CVisL+CVisS
	CTstSpUH		<- CSmr+CVisL+CVisS
	CTstSpRH		<- CSmr+CVisL+CVisS
	CTstSnUT		<- CSmr+CVisL+CVisS		
	CTstSnRT		<- CSmr+CVisL+CVisS	
	CTstSpUT		<- CSmr+CVisL+CVisS
	CTstSpRT		<- CSmr+CVisL+CVisS
	VTestCostND		<- c(c(CTstSnU,rep(c(CTstSnU,CTstSnU,CTstSpU,rep(0,4)),5)),c(CTstSnR,rep(c(CTstSnR,CTstSnR,CTstSpR,rep(0,4)),5)),
					rep(c(c(CTstSnUH,rep(c(CTstSnUH,CTstSnUH,CTstSpUH,rep(0,4)),5)),c(CTstSnRH,rep(c(CTstSnRH,CTstSnRH,CTstSpRH,rep(0,4)),5)),
					c(CTstSnUT,rep(c(CTstSnUT,CTstSnUT,CTstSpUT,rep(0,4)),5)),c(CTstSnRT,rep(c(CTstSnRT,CTstSnRT,CTstSpRT,rep(0,4)),5))),3))

# Early DOTS Algorithm (all DOTS onto 1st-line regimen)
	RegMat1 <- RegMat0
	RegMat1[,c(Vtemp1+5,Vtemp1+7)] <- c(1,0,0,0,0,0)

	TxMatAlg1	<- TxMatAlg0
	TxMatAlg1[1,]	<- VDur%*%RegMat1
	for(i in 0:4) { TxMatAlg1[2,Vtemp2+7*i]  <-  TxEfMat[,i+1]%*%RegMat1[,Vtemp2+7*i]  }
	for(i in 1:6) { TxMatAlg1[i+2,Vtemp2+7*c(0,0,0,1,2,3)[i]]  <-	ARMat[,i]%*%RegMat1[,Vtemp2+7*c(0,0,0,1,2,3)[i]]  }
	for(i in 1:6) { TxMatAlg1[i+8,]  <-  OthMat[,i]%*%RegMat1  }

	TruPosDAlg1		<- rep(c(rep(c(SensSpIn*(1-pLtfuS),SensSpIp*(1-pLtfuS)),5),rep(c(SensSpIn*(1-pLtfuS),SensSpIp*(1-pLtfuS)),5)),7)
										# True positive diagnoses starting treatment
	TruPosDAlgB1	<- rep(c(rep(c(SensSpIn,SensSpIp),5),rep(c(SensSpIn,SensSpIp),5)),7)
										# True positive diagnoses identified (whether started treatment or not)
	FalsPosDAlg1	<- rep(c(rep((1-SpecSpSL)*(1-pLtfuS),6),rep((1-SpecSpSL)*(1-pLtfuS),6)),7) 
										# False positive diagnoses starting treatment
	FalsPosDAlgB1	<- rep(c(rep((1-SpecSpSL),6),rep((1-SpecSpSL),6)),7) 
										# False positive diagnoses identified (whether started treatment or not)
	# Vector of all states
	GetXpt1		<- rep(0,504)
	VTxCost1		<- t(t(TxMatAlg1[10:14,])%*%c(CVisS,CSmr,CCltr,CXry,CIpd))
	CTstSnU		<- CSmr+CVisL+CVisS
	CTstSnR		<- CSmr+CVisL+CVisS
	CTstSpU		<- CSmr+CVisL+CVisS		
	CTstSpR		<- CSmr+CVisL+CVisS	
	CTstSnUH		<- CSmr+CVisL+CVisS		
	CTstSnRH		<- CSmr+CVisL+CVisS	
	CTstSpUH		<- CSmr+CVisL+CVisS	
	CTstSpRH		<- CSmr+CVisL+CVisS
	CTstSnUT		<- CSmr+CVisL+CVisS		
	CTstSnRT		<- CSmr+CVisL+CVisS
	CTstSpUT		<- CSmr+CVisL+CVisS	
	CTstSpRT		<- CSmr+CVisL+CVisS
	VTestCostD1		<- c(c(CTstSnU,rep(c(CTstSnU,CTstSnU,CTstSpU,rep(0,4)),5)),c(CTstSnR,rep(c(CTstSnR,CTstSnR,CTstSpR,rep(0,4)),5)),
					rep(c(c(CTstSnUH,rep(c(CTstSnUH,CTstSnUH,CTstSpUH,rep(0,4)),5)),c(CTstSnRH,rep(c(CTstSnRH,CTstSnRH,CTstSpRH,rep(0,4)),5)),
					c(CTstSnUT,rep(c(CTstSnUT,CTstSnUT,CTstSpUT,rep(0,4)),5)),c(CTstSnRT,rep(c(CTstSnRT,CTstSnRT,CTstSpRT,rep(0,4)),5))),3))

# Standard DOTS Algorithm (all tx naive onto 1st-line regimen, pDST of all retreatment get DST, on to appropriate regimen)
	RegMat2 <- RegMat0
	RegMat2[,c(Vtemp4+5,Vtemp4+7)] <- c(1,0,0,0,0,0)						# Tx naive onto 1st line regimen
	RegMat2[,c(0:6*72,0:6*72+2)+41] <- c(1,0,0,0,0,0)*(1-pDstR)+ c(1,0,0,0,0,0)*pDstR	# Tx expcd, pansen
	RegMat2[,c(0:6*72,0:6*72+2)+48] <- c(1,0,0,0,0,0)*(1-pDstR)+ c(0,1,0,0,0,0)*pDstR	# Tx expcd, mono-INH
	RegMat2[,c(0:6*72,0:6*72+2)+55] <- c(1,0,0,0,0,0)*(1-pDstR)+ c(0,0,1,0,0,0)*pDstR	# Tx expcd, mono-RIF
	RegMat2[,c(0:6*72,0:6*72+2)+62] <- c(1,0,0,0,0,0)*(1-pDstR)+ c(0,0,0,1,0,0)*pDstR	# Tx expcd, MDR
	RegMat2[,c(0:6*72,0:6*72+2)+69] <- c(1,0,0,0,0,0)*(1-pDstR)+ c(0,0,0,1,0,0)*pDstR	# Tx expcd, XDR, but put on MDR treatment

	TxMatAlg2	<- TxMatAlg0
	TxMatAlg2[1,]	<- VDur%*%RegMat2
	for(i in 0:4) { TxMatAlg2[2,Vtemp2+7*i]  <-  TxEfMat[,i+1]%*%RegMat2[,Vtemp2+7*i]  }
	for(i in 1:6) { TxMatAlg2[i+2,Vtemp2+7*c(0,0,0,1,2,3)[i]]  <-	ARMat[,i]%*%RegMat2[,Vtemp2+7*c(0,0,0,1,2,3)[i]]  }
	for(i in 1:6) { TxMatAlg2[i+8,]  <-  OthMat[,i]%*%RegMat2  }

	TruPosDAlg2		<- rep(c(rep(c(SensSpIn*(1-pLtfuS)+pCulU*(1-SensSpIn)*SensCuIn*(1-pLtfuL),SensSpIp*(1-pLtfuS)+pCulU*(1-SensSpIp)*SensCuIp*(1-pLtfuL)),5),
					rep(c(SensSpIn*(1-pLtfuS)+pCulR*(1-SensSpIn)*SensCuIn*(1-pLtfuL) ,SensSpIp*(1-pLtfuS)+pCulR*(1-SensSpIp)*SensCuIp*(1-pLtfuL)),5)),7)
										# True positive diagnoses starting treatment
	TruPosDAlgB2	<- rep(c(rep(c(SensSpIn+pCulU*(1-SensSpIn)*SensCuIn,SensSpIp+pCulU*(1-SensSpIp)*SensCuIp),5),rep(c(SensSpIn+pCulR*(1-SensSpIn)*SensCuIn,SensSpIp+pCulR*(1-SensSpIp)*SensCuIp),5)),7)
										# True positive diagnoses identified (whether started treatment or not)
	FalsPosDAlg2	<- rep(c(rep((1-SpecSpSL)*(1-pLtfuS)+pCulU*SpecSpSL*(1-SpecCuSL)*(1-pLtfuL),6),rep((1-SpecSpSL)*(1-pLtfuS)+pCulR*SpecSpSL*(1-SpecCuSL)*(1-pLtfuL),6)),7) 
										# False positive diagnoses starting treatment
	FalsPosDAlgB2	<- rep(c(rep((1-SpecSpSL)+pCulU*SpecSpSL*(1-SpecCuSL),6),rep((1-SpecSpSL)+pCulR*SpecSpSL*(1-SpecCuSL),6)),7) 
										# False positive diagnoses identified (whether started treatment or not)
	# Vector of all states
	GetXpt2		<- rep(0,504)
	VTxCost2		<- t(t(TxMatAlg2[10:14,])%*%c(CVisS,CSmr,CCltr,CXry,CIpd))
	CTstSnU		<- CSmr+CVisL+CVisS + pCulU*(CCltr+CVisL)
	CTstSnR		<- CSmr+CVisL+CVisS + pCulR*(CCltr+CVisL)
	CTstSpU		<- CSmr+CVisL+CVisS + pCulU*(CCltr+CVisL)		
	CTstSpR		<- CSmr+CVisL+CVisS + pCulR*(CCltr+CVisL)	
	CTstSnUH		<- CSmr+CVisL+CVisS + pCulU*(CCltr+CVisL)		
	CTstSnRH		<- CSmr+CVisL+CVisS + pCulR*(CCltr+CVisL)	
	CTstSpUH		<- CSmr+CVisL+CVisS + pCulU*(CCltr+CVisL)	
	CTstSpRH		<- CSmr+CVisL+CVisS + pCulR*(CCltr+CVisL)
	CTstSnUT		<- CSmr+CVisL+CVisS + pCulU*(CCltr+CVisL)		
	CTstSnRT		<- CSmr+CVisL+CVisS + pCulR*(CCltr+CVisL)
	CTstSpUT		<- CSmr+CVisL+CVisS + pCulU*(CCltr+CVisL)	
	CTstSpRT		<- CSmr+CVisL+CVisS + pCulR*(CCltr+CVisL)

	VTestCostD2		<- c(c(CTstSnU,rep(c(CTstSnU,CTstSnU,CTstSpU,rep(0,4)),5)),c(CTstSnR,rep(c(CTstSnR,CTstSnR,CTstSpR,rep(0,4)),5)),
					rep(c(c(CTstSnUH,rep(c(CTstSnUH,CTstSnUH,CTstSpUH,rep(0,4)),5)),c(CTstSnRH,rep(c(CTstSnRH,CTstSnRH,CTstSpRH,rep(0,4)),5)),
					c(CTstSnUT,rep(c(CTstSnUT,CTstSnUT,CTstSpUT,rep(0,4)),5)),c(CTstSnRT,rep(c(CTstSnRT,CTstSnRT,CTstSpRT,rep(0,4)),5))),3))
	VTestCostD2[Vtemp9[rep(1:10,7)+rep(0:6*20,each=10)]] 	<- VTestCostD2[Vtemp9[rep(1:10,7)+rep(0:6*20,each=10)]] + TruPosDAlgB2[rep(1:10,7)+rep(0:6*20,each=10)]*pDstU*(CDst+CVisL)
	VTestCostD2[Vtemp9[rep(11:20,7)+rep(0:6*20,each=10)]] <- VTestCostD2[Vtemp9[rep(11:20,7)+rep(0:6*20,each=10)]] + TruPosDAlgB2[rep(11:20,7)+rep(0:6*20,each=10)]*pDstR*(CDst+CVisL)
	VTestCostD2[Vtemp8[rep(1:6,7)+rep(0:6*12,each=6)]] 	<- VTestCostD2[Vtemp8[rep(1:6,7)+rep(0:6*12,each=6)]] + FalsPosDAlgB2[rep(1:6,7)+rep(0:6*12,each=6)]*pDstU*(CDst+CVisL)
	VTestCostD2[Vtemp8[rep(7:12,7)+rep(0:6*12,each=6)]] 	<- VTestCostD2[Vtemp8[rep(7:12,7)+rep(0:6*12,each=6)]] + FalsPosDAlgB2[rep(7:12,7)+rep(0:6*12,each=6)]*pDstR*(CDst+CVisL)

# XPERT Algorithm 1 (all get Xpert, those testing RIF-Neg to 1st line)
	RegMat3 <- RegMat0
	RegMat3[,Vtemp2[c(0:13*7+4,0:13*7+6)]+0]  <- c(1,0,0,0,0,0)	# PanSens on to 1st line (impfct Xprt Spec made irrelevant by DST) 
	RegMat3[,Vtemp2[c(0:13*7+4,0:13*7+6)]+7]  <- c(1,0,0,0,0,0)*SpecXpRIF + c(0,1,0,0,0,0)*(1-SpecXpRIF)	# MonoINH 
	RegMat3[,Vtemp2[c(0:13*7+4,0:13*7+6)]+14] <- c(1,0,0,0,0,0)*(1-SensXpRIF)+c(0,0,1,0,0,0)*SensXpRIF	# MonoRIF
	RegMat3[,Vtemp2[c(0:13*7+4,0:13*7+6)]+21] <- c(1,0,0,0,0,0)*(1-SensXpRIF)+c(0,0,0,1,0,0)*SensXpRIF	# MDR
	RegMat3[,Vtemp2[c(0:13*7+4,0:13*7+6)]+28] <- c(1,0,0,0,0,0)*(1-SensXpRIF)+c(0,0,0,1,0,0)*SensXpRIF	# XDR, but put on MDR treatment

	TxMatAlg3	<- TxMatAlg0
	TxMatAlg3[1,]	<- VDur%*%RegMat3
	for(i in 0:4) { TxMatAlg3[2,Vtemp2+7*i]  <-  TxEfMat[,i+1]%*%RegMat3[,Vtemp2+7*i]  }
	for(i in 1:6) { TxMatAlg3[i+2,Vtemp2+7*c(0,0,0,1,2,3)[i]]  <-	ARMat[,i]%*%RegMat3[,Vtemp2+7*c(0,0,0,1,2,3)[i]]  }
	for(i in 1:6) { TxMatAlg3[i+8,]  <-  OthMat[,i]%*%RegMat3  }

	TruPosDAlg3		<- rep(c(rTstIn*SensXpIn,SensXpIp),70)*(1-pLtfuS) 	# True positive diagnoses starting treatment
	TruPosDAlgB3	<- rep(c(rTstIn*SensXpIn,SensXpIp),70) 			# True positive diagnoses identified (whether started treatment or not)
	FalsPosDAlg3	<- rep(c(rep((1-SpecXpSL)*(1-pLtfuS),6),rep((1-SpecXpSL)*(1-pLtfuS),6)),7) 
	FalsPosDAlgB3	<- rep(c(rep((1-SpecXpSL),6),rep((1-SpecXpSL),6)),7) 
	GetXpt3		<- rep(1,504)
	VTxCost3		<- t(t(TxMatAlg3[10:14,])%*%c(CVisS,CSmr,CCltr,CXry,CIpd))
	CTstSnU		<- CXpt+CVisL+CVisS
	CTstSnR		<- CXpt+CVisL+CVisS
	CTstSpU		<- CXpt+CVisL+CVisS		
	CTstSpR		<- CXpt+CVisL+CVisS		
	CTstSnUH		<- CXpt+CVisL+CVisS		
	CTstSnRH		<- CXpt+CVisL+CVisS	
	CTstSpUH		<- CXpt+CVisL+CVisS	
	CTstSpRH		<- CXpt+CVisL+CVisS
	CTstSnUT		<- CXpt+CVisL+CVisS		
	CTstSnRT		<- CXpt+CVisL+CVisS	
	CTstSpUT		<- CXpt+CVisL+CVisS	
	CTstSpRT		<- CXpt+CVisL+CVisS
	VTestCostD3		<- c(c(CTstSnU,rep(c(CTstSnU,CTstSnU,CTstSpU,rep(0,4)),5)),c(CTstSnR,rep(c(CTstSnR,CTstSnR,CTstSpR,rep(0,4)),5)),
					rep(c(c(CTstSnUH,rep(c(CTstSnUH,CTstSnUH,CTstSpUH,rep(0,4)),5)),c(CTstSnRH,rep(c(CTstSnRH,CTstSnRH,CTstSpRH,rep(0,4)),5)),
					c(CTstSnUT,rep(c(CTstSnUT,CTstSnUT,CTstSpUT,rep(0,4)),5)),c(CTstSnRT,rep(c(CTstSnRT,CTstSnRT,CTstSpRT,rep(0,4)),5))),3))
	VTestCostD3		<- VTestCostD3 + (CDst+CVisL)*rep(c((1-SpecXpSL)*(1-SpecXpRIF),rep(c((1-SpecXpSL)*(1-SpecXpRIF),SensXpIn*(1-SpecXpRIF),SensXpIp*
					(1-SpecXpRIF),0,0,0,0),2),rep(c((1-SpecXpSL)*(1-SpecXpRIF),SensXpIn*SensXpRIF,SensXpIp*SensXpRIF,0,0,0,0),3)),14)

### A. SETTING UP STATE VECTOR AND OUTPUTS VECTOR

	Vcurrent 	<- rep(NA, 505)
 	StatNam <- c("Ls","In","Ip","N1","N2","P1","P2"); temp <- StatNam
		for (i in 1:5) {StatNam <- c(StatNam,paste(i,temp,sep=""))}; StatNam <- c("Su",StatNam[-(1:7)])
	StatNam <- c(paste("U",StatNam,sep=""),paste("R",StatNam,sep=""))
	StatNam <- c(paste("N0",StatNam,sep=""),paste("H1",StatNam,sep=""),paste("T1",StatNam,sep=""),paste("H2",StatNam,sep=""),
		paste("T2",StatNam,sep=""),paste("H3",StatNam,sep=""),paste("T3",StatNam,sep=""),"DEAD")

	names(Vcurrent) <- StatNam
	
	Vout			<- rep(NA,92)
	names(Vout)		<- c("NAll","Ndaly","NAnyTb","NActDis","NUnTx","NUnTxH","NSmP","NStr1n","NStr2n","NStr3n","NStr4n","NStr5n","NStr1e",
					"NStr2e","NStr3e","NStr4e","NStr5e","NTxD","NTxND","NHiv","NHiv350","NArt","NTbH","NMort","NHivMort",
					"NTbMort","NSmPMort","NTbHMort","NInf","NCase","NCaseNF","NCaseNS","NCaseHF","NCaseHS","NCaseIp","NCaseIpHiv",
					"NCdIpD","NCdInD","NCdIpND","NCdInND","NTxResU","NTxMdrU","NTxXdrU","NTxResR","NTxMdrR","NTxXdrR",
					"NTxSp","SuspctD","SuspctDTB","SuspctND","CostTxD","CostTxND","CostRegD","CostRegND","CostART","CostTestD",
					"CostTestND","CostFalsTxD","CostFalsTxND","CostFalsRegD","CostFalsRegND","NCdFalsD","NCdFalsND","Check1",
					"DurInfSn","DurInfSp","DurInfAll","PfailDtx","PcureDtx","PdfltDtx","PmortDtx","EffContRate","NotifD", "NotifTBD",
					"NotifND","PPVTb","NPVTb","PPVRif","NPVRif","PDst","ExTbC","ExTbT","ExTbD","ExTbCH","ExTbTH","ExTbDH","GetXpt",
					"TC","NMDR","ArtCov","ArtNdCov","Art200Cov")
	OutMat 		<- matrix(0,nrow=(150*12),ncol=length(Vout))
	colnames(OutMat)	<- names(Vout)




