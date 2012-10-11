##########   TB XPERT DIAGNOSTIC MODEL 2012  ##########   
##########           TIMESTEP FILE           ##########  
## NOTES
# a.	This is the workhorse function, calculating each new month of the model.
# b.	This is a single function, taking as its arguments (i) a state vector and (ii) a subset of variables (all the time varying ones), and (iii) a partially filled model matrix.
# c.	New entrants are added to state vector, deaths removed
# d.	Variables are retransformed and dynamic functions calculated (TB force of infection)
# e.	Static model matrix is updated with time-varying values, diagonal elements calculated so that rowsums = 1  
# f.	Matrix multiplication to update state vector
# g.	Results calculated
# h.	State vector and results vector outputed 

	aRateMatStat <- matrix(0,nrow=505,ncol=505); rownames(aRateMatStat) <- StatNam; colnames(aRateMatStat) <- StatNam
	


#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%
#%#%
	timestep	<- function(Vcurrent,t)  {   # Start of function!
#%#%
#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%

### C1. ADD NEW ENTRANTS TO STATE VECTOR
	Vnext 		<- Vcurrent		# Initializes new vector
	Vnext[505] 		<- 0			# Clears out deaths
	Vnext[1]		<- Vnext[1]	+ NewEntt[t]*1000000	  # Adds new entrants to NU0Su based on birth rate

### B. SETTING UP STATIC PARTS OF RATE MATRIX (FOR SPEED OF TIMESTEP FUNCTION)
	RateMatStat			<- aRateMatStat

### B1. BREAKDOWN TO ACTIVE DISEASE (Stay in HIV / Resistance / Treatment subdivisions) 
	j <- rep(0:6,each=10)
	i <- Vtemp1[rep(1:10,7)+j*10]+2

	RateMatStat[cbind(i,i+1)]  <- VrBreakD[j+1]*VpToIn[j+1]
	RateMatStat[cbind(i,i+2)]  <- VrBreakD[j+1]*(1-VpToIn[j+1])

#	for(j in 0:6) {   for (i in Vtemp1[1:10+j*10]+2)  {
#	RateMatStat[i,i+1:2]  <- VrBreakD[j+1]*c(VpToIn[j+1],(1-VpToIn[j+1]))	} }
 
### B2. SMEAR NEG CONVERT TO SMEAR POS (Stay in HIV / Resistance / Treatment subdivisions) 
	RateMatStat[cbind(i+1,i+2)]  <- rNtoP 

#	for(j in 0:6) {   for (i in Vtemp1[1:10+j*10]+3)  {  RateMatStat[i,i+1]  <- rNtoP  }  }

### B3. SPONTANEOUS CURE (Stay in HIV / Resistance / Treatment subdivisions) 
	RateMatStat[cbind(i+1,i)]  <- VrIToLs[j+1]
	RateMatStat[cbind(i+2,i)]  <- VrIToLs[j+1]

#	for(j in 0:6) {   for (i in Vtemp1[1:10+j*10]+2)  {  RateMatStat[i+1:2,i]  <- VrIToLs[j+1]  }  }

### B4. HIV Progression 
	RateMatStat[cbind(1:72+72,1:72+216)] <- H1toH2	
	RateMatStat[cbind(1:72+216,1:72+360)] <- H2toH3  

#	for (i in 73:144){  RateMatStat[i,i+144] <- H1toH2;	RateMatStat[i+144,i+288] <- H2toH3  }
 
### TREATMENT TRANSITIONS, BY ALGORITHM
	for (z in 1:Nalg) {  # this needs to be for all algorithms z
	RateMatStatz	<- RateMatStat
	TxMat			<- get(paste("TxMatAlg",z,sep=""))

### B5. POPULATE MORTALITY RATES
	RateMatStatz[cbind(73:504,rep(505,432))] <- rep(VmuHIV,each=72)
#	for (i in 1:6) { RateMatStatz[1:72+i*72,505] <- rep(VmuHIV[i],72)	}		# HIV mortality
	RateMatStatz[Vtemp1+3,505]		<- RateMatStatz[Vtemp1+3,505] + muIn	# Untreated Smear-neg TB mortality
	RateMatStatz[Vtemp1+4,505]		<- RateMatStatz[Vtemp1+4,505] + muIp	# Untreated Smear-pos TB mortality
	RateMatStatz[Vtemp6,505]		<- RateMatStatz[Vtemp6,505]+(1-TxMat[2,Vtemp6]*TxEft[t])*rep(c(muIn,muIn,muIp,muIp),70)
														# Treatment TB mortality
	RateMatStatz[Vtemp9[c(61:80,100:120)],505]	<- muTbH				# TB-HIV mortality for CD4 350
	RateMatStatz[Vtemp6,505]		<- RateMatStatz[Vtemp6,505]+TunTxMort	# Tuning tx mortality to observed levels

# Vector of current contact rates...
	VTrStatz			<- rep(0,504) 				# Creates a vector for contact rates
	VTrStatz[Vtemp1+3]	<- rep(TrIn*RelFit,14) 			# Contact rates for smear neg 
	VTrStatz[Vtemp1+4]	<- rep(RelFit,14) 			# Contact rates for smear pos 
	VTrStatz[Vtemp6]		<- (1-TxMat[2,Vtemp6]*TxEft[t])*rep(rep(c(TrIn,TrIn,1,1),5)*rep(RelFit,each=4),14)*TunTxEff
											# Contact rates for individuals on treatment
	assign(paste('VTrStat',z,sep=""),VTrStatz)

### B5. TREATMENT OUTCOMES (Stay in HIV / Resistance subdivisions) 

 # Cures back to Ls state, treatment experienced subdivision:
	RateMatStatz[cbind(Vtemp4[1:35]+5,Vtemp4[1:35]+2+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+5,Vtemp4[1:35]+2+36)] + 12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+5)]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+5)]*TxEft[t]
	RateMatStatz[cbind(Vtemp4[1:35]+6,Vtemp4[1:35]+2+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+6,Vtemp4[1:35]+2+36)] + 12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+6)]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+6)]*TxEft[t]
	RateMatStatz[cbind(Vtemp4[1:35]+7,Vtemp4[1:35]+2+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+7,Vtemp4[1:35]+2+36)] + 12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+7)]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+7)]*TxEft[t]
	RateMatStatz[cbind(Vtemp4[1:35]+8,Vtemp4[1:35]+2+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+8,Vtemp4[1:35]+2+36)] + 12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+8)]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+8)]*TxEft[t]

	RateMatStatz[cbind(Vtemp4[1:35]+5+36,Vtemp4[1:35]+2+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+5+36,Vtemp4[1:35]+2+36)] + 12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+5+36)]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+5+36)]*TxEft[t]
	RateMatStatz[cbind(Vtemp4[1:35]+6+36,Vtemp4[1:35]+2+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+6+36,Vtemp4[1:35]+2+36)] + 12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+6+36)]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+6+36)]*TxEft[t]
	RateMatStatz[cbind(Vtemp4[1:35]+7+36,Vtemp4[1:35]+2+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+7+36,Vtemp4[1:35]+2+36)] + 12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+7+36)]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+7+36)]*TxEft[t]
	RateMatStatz[cbind(Vtemp4[1:35]+8+36,Vtemp4[1:35]+2+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+8+36,Vtemp4[1:35]+2+36)] + 12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+8+36)]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+8+36)]*TxEft[t]

	# for(i in 1:35) {
	#	RateMatStatz[Vtemp4[i]+5,Vtemp4[i]+2+36] <-  RateMatStatz[Vtemp4[i]+5,Vtemp4[i]+2+36] + 12/TxMat[1,Vtemp4[i]+5]*TxMat[2,Vtemp4[i]+5]*TxEft[t]
	#	RateMatStatz[Vtemp4[i]+6,Vtemp4[i]+2+36] <-  RateMatStatz[Vtemp4[i]+6,Vtemp4[i]+2+36] + 12/TxMat[1,Vtemp4[i]+6]*TxMat[2,Vtemp4[i]+6]*TxEft[t]
	#	RateMatStatz[Vtemp4[i]+7,Vtemp4[i]+2+36] <-  RateMatStatz[Vtemp4[i]+7,Vtemp4[i]+2+36] + 12/TxMat[1,Vtemp4[i]+7]*TxMat[2,Vtemp4[i]+7]*TxEft[t]
	#	RateMatStatz[Vtemp4[i]+8,Vtemp4[i]+2+36] <-  RateMatStatz[Vtemp4[i]+8,Vtemp4[i]+2+36] + 12/TxMat[1,Vtemp4[i]+8]*TxMat[2,Vtemp4[i]+8]*TxEft[t]

	#	RateMatStatz[Vtemp4[i]+5+36,Vtemp4[i]+2+36] <-  RateMatStatz[Vtemp4[i]+5+36,Vtemp4[i]+2+36] + 12/TxMat[1,Vtemp4[i]+5+36]*TxMat[2,Vtemp4[i]+5+36]*TxEft[t]
	#	RateMatStatz[Vtemp4[i]+6+36,Vtemp4[i]+2+36] <-  RateMatStatz[Vtemp4[i]+6+36,Vtemp4[i]+2+36] + 12/TxMat[1,Vtemp4[i]+6+36]*TxMat[2,Vtemp4[i]+6+36]*TxEft[t]
	#	RateMatStatz[Vtemp4[i]+7+36,Vtemp4[i]+2+36] <-  RateMatStatz[Vtemp4[i]+7+36,Vtemp4[i]+2+36] + 12/TxMat[1,Vtemp4[i]+7+36]*TxMat[2,Vtemp4[i]+7+36]*TxEft[t]
	#	RateMatStatz[Vtemp4[i]+8+36,Vtemp4[i]+2+36] <-  RateMatStatz[Vtemp4[i]+8+36,Vtemp4[i]+2+36] + 12/TxMat[1,Vtemp4[i]+8+36]*TxMat[2,Vtemp4[i]+8+36]*TxEft[t]
	#	}

 # Failures identified and reinitiated on treatment, treatment experienced subdivision:

	RateMatStatz[cbind(Vtemp4[1:35]+5,Vtemp4[1:35]+5+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+5,Vtemp4[1:35]+5+36)] + 12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+5)]*(1-TxEft[t]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+5)])*pReTx
	RateMatStatz[cbind(Vtemp4[1:35]+6,Vtemp4[1:35]+6+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+6,Vtemp4[1:35]+6+36)] + 12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+6)]*(1-TxEft[t]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+6)])*pReTx
	RateMatStatz[cbind(Vtemp4[1:35]+7,Vtemp4[1:35]+7+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+7,Vtemp4[1:35]+7+36)] + 12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+7)]*(1-TxEft[t]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+7)])*pReTx
	RateMatStatz[cbind(Vtemp4[1:35]+8,Vtemp4[1:35]+8+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+8,Vtemp4[1:35]+8+36)] + 12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+8)]*(1-TxEft[t]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+8)])*pReTx

	RateMatStatz[cbind(Vtemp4[1:35]+5+36,Vtemp4[1:35]+5+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+5+36,Vtemp4[1:35]+5+36)] + 12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+5+36)]*(1-TxEft[t]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+5+36)])*pReTx
	RateMatStatz[cbind(Vtemp4[1:35]+6+36,Vtemp4[1:35]+6+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+6+36,Vtemp4[1:35]+6+36)] + 12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+6+36)]*(1-TxEft[t]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+6+36)])*pReTx
	RateMatStatz[cbind(Vtemp4[1:35]+7+36,Vtemp4[1:35]+7+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+7+36,Vtemp4[1:35]+7+36)] + 12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+7+36)]*(1-TxEft[t]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+7+36)])*pReTx
	RateMatStatz[cbind(Vtemp4[1:35]+8+36,Vtemp4[1:35]+8+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+8+36,Vtemp4[1:35]+8+36)] + 12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+8+36)]*(1-TxEft[t]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+8+36)])*pReTx
		

	# for(i in 1:35) {
	#	RateMatStatz[Vtemp4[i]+5,Vtemp4[i]+5+36] <-  RateMatStatz[Vtemp4[i]+5,Vtemp4[i]+5+36] + 12/TxMat[1,Vtemp4[i]+5]*(1-TxEft[t]*TxMat[2,Vtemp4[i]+5])*pReTx
	#	RateMatStatz[Vtemp4[i]+6,Vtemp4[i]+6+36] <-  RateMatStatz[Vtemp4[i]+6,Vtemp4[i]+6+36] + 12/TxMat[1,Vtemp4[i]+6]*(1-TxEft[t]*TxMat[2,Vtemp4[i]+6])*pReTx
	#	RateMatStatz[Vtemp4[i]+7,Vtemp4[i]+7+36] <-  RateMatStatz[Vtemp4[i]+7,Vtemp4[i]+7+36] + 12/TxMat[1,Vtemp4[i]+7]*(1-TxEft[t]*TxMat[2,Vtemp4[i]+7])*pReTx
	#	RateMatStatz[Vtemp4[i]+8,Vtemp4[i]+8+36] <-  RateMatStatz[Vtemp4[i]+8,Vtemp4[i]+8+36] + 12/TxMat[1,Vtemp4[i]+8]*(1-TxEft[t]*TxMat[2,Vtemp4[i]+8])*pReTx

	#	RateMatStatz[Vtemp4[i]+5+36,Vtemp4[i]+5+36] <-  RateMatStatz[Vtemp4[i]+5+36,Vtemp4[i]+5+36] + 12/TxMat[1,Vtemp4[i]+5+36]*(1-TxEft[t]*TxMat[2,Vtemp4[i]+5+36])*pReTx
	#	RateMatStatz[Vtemp4[i]+6+36,Vtemp4[i]+6+36] <-  RateMatStatz[Vtemp4[i]+6+36,Vtemp4[i]+6+36] + 12/TxMat[1,Vtemp4[i]+6+36]*(1-TxEft[t]*TxMat[2,Vtemp4[i]+6+36])*pReTx
	#	RateMatStatz[Vtemp4[i]+7+36,Vtemp4[i]+7+36] <-  RateMatStatz[Vtemp4[i]+7+36,Vtemp4[i]+7+36] + 12/TxMat[1,Vtemp4[i]+7+36]*(1-TxEft[t]*TxMat[2,Vtemp4[i]+7+36])*pReTx
	#	RateMatStatz[Vtemp4[i]+8+36,Vtemp4[i]+8+36] <-  RateMatStatz[Vtemp4[i]+8+36,Vtemp4[i]+8+36] + 12/TxMat[1,Vtemp4[i]+8+36]*(1-TxEft[t]*TxMat[2,Vtemp4[i]+8+36])*pReTx
	#	}

 # Defaulters and failures to active disease:

		RateMatStatz[cbind(Vtemp4[1:35]+5,Vtemp4[1:35]+3+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+5,Vtemp4[1:35]+3+36)] + apply(cbind(rep(0,35),12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+5)]*(pDeft[t]+(1-TxEft[t]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+5)])*(1-pReTx))-colSums(TxMat[3:8,Vtemp4[1:35]+5])),1,max)
		RateMatStatz[cbind(Vtemp4[1:35]+6,Vtemp4[1:35]+3+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+6,Vtemp4[1:35]+3+36)] + apply(cbind(rep(0,35),12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+6)]*(pDefND+(1-TxEft[t]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+6)])*(1-pReTx))-colSums(TxMat[3:8,Vtemp4[1:35]+6])),1,max)
		RateMatStatz[cbind(Vtemp4[1:35]+7,Vtemp4[1:35]+4+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+7,Vtemp4[1:35]+4+36)] + apply(cbind(rep(0,35),12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+7)]*(pDeft[t]+(1-TxEft[t]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+7)])*(1-pReTx))-colSums(TxMat[3:8,Vtemp4[1:35]+7])),1,max)
		RateMatStatz[cbind(Vtemp4[1:35]+8,Vtemp4[1:35]+4+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+8,Vtemp4[1:35]+4+36)] + apply(cbind(rep(0,35),12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+8)]*(pDefND+(1-TxEft[t]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+8)])*(1-pReTx))-colSums(TxMat[3:8,Vtemp4[1:35]+8])),1,max)

		RateMatStatz[cbind(Vtemp4[1:35]+5+36,Vtemp4[1:35]+3+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+5+36,Vtemp4[1:35]+3+36)] + apply(cbind(rep(0,35),12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+5+36)]*(pDeft[t]+(1-TxEft[t]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+5+36)])*(1-pReTx))-colSums(TxMat[3:8,Vtemp4[1:35]+5+36])),1,max)
		RateMatStatz[cbind(Vtemp4[1:35]+6+36,Vtemp4[1:35]+3+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+6+36,Vtemp4[1:35]+3+36)] + apply(cbind(rep(0,35),12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+6+36)]*(pDefND+(1-TxEft[t]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+6+36)])*(1-pReTx))-colSums(TxMat[3:8,Vtemp4[1:35]+6+36])),1,max)
		RateMatStatz[cbind(Vtemp4[1:35]+7+36,Vtemp4[1:35]+4+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+7+36,Vtemp4[1:35]+4+36)] + apply(cbind(rep(0,35),12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+7+36)]*(pDeft[t]+(1-TxEft[t]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+7+36)])*(1-pReTx))-colSums(TxMat[3:8,Vtemp4[1:35]+7+36])),1,max)
		RateMatStatz[cbind(Vtemp4[1:35]+8+36,Vtemp4[1:35]+4+36)] <-  RateMatStatz[cbind(Vtemp4[1:35]+8+36,Vtemp4[1:35]+4+36)] + apply(cbind(rep(0,35),12/TxMat[cbind(rep(1,35),Vtemp4[1:35]+8+36)]*(pDefND+(1-TxEft[t]*TxMat[cbind(rep(2,35),Vtemp4[1:35]+8+36)])*(1-pReTx))-colSums(TxMat[3:8,Vtemp4[1:35]+8+36])),1,max)

	# for(i in 1:35) {
	#	RateMatStatz[Vtemp4[i]+5,Vtemp4[i]+3+36] <-  RateMatStatz[Vtemp4[i]+5,Vtemp4[i]+3+36] + max(0,12/TxMat[1,Vtemp4[i]+5]*(pDeft[t]+(1-TxEft[t]*TxMat[2,Vtemp4[i]+5])*(1-pReTx))-sum(TxMat[3:8,Vtemp4[i]+5]))
	#	RateMatStatz[Vtemp4[i]+6,Vtemp4[i]+3+36] <-  RateMatStatz[Vtemp4[i]+6,Vtemp4[i]+3+36] + max(0,12/TxMat[1,Vtemp4[i]+6]*(pDefND+(1-TxEft[t]*TxMat[2,Vtemp4[i]+6])*(1-pReTx))-sum(TxMat[3:8,Vtemp4[i]+6]))
	#	RateMatStatz[Vtemp4[i]+7,Vtemp4[i]+4+36] <-  RateMatStatz[Vtemp4[i]+7,Vtemp4[i]+4+36] + max(0,12/TxMat[1,Vtemp4[i]+7]*(pDeft[t]+(1-TxEft[t]*TxMat[2,Vtemp4[i]+7])*(1-pReTx))-sum(TxMat[3:8,Vtemp4[i]+7]))
	#	RateMatStatz[Vtemp4[i]+8,Vtemp4[i]+4+36] <-  RateMatStatz[Vtemp4[i]+8,Vtemp4[i]+4+36] + max(0,12/TxMat[1,Vtemp4[i]+8]*(pDefND+(1-TxEft[t]*TxMat[2,Vtemp4[i]+8])*(1-pReTx))-sum(TxMat[3:8,Vtemp4[i]+8]))

	#	RateMatStatz[Vtemp4[i]+5+36,Vtemp4[i]+3+36] <-  RateMatStatz[Vtemp4[i]+5+36,Vtemp4[i]+3+36] + max(0,12/TxMat[1,Vtemp4[i]+5+36]*(pDeft[t]+(1-TxEft[t]*TxMat[2,Vtemp4[i]+5+36])*(1-pReTx))-sum(TxMat[3:8,Vtemp4[i]+5+36]))
	#	RateMatStatz[Vtemp4[i]+6+36,Vtemp4[i]+3+36] <-  RateMatStatz[Vtemp4[i]+6+36,Vtemp4[i]+3+36] + max(0,12/TxMat[1,Vtemp4[i]+6+36]*(pDefND+(1-TxEft[t]*TxMat[2,Vtemp4[i]+6+36])*(1-pReTx))-sum(TxMat[3:8,Vtemp4[i]+6+36]))
	#	RateMatStatz[Vtemp4[i]+7+36,Vtemp4[i]+4+36] <-  RateMatStatz[Vtemp4[i]+7+36,Vtemp4[i]+4+36] + max(0,12/TxMat[1,Vtemp4[i]+7+36]*(pDeft[t]+(1-TxEft[t]*TxMat[2,Vtemp4[i]+7+36])*(1-pReTx))-sum(TxMat[3:8,Vtemp4[i]+7+36]))
	#	RateMatStatz[Vtemp4[i]+8+36,Vtemp4[i]+4+36] <-  RateMatStatz[Vtemp4[i]+8+36,Vtemp4[i]+4+36] + max(0,12/TxMat[1,Vtemp4[i]+8+36]*(pDefND+(1-TxEft[t]*TxMat[2,Vtemp4[i]+8+36])*(1-pReTx))-sum(TxMat[3:8,Vtemp4[i]+8+36]))
	#	}

 # Defaulters and failures to active disease with Acquired Resistance:
	j <- rep(rep(rep(0:6*5,each=5),6),each=8)
	i <- Vtemp4[rep(rep(1:5,7*6),each=8)+j]+rep(c(5:8,5:8+36),5*6*7)
	m <- rep(c(3,3,4,4,3,3,4,4)+36,5*6*7)

	 for(k in 1:6)   { RateMatStatz[cbind(i,Vtemp4[c(2,3,4,4,4,5)[k]+j]+m)] <-  RateMatStatz[cbind(i,Vtemp4[c(2,3,4,4,4,5)[k]+j]+m)]+TxMat[k+2,i] }

	# for(k in 1:6)   {
	# for(j in 0:6*5) {
	# for(i in 1:5)   {
	#	RateMatStatz[Vtemp4[i+j]+5,Vtemp4[c(2,3,4,4,4,5)[k]+j]+3+36] <-  RateMatStatz[Vtemp4[i+j]+5,Vtemp4[c(2,3,4,4,4,5)[k]+j]+3+36] + TxMat[k+2,Vtemp4[i+j]+5]
	#	RateMatStatz[Vtemp4[i+j]+6,Vtemp4[c(2,3,4,4,4,5)[k]+j]+3+36] <-  RateMatStatz[Vtemp4[i+j]+6,Vtemp4[c(2,3,4,4,4,5)[k]+j]+3+36] + TxMat[k+2,Vtemp4[i+j]+6]
	#	RateMatStatz[Vtemp4[i+j]+7,Vtemp4[c(2,3,4,4,4,5)[k]+j]+4+36] <-  RateMatStatz[Vtemp4[i+j]+7,Vtemp4[c(2,3,4,4,4,5)[k]+j]+4+36] + TxMat[k+2,Vtemp4[i+j]+7]
	#	RateMatStatz[Vtemp4[i+j]+8,Vtemp4[c(2,3,4,4,4,5)[k]+j]+4+36] <-  RateMatStatz[Vtemp4[i+j]+8,Vtemp4[c(2,3,4,4,4,5)[k]+j]+4+36] + TxMat[k+2,Vtemp4[i+j]+8]
	
	#	RateMatStatz[Vtemp4[i+j]+5+36,Vtemp4[c(2,3,4,4,4,5)[k]+j]+3+36] <-  RateMatStatz[Vtemp4[i+j]+5+36,Vtemp4[c(2,3,4,4,4,5)[k]+j]+3+36] + TxMat[k+2,Vtemp4[i+j]+5+36]
	#	RateMatStatz[Vtemp4[i+j]+6+36,Vtemp4[c(2,3,4,4,4,5)[k]+j]+3+36] <-  RateMatStatz[Vtemp4[i+j]+6+36,Vtemp4[c(2,3,4,4,4,5)[k]+j]+3+36] + TxMat[k+2,Vtemp4[i+j]+6+36]
	#	RateMatStatz[Vtemp4[i+j]+7+36,Vtemp4[c(2,3,4,4,4,5)[k]+j]+4+36] <-  RateMatStatz[Vtemp4[i+j]+7+36,Vtemp4[c(2,3,4,4,4,5)[k]+j]+4+36] + TxMat[k+2,Vtemp4[i+j]+7+36]
	#	RateMatStatz[Vtemp4[i+j]+8+36,Vtemp4[c(2,3,4,4,4,5)[k]+j]+4+36] <-  RateMatStatz[Vtemp4[i+j]+8+36,Vtemp4[c(2,3,4,4,4,5)[k]+j]+4+36] + TxMat[k+2,Vtemp4[i+j]+8+36]
	#	}  }  }

	assign(paste("RateMatStat",z,sep=""),RateMatStatz)	
	} 

### C2. CREATE RATE MATRIX AND TRANSITION PARAMETER VECTOR
	RateMat 		<- (RateMatStat1*(1-PhaseIn1[t])+RateMatStat2*PhaseIn1[t])*(1-PhaseIn2[t]) +
					get(paste("RateMatStat",DIAG,sep=""))*PhaseIn2[t]
	VTrStat		<- (VTrStat1*(1-PhaseIn1[t])+VTrStat2*PhaseIn1[t])*(1-PhaseIn2[t]) +
					get(paste("VTrStat",DIAG,sep=""))*PhaseIn2[t]

### C3. UPDATE MORTALITY RATES WITH BACKGROUND MORTALITY
	RateMat[-505,505]	<- RateMat[-505,505]+mubt[t]

### C4. TB INCIDENCE (Can change strain subdivision, stay in HIV / treatment subd.) 
	VInf 		<- Vnext[1:504]/sum(Vnext[1:504])*VTrStat*CRt[t]			# P(meet carrier)*CR|carrier, homogeneous mixing

	k <- rep(0:6,each=10)
	j <- rep(rep(c(1,37),each=5),7)+k*72
	i <- rep(0:4,14)
	m <- c(sum(VInf[Vtemp2+0*7]),sum(VInf[Vtemp2+1*7]),sum(VInf[Vtemp2+2*7]),sum(VInf[Vtemp2+3*7]),sum(VInf[Vtemp2+4*7])); m <- rep(m,14)

		RateMat[cbind(j,j+1+i*7)]	<- RateMat[cbind(j,j+1+i*7)]+m*(1-Vpfast[k+1])
		RateMat[cbind(j,j+2+i*7)]	<- RateMat[cbind(j,j+2+i*7)]+m*Vpfast[k+1]*VpToIn[k+1]
		RateMat[cbind(j,j+3+i*7)]	<- RateMat[cbind(j,j+3+i*7)]+m*Vpfast[k+1]*(1-VpToIn[k+1])

	# for(k in 0:6) {
	# for(j in c(1,37)+k*72) {
	# for(i in 0:4) {
	#	RateMat[j,j+1+i*7]	<- RateMat[j,j+1+i*7]+sum(VInf[Vtemp2+i*7])*(1-Vpfast[k+1])
	#	RateMat[j,j+2+i*7]	<- RateMat[j,j+2+i*7]+sum(VInf[Vtemp2+i*7])*Vpfast[k+1]*VpToIn[k+1]
	#	RateMat[j,j+3+i*7]	<- RateMat[j,j+3+i*7]+sum(VInf[Vtemp2+i*7])*Vpfast[k+1]*(1-VpToIn[k+1])
	#	} } }	# Contact rates for HIV neg, applied to Su states

### C5. SUPERINFECTION (Can change strain subdivision, stay in HIV / treatment subd.) 
	VSupInf 		<- VInf*(1-rep(VPartIm,each=72))  # As above, with partial immunity, homogeneous mixing

	k <- rep(0:6,each=2*5*5)
	j <- rep(rep(c(2,38),each=5*5),7)+k*72
	i <- rep(rep(0:4,each=5),7*2)
	m <- rep(0:4*7,7*2*5)
	v <- c(sum(VSupInf[Vtemp2+0*7]),sum(VSupInf[Vtemp2+1*7]),sum(VSupInf[Vtemp2+2*7]),sum(VSupInf[Vtemp2+3*7]),sum(VSupInf[Vtemp2+4*7])); v <- v[i+1]
 
		RateMat[cbind(j+m,j+0+i*7)]	<- RateMat[cbind(j+m,j+0+i*7)]+v*(1-Vpfast[k+1])
		RateMat[cbind(j+m,j+1+i*7)]	<- RateMat[cbind(j+m,j+1+i*7)]+v*Vpfast[k+1]*VpToIn[k+1]
		RateMat[cbind(j+m,j+2+i*7)]	<- RateMat[cbind(j+m,j+2+i*7)]+v*Vpfast[k+1]*(1-VpToIn[k+1])


	# for(k in 0:6) {
	# for(j in c(2,38)+k*72) {
	# for(i in 0:4) {
	#  for(m in 0:4*7) {
	#	RateMat[j+m,j+0+i*7]	<- RateMat[j+m,j+0+i*7]+sum(VSupInf[Vtemp2+i*7])*(1-Vpfast[k+1])/2
	#	RateMat[j+m,j+1+i*7]	<- RateMat[j+m,j+1+i*7]+sum(VSupInf[Vtemp2+i*7])*Vpfast[k+1]*VpToIn[k+1]
	#	RateMat[j+m,j+2+i*7]	<- RateMat[j+m,j+2+i*7]+sum(VSupInf[Vtemp2+i*7])*Vpfast[k+1]*(1-VpToIn[k+1])
	#	}}
	#for(i in 0:4) {	
	#	RateMat[j+0+i*7,j+0+i*7]	<- RateMat[j+0+i*7,j+0+i*7]+sum(RateMat[j+0+i*7,j+c(0:4*7)]) 
	# 	} } }  	# Contact rates for HIV neg,with partial immunity, applied to Ls states


### C6. DIAGNOSIS AND TREATMENT STRATEGY (Stay in HIV / Resistance / Treatment subdivisions) 

#  C6a. Specifying diagnosis and treatment as a result of algorithm
	TxMat		<- (TxMatAlg1*(1-PhaseIn1[t])+TxMatAlg2*PhaseIn1[t])*(1-PhaseIn2[t]) +
				get(paste("TxMatAlg",DIAG,sep=""))*PhaseIn2[t]
	TruPosD	<- (TruPosDAlg1*(1-PhaseIn1[t])+TruPosDAlg2*PhaseIn1[t])*(1-PhaseIn2[t]) +
				get(paste("TruPosDAlg",DIAG,sep=""))*PhaseIn2[t]
	FalsPosD	<- (FalsPosDAlg1*(1-PhaseIn1[t])+FalsPosDAlg2*PhaseIn1[t])*(1-PhaseIn2[t]) +
				get(paste("FalsPosDAlg",DIAG,sep=""))*PhaseIn2[t]
	TruPosDB	<- (TruPosDAlgB1*(1-PhaseIn1[t])+TruPosDAlgB2*PhaseIn1[t])*(1-PhaseIn2[t]) +
				get(paste("TruPosDAlgB",DIAG,sep=""))*PhaseIn2[t]
	FalsPosDB	<- (FalsPosDAlgB1*(1-PhaseIn1[t])+FalsPosDAlgB2*PhaseIn1[t])*(1-PhaseIn2[t]) +
				get(paste("FalsPosDAlgB",DIAG,sep=""))*PhaseIn2[t]
	VTestCostD	<- (VTestCostD1*(1-PhaseIn1[t])+VTestCostD2*PhaseIn1[t])*(1-PhaseIn2[t]) +
				get(paste("VTestCostD",DIAG,sep=""))*PhaseIn2[t]
	VTxCost	<- (VTxCost1*(1-PhaseIn1[t])+VTxCost2*PhaseIn1[t])*(1-PhaseIn2[t]) +
				get(paste("VTxCost",DIAG,sep=""))*PhaseIn2[t]
	GetXpt	<- (GetXpt1*(1-PhaseIn1[t])+GetXpt2*PhaseIn1[t])*(1-PhaseIn2[t]) +
				get(paste("GetXpt",DIAG,sep=""))*PhaseIn2[t]

#  C6b. Diagnosis and tx initiation

		RateMat[cbind(Vtemp1[1:70]+3,Vtemp1[1:70]+5)]	<-	DTestt[t]*TruPosD[1:70*2-1]*rTstIn 		# From In to Tn1
		RateMat[cbind(Vtemp1[1:70]+3,Vtemp1[1:70]+6)]	<-	NDTestt[t]*TruPosND[1:70*2-1]*rTstIn  	# From In to Tn2
		RateMat[cbind(Vtemp1[1:70]+4,Vtemp1[1:70]+7)]	<-	DTestt[t]*TruPosD[1:70*2]   			# From Ip to Tp1
		RateMat[cbind(Vtemp1[1:70]+4,Vtemp1[1:70]+8)]	<-	NDTestt[t]*TruPosND[1:70*2]   		# From Ip to Tp2

	# for (i in 1:70)  {
	#	RateMat[Vtemp1[i]+3,Vtemp1[i]+5]	<-	DTestt[t]*TruPosD[i*2-1]*rTstIn 	# From In to Tn1
	#	RateMat[Vtemp1[i]+3,Vtemp1[i]+6]	<-	NDTestt[t]*TruPosND[i*2-1]*rTstIn  	# From In to Tn2
	#	RateMat[Vtemp1[i]+4,Vtemp1[i]+7]	<-	DTestt[t]*TruPosD[i*2]   		# From Ip to Tp1
	#	RateMat[Vtemp1[i]+4,Vtemp1[i]+8]	<-	NDTestt[t]*TruPosND[i*2]   		# From Ip to Tp2
	#	}

### C8. HIV INCIDENCE and ART ENROLLMENT
	# HIV incidence
	RateMat[cbind(1:72,1:72+72)] 	<- rHIVt[t]
	
	# for (i in 1:72){	RateMat[i,i+72] 	<- rHIVt[t]  }

	OnTx <- sum(Vnext[c(145:216,289:360,433:504)]%*%(1-RateMat[c(145:216,289:360,433:504),-c(145:216,289:360,433:504)]/12))
	TxNeed200 <- sum(Vnext[361:432]%*%(1-(RateMat-diag(RateMat))[361:432,]/12)) 
	TxNeed350 <- sum(Vnext[217:288]%*%(1-(RateMat-diag(RateMat))[217:288,]/12))

	# ART Enrollment up to end 2011
	if(t<12*61+1) {
	# Below assumes preferential uptake from CD4<200
	VH3toT3A <- max(0,min(1,(ArtHistt[t]-OnTx)/(TxNeed200+10^-6)))*(1-rowSums((RateMat-diag(RateMat))[361:432,]/12))*12
	VH2toT2A <- max(0,min(1,(ArtHistt[t]-OnTx-TxNeed200)/(TxNeed350+10^-6)))*(1-rowSums((RateMat-diag(RateMat))[217:288,]/12))*12

	# Below assumes equal probability of uptake from CD4<200 and 200-350
	VH3toT3B <- max(0,min(1,(ArtHistt[t]-OnTx)/(TxNeed200+TxNeed350+10^-6)))*(1-rowSums((RateMat-diag(RateMat))[361:432,]/12))*12
	VH2toT2B <- max(0,min(1,(ArtHistt[t]-OnTx)/(TxNeed200+TxNeed350+10^-6)))*(1-rowSums((RateMat-diag(RateMat))[217:288,]/12))*12	 } else { 

	# ART Enrollment post 2011
	
	# ART enrollment under demand constraint
	if(ARTConstr==1) {
	# Below assumes preferential uptake from CD4<200
	VH3toT3A <- max(0,min(1,(ARTVolt[t-732]-OnTx)/(TxNeed200+10^-6)))*(1-rowSums((RateMat-diag(RateMat))[361:432,]/12))*12
	VH2toT2A <- max(0,min(1,(ARTVolt[t-732]-OnTx-TxNeed200)/(TxNeed350+10^-6)))*(1-rowSums((RateMat-diag(RateMat))[217:288,]/12))*12

	# Below assumes equal probability of uptake from CD4<200 and 200-350
	VH3toT3B <- max(0,min(1,(ARTVolt[t-732]-OnTx)/(TxNeed200+TxNeed350+10^-6)))*(1-rowSums((RateMat-diag(RateMat))[361:432,]/12))*12
	VH2toT2B <- max(0,min(1,(ARTVolt[t-732]-OnTx)/(TxNeed200+TxNeed350+10^-6)))*(1-rowSums((RateMat-diag(RateMat))[217:288,]/12))*12  }	else { 
		
	# ART enrollment without demand constraint

	PctCov <- c(seq(ArtCov11,0.8,length.out=10*12),rep(0.8,21*12))[t]
	VH3toT3A <- max(0,min(1,(PctCov*sum(Vnext[73:504])-OnTx)/(TxNeed200+10^-6)))*(1-rowSums((RateMat-diag(RateMat))[361:432,]/12))*12
	VH2toT2A <- max(0,min(1,(PctCov*sum(Vnext[73:504])-OnTx-TxNeed200)/(TxNeed350+10^-6)))*(1-rowSums((RateMat-diag(RateMat))[217:288,]/12))*12

	# Below assumes equal probability of uptake from CD4<200 and 200-350
	VH3toT3B <- max(0,min(1,(PctCov*sum(Vnext[73:504])-OnTx)/(TxNeed200+TxNeed350+10^-6)))*(1-rowSums((RateMat-diag(RateMat))[361:432,]/12))*12
	VH2toT2B <- max(0,min(1,(PctCov*sum(Vnext[73:504])-OnTx)/(TxNeed200+TxNeed350+10^-6)))*(1-rowSums((RateMat-diag(RateMat))[217:288,]/12))*12  }  }

	RateMat[cbind(1:72+360,1:72+432)] <- VH3toT3A[1:72]*PriCD4200t[t] + VH3toT3B[1:72]*(1-PriCD4200t[t])
	RateMat[cbind(1:72+216,1:72+288)] <- VH2toT2A[1:72]*PriCD4200t[t] + VH2toT2B[1:72]*(1-PriCD4200t[t])  

	# for(i in 1:72) { RateMat[i+360,i+432] <- VH3toT3A[i]*PriCD4200t[t] + VH3toT3B[i]*(1-PriCD4200t[t]) }
	# for(i in 1:72) { RateMat[i+216,i+288] <- VH2toT2A[i]*PriCD4200t[t] + VH2toT2B[i]*(1-PriCD4200t[t]) }  
	
# C9. CONSTRUCT TRANSITION MATRIX
	TransMat	<- RateMat/12  # uses the rates to approximate the probabilities (means that probabilities are independent) 

	#  for (i in 1:504) { TransMat[i,]	<- (1-exp(-sum(RateMat[i,])/12))*RateMat[i,]/sum(RateMat[i,]) }	#  Calculates probabilities as cumulative incidence

	diag(TransMat) <- 1-(rowSums(TransMat)-diag(TransMat))

# C10. Calculate costs etc
	CostTxD	<- sum(Vnext[Vtemp1+5]*VTxCost[Vtemp1+5]+Vnext[Vtemp1+7]*VTxCost[Vtemp1+7])
	CostTxND	<- sum(Vnext[Vtemp1+6]*VTxCost[Vtemp1+6]+Vnext[Vtemp1+8]*VTxCost[Vtemp1+8])
	CostRegD	<- sum(Vnext[Vtemp1+5]*TxMat[9,Vtemp1+5]+Vnext[Vtemp1+7]*TxMat[9,Vtemp1+7])
	CostRegND	<- sum(Vnext[Vtemp1+6]*TxMat[9,Vtemp1+6]+Vnext[Vtemp1+8]*TxMat[9,Vtemp1+8])

	CostFalsTxD   <- sum(Vnext[Vtemp8]*rTstSL*DTestt[t]/12*FalsPosD*VTxCost[5]*1/(2+pDeft[t])*12)
	CostFalsRegD  <- sum(Vnext[Vtemp8]*rTstSL*DTestt[t]/12*FalsPosD*TxMat[9,5]*1/(2+pDeft[t])*12)
	CostFalsTxND  <- sum(Vnext[Vtemp8]*rTstSL*NDTestt[t]/12*FalsPosND*VTxCost[6]*1/(2+pDefND)*12)
	CostFalsRegND <- sum(Vnext[Vtemp8]*rTstSL*NDTestt[t]/12*FalsPosND*TxMat[9,6]*1/(2+pDefND)*12)

	CostART	<- sum(Vnext[c(145:216,289:360,433:504)])*CArt

	CostTestD  	<- sum(Vnext[-505]*Vtestfreq*DTestt[t]/12*VTestCostD)
	CostTestND 	<- sum(Vnext[-505]*Vtestfreq*NDTestt[t]/12*VTestCostND)

# C10. MATRIX MULTIPLY TO UPDATE STATE VECTOR
	Vnext	<- Vnext%*%TransMat

# C11. OUTPUTS

# State Membership
		
	Vout["NAll"]	<- sum(Vnext[-505])		 							# Total N
	Vout["Ndaly"]	<- Vnext[-505]%*%(1-VDwt)								# Total N, adjusted for YLD from HIV and TB
	Vout["NAnyTb"]	<- sum(Vnext[-505])- sum(Vnext[1+0:13*36])	 				# Any TB, incl latent infection and on treatment
	Vout["NActDis"]	<- sum(Vnext[Vtemp7])									# Active TB, incl those on treatment
	Vout["NUnTx"]	<- sum(Vnext[Vtemp9])									# Active TB, excl those on treatment
	Vout["NUnTxH"]	<- sum(Vnext[Vtemp9[21:140]])								# Active TB with HIV, excl those on treatment
	Vout["NSmP"]	<- sum(Vnext[Vtemp1+4])+sum(Vnext[Vtemp1+7])+sum(Vnext[Vtemp1+8])		# Smear positive, incl those on treatment
 
	Vout["NStr1n"]	<- sum(Vnext[rep(3:4,7)+rep(0:6*72,each=2)])					# Active TB not on tx, Pansensitive strain, tx naive
	Vout["NStr2n"]	<- sum(Vnext[rep(3:4,7)+rep(0:6*72,each=2)+7])					# Active TB not on tx, INH monores strain, tx naive
	Vout["NStr3n"]	<- sum(Vnext[rep(3:4,7)+rep(0:6*72,each=2)+14])					# Active TB not on tx, RIF monores strain, tx naive
	Vout["NStr4n"]	<- sum(Vnext[rep(3:4,7)+rep(0:6*72,each=2)+21])					# Active TB not on tx, MDR-TB strain, tx naive
	Vout["NStr5n"]	<- sum(Vnext[rep(3:4,7)+rep(0:6*72,each=2)+28])					# Active TB not on tx, MDR+ / XDR-TB strain, tx naive
	Vout["NStr1e"]	<- sum(Vnext[rep(39:40,7)+rep(0:6*72,each=2)])					# Active TB not on tx, Pansensitive strain, tx experienced
	Vout["NStr2e"]	<- sum(Vnext[rep(39:40,7)+rep(0:6*72,each=2)+7])				# Active TB not on tx, INH monores strain, tx experienced
	Vout["NStr3e"]	<- sum(Vnext[rep(39:40,7)+rep(0:6*72,each=2)+14])				# Active TB not on tx, RIF monores strain, tx experienced
	Vout["NStr4e"]	<- sum(Vnext[rep(39:40,7)+rep(0:6*72,each=2)+21])				# Active TB not on tx, MDR-TB strain, tx experienced
	Vout["NStr5e"]	<- sum(Vnext[rep(39:40,7)+rep(0:6*72,each=2)+28])				# Active TB not on tx, MDR+ / XDR-TB strain, tx experienced
	Vout["NTxD"]	<- sum(Vnext[Vtemp1+5])+sum(Vnext[Vtemp1+7])					# DOTS Treatment
	Vout["NTxND"]	<- sum(Vnext[Vtemp1+6])+sum(Vnext[Vtemp1+8])					# Non-DOTS Treatment
	Vout["NHiv"]	<- sum(Vnext[-(1:72)])									# HIV CD4 <500
	Vout["NArt"]	<- sum(Vnext[c(145:216,289:360,433,505)])						# On HAART 
	Vout["NTbH"]	<- sum(Vnext[Vtemp7[-(1:60)]])							# TB-HIV (HIV) incl those on treatment
	Vout["NTxSp"]	<- sum(Vnext[Vtemp1+7])+sum(Vnext[Vtemp1+8])					# Smear Positive on Treatment
	Vout["NMDR"]	<- sum(Vnext[rep(c(24,25,31,32),7)+rep(0:13*36,each=4)])			# MDR, Active disease
# State Transitions
	Vout["NMort"]	<- Vnext[505]										# All cause mortality
	Vout["NHivMort"]	<- as.vector(Vnext[73:504]%*%TransMat[73:504,505])				# Mortality in HIV +ve
	Vout["NTbMort"]	<- as.vector(Vnext[Vtemp7]%*%TransMat[Vtemp7,505])				# Mortality in Active TB / on treatment
	Vout["NSmPMort"]	<- as.vector(Vnext[c(Vtemp1+4,Vtemp1+7,Vtemp1+8)]%*%TransMat[c(Vtemp1+4,Vtemp1+7,Vtemp1+8),505])				
																# Mortality in Sm Pos Active TB / on treatment
	Vout["NTbHMort"]	<- as.vector(Vnext[Vtemp7[61:420]]%*%TransMat[Vtemp7[61:420],505]) 	# Mortality in TB-HIV

	Vout["NInf"]	<- sum(Vnext[rep(0:13*36,each=35)+1]*TransMat[cbind(rep(0:13*36,each=35)+1,rep(0:13*36,each=35)+rep(2:36,14))])
	# Vout["NInf"]	<- 0; for(i in 0:13*36){ Vout["NInf"]<- Vout["NInf"]+Vnext[i+1]*sum(TransMat[i+1,i+2:36])}
																# New infections (ignores superinfection)
	m <- 	rep(0:13*36,each=30*6)+rep(rep(c(1,2,9,16,23,30),each=30),14)
	Vout["NCase"]	<- sum(Vnext[m]*TransMat[cbind(m,rep(Vtemp7[1:30],6*14)+rep(0:13*36,each=30*6))])  

	# Vout["NCase"]	<- 0; for(i in 0:13*36){ for(j in c(1,2,9,16,23,30)) { Vout["NCase"]<- Vout["NCase"]+Vnext[i+j]*sum(TransMat[i+j,Vtemp7[1:30]+i])  }}
																# New TB Cases (active disease)
	Vout["NCaseNF"]	<- sum(Vnext[rep(0:1*36,each=30)+1]*TransMat[cbind(rep(0:1*36,each=30)+1,rep(Vtemp7[1:30],2)+rep(0:1*36,each=30))]) 

	# Vout["NCaseNF"]	<- 0; for(i in 0:1*36){  Vout["NCaseNF"]<- Vout["NCaseNF"]+Vnext[i+1]*sum(TransMat[i+1,Vtemp7[1:30]+i])  }
																# New TB Cases, HIV-Neg, Fast (active disease)
	Vout["NCaseHF"]<- sum(Vnext[rep(2:13*36,each=30)+1]*TransMat[cbind(rep(2:13*36,each=30)+1,rep(Vtemp7[1:30],12)+rep(2:13*36,each=30))])  

	# Vout["NCaseHF"]	<- 0; for(i in 2:13*36){  Vout["NCaseHF"]<- Vout["NCaseHF"]+Vnext[i+1]*sum(TransMat[i+1,Vtemp7[1:30]+i])  }
																# New TB Cases, HIV-Pos, Fast (active disease)
	m <- 	rep(0:1*36,each=30*5)+rep(rep(c(2,9,16,23,30),each=30),2)
	Vout["NCaseNS"]<- sum(Vnext[m]*TransMat[cbind(m,rep(Vtemp7[1:30],2*5)+rep(0:1*36,each=30*5))]) 

	# Vout["NCaseNS"]	<- 0; for(i in 0:1*36){ for(j in c(2,9,16,23,30)) { Vout["NCaseNS"]<- Vout["NCaseNS"]+Vnext[i+j]*sum(TransMat[i+j,Vtemp7[1:30]+i])  }}
																# New TB Cases, HIV-Neg, Slow (active disease)
	m <- 	rep(2:13*36,each=30*5)+rep(rep(c(2,9,16,23,30),each=30),12)
	Vout["NCaseHS"]	<- sum(Vnext[m]*TransMat[cbind(m,rep(Vtemp7[1:30],12*5)+rep(2:13*36,each=30*5))]) 

	# Vout["NCaseHS"]	<- 0; for(i in 2:13*36){ for(j in c(2,9,16,23,30)) { Vout["NCaseHS"]<- Vout["NCaseHS"]+Vnext[i+j]*sum(TransMat[i+j,Vtemp7[1:30]+i])  }}
																# New TB Cases, HIV-Pos, Slow (active disease)
	m <- 	rep(0:13*36,each=11*5)+rep(rep(c(1,2,3,9,10,16,17,23,24,30,31),each=5),14)
	Vout["NCaseIp"]<- sum(Vnext[m]*TransMat[cbind(m,rep(c(4,11,18,25,32),11*14)+rep(0:13*36,each=11*5))])  

	# Vout["NCaseIp"]	<- 0; for(i in 0:13*36){ for(j in c(1,2,3,9,10,16,17,23,24,30,31)) { Vout["NCaseIp"]<- Vout["NCaseIp"]+Vnext[i+j]*sum(TransMat[i+j,c(4,11,18,25,32)+i])  }}

	m <- 	rep(2:13*36,each=11*5)+rep(rep(c(1,2,3,9,10,16,17,23,24,30,31),each=5),12)
	Vout["NCaseIpHiv"]<- sum(Vnext[m]*TransMat[cbind(m,rep(c(4,11,18,25,32),11*12)+rep(2:13*36,each=11*5))])  			 													# New Smear-positive TB Cases (from Su,Ls and In)

	# Vout["NCaseIpHiv"] <- 0; for(i in 2:13*36){ for(j in c(1,2,3,9,10,16,17,23,24,30,31)) { Vout["NCaseIpHiv"]<- Vout["NCaseIpHiv"]+Vnext[i+j]*sum(TransMat[i+j,c(4,11,18,25,32)+i])  }}
			 													# New Smear-positive TB Cases in HIV CD4<500 (from Su,Ls and In)
	Vout["SuspctD"]	<- sum(Vnext[-505]*Vtestfreq*DTestt[t]/12)					# No suspects, DOTS programs
	Vout["SuspctND"]	<- sum(Vnext[-505]*Vtestfreq*NDTestt[t]/12)					# No suspects, Non-DOTS programs

	m <- rep(0:13*36,each=5)+rep(c(4,11,18,25,32),14)
	Vout["NCdIpD"]<- sum(Vnext[m]*TransMat[cbind(m,m+3)])

	# Vout["NCdIpD"]	<- 0; for(i in 0:13*36){ for(j in c(4,11,18,25,32)) { Vout["NCdIpD"]<- Vout["NCdIpD"]+Vnext[i+j]*sum(TransMat[i+j,j+3+i])  }}
			 													# TB Case detections, Smear Pos, DOTS,(minus losses before tx init)
	m <- rep(0:13*36,each=5)+rep(c(3,10,17,24,31),14)
	Vout["NCdInD"]<- sum(Vnext[m]*TransMat[cbind(m,m+2)])

	# Vout["NCdInD"]	<- 0; for(i in 0:13*36){ for(j in c(3,10,17,24,31)) { Vout["NCdInD"]<- Vout["NCdInD"]+Vnext[i+j]*sum(TransMat[i+j,j+2+i])  }}
			 													# TB Case detections, Smear Neg, DOTS,(minus losses before tx init)
	m <- rep(0:13*36,each=5)+rep(c(4,11,18,25,32),14)
	Vout["NCdIpND"]<- sum(Vnext[m]*TransMat[cbind(m,m+4)])

	# Vout["NCdIpND"]	<- 0; for(i in 0:13*36){ for(j in c(4,11,18,25,32)) { Vout["NCdIpND"]<- Vout["NCdIpND"]+Vnext[i+j]*sum(TransMat[i+j,j+4+i])  }}
			 													# TB Case detections, Smear Pos, NonDOTS,(minus losses before tx init)
	m <- rep(0:13*36,each=5)+rep(c(3,10,17,24,31),14)
	Vout["NCdInND"]<- sum(Vnext[m]*TransMat[cbind(m,m+3)])

	# Vout["NCdInND"]	<- 0; for(i in 0:13*36){ for(j in c(3,10,17,24,31)) { Vout["NCdInND"]<- Vout["NCdInND"]+Vnext[i+j]*sum(TransMat[i+j,j+3+i])  }}
			 													# TB Case detections, Smear Neg, NonDOTS,(minus losses before tx init)
	Vout["NCdFalsD"] 	<- sum(Vnext[Vtemp8]*rTstSL*DTestt[t]/12*FalsPosD)				# False-positive diagnoses, DOTS			
	Vout["NCdFalsND"]	<- sum(Vnext[Vtemp8]*rTstSL*NDTestt[t]/12*FalsPosND)				# False-positive diagnoses, Non-DOTS	


	m <- rep(0:6*72,each=8*16)+rep(rep(c(10:11,17:18,24:25,31:32),each=16),7)
	Vout["NTxResU"]<- sum(Vnext[m]*TransMat[cbind(m,rep(c(12:15,19:22,26:29,33:36),8*7)+rep(0:6*72,each=8*16))])

	# Vout["NTxResU"]	<- 0; for(i in 0:6*72){ for(j in c(10:11,17:18,24:25,31:32)) { Vout["NTxResU"]<- Vout["NTxResU"]+Vnext[i+j]*sum(TransMat[i+j,c(12:15,19:22,26:29,33:36)+i])  }}
			 													# Any resistance starting treatment (tx naive)
	m <- rep(0:6*72,each=4*8)+rep(rep(c(24:25,31:32),each=8),7)
	Vout["NTxMdrU"]<- sum(Vnext[m]*TransMat[cbind(m,rep(c(26:29,33:36),4*7)+rep(0:6*72,each=4*8))])

	# Vout["NTxMdrU"]	<- 0; for(i in 0:6*72){ for(j in c(24:25,31:32)) { Vout["NTxMdrU"]<- Vout["NTxMdrU"]+Vnext[i+j]*sum(TransMat[i+j,c(26:29,33:36)+i])  }}
			 													# MDR starting treatment (incl XDR)  (tx naive)
	m <- rep(0:6*72,each=2*4)+rep(rep(31:32,each=4),7)
	Vout["NTxXdrU"]<- sum(Vnext[m]*TransMat[cbind(m,rep(33:36,2*7)+rep(0:6*72,each=2*4))])

	# Vout["NTxXdrU"]	<- 0; for(i in 0:6*72){ for(j in c(31:32)) { Vout["NTxXdrU"]<- Vout["NTxXdrU"]+Vnext[i+j]*sum(TransMat[i+j,c(33:36)+i])  }}
			 													# MDR+/XDR starting treatment (tx naive)
	m <- rep(0:6*72+36,each=8*16)+rep(rep(c(10:11,17:18,24:25,31:32),each=16),7)
	Vout["NTxResR"]<- sum(Vnext[m]*TransMat[cbind(m,rep(c(12:15,19:22,26:29,33:36),8*7)+rep(0:6*72+36,each=8*16))])

	# Vout["NTxResR"]	<- 0; for(i in 0:6*72+36){ for(j in c(10:11,17:18,24:25,31:32)) { Vout["NTxResR"]<- Vout["NTxResR"]+Vnext[i+j]*sum(TransMat[i+j,c(12:15,19:22,26:29,33:36)+i])  }}
			 													# Any resistance starting treatment (tx experienced)
	m <- rep(0:6*72+36,each=4*8)+rep(rep(c(24:25,31:32),each=8),7)
	Vout["NTxMdrR"]<- sum(Vnext[m]*TransMat[cbind(m,rep(c(26:29,33:36),4*7)+rep(0:6*72+36,each=4*8))])

	# Vout["NTxMdrR"]	<- 0; for(i in 0:6*72+36){ for(j in c(24:25,31:32)) { Vout["NTxMdrR"]<- Vout["NTxMdrR"]+Vnext[i+j]*sum(TransMat[i+j,c(26:29,33:36)+i])  }}
			 													# MDR starting treatment (incl XDR) (tx experienced)
	m <- rep(0:6*72+36,each=2*4)+rep(rep(31:32,each=4),7)
	Vout["NTxXdrR"]<- sum(Vnext[m]*TransMat[cbind(m,rep(33:36,2*7)+rep(0:6*72+36,each=2*4))])

	# Vout["NTxXdrR"]	<- 0; for(i in 0:6*72+36){ for(j in c(10:11,17:18,24:25,31:32)) { Vout["NTxXdrR"]<- Vout["NTxXdrR"]+Vnext[i+j]*sum(TransMat[i+j,c(33:36)+i])  }}
			 													# MDR+/XDR starting treatment (tx experienced)

# Costs
	Vout["CostTxD"]	<- CostTxD											# Non-drug cost for treatment in DOTS programs
	Vout["CostTxND"]	<- CostTxND											# Non-drug cost for treatment in Non-DOTS programs
	Vout["CostRegD"]	<- CostRegD											# Drug cost for treatment in DOTS programs
	Vout["CostRegND"]	<- CostRegND										# Drug cost for treatment in Non-DOTS programs
	Vout["CostFalsTxD"]   <- CostFalsTxD									# Non-drug cost for treatment in DOTS programs
	Vout["CostFalsTxND"]  <- CostFalsTxND									# Non-drug cost for treatment in Non-DOTS programs
	Vout["CostFalsRegD"]  <- CostFalsRegD									# Drug cost for treatment in DOTS programs
	Vout["CostFalsRegND"] <- CostFalsRegND									# Drug cost for treatment in Non-DOTS programs
	Vout["CostART"]	<- CostART											# HAART costs
	Vout["CostTestD"]	<- CostTestD										# Diagnosis costs in DOTS programs
	Vout["CostTestND"]<- CostTestND 										# Diagnosis costs in Non-DOTS programs
	Vout["TC"]		<- CostTxD+CostTxND+CostRegD+CostRegND+CostFalsTxD+CostFalsTxND+CostFalsRegD+CostFalsRegND+CostART+CostTestD+CostTestND 
																# Total Costs
# Additional outcomes

	Vout["Check1"]	<- min(diag(TransMat))									# Check to see p(stay in state) doesnt become negative
	Vout["PfailDtx"]<- sum(Vnext[Vtemp6[1:140*2-1]]*(12/TxMat[1,Vtemp6[1:140*2-1]]*(1-TxEft[t]*TxMat[2,Vtemp6[1:140*2-1]])))/
					(sum(Vnext[Vtemp6[1:140*2-1]]*(12/TxMat[1,Vtemp6[1:140*2-1]]+12/TxMat[1,Vtemp6[1:140*2-1]]*pDeft[t]+RateMat[Vtemp6[1:140*2-1],505]))+10^-6)
																# Average failure probability in DOTS programs
	Vout["PcureDtx"] <- sum(Vnext[Vtemp6[1:140*2-1]]*(12/TxMat[1,Vtemp6[1:140*2-1]]*TxEft[t]*TxMat[2,Vtemp6[1:140*2-1]]))/
					(sum(Vnext[Vtemp6[1:140*2-1]]*(12/TxMat[1,Vtemp6[1:140*2-1]]+12/TxMat[1,Vtemp6[1:140*2-1]]*pDeft[t]+RateMat[Vtemp6[1:140*2-1],505]))+10^-6)
																# Average cure probability in DOTS programs
	Vout["PdfltDtx"] <- sum(Vnext[Vtemp6[1:140*2-1]]* 12/TxMat[1,Vtemp6[1:140*2-1]]*pDeft[t])/
					(sum(Vnext[Vtemp6[1:140*2-1]]*(12/TxMat[1,Vtemp6[1:140*2-1]]+ 12/TxMat[1,Vtemp6[1:140*2-1]]*pDeft[t]+RateMat[Vtemp6[1:140*2-1],505]))+10^-6)
																# Average default probability in DOTS programs
	Vout["PmortDtx"] <- sum(Vnext[Vtemp6[1:140*2-1]]*RateMat[Vtemp6[1:140*2-1],505])/
					(sum(Vnext[Vtemp6[1:140*2-1]]*(12/TxMat[1,Vtemp6[1:140*2-1]]+ 12/TxMat[1,Vtemp6[1:140*2-1]]*pDeft[t]+RateMat[Vtemp6[1:140*2-1],505]))+10^-6)
																# Average mortality probability in DOTS programs
	Vout["DurInfSn"]	<- 1/((sum(Vnext[Vtemp1+3]*apply(TransMat[Vtemp1+3,c(Vtemp1+2,Vtemp6,505)],1,sum))+10^-6)/(sum(Vnext[Vtemp1+3])+10^-6))/12
																# Duration of infectiousness smear negative
 	Vout["DurInfSp"]	<-1/((sum(Vnext[Vtemp1+4]*apply(TransMat[Vtemp1+4,c(Vtemp1+2,Vtemp6,505)],1,sum))+10^-6)/(sum(Vnext[Vtemp1+4])+10^-6))/12
																# Duration of infectiousness smear positive
 	Vout["DurInfAll"]	<-1/((sum(Vnext[c(Vtemp1+3,Vtemp1+4)]*apply(TransMat[c(Vtemp1+3,Vtemp1+4),c(Vtemp1+2,Vtemp6,505)],1,sum))+10^-6)/(sum(Vnext[c(Vtemp1+3,Vtemp1+4)])+10^-6))/12
																# Duration of infectiousness, all
	Vout["EffContRate"] <- sum(Vnext[c(Vtemp1+3,Vtemp1+4)]*VTrStat[c(Vtemp1+3,Vtemp1+4)])/(sum(Vnext[c(Vtemp1+3,Vtemp1+4)])+10^-6)*CRt[t]
																# Effective contact rate, untreated active disease

	Vout["ExTbC"]	<- sum(Vnext[Vtemp9[1:20]]*apply(TransMat[Vtemp9[1:20],Vtemp1+2],1,sum))
	Vout["ExTbT"]	<- sum(Vnext[Vtemp9[1:20]]*apply(TransMat[Vtemp9[1:20],Vtemp6],1,sum))
	Vout["ExTbD"]	<- sum(Vnext[Vtemp9[1:20]]*TransMat[Vtemp9[1:20],505])			# Non-HIV Exits from treatment, self-cure/treatment/death 
	Vout["ExTbCH"]	<- sum(Vnext[Vtemp9[21:140]]*apply(TransMat[Vtemp9[21:140],Vtemp1+2],1,sum))
	Vout["ExTbTH"]	<- sum(Vnext[Vtemp9[21:140]]*apply(TransMat[Vtemp9[21:140],Vtemp6],1,sum))
	Vout["ExTbDH"]	<- sum(Vnext[Vtemp9[21:140]]*TransMat[Vtemp9[21:140],505])			# HIV Exits from treatment, self-cure/treatment/death 

# Test Characteristics
	Vout["NotifD"]	<- sum(Vnext[Vtemp9]*Vtestfreq[Vtemp9]*DTestt[t]/12*TruPosDB)  +
					sum(Vnext[Vtemp8]*Vtestfreq[Vtemp8]*DTestt[t]/12*FalsPosDB)		# DOTS Notifications (true and false positive, ignoring LTFU)
	Vout["NotifND"]	<- sum(Vnext[Vtemp9]*Vtestfreq[Vtemp9]*NDTestt[t]/12*TruPosNDB)  +
					sum(Vnext[Vtemp8]*Vtestfreq[Vtemp8]*NDTestt[t]/12*FalsPosNDB)	# Non-DOTS Notifications (true and false positive, ignoring LTFU)
	Vout["PPVTb"] 	<- (sum(Vnext[Vtemp1+3]*Vtestfreq[Vtemp1+3]*DTestt[t]/12*TruPosDB[1:70*2-1]) + 
					sum(Vnext[Vtemp1+4]*Vtestfreq[Vtemp1+4]*DTestt[t]/12*TruPosDB[1:70*2]))/(Vout["NotifD"]+10^-6)
																# Positive predictive value, DOTS TB diagnosis
	Vout["NPVTb"]	<- sum(Vnext[Vtemp8]*Vtestfreq[Vtemp8]*DTestt[t]/12*(1-FalsPosDB))/((Vout["SuspctD"]-Vout["NotifD"])+10^-6)
																# Negative predictive value, TB diagnosis
	Vout["PPVRif"]	<- (sum(Vnext[Vtemp1+3]*Vtestfreq[Vtemp1+3]*DTestt[t]/12*TruPosDB[1:70*2-1]*GetXpt[Vtemp1+3]*rep(c(rep(0,2),rep(SensXpRIF,3)),14)) + 
					sum(Vnext[Vtemp1+4]*Vtestfreq[Vtemp1+4]*DTestt[t]/12*TruPosDB[1:70*2]*GetXpt[Vtemp1+4]*rep(c(rep(0,2),rep(SensXpRIF,3)),14)))/		
					((sum(Vnext[Vtemp1+3]*Vtestfreq[Vtemp1+3]*DTestt[t]/12*TruPosDB[1:70*2-1]*GetXpt[Vtemp1+3]*rep(c(rep((1-SpecXpRIF),2),rep(SensXpRIF,3)),14)) + 
					sum(Vnext[Vtemp1+4]*Vtestfreq[Vtemp1+4]*DTestt[t]/12*TruPosDB[1:70*2]*GetXpt[Vtemp1+4]*rep(c(rep((1-SpecXpRIF),2),rep(SensXpRIF,3)),14))+
					sum(Vnext[Vtemp8]*Vtestfreq[Vtemp8]*DTestt[t]/12*FalsPosDB*GetXpt[Vtemp8])*(1-SpecXpRIF))+10^-6)
																# Positive predictive value, RIF resistant diagnosis (with Xpert for all scenario)
	Vout["NPVRif"]	<- (sum(Vnext[Vtemp1+3]*Vtestfreq[Vtemp1+3]*DTestt[t]/12*TruPosDB[1:70*2-1]*GetXpt[Vtemp1+3]*rep(c(rep(SpecXpRIF,2),rep(0,3)),14)) + 
					sum(Vnext[Vtemp1+4]*Vtestfreq[Vtemp1+4]*DTestt[t]/12*TruPosDB[1:70*2]*GetXpt[Vtemp1+4]*rep(c(rep(SpecXpRIF,2),rep(0,3)),14))+
					sum(Vnext[Vtemp8]*Vtestfreq[Vtemp8]*DTestt[t]/12*FalsPosDB*GetXpt[Vtemp8])*SpecXpRIF)/		
					((sum(Vnext[Vtemp1+3]*Vtestfreq[Vtemp1+3]*DTestt[t]/12*TruPosDB[1:70*2-1]*GetXpt[Vtemp1+3]*rep(c(rep(SpecXpRIF,2),rep((1-SensXpRIF),3)),14)) + 
					sum(Vnext[Vtemp1+4]*Vtestfreq[Vtemp1+4]*DTestt[t]/12*TruPosDB[1:70*2]*GetXpt[Vtemp1+4]*rep(c(rep(SpecXpRIF,2),rep((1-SensXpRIF),3)),14))+
					sum(Vnext[Vtemp8]*Vtestfreq[Vtemp8]*DTestt[t]/12*FalsPosDB*GetXpt[Vtemp8])*SpecXpRIF)+10^-6)
																# Negative predictive value, RIF resistant diagnosis (with Xpert for all scenario)

	Vout["PDst"]	<-  (sum(Vnext[Vtemp9]*Vtestfreq[Vtemp9]*DTestt[t]/12*TruPosDB*rep(c(rep(pDstU*PhaseIn1[t],10),rep(pDstR*PhaseIn1[t],10)),7))+
					sum(Vnext[Vtemp8]*Vtestfreq[Vtemp8]*DTestt[t]/12*FalsPosDB*rep(c(rep(pDstU*PhaseIn1[t],6),rep(pDstR*PhaseIn1[t],6)),7)))/(Vout["NotifD"]+10^-6)
																# No.  getting a DST, UNDER BASECASE ALGORITHM

	Vout["GetXpt"]	<- sum(Vnext[-505]*Vtestfreq*DTestt[t]/12*GetXpt)

	Vout["ArtCov"]	<- sum(Vnext[c(145:216,289:360,433:504)])/(sum(Vnext[73:504])+10^-6)

#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%
#%#%
	return(list(Vnext=Vnext,Vout=Vout))
	} 	# End of function!
#%#%
#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%


# sum(RateMateInit-RateMat);sum(VoutInit-Vout);sum(VnextInit-Vnext)

