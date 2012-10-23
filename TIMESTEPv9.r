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

# NOTE: The difference between v6 and v5 is that the former does not have a local RateMat inside the timestep function, instead modifying the global RateMatStat.
# NOTE: The difference between v7 and v6 is that the former does not have a loop over z inside the timestep function, instead using explicit linear combinations.
# NOTE: The difference between v8 and v7 is that the former does not have an explicitly constructed temporary copy of the rate matrix with a zeroed-out diagonal.

######## EXTRA INDICES ##########
	i1 <- j1 <- rep(0,0); for(j in 0:6) { for (i in Vtemp1[1:10+j*10]+2) { i1[i]<-i; j1[i]<-j } }
	i1 <- as.numeric(na.omit(i1));  j1 <- as.numeric(na.omit(j1))+1
	a1 <- cbind(i1,i1+1); a2 <- cbind(i1,i1+2); a3 <- cbind(i1+1,i1+2)
	a4 <- cbind(i1+1,i1); a5 <- cbind(i1+2,i1); a6 <- cbind(1:72+72,1:72+216)
	a7 <- cbind(1:72+216,1:72+360); a8 <- cbind(73:504,rep(505,432))
	i2 <- rep(Vtemp4[1:35],8)+rep(c(5:8,5:8+36),each=35)
	a9 <- cbind(i2,rep(Vtemp4[1:35],8)+2+36)
	a10 <- cbind(1,i2); a11 <- cbind(2,i2)
	a12 <- cbind(i2,rep(Vtemp4[1:35],8)+rep(c(5:8+36,5:8+36),each=35))
	a13 <- cbind(i2,rep(Vtemp4[1:35],8)+rep(c(3,3,4,4,3,3,4,4)+36,each=35))
	a14 <- matrix(NA,7*6*8*5,7)
	i3 <- Vtemp4[rep(rep(1:5,7*6),each=8)+rep(rep(rep(0:6*5,each=5),6),each=8)]+rep(c(5:8,5:8+36),5*6*7)
	a14[,1] <- Vtemp4[rep(rep(1:5,7*6),each=8)+rep(rep(rep(0:6*5,each=5),6),each=8)]+rep(c(5:8,5:8+36),5*6*7)
	a14[,2] <- Vtemp4[2+rep(rep(rep(0:6*5,each=5),6),each=8)]+rep(c(3,3,4,4,3,3,4,4)+36,5*6*7)
	a14[,3] <- Vtemp4[3+rep(rep(rep(0:6*5,each=5),6),each=8)]+rep(c(3,3,4,4,3,3,4,4)+36,5*6*7)
	a14[,4] <- Vtemp4[4+rep(rep(rep(0:6*5,each=5),6),each=8)]+rep(c(3,3,4,4,3,3,4,4)+36,5*6*7)
	a14[,5] <- Vtemp4[4+rep(rep(rep(0:6*5,each=5),6),each=8)]+rep(c(3,3,4,4,3,3,4,4)+36,5*6*7)
	a14[,6] <- Vtemp4[4+rep(rep(rep(0:6*5,each=5),6),each=8)]+rep(c(3,3,4,4,3,3,4,4)+36,5*6*7)
	a14[,7] <- Vtemp4[5+rep(rep(rep(0:6*5,each=5),6),each=8)]+rep(c(3,3,4,4,3,3,4,4)+36,5*6*7)
	i4 <- rep(0:6,each=10)+1; i5 <- rep(0:6,each=2*5*5)+1
	a15 <- cbind(rep(rep(c(1,37),each=5),7)+(i4-1)*72,rep(rep(c(1,37),each=5),7)+(i4-1)*72+1+rep(0:4,14)*7)
	a16 <- cbind(rep(rep(c(1,37),each=5),7)+(i4-1)*72,rep(rep(c(1,37),each=5),7)+(i4-1)*72+2+rep(0:4,14)*7)
	a17 <- cbind(rep(rep(c(1,37),each=5),7)+(i4-1)*72,rep(rep(c(1,37),each=5),7)+(i4-1)*72+3+rep(0:4,14)*7)
	a18 <- cbind(rep(rep(c(2,38),each=5*5),7)+(i5-1)*72+rep(0:4*7,7*2*5),rep(rep(c(2,38),each=5*5),7)+(i5-1)*72+0+rep(rep(0:4,each=5),7*2)*7)
	a19 <- cbind(rep(rep(c(2,38),each=5*5),7)+(i5-1)*72+rep(0:4*7,7*2*5),rep(rep(c(2,38),each=5*5),7)+(i5-1)*72+1+rep(rep(0:4,each=5),7*2)*7)
	a20 <- cbind(rep(rep(c(2,38),each=5*5),7)+(i5-1)*72+rep(0:4*7,7*2*5),rep(rep(c(2,38),each=5*5),7)+(i5-1)*72+2+rep(rep(0:4,each=5),7*2)*7)
	a21 <- cbind(Vtemp1[1:70]+3,Vtemp1[1:70]+5); a22 <- cbind(Vtemp1[1:70]+3,Vtemp1[1:70]+6)
	a23 <- cbind(Vtemp1[1:70]+4,Vtemp1[1:70]+7); a24 <- cbind(Vtemp1[1:70]+4,Vtemp1[1:70]+8)
	a25 <- cbind(1:72,1:72+72); a26 <- cbind(1:72+360,1:72+432); a27 <- cbind(1:72+216,1:72+288)
	i6 <- rep(0:13*36,each=35)+1; a28 <- cbind(rep(0:13*36,each=35)+1,rep(0:13*36,each=35)+rep(2:36,14))
	i7 <- rep(0:13*36,each=30*6)+rep(rep(c(1,2,9,16,23,30),each=30),14); a29 <- cbind(rep(0:13*36,each=30*6)+rep(rep(c(1,2,9,16,23,30),each=30),14),rep(Vtemp7[1:30],6*14)+rep(0:13*36,each=30*6))
	i8 <- rep(0:1*36,each=30)+1; a30 <- cbind(rep(0:1*36,each=30)+1,rep(Vtemp7[1:30],2)+rep(0:1*36,each=30)) 
	i9 <- rep(2:13*36,each=30)+1; a31 <- cbind(rep(2:13*36,each=30)+1,rep(Vtemp7[1:30],12)+rep(2:13*36,each=30))	
	i10 <- rep(0:1*36,each=30*5)+rep(rep(c(2,9,16,23,30),each=30),2); a32 <- cbind(rep(0:1*36,each=30*5)+rep(rep(c(2,9,16,23,30),each=30),2),rep(Vtemp7[1:30],2*5)+rep(0:1*36,each=30*5))
	i11 <- rep(2:13*36,each=30*5)+rep(rep(c(2,9,16,23,30),each=30),12); a33 <- cbind(rep(2:13*36,each=30*5)+rep(rep(c(2,9,16,23,30),each=30),12),rep(Vtemp7[1:30],12*5)+rep(2:13*36,each=30*5))
	i12 <- rep(0:13*36,each=11*5)+rep(rep(c(1,2,3,9,10,16,17,23,24,30,31),each=5),14); a34 <- cbind(rep(0:13*36,each=11*5)+rep(rep(c(1,2,3,9,10,16,17,23,24,30,31),each=5),14),rep(c(4,11,18,25,32),11*14)+rep(0:13*36,each=11*5))
	i13 <- rep(2:13*36,each=11*5)+rep(rep(c(1,2,3,9,10,16,17,23,24,30,31),each=5),12); a35 <- cbind(rep(2:13*36,each=11*5)+rep(rep(c(1,2,3,9,10,16,17,23,24,30,31),each=5),12),rep(c(4,11,18,25,32),11*12)+rep(2:13*36,each=11*5))
	i14 <- rep(0:13*36,each=5)+rep(c(4,11,18,25,32),14); a36 <- cbind(rep(0:13*36,each=5)+rep(c(4,11,18,25,32),14),rep(0:13*36,each=5)+rep(c(4,11,18,25,32),14)+3)
	i15 <- rep(0:13*36,each=5)+rep(c(3,10,17,24,31),14); a37 <- cbind(rep(0:13*36,each=5)+rep(c(3,10,17,24,31),14),rep(0:13*36,each=5)+rep(c(3,10,17,24,31),14)+2)
	i16 <- rep(0:13*36,each=5)+rep(c(4,11,18,25,32),14); a38 <- cbind(rep(0:13*36,each=5)+rep(c(4,11,18,25,32),14),rep(0:13*36,each=5)+rep(c(4,11,18,25,32),14)+4)
	i17 <- rep(0:13*36,each=5)+rep(c(3,10,17,24,31),14); a39 <- cbind(rep(0:13*36,each=5)+rep(c(3,10,17,24,31),14),rep(0:13*36,each=5)+rep(c(3,10,17,24,31),14)+3)
	i18 <- rep(0:6*72,each=8*16)+rep(rep(c(10:11,17:18,24:25,31:32),each=16),7); a40 <- cbind(rep(0:6*72,each=8*16)+rep(rep(c(10:11,17:18,24:25,31:32),each=16),7),rep(c(12:15,19:22,26:29,33:36),8*7)+rep(0:6*72,each=8*16))
	i19 <- rep(0:6*72,each=4*8)+rep(rep(c(24:25,31:32),each=8),7); a41 <- cbind(rep(0:6*72,each=4*8)+rep(rep(c(24:25,31:32),each=8),7),rep(c(26:29,33:36),4*7)+rep(0:6*72,each=4*8))
	i20 <- rep(0:6*72,each=2*4)+rep(rep(31:32,each=4),7); a42 <- cbind(rep(0:6*72,each=2*4)+rep(rep(31:32,each=4),7),rep(33:36,2*7)+rep(0:6*72,each=2*4))
	i21 <- rep(0:6*72+36,each=8*16)+rep(rep(c(10:11,17:18,24:25,31:32),each=16),7); a43 <- cbind(rep(0:6*72+36,each=8*16)+rep(rep(c(10:11,17:18,24:25,31:32),each=16),7),rep(c(12:15,19:22,26:29,33:36),8*7)+rep(0:6*72+36,each=8*16))
	i22 <- rep(0:6*72+36,each=4*8)+rep(rep(c(24:25,31:32),each=8),7); a44 <- cbind(rep(0:6*72+36,each=4*8)+rep(rep(c(24:25,31:32),each=8),7),rep(c(26:29,33:36),4*7)+rep(0:6*72+36,each=4*8))
	i23 <- rep(0:6*72+36,each=2*4)+rep(rep(31:32,each=4),7); a45 <- cbind(rep(0:6*72+36,each=2*4)+rep(rep(31:32,each=4),7),rep(33:36,2*7)+rep(0:6*72+36,each=2*4))
	a46 <- cbind(1:505,1:505)

########################################################
########################################################

### B. SETTING UP STATIC PARTS OF RATE MATRIX (FOR SPEED OF TIMESTEP FUNCTION)
	RateMatStat<- matrix(0,nrow=505,ncol=505); rownames(RateMatStat) <- StatNam; colnames(RateMatStat) <- StatNam

### B1. BREAKDOWN TO ACTIVE DISEASE (Stay in HIV / Resistance / Treatment subdivisions) 
 	RateMatStat[a1] <- VrBreakD[j1]*(1-VpToIp[j1]);  RateMatStat[a2] <- VrBreakD[j1]*VpToIp[j1]

### B2. SMEAR NEG CONVERT TO SMEAR POS (Stay in HIV / Resistance / Treatment subdivisions) 
	RateMatStat[a3]  <- rNtoP 

### B3. SPONTANEOUS CURE (Stay in HIV / Resistance / Treatment subdivisions) 
	RateMatStat[a4] <- VrIToLs[j1];  RateMatStat[a5] <- VrIToLs[j1]

### B4. HIV Progression 
	RateMatStat[a6] <- H1toH2;  RateMatStat[a7] <- H2toH3  

### B5. POPULATE MORTALITY RATES
	RateMatStat[a8] 			<- rep(VmuHIV,each=72)			# HIV mortality
	Vtemp1shift3 = Vtemp1 + 3
	Vtemp1shift4 = Vtemp1 + 4
	RateMatStat[Vtemp1shift3,505]	<- RateMatStat[Vtemp1shift3,505] + muIn	# Untreated Smear-neg TB mortality
	RateMatStat[Vtemp1shift4,505]	<- RateMatStat[Vtemp1shift4,505] + muIp	# Untreated Smear-pos TB mortality

	VTrStatz			<- rep(0,504) 				# Creates a vector for contact rates
	VTrStatz[Vtemp1shift3]	<- rep(TrIn*RelFit,14) 	# Contact rates for smear neg 
	VTrStatz[Vtemp1shift4]	<- rep(RelFit,14) 			# Contact rates for smear pos 
	
#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#% 
#%#%
	timestep	<- function(Vcurrent,t,ArtNdCov11,DIAG,OutMat1)  {   # Start of function!
#%#%
#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%

### C1. ADD NEW ENTRANTS TO STATE VECTOR
	Vnext 		<- Vcurrent		# Initializes new vector
	Vnext[505] 		<- 0			# Clears out deaths
	Vnext[1]		<- Vnext[1]	+ NewEntt[t]*1000000	  # Adds new entrants to NU0Su based on birth rate
	
  ### Compute the relative weights of each diacgnostic
  coeff1 = (1 - PhaseIn1[t]) * (1 - PhaseIn2[t])
  coeff2 =  PhaseIn1[t] * (1 - PhaseIn2[t])
  coeff3 = 0
  if (DIAG == 1) {
    coeff1 = coeff1 + PhaseIn2[t]
  }
  else if (DIAG == 2) {
    coeff2 = coeff2 + PhaseIn2[t]
  }
  else { # DIAG = 3
    coeff3 = coeff3 + PhaseIn2[t]
  }

### TREATMENT TRANSITIONS, BY ALGORITHM  

### B5. POPULATE MORTALITY RATES 
	RateMatStat[Vtemp6,505] <- RateMatStat[Vtemp6,505] + (coeff1/TxMatAlg1[1,Vtemp6] + coeff2/TxMatAlg2[1,Vtemp6] + coeff3/TxMatAlg3[1,Vtemp6]) *2*rep(c(muIn,muIn,muIp,muIp),70)*TunTxMort	# Treatment TB mortality
	RateMatStat[Vtemp9[101:120],505]	<- muTBH			# TB-HIV mortality for CD4 350

# Vector of current contact rates...
	VTrStatz[Vtemp6] = (1 - (coeff1 * TxMatAlg1[2,Vtemp6] + coeff2 * TxMatAlg2[2,Vtemp6] + coeff3 * TxMatAlg3[2,Vtemp6]) * TxEft[t]) * rep(rep(c(TrIn, TrIn, 1, 1), 5) * rep(RelFit, each=4), 14)
											# Contact rates for individuals on treatment

### B5. TREATMENT OUTCOMES (Stay in HIV / Resistance subdivisions) 
	RateMatStat[a9]  <- RateMatStat[a9] + (coeff1*TxMatAlg1[a11]/TxMatAlg1[a10] + coeff2*TxMatAlg2[a11]/TxMatAlg2[a10] + coeff3*TxMatAlg3[a11]/TxMatAlg3[a10])*(12*TxEft[t])   # Cures back to Ls state, treatment experienced subdivision
	RateMatStat[a12] <-  RateMatStat[a12] + (coeff1*(1-TxEft[t]*TxMatAlg1[a11])/TxMatAlg1[a10] +coeff2*(1-TxEft[t]*TxMatAlg2[a11])/TxMatAlg2[a10] + coeff3*(1-TxEft[t]*TxMatAlg3[a11]) / TxMatAlg3[a10]) * (12*pReTx)   # Failures identified and reinitiated on treatment, treatment experienced subdivision
  vec1 = apply(cbind(0, 12 / TxMatAlg1[a10] * (rep(rep(c(pDeft[t], pDefND), each=35), 4) + (1 - TxEft[t] * TxMatAlg1[a11]) * (1-pReTx)) - colSums(TxMatAlg1[3:8, i2])), 1, max)
  vec2 = apply(cbind(0, 12 / TxMatAlg2[a10] * (rep(rep(c(pDeft[t], pDefND), each=35), 4) + (1 - TxEft[t] * TxMatAlg2[a11]) * (1-pReTx)) - colSums(TxMatAlg2[3:8, i2])), 1, max)
  vec3 = apply(cbind(0, 12 / TxMatAlg3[a10] * (rep(rep(c(pDeft[t], pDefND), each=35), 4) + (1 - TxEft[t] * TxMatAlg3[a11]) * (1-pReTx)) - colSums(TxMatAlg3[3:8, i2])), 1, max)
	RateMatStat[a13] <-  RateMatStat[a13] + (coeff1 * vec1 + coeff2 * vec2 + coeff3 * vec3)   # Defaulters and failures to active disease
 	for(k in 1:6)  {
   RateMatStat[a14[,c(1,k+1)]] <- RateMatStat[a14[,c(1,k+1)]] + coeff1 * TxMatAlg1[k+2,i3] + coeff2 * TxMatAlg2[k+2,i3] + coeff3 * TxMatAlg3[k+2,i3]
  }  # Defaulters and failures to active disease with Acquired Resistance

### C3. UPDATE MORTALITY RATES WITH BACKGROUND MORTALITY
	RateMatStat[-505,505]	<- RateMatStat[-505,505]+mubt[t]

### C4. TB INCIDENCE (Can change strain subdivision, stay in HIV / treatment subd.) 
	VInf 		<- Vnext[1:504]/sum(Vnext[1:504])*VTrStatz*CRt[t]			# P(meet carrier)*CR|carrier, homogeneous mixing
	m <- c(sum(VInf[Vtemp2+0*7]),sum(VInf[Vtemp2+1*7]),sum(VInf[Vtemp2+2*7]),sum(VInf[Vtemp2+3*7]),sum(VInf[Vtemp2+4*7])); m <- rep(m,14)

		RateMatStat[a15]	<- RateMatStat[a15]+m*(1-Vpfast[i4])
		RateMatStat[a16]	<- RateMatStat[a16]+m*Vpfast[i4]*(1-VpToIp[i4])
		RateMatStat[a17]	<- RateMatStat[a17]+m*Vpfast[i4]*VpToIp[i4]

### C5. SUPERINFECTION (Can change strain subdivision, stay in HIV / treatment subd.) 
	VSupInf 		<- VInf*(1-rep(VPartIm,each=72))  # As above, with partial immunity, homogeneous mixing
	v <- c(sum(VSupInf[Vtemp2+0*7]),sum(VSupInf[Vtemp2+1*7]),sum(VSupInf[Vtemp2+2*7]),sum(VSupInf[Vtemp2+3*7]),sum(VSupInf[Vtemp2+4*7])); v <- v[rep(rep(1:5,each=5),7*2)]
 
		RateMatStat[a18]	<- RateMatStat[a18]+v*(1-Vpfast[i5])
		RateMatStat[a19]	<- RateMatStat[a19]+v*Vpfast[i5]*(1-VpToIp[i5])
		RateMatStat[a20]	<- RateMatStat[a20]+v*Vpfast[i5]*VpToIp[i5]

### C6. DIAGNOSIS AND TREATMENT STRATEGY (Stay in HIV / Resistance / Treatment subdivisions) 
#  C6a. Specifying diagnosis and treatment as a result of algorithm
	TxMat		  = coeff1 * TxMatAlg1 + coeff2 * TxMatAlg2 + coeff3 * TxMatAlg3
	TruPosD	  = coeff1 * TruPosDAlg1 + coeff2 * TruPosDAlg2 + coeff3 * TruPosDAlg3
	FalsPosD	= coeff1 * FalsPosDAlg1 + coeff2 * FalsPosDAlg2 + coeff3 * FalsPosDAlg3
	TruPosDB	= coeff1 * TruPosDAlgB1 + coeff2 * TruPosDAlgB2 + coeff3 * TruPosDAlgB3
	FalsPosDB	= coeff1 * FalsPosDAlgB1 + coeff2 * FalsPosDAlgB2 + coeff3 * FalsPosDAlgB3
	VTestCostD= coeff1 * VTestCostD1 + coeff2 * VTestCostD2 + coeff3 * VTestCostD3
	VTxCost   = coeff1 * VTxCost1 + coeff2 * VTxCost2 + coeff3 * VTxCost3
	GetXpt	  = coeff1 * GetXpt1 + coeff2 * GetXpt2 + coeff3 * GetXpt3

#  C6b. Diagnosis and tx initiation
	RateMatStat[a21]	<-	DTestt[t]*TruPosD[1:70*2-1]*rTstIn 		# From In to Tn1
	RateMatStat[a22]	<-	NDTestt[t]*TruPosND[1:70*2-1]*rTstIn  	# From In to Tn2
	RateMatStat[a23]	<-	DTestt[t]*TruPosD[1:70*2]   			# From Ip to Tp1
	RateMatStat[a24]	<-	NDTestt[t]*TruPosND[1:70*2]   		# From Ip to Tp2

### C8. HIV INCIDENCE and ART ENROLLMENT
	# HIV incidence
	RateMatStat[a25] 	<- rHIVt[t]
	
	firstInds  = 217:288
	secondInds = 361:432
	RMDiag  = RateMatStat[a46]
	RMDiag1 = RMDiag[firstInds]
	RMtemp1 = RateMatStat[firstInds,  ]
	RMDiag2 = RMDiag[secondInds]
  RMtemp2 = RateMatStat[secondInds, ]

	RMrowsum = c(rowSums(RMtemp1), rowSums(RMtemp2)) - c(RMDiag1, RMDiag2)
	HAARTInds = c(145:216, 289:360, 433:504)
	OnTx = sum(Vnext[HAARTInds]) - sum(Vnext[HAARTInds] %*% RateMatStat[HAARTInds, -HAARTInds]/12)
	firstNext = Vnext[firstInds]
	TxNeed350 = sum(firstNext)  - (sum(firstNext %*% RMtemp1)  - sum(firstNext  * RMDiag1)) / 12
	secondNext = Vnext[secondInds]
	TxNeed200 = sum(secondNext) - (sum(secondNext %*% RMtemp2) - sum(secondNext * RMDiag2)) / 12

  ##### ART Enrollment up to end 2011 
	if(t<(12*61+1)) {
	# Below assumes preferential uptake from CD4<200
	VH3toT3A <- max(0,min(1,(ArtHistt[t]-OnTx)/(TxNeed200+10^-6)))*(12-RMrowsum[73:144])
	VH2toT2A <- max(0,min(1,(ArtHistt[t]-OnTx-TxNeed200)/(TxNeed350+10^-6)))*(12-RMrowsum[1:72])

	# Below assumes equal probability of uptake from CD4<200 and 200-350
	VH3toT3B <- max(0,min(1,(ArtHistt[t]-OnTx)/(TxNeed200+TxNeed350+10^-6)))*(12-RMrowsum[73:144])
	VH2toT2B <- max(0,min(1,(ArtHistt[t]-OnTx)/(TxNeed200+TxNeed350+10^-6)))*(12-RMrowsum[1:72])	   } else {

  ##### ART Enrollment post 2011 
	# ART enrollment under demand constraint
	if(ARTConstr==1) {
	# Below assumes preferential uptake from CD4<200
	VH3toT3A <- max(0,min(1,(ARTVolt[t-732]-OnTx)/(TxNeed200+10^-6)))*(12-RMrowsum[73:144])
	VH2toT2A <- max(0,min(1,(ARTVolt[t-732]-OnTx-TxNeed200)/(TxNeed350+10^-6)))*(12-RMrowsum[1:72])

	# Below assumes equal probability of uptake from CD4<200 and 200-350
	VH3toT3B <- max(0,min(1,(ARTVolt[t-732]-OnTx)/(TxNeed200+TxNeed350+10^-6)))*(12-RMrowsum[73:144])
	VH2toT2B <- max(0,min(1,(ARTVolt[t-732]-OnTx)/(TxNeed200+TxNeed350+10^-6)))*(12-RMrowsum[1:72])    } 
		
  ##### ART enrollment without demand constraint 
	if(ARTConstr==0|(ARTConstr==2&DIAG==2)) { 
	PctCov <- c(seq(ArtNdCov11,ArtFutCov,length.out=10*12),rep(ArtFutCov,21*12))[t-732]
	VH3toT3A <- max(0,min(1,(PctCov*sum(Vnext[217:504])-OnTx)/(TxNeed200+10^-6)))*(12-RMrowsum[73:144])
	VH2toT2A <- max(0,min(1,(PctCov*sum(Vnext[217:504])-OnTx-TxNeed200)/(TxNeed350+10^-6)))*(12-RMrowsum[1:72])

	# Below assumes equal probability of uptake from CD4<200 and 200-350
	VH3toT3B <- max(0,min(1,(PctCov*sum(Vnext[217:504])-OnTx)/(TxNeed200+TxNeed350+10^-6)))*(12-RMrowsum[73:144])
	VH2toT2B <- max(0,min(1,(PctCov*sum(Vnext[217:504])-OnTx)/(TxNeed200+TxNeed350+10^-6)))*(12-RMrowsum[1:72])   } 

 ##### ART enrollment with constraint formed by Status Quo
	if(ARTConstr==2&DIAG==3) {

	xxx <- OutMat1[,"NArt"]
	VH3toT3A <- max(0,min(1,(xxx[t]-OnTx)/(TxNeed200+10^-6)))*(12-RMrowsum[73:144])
	VH2toT2A <- max(0,min(1,(xxx[t]-OnTx-TxNeed200)/(TxNeed350+10^-6)))*(12-RMrowsum[1:72])

	# Below assumes equal probability of uptake from CD4<200 and 200-350
	VH3toT3B <- max(0,min(1,(xxx[t]-OnTx)/(TxNeed200+TxNeed350+10^-6)))*(12-RMrowsum[73:144])
	VH2toT2B <- max(0,min(1,(xxx[t]-OnTx)/(TxNeed200+TxNeed350+10^-6)))*(12-RMrowsum[1:72])    }   }

	RateMatStat[a26] <- VH3toT3A*PriCD4200t[t] + VH3toT3B*(1-PriCD4200t[t])
	RateMatStat[a27] <- VH2toT2A*PriCD4200t[t] + VH2toT2B*(1-PriCD4200t[t])  
 
# C9. CONSTRUCT TRANSITION MATRIX 
	TransMat	<- RateMatStat/12  # uses the rates to approximate the probabilities (means that probabilities are independent) 
	TransMat[a46] <- 1-(rowSums(TransMat)-TransMat[a46])

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
  Dead = Vnext[505]
  Alive = Vnext[-505]
  NAlive = sum(Alive)
  
	Vout["NAll"]	<- NAlive		 							                                # Total N
	Vout["Ndaly"]	<- Alive%*%(1-VDwt)								                        # Total N, adjusted for YLD from HIV and TB
	Vout["NAnyTb"]	<- NAlive - sum(Vnext[1+0:13*36])	 				              # Any TB, incl latent infection and on treatment
	Vout["NActDis"]	<- sum(Vnext[Vtemp7])									                  # Active TB, incl those on treatment
	Vout["NUnTx"]	<- sum(Vnext[Vtemp9])									                    # Active TB, excl those on treatment
	Vout["NUnTxH"]	<- sum(Vnext[Vtemp9[21:140]])								            # Active TB with HIV, excl those on treatment
	
	smearPositiveInds = c(Vtemp1+4,Vtemp1+7,Vtemp1+8)
	Vout["NSmP"]	<- sum(Vnext[smearPositiveInds])		                      # Smear positive, incl those on treatment

  indexVector1 = rep(3:4,7)+rep(0:6*72,each=2)
	Vout["NStr1n"]	<- sum(Vnext[indexVector1])					                    # Active TB not on tx, Pansensitive strain, tx naive
	Vout["NStr2n"]	<- sum(Vnext[indexVector1+7])					                  # Active TB not on tx, INH monores strain, tx naive
	Vout["NStr3n"]	<- sum(Vnext[indexVector1+14])					                # Active TB not on tx, RIF monores strain, tx naive
	Vout["NStr4n"]	<- sum(Vnext[indexVector1+21])					                # Active TB not on tx, MDR-TB strain, tx naive
	Vout["NStr5n"]	<- sum(Vnext[indexVector1+28])					                # Active TB not on tx, MDR+ / XDR-TB strain, tx naive
	indexVector2 = rep(39:40,7)+rep(0:6*72,each=2)
	Vout["NStr1e"]	<- sum(Vnext[indexVector2])					                    # Active TB not on tx, Pansensitive strain, tx experienced
	Vout["NStr2e"]	<- sum(Vnext[indexVector2+7])				                    # Active TB not on tx, INH monores strain, tx experienced
	Vout["NStr3e"]	<- sum(Vnext[indexVector2+14])				                  # Active TB not on tx, RIF monores strain, tx experienced
	Vout["NStr4e"]	<- sum(Vnext[indexVector2+21])				                  # Active TB not on tx, MDR-TB strain, tx experienced
	Vout["NStr5e"]	<- sum(Vnext[indexVector2+28])				                  # Active TB not on tx, MDR+ / XDR-TB strain, tx experienced
	Vout["NTxD"]		<- sum(Vnext[Vtemp1+5])+sum(Vnext[Vtemp1+7])					  # DOTS Treatment
	Vout["NTxND"]		<- sum(Vnext[Vtemp1+6])+sum(Vnext[Vtemp1+8])					  # Non-DOTS Treatment
	Vout["NHiv"]		<- sum(Vnext[73:504])									                  # HIV 
	Vout["NHiv350"]	<- sum(Vnext[217:504])									                # HIV CD4 <350
	Vout["NArt"]		<- sum(Vnext[HAARTInds])						                    # On HAART 
	Vout["NTbH"]		<- sum(Vnext[Vtemp7[-(1:60)]])							            # TB-HIV (HIV) incl those on treatment
	Vout["NTxSp"]		<- sum(Vnext[Vtemp1+7])+sum(Vnext[Vtemp1+8])					  # Smear Positive on Treatment
	Vout["NMDR"]		<- sum(Vnext[rep(c(24,25,31,32),7)+rep(0:13*36,each=4)])# MDR, Active disease, not on treatment
# State Transitions 
	Vout["NMort"]		<- Dead										                              # All cause mortality
	allMortality = Vnext * TransMat[,505]
	Vout["NHivMort"]	<- sum(allMortality[73:504])		                      # Mortality in HIV +ve
	Vout["NTbMort"]	  <- sum(allMortality[Vtemp7])		                      # Mortality in Active TB / on treatment
	Vout["NSmPMort"]	<- sum(allMortality[smearPositiveInds])		            # Mortality in Sm Pos Active TB / on treatment
	Vout["NTbHMort"]	<- sum(allMortality[Vtemp7[61:420]]) 	                # Mortality in TB-HIV

	Vout["NInf"]		<- sum(Vnext[i6]*TransMat[a28])                 # New infections (ignores superinfection)
	Vout["NCase"]		<- sum(Vnext[i7]*TransMat[a29])  								# New TB Cases (active disease)
	Vout["NCaseNF"]	<- sum(Vnext[i8]*TransMat[a30]) 								# New TB Cases, HIV-Neg, Fast (active disease)
	Vout["NCaseHF"]	<- sum(Vnext[i9]*TransMat[a31])  								# New TB Cases, HIV-Pos, Fast (active disease)
	Vout["NCaseNS"]	<- sum(Vnext[i10]*TransMat[a32]) 								# New TB Cases, HIV-Neg, Slow (active disease)
	Vout["NCaseHS"]	<- sum(Vnext[i11]*TransMat[a33]) 								# New TB Cases, HIV-Pos, Slow (active disease)
	Vout["NCaseIp"]	<- sum(Vnext[i12]*TransMat[a34]) 								# New Smear-positive TB Cases (from Su,Ls and In)
	Vout["NCaseIpHiv"]	<- sum(Vnext[i13]*TransMat[a35])						# New Smear-positive TB Cases in HIV CD4<500 (from Su,Ls and In)
	Vout["SuspctD"]	<- sum(Alive*Vtestfreq*DTestt[t]/12)						# No suspects, DOTS programs
	Vout["SuspctDTB"]	<- sum((Alive*Vtestfreq*DTestt[t]/12)[Vtemp9])# No suspects, DOTS programs
	Vout["SuspctND"]	<- sum(Alive*Vtestfreq*NDTestt[t]/12)					# No suspects, Non-DOTS programs
	Vout["NCdIpD"]	<- sum(Vnext[i14]*TransMat[a36])								# TB Case detections, Smear Pos, DOTS,(minus losses before tx init)
	Vout["NCdInD"]	<- sum(Vnext[i15]*TransMat[a37])								# TB Case detections, Smear Neg, DOTS,(minus losses before tx init)
	Vout["NCdIpND"]	<- sum(Vnext[i16]*TransMat[a38])								# TB Case detections, Smear Pos, NonDOTS,(minus losses before tx init)
	Vout["NCdInND"]	<- sum(Vnext[i17]*TransMat[a39])								# TB Case detections, Smear Neg, NonDOTS,(minus losses before tx init)
	Vout["NCdFalsD"] 	<- sum(Vnext[Vtemp8]*rTstSL*DTestt[t]/12*FalsPosD)					# False-positive diagnoses, DOTS ,(minus losses before tx init)		
	Vout["NCdFalsND"]	<- sum(Vnext[Vtemp8]*rTstSL*NDTestt[t]/12*FalsPosND)					# False-positive diagnoses, Non-DOTS ,(minus losses before tx init)
	Vout["NTxResU"]	<- sum(Vnext[i18]*TransMat[a40])								# Any resistance starting treatment (tx naive)
	Vout["NTxMdrU"]	<- sum(Vnext[i19]*TransMat[a41])								# MDR starting treatment (incl XDR)  (tx naive)
	Vout["NTxXdrU"]	<- sum(Vnext[i20]*TransMat[a42])								# MDR+/XDR starting treatment (tx naive)
	Vout["NTxResR"]	<- sum(Vnext[i21]*TransMat[a43])								# Any resistance starting treatment (tx experienced)
	Vout["NTxMdrR"]	<- sum(Vnext[i22]*TransMat[a44])								# MDR starting treatment (incl XDR) (tx experienced)
	Vout["NTxXdrR"]	<- sum(Vnext[i23]*TransMat[a45])								# MDR+/XDR starting treatment (tx experienced)

# Costs
	Vout["CostTxD"]	<- CostTxD											        # Non-drug cost for treatment in DOTS programs
	Vout["CostTxND"]	<- CostTxND											      # Non-drug cost for treatment in Non-DOTS programs
	Vout["CostRegD"]	<- CostRegD											      # Drug cost for treatment in DOTS programs
	Vout["CostRegND"]	<- CostRegND										      # Drug cost for treatment in Non-DOTS programs
	Vout["CostFalsTxD"]   <- CostFalsTxD									  # Non-drug cost for treatment in DOTS programs
	Vout["CostFalsTxND"]  <- CostFalsTxND									  # Non-drug cost for treatment in Non-DOTS programs
	Vout["CostFalsRegD"]  <- CostFalsRegD									  # Drug cost for treatment in DOTS programs
	Vout["CostFalsRegND"] <- CostFalsRegND									# Drug cost for treatment in Non-DOTS programs
	Vout["CostART"]	<- CostART											        # HAART costs
	Vout["CostTestD"]	<- CostTestD										      # Diagnosis costs in DOTS programs
	Vout["CostTestND"]<- CostTestND 										    # Diagnosis costs in Non-DOTS programs
	Vout["TC"]		<- CostTxD+CostTxND+CostRegD+CostRegND+CostFalsTxD+CostFalsTxND+CostFalsRegD+CostFalsRegND+CostART+CostTestD+CostTestND        # Vectorize?
																# Total Costs
# Additional outcomes
  eps = 10^-6
	Vout["Check1"]	<- min(TransMat[a46])									# Check to see p(stay in state) doesnt become negative
	specialInds = Vtemp6[1:140*2-1]
	VSpecialInds = Vnext[specialInds]
	term0 = sum(VSpecialInds * (12/TxMat[1,specialInds] * TxEft[t] * TxMat[2,specialInds]))
	term1 = sum(VSpecialInds * RateMatStat[specialInds,505])
	term2 = sum(VSpecialInds * 12/TxMat[1,specialInds])
	denominator = term1 + term2 * (1 + pDeft[t]) + eps
	Vout["PfailDtx"] <- (term2 - term0) / denominator    # Average failure probability in DOTS programs
	Vout["PcureDtx"] <- term0 / denominator	             # Average cure probability in DOTS programs
	Vout["PdfltDtx"] <- term2 * pDeft[t] / denominator   # Average default probability in DOTS programs
	Vout["PmortDtx"] <- term1 / denominator						   # Average mortality probability in DOTS programs
																
  Vtemp1shift3and4 = c(Vtemp1+3,Vtemp1+4)
  rowSumProducts = Vnext * rowSums(TransMat[,c(Vtemp1+2,Vtemp6,505)])
	Vout["DurInfSn"]	<- (1/12)/((sum(rowSumProducts[Vtemp1shift3])    + eps) /(sum(Vnext[Vtemp1shift3])     + eps))                        # Duration of infectiousness smear negative
 	Vout["DurInfSp"]	<- (1/12)/((sum(rowSumProducts[Vtemp1shift4])    + eps) /(sum(Vnext[Vtemp1shift4])     + eps ))                       # Duration of infectiousness smear positive
 	Vout["DurInfAll"]	<- (1/12)/((sum(rowSumProducts[Vtemp1shift3and4]) + eps)/(sum(Vnext[Vtemp1shift3and4]) + eps))                        # Duration of infectiousness, all
	Vout["EffContRate"] <- sum(Vnext[Vtemp1shift3and4]*VTrStatz[Vtemp1shift3and4])/(sum(Vnext[Vtemp1shift3and4])+eps)*CRt[t]    					  # Effective contact rate, untreated active disease

  Vtemp9prefix = Vtemp9[1:20]
  subMat2 = TransMat[Vtemp9prefix,]
  Vnext2  = Vnext[Vtemp9prefix]
	Vout["ExTbC"]	<- sum(Vnext2 * rowSums(subMat2[,Vtemp1+2]))
	Vout["ExTbT"]	<- sum(Vnext2 * rowSums(subMat2[,Vtemp6]))
	Vout["ExTbD"]	<- sum(Vnext2 * subMat2[,505])			# Non-HIV Exits from active TB, self-cure/treatment/death
  
  Vtemp9suffix = Vtemp9[21:140]
  subMat3 = TransMat[Vtemp9suffix,]
  Vnext3  = Vnext[Vtemp9suffix]
	Vout["ExTbCH"]	<- sum(Vnext3 * rowSums(subMat3[,Vtemp1+2]))
	Vout["ExTbTH"]	<- sum(Vnext3 * rowSums(subMat3[,Vtemp6]))
	Vout["ExTbDH"]	<- sum(Vnext3 * subMat3[,505])			# HIV Exits from active TB, self-cure/treatment/death 

# Test Characteristics
  NotifTBD = sum(Vnext[Vtemp9]*Vtestfreq[Vtemp9]*DTestt[t]/12*TruPosDB)
	Vout["NotifD"]	<- NotifTBD + sum(Vnext[Vtemp8]*Vtestfreq[Vtemp8]*DTestt[t]/12*FalsPosDB)		# DOTS Notifications (true and false positive, ignoring LTFU)
	Vout["NotifTBD"]	<- 	NotifTBD	# DOTS Notifications (true positive, ignoring LTFU)
	Vout["NotifND"]	<- sum(Vnext[Vtemp9]*Vtestfreq[Vtemp9]*NDTestt[t]/12*TruPosNDB)  +
					sum(Vnext[Vtemp8]*Vtestfreq[Vtemp8]*NDTestt[t]/12*FalsPosNDB)	# Non-DOTS Notifications (true and false positive, ignoring LTFU)
	Vout["PPVTb"] 	<- (sum(Vnext[Vtemp1+3]*Vtestfreq[Vtemp1+3]*DTestt[t]/12*TruPosDB[1:70*2-1]) + 
					sum(Vnext[Vtemp1+4]*Vtestfreq[Vtemp1+4]*DTestt[t]/12*TruPosDB[1:70*2]))/(Vout["NotifD"]+eps)
																# Positive predictive value, DOTS TB diagnosis
	Vout["NPVTb"]	<- sum(Vnext[Vtemp8]*Vtestfreq[Vtemp8]*DTestt[t]/12*(1-FalsPosDB))/((Vout["SuspctD"]-Vout["NotifD"])+eps)
																# Negative predictive value, TB diagnosis
	Vout["PPVRif"]	<- (sum(Vnext[Vtemp1+3]*Vtestfreq[Vtemp1+3]*DTestt[t]/12*TruPosDB[1:70*2-1]*GetXpt[Vtemp1+3]*rep(c(rep(0,2),rep(SensXpRIF,3)),14)) + 
					sum(Vnext[Vtemp1+4]*Vtestfreq[Vtemp1+4]*DTestt[t]/12*TruPosDB[1:70*2]*GetXpt[Vtemp1+4]*rep(c(rep(0,2),rep(SensXpRIF,3)),14)))/		
					((sum(Vnext[Vtemp1+3]*Vtestfreq[Vtemp1+3]*DTestt[t]/12*TruPosDB[1:70*2-1]*GetXpt[Vtemp1+3]*rep(c(rep((1-SpecXpRIF),2),rep(SensXpRIF,3)),14)) + 
					sum(Vnext[Vtemp1+4]*Vtestfreq[Vtemp1+4]*DTestt[t]/12*TruPosDB[1:70*2]*GetXpt[Vtemp1+4]*rep(c(rep((1-SpecXpRIF),2),rep(SensXpRIF,3)),14))+
					sum(Vnext[Vtemp8]*Vtestfreq[Vtemp8]*DTestt[t]/12*FalsPosDB*GetXpt[Vtemp8])*(1-SpecXpRIF))+eps)
																# Positive predictive value, RIF resistant diagnosis (with Xpert for all scenario)
	Vout["NPVRif"]	<- (sum(Vnext[Vtemp1+3]*Vtestfreq[Vtemp1+3]*DTestt[t]/12*TruPosDB[1:70*2-1]*GetXpt[Vtemp1+3]*rep(c(rep(SpecXpRIF,2),rep(0,3)),14)) + 
					sum(Vnext[Vtemp1+4]*Vtestfreq[Vtemp1+4]*DTestt[t]/12*TruPosDB[1:70*2]*GetXpt[Vtemp1+4]*rep(c(rep(SpecXpRIF,2),rep(0,3)),14))+
					sum(Vnext[Vtemp8]*Vtestfreq[Vtemp8]*DTestt[t]/12*FalsPosDB*GetXpt[Vtemp8])*SpecXpRIF)/		
					((sum(Vnext[Vtemp1+3]*Vtestfreq[Vtemp1+3]*DTestt[t]/12*TruPosDB[1:70*2-1]*GetXpt[Vtemp1+3]*rep(c(rep(SpecXpRIF,2),rep((1-SensXpRIF),3)),14)) + 
					sum(Vnext[Vtemp1+4]*Vtestfreq[Vtemp1+4]*DTestt[t]/12*TruPosDB[1:70*2]*GetXpt[Vtemp1+4]*rep(c(rep(SpecXpRIF,2),rep((1-SensXpRIF),3)),14))+
					sum(Vnext[Vtemp8]*Vtestfreq[Vtemp8]*DTestt[t]/12*FalsPosDB*GetXpt[Vtemp8])*SpecXpRIF)+eps)
																# Negative predictive value, RIF resistant diagnosis (with Xpert for all scenario)

	Vout["PDst"]	<-  (sum(Vnext[Vtemp9]*Vtestfreq[Vtemp9]*DTestt[t]/12*TruPosDB*rep(c(rep(pDstU*PhaseIn1[t],10),rep(pDstR*PhaseIn1[t],10)),7))+
					sum(Vnext[Vtemp8]*Vtestfreq[Vtemp8]*DTestt[t]/12*FalsPosDB*rep(c(rep(pDstU*PhaseIn1[t],6),rep(pDstR*PhaseIn1[t],6)),7)))/(Vout["NotifD"]+eps)
																# No.  getting a DST, UNDER BASECASE ALGORITHM

	Vout["GetXpt"]	  <- sum(Alive*Vtestfreq*DTestt[t]/12*GetXpt)                             
	Vout["ArtCov"]	  <- sum(Vnext[HAARTInds])/sum(Vnext[73:504])
	Vout["ArtNdCov"]	<- sum(Vnext[c(289:360,433:504)])/sum(Vnext[217:504])
	Vout["Art200Cov"]	<- sum(Vnext[433:504])/sum(Vnext[361:504])

#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#% 
#%#%
	return(list(Vnext=Vnext,Vout=Vout))
	} 	# End of function!
#%#%
#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%


# sum(RateMateInit-RateMat);sum(VoutInit-Vout);sum(VnextInit-Vnext)

