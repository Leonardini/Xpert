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
	RateMatStat[Vtemp1+3,505]	<- RateMatStat[Vtemp1+3,505] + muIn	# Untreated Smear-neg TB mortality
	RateMatStat[Vtemp1+4,505]	<- RateMatStat[Vtemp1+4,505] + muIp	# Untreated Smear-pos TB mortality

	VTrStatz			<- rep(0,504) 				# Creates a vector for contact rates
	VTrStatz[Vtemp1+3]	<- rep(TrIn*RelFit,14) 			# Contact rates for smear neg 
	VTrStatz[Vtemp1+4]	<- rep(RelFit,14) 			# Contact rates for smear pos 
	
#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#% 
#%#%
	timestep	<- function(Vcurrent,t,ArtNdCov11,DIAG,OutMat1)  {   # Start of function!
#%#%
#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%

### C1. ADD NEW ENTRANTS TO STATE VECTOR
	Vnext 		<- Vcurrent		# Initializes new vector
	Vnext[505] 		<- 0			# Clears out deaths
	Vnext[1]		<- Vnext[1]	+ NewEntt[t]*1000000	  # Adds new entrants to NU0Su based on birth rate

### TREATMENT TRANSITIONS, BY ALGORITHM  
	for (z in 1:Nalg) {  # this needs to be for all algorithms z
	RateMatStatz	<- RateMatStat
	TxMat			<- get(paste("TxMatAlg",z,sep=""))

### B5. POPULATE MORTALITY RATES 
	RateMatStatz[Vtemp6,505] <- RateMatStatz[Vtemp6,505] + 2/TxMat[1,Vtemp6]*rep(c(muIn,muIn,muIp,muIp),70)*TunTxMort	# Treatment TB mortality
  	RateMatStatz[Vtemp9[101:120],505]	<- muTBH			# TB-HIV mortality for CD4 350

# Vector of current contact rates...
	VTrStatz[Vtemp6]		<- (1-TxMat[2,Vtemp6]*TxEft[t])*rep(rep(c(TrIn,TrIn,1,1),5)*rep(RelFit,each=4),14)
											# Contact rates for individuals on treatment
	assign(paste('VTrStat',z,sep=""),VTrStatz)

### B5. TREATMENT OUTCOMES (Stay in HIV / Resistance subdivisions) 
	RateMatStatz[a9]  <-  RateMatStatz[a9] + 12/TxMat[a10]*TxMat[a11]*TxEft[t]   # Cures back to Ls state, treatment experienced subdivision
	RateMatStatz[a12] <-  RateMatStatz[a12] + 12/TxMat[a10]*(1-TxEft[t]*TxMat[a11])*pReTx   # Failures identified and reinitiated on treatment, treatment experienced subdivision
	RateMatStatz[a13] <-  RateMatStatz[a13] + apply(cbind(0,12/TxMat[a10]*(rep(rep(c(pDeft[t],pDefND),each=35),4)+(1-TxEft[t]*TxMat[a11])*(1-pReTx))-colSums(TxMat[3:8,i2])),1,max)   # Defaulters and failures to active disease
 	for(k in 1:6)  { RateMatStatz[a14[,c(1,k+1)]] <- RateMatStatz[a14[,c(1,k+1)]]+TxMat[k+2,i3] }  # Defaulters and failures to active disease with Acquired Resistance
	assign(paste("RateMatStat",z,sep=""),RateMatStatz)	} 

### C2. CREATE RATE MATRIX AND TRANSITION PARAMETER VECTOR
	RateMatStat 	<- (RateMatStat1*(1-PhaseIn1[t])+RateMatStat2*PhaseIn1[t])*(1-PhaseIn2[t]) +
				get(paste("RateMatStat",DIAG,sep=""))*PhaseIn2[t]
	VTrStat	<- (VTrStat1*(1-PhaseIn1[t])+VTrStat2*PhaseIn1[t])*(1-PhaseIn2[t]) +
				get(paste("VTrStat",DIAG,sep=""))*PhaseIn2[t]

### C3. UPDATE MORTALITY RATES WITH BACKGROUND MORTALITY
	RateMatStat[-505,505]	<- RateMatStat[-505,505]+mubt[t]

### C4. TB INCIDENCE (Can change strain subdivision, stay in HIV / treatment subd.) 
	VInf 		<- Vnext[1:504]/sum(Vnext[1:504])*VTrStat*CRt[t]			# P(meet carrier)*CR|carrier, homogeneous mixing
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
	RateMatStat[a21]	<-	DTestt[t]*TruPosD[1:70*2-1]*rTstIn 		# From In to Tn1
	RateMatStat[a22]	<-	NDTestt[t]*TruPosND[1:70*2-1]*rTstIn  	# From In to Tn2
	RateMatStat[a23]	<-	DTestt[t]*TruPosD[1:70*2]   			# From Ip to Tp1
	RateMatStat[a24]	<-	NDTestt[t]*TruPosND[1:70*2]   		# From Ip to Tp2

### C8. HIV INCIDENCE and ART ENROLLMENT
	# HIV incidence
	RateMatStat[a25] 	<- rHIVt[t]
	RMtemp <- RateMatStat; diag(RMtemp) <- 0
	RMrowsum <- rowSums(RMtemp[c(217:288,361:432),])
	OnTx <- sum(Vnext[c(145:216,289:360,433:504)])-sum(Vnext[c(145:216,289:360,433:504)]%*%(RateMatStat[c(145:216,289:360,433:504),-c(145:216,289:360,433:504)]/12))
	TxNeed200 <- sum(Vnext[361:432])-sum(Vnext[361:432]%*%(RMtemp[361:432,]/12)) 
	TxNeed350 <- sum(Vnext[217:288])-sum(Vnext[217:288]%*%(RMtemp[217:288,]/12))

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
	Vout["NTxD"]		<- sum(Vnext[Vtemp1+5])+sum(Vnext[Vtemp1+7])					# DOTS Treatment
	Vout["NTxND"]		<- sum(Vnext[Vtemp1+6])+sum(Vnext[Vtemp1+8])					# Non-DOTS Treatment
	Vout["NHiv"]		<- sum(Vnext[73:504])									# HIV 
	Vout["NHiv350"]		<- sum(Vnext[217:504])									# HIV CD4 <350
	Vout["NArt"]		<- sum(Vnext[c(145:216,289:360,433:504)])						# On HAART 
	Vout["NTbH"]		<- sum(Vnext[Vtemp7[-(1:60)]])							# TB-HIV (HIV) incl those on treatment
	Vout["NTxSp"]		<- sum(Vnext[Vtemp1+7])+sum(Vnext[Vtemp1+8])					# Smear Positive on Treatment
	Vout["NMDR"]		<- sum(Vnext[rep(c(24,25,31,32),7)+rep(0:13*36,each=4)])			# MDR, Active disease, not on treatment
# State Transitions 
	Vout["NMort"]		<- Vnext[505]										# All cause mortality
	Vout["NHivMort"]	<- as.vector(Vnext[73:504]%*%TransMat[73:504,505])				# Mortality in HIV +ve
	Vout["NTbMort"]	<- as.vector(Vnext[Vtemp7]%*%TransMat[Vtemp7,505])				# Mortality in Active TB / on treatment
	Vout["NSmPMort"]	<- as.vector(Vnext[c(Vtemp1+4,Vtemp1+7,Vtemp1+8)]%*%TransMat[c(Vtemp1+4,Vtemp1+7,Vtemp1+8),505])				
																# Mortality in Sm Pos Active TB / on treatment
	Vout["NTbHMort"]	<- as.vector(Vnext[Vtemp7[61:420]]%*%TransMat[Vtemp7[61:420],505]) 	# Mortality in TB-HIV

	Vout["NInf"]		<- sum(Vnext[i6]*TransMat[a28])
																# New infections (ignores superinfection)
	Vout["NCase"]		<- sum(Vnext[i7]*TransMat[a29])  								# New TB Cases (active disease)
	Vout["NCaseNF"]	<- sum(Vnext[i8]*TransMat[a30]) 								# New TB Cases, HIV-Neg, Fast (active disease)
	Vout["NCaseHF"]	<- sum(Vnext[i9]*TransMat[a31])  								# New TB Cases, HIV-Pos, Fast (active disease)
	Vout["NCaseNS"]	<- sum(Vnext[i10]*TransMat[a32]) 								# New TB Cases, HIV-Neg, Slow (active disease)
	Vout["NCaseHS"]	<- sum(Vnext[i11]*TransMat[a33]) 								# New TB Cases, HIV-Pos, Slow (active disease)
	Vout["NCaseIp"]	<- sum(Vnext[i12]*TransMat[a34]) 								# New Smear-positive TB Cases (from Su,Ls and In)
	Vout["NCaseIpHiv"]	<- sum(Vnext[i13]*TransMat[a35])							# New Smear-positive TB Cases in HIV CD4<500 (from Su,Ls and In)
	Vout["SuspctD"]	<- sum(Vnext[-505]*Vtestfreq*DTestt[t]/12)						# No suspects, DOTS programs
	Vout["SuspctDTB"]	<- sum((Vnext[-505]*Vtestfreq*DTestt[t]/12)[Vtemp9])					# No suspects, DOTS programs

	Vout["SuspctND"]	<- sum(Vnext[-505]*Vtestfreq*NDTestt[t]/12)						# No suspects, Non-DOTS programs
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

	Vout["Check1"]	<- min(TransMat[a46])									# Check to see p(stay in state) doesnt become negative
	Vout["PfailDtx"]<- sum(Vnext[Vtemp6[1:140*2-1]]*(12/TxMat[1,Vtemp6[1:140*2-1]]*(1-TxEft[t]*TxMat[2,Vtemp6[1:140*2-1]])))/
					(sum(Vnext[Vtemp6[1:140*2-1]]*(12/TxMat[1,Vtemp6[1:140*2-1]]+12/TxMat[1,Vtemp6[1:140*2-1]]*pDeft[t]+RateMatStat[Vtemp6[1:140*2-1],505]))+10^-6)
																# Average failure probability in DOTS programs
	Vout["PcureDtx"] <- sum(Vnext[Vtemp6[1:140*2-1]]*(12/TxMat[1,Vtemp6[1:140*2-1]]*TxEft[t]*TxMat[2,Vtemp6[1:140*2-1]]))/
					(sum(Vnext[Vtemp6[1:140*2-1]]*(12/TxMat[1,Vtemp6[1:140*2-1]]+12/TxMat[1,Vtemp6[1:140*2-1]]*pDeft[t]+RateMatStat[Vtemp6[1:140*2-1],505]))+10^-6)
																# Average cure probability in DOTS programs
	Vout["PdfltDtx"] <- sum(Vnext[Vtemp6[1:140*2-1]]* 12/TxMat[1,Vtemp6[1:140*2-1]]*pDeft[t])/
					(sum(Vnext[Vtemp6[1:140*2-1]]*(12/TxMat[1,Vtemp6[1:140*2-1]]+ 12/TxMat[1,Vtemp6[1:140*2-1]]*pDeft[t]+RateMatStat[Vtemp6[1:140*2-1],505]))+10^-6)
																# Average default probability in DOTS programs
	Vout["PmortDtx"] <- sum(Vnext[Vtemp6[1:140*2-1]]*RateMatStat[Vtemp6[1:140*2-1],505])/
					(sum(Vnext[Vtemp6[1:140*2-1]]*(12/TxMat[1,Vtemp6[1:140*2-1]]+ 12/TxMat[1,Vtemp6[1:140*2-1]]*pDeft[t]+RateMatStat[Vtemp6[1:140*2-1],505]))+10^-6)
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
	Vout["ExTbD"]	<- sum(Vnext[Vtemp9[1:20]]*TransMat[Vtemp9[1:20],505])			# Non-HIV Exits from active TB, self-cure/treatment/death 
	Vout["ExTbCH"]	<- sum(Vnext[Vtemp9[21:140]]*apply(TransMat[Vtemp9[21:140],Vtemp1+2],1,sum))
	Vout["ExTbTH"]	<- sum(Vnext[Vtemp9[21:140]]*apply(TransMat[Vtemp9[21:140],Vtemp6],1,sum))
	Vout["ExTbDH"]	<- sum(Vnext[Vtemp9[21:140]]*TransMat[Vtemp9[21:140],505])			# HIV Exits from active TB, self-cure/treatment/death 

# Test Characteristics
	Vout["NotifD"]	<- sum(Vnext[Vtemp9]*Vtestfreq[Vtemp9]*DTestt[t]/12*TruPosDB)  +
					sum(Vnext[Vtemp8]*Vtestfreq[Vtemp8]*DTestt[t]/12*FalsPosDB)		# DOTS Notifications (true and false positive, ignoring LTFU)
	Vout["NotifTBD"]	<- sum(Vnext[Vtemp9]*Vtestfreq[Vtemp9]*DTestt[t]/12*TruPosDB)		# DOTS Notifications (true positive, ignoring LTFU)


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

	Vout["ArtCov"]	<- sum(Vnext[c(145:216,289:360,433:504)])/sum(Vnext[73:504])
	Vout["ArtNdCov"]	<- sum(Vnext[c(289:360,433:504)])/sum(Vnext[217:504])
	Vout["Art200Cov"]	<- sum(Vnext[433:504])/sum(Vnext[361:504])

#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#% 
#%#%
	return(list(Vnext=Vnext,Vout=Vout))
	} 	# End of function!
#%#%
#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%


# sum(RateMateInit-RateMat);sum(VoutInit-Vout);sum(VnextInit-Vnext)

