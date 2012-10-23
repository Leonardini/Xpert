sourceA=##########   TB XPERT DIAGNOSTIC MODEL 2012  ##########
##########            ONESIM FILE            ##########
## NOTES
# a.	This is another key one. It runs one simulation, summarizes results
# b.	Needs PARAM and TIMESTEP to have been run
# c.	Takes starting vector, loops over TIMESTEP for 112 years (i.e. till the very end of 2011). At each cycle it inserts results and markov trace vectors into corresponding matrices.
# d.	For each scenario, continues simulation for 30 years from  2012 to end 2041.
# e.	Summarizes results (annualizes full results and markov trace) and calculates key quantities (10/30 year costs etc)
# f.	Outputs results for all matrices

################### START OF FUNCTION ######################
	OneSim <- function(Vparam,Setting)  {

	source("PARAMv5.r")
	source("TIMESTEPOPTv1.r")

### RUN BASECASE TO END 2012
	DIAG <- 2
  OutMat1 <- OutMatInit <- OutMat
	Vcurrent <- V1950/sum(V1950)*InitPop[Setting]*10^6
	for(t in 1:(62*12)) {
    # print(paste("Processing timestep", t))
    Out <- timestep(Vcurrent,t,ArtNdCov11,DIAG,OutMat1)
		Vcurrent <- Out$Vnext
    OutMatInit[t,] <- Out$Vout
		if (t == 732) {
      ArtNdCov11 <- Out$Vout["ArtNdCov"]
      # break
    }
  }
	V2012 <- Vcurrent

### RUN SCENARIO 2012 TO 2042
	for(j in 2:3) {
		DIAG <- j; OutMati <- OutMatInit
		Vcurrent <- V2012
	for(t in (62*12+1):(92*12)) { Out <- timestep(Vcurrent,t,ArtNdCov11,DIAG,OutMat1);Vcurrent <- Out$Vnext; OutMati[t,] <- Out$Vout  }
	assign(paste("x",j-1,sep=""),OutMati)	}

### Calculate some results
	Drt <- rep(1,30*12); for(i in 1:length(Drt)) { Drt[i] <- (1+0.03)^(-(i-1)/12) }

	BcPre12 <-  x1[744,"NActDis"]/x1[744,"NAll"]*100000
	BcInc12 <-  x1[744,"NCase"]/x1[744,"NAll"]*100000*12
	BcMdr12 <-  x1[744,"NMDR"]/x1[744,"NActDis"]*100
	BcTbM12 <- x1[744,"NTbMort"]/x1[744,"NAll"]*100000*12
	InPre22 <-  (x2[864,"NActDis"]/x2[864,"NAll"]-x1[864,"NActDis"]/x1[864,"NAll"])*100000
	InInc22 <-  (x2[864,"NCase"]/x2[864,"NAll"]-x1[864,"NCase"]/x1[864,"NAll"])*100000*12
	InMdr22 <-  (x2[864,"NMDR"]/x2[864,"NActDis"]-x1[864,"NMDR"]/x1[864,"NActDis"])*100
	InTbM22 <-  (x2[864,"NTbMort"]/x2[864,"NAll"]-x1[864,"NTbMort"]/x1[864,"NAll"])*100000*12
	InDal10d <- sum((x2[745:864,"Ndaly"]-x1[745:864,"Ndaly"])/12*Drt[1:(12*10)])
	InDal30d <- sum((x2[745:1104,"Ndaly"]-x1[745:1104,"Ndaly"])/12*Drt[1:(12*30)])
	InLY10d <- sum((x2[745:864,"NAll"]-x1[745:864,"NAll"])/12*Drt[1:(12*10)])
	InLY30d <- sum((x2[745:1104,"NAll"]-x1[745:1104,"NAll"])/12*Drt[1:(12*30)])
	InCst10d <- sum((x2[745:864,"TC"]-x1[745:864,"TC"])*Drt[1:(12*10)])
	InCst30d <- sum((x2[745:1104,"TC"]-x1[745:1104,"TC"])*Drt[1:(12*30)])
	InXpt10d <- sum((x2[745:864,"GetXpt"]-x1[745:864,"GetXpt"])*Drt[1:(12*10)])
	InXpt30d <- sum((x2[745:1104,"GetXpt"]-x1[745:1104,"GetXpt"])*Drt[1:(12*30)])

	SumRes <- c(BcPre12,BcInc12,BcMdr12,BcTbM12,InPre22,InInc22,InMdr22,InTbM22,InDal10d,
			InDal30d,InLY10d,InLY30d,InCst10d,InCst30d,InXpt10d,InXpt30d)
	names(SumRes) <- c("BcPre12","BcInc12","BcMdr12","BcTbM12","InPre22","InInc22","InMdr22",
			"InTbM22","InDal10d","InDal30d","InLY10d","InLY30d","InCst10d","InCst30d","InXpt10d","InXpt30d")

### Output results
	BcRes <- x1[c(1,1:92*12),]
	AltRes <- x2[c(1,1:92*12),]
	return(list(z1=SumRes, z2=BcRes, z3=AltRes))	}

################### END OF FUNCTION ######################