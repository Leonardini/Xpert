### Calibrating each country
	rm(list=ls())

	library(abind)
	setwd("C:/Users/nick/Documents/Harvard/TB Diagnostics/ANALYSIS")
	source("LIKELIHOODv3.r")
	Setting <- "Swaziland"
	setwd(paste("C:/Users/nick/Documents/Harvard/TB Diagnostics/ANALYSIS/SIR Results 3-20-2012/",Setting,sep=""))

	BcAll <- AltAll <- array(dim=c(93,92,20000)); SumAll <- matrix(NA,nrow=20000,ncol=16)
	colnames(BcAll) <- colnames(AltAll)	<- c("NAll","Ndaly","NAnyTb","NActDis","NUnTx","NUnTxH","NSmP","NStr1n","NStr2n","NStr3n","NStr4n","NStr5n","NStr1e",
					"NStr2e","NStr3e","NStr4e","NStr5e","NTxD","NTxND","NHiv","NHiv350","NArt","NTbH","NMort","NHivMort",
					"NTbMort","NSmPMort","NTbHMort","NInf","NCase","NCaseNF","NCaseNS","NCaseHF","NCaseHS","NCaseIp","NCaseIpHiv",
					"NCdIpD","NCdInD","NCdIpND","NCdInND","NTxResU","NTxMdrU","NTxXdrU","NTxResR","NTxMdrR","NTxXdrR",
					"NTxSp","SuspctD","SuspctDTB","SuspctND","CostTxD","CostTxND","CostRegD","CostRegND","CostART","CostTestD",
					"CostTestND","CostFalsTxD","CostFalsTxND","CostFalsRegD","CostFalsRegND","NCdFalsD","NCdFalsND","Check1",
					"DurInfSn","DurInfSp","DurInfAll","PfailDtx","PcureDtx","PdfltDtx","PmortDtx","EffContRate","NotifD", "NotifTBD",
					"NotifND","PPVTb","NPVTb","PPVRif","NPVRif","PDst","ExTbC","ExTbT","ExTbD","ExTbCH","ExTbTH","ExTbDH","GetXpt",
					"TC","NMDR","ArtCov","ArtNdCov","Art200Cov")

	for (i in 1:20) { 
		load(paste("BcOut",Setting,i,".rData",sep=""));  BcAll[,,1:1000+1000*(i-1)]  <- BcOut
		load(paste("AltOut",Setting,i,".rData",sep="")); AltAll[,,1:1000+1000*(i-1)] <- AltOut
		load(paste("SumOut",Setting,i,".rData",sep="")); SumAll[1:1000+1000*(i-1),]  <- SumOut 
		print(paste(i," of ",20," completed",sep=""));  flush.console()  }

	LLPrev2 <- LLIncid2 <- LLMdr <- rep(NA,20000)

	for (i in 1:20000) {
		LLPrev2[i] <- PrevLnLik(BcAll[,"NUnTx",i]/BcAll[,"NAll",i]*100000,Setting,2) 
		LLIncid2[i] <- IncidLnLik(BcAll[,"NCase",i]/BcAll[,"NAll",i]*100000*12,Setting,2) 

		MdrN   <- (BcAll[,"NStr4n",i]+BcAll[,"NStr5n",i])/(BcAll[,"NStr1n",i]+BcAll[,"NStr2n",i]+BcAll[,"NStr3n",i]+BcAll[,"NStr4n",i]+BcAll[,"NStr5n",i])
		MdrR   <- (BcAll[,"NStr4e",i]+BcAll[,"NStr5e",i])/(BcAll[,"NStr1e",i]+BcAll[,"NStr2e",i]+BcAll[,"NStr3e",i]+BcAll[,"NStr4e",i]+BcAll[,"NStr5e",i])

		LLMdr[i] <- MdrLnLik(MdrN,MdrR,Setting,2) 

		if((i/1000)==round(i/1000)) { print(paste(i," of ",20000," completed",sep=""));  flush.console() }  }

	Score1 <- exp(LLPrev2+LLIncid2+LLMdr)   	## P/I/MDR

## Resample
	seed <- 51:55; names(seed) <- c("Botswana","Lesotho","Namibia","SouthAfrica","Swaziland")
	set.seed(seed[Setting])
	resampID <- sample(1:20000,100000,replace=T,prob=Score1)
		zz1 <- unique(resampID); length(zz1)  			# No. unique resamples B,L,N,S,S = 2534, 1900, 1628, 1057, 2601
		max(table(resampID))/length(resampID)*100 		# Pct from single resample B,L,N,S,S = 1.7%, 1.9%, 3.0%, 3.2%, 3.0%
		sum(table(resampID))^2/sum(table(resampID)^2)		# ESS for estmating mean B,L,N,S,S = 269, 236, 202, 111, 178



## Saving resampID
	save(resampID,file=paste("resampID_",Setting,"NEWER1.rData",sep=""))

## Creating and saving resampID
	resampID2<- resampID
	for(i in 1:length(zz1))  { resampID2[resampID==zz1[i]] <- i }

	save(resampID2 , file=paste("resampID2_",Setting,"NEWER1.rData",sep=""))

############ Basecase
	BcSmall <- BcAll[,,zz1]
	save(BcSmall,file=paste("BcSmall",Setting,"NEWER1.Rdata",sep=""))

	BcMean <- matrix (NA,93,92)
	for(i in 1:93) { for(j in 1:92) { BcMean[i,j] <- mean(BcSmall[i,j,resampID2],na.rm=T) } }
	save(BcMean,file=paste("BcMean",Setting,"NEWER1.Rdata",sep=""))
	
	rm(BcAll,BcSmall,BcMean)

############ Alt
	AltSmall <- AltAll[,,zz1]
	save(AltSmall ,file=paste("AltSmall",Setting,"NEWER1.Rdata",sep=""))

	AltMean <- matrix (NA,93,92)
	for(i in 1:93) { for(j in 1:92) { AltMean[i,j] <- mean(AltSmall[i,j,resampID2],na.rm=T) } }
	save(AltMean ,file=paste("AltMean",Setting,"NEWER1.Rdata",sep=""))
	rm(AltAll,AltSmall,AltMean)

############ Sum
	SumSmall <- SumAll[zz1,]
	save(SumSmall ,file=paste("SumSmall",Setting,"NEWER1.Rdata",sep=""))

# Mean of resample
	SumMean <- apply(SumSmall[resampID2,],2,mean)
	save(SumMean,file=paste("SumMean",Setting,"NEWER1.Rdata",sep=""))
	rm(SumAll,SumSmall,SumMean)

############ 
############ Calculating Mean of Posterior ############ 
############ 
	paste(load("C:/Users/nick/Documents/Harvard/TB Diagnostics/ANALYSIS/PriorDrawsU3-19.rData")) #	PriorDrawsU
	ParamInit <- as.data.frame(read.csv("C:/Users/nick/Documents/Harvard/TB Diagnostics/ANALYSIS/ParamInit4.csv")[,2:6])

	MeanPost 	<- rep(NA,166); names(MeanPost) <- dimnames(PriorDrawsU)[[2]]; MeanPostU <- MeanPost

	MeanPostU 	<- apply(PriorDrawsU[resampID,],2,mean)

	for (i in 1:nrow(ParamInit))  { 
		if(ParamInit[i,5]==1) 	{ MeanPost[i] <- log(MeanPostU[i]/(1-MeanPostU[i])) }
			else 			{ MeanPost[i] <- log(MeanPostU[i])  }
		if(MeanPost[i]==Inf)  	{ MeanPost[i] <-  100   }
		if(MeanPost[i]==-Inf) 	{ MeanPost[i] <- -100   }  }

	write.csv(MeanPost,file=paste("MeanPost",Setting,"NEWER1.rData",sep=""))
	write.csv(MeanPostU,file=paste("MeanPostU",Setting,"NEWER1.rData",sep=""))


