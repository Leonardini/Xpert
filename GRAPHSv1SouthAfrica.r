##########   TB XPERT DIAGNOSTIC MODEL 2012  ##########   
##########            GRAPHS FILE            ########## 
	# source("BURNINv2.r")

######### SIMULATIONS ######### 
	rm(list=ls())
	setwd("C:/Users/nick/Documents/Harvard/TB Diagnostics/ANALYSIS")
	MeanParam <- read.csv(file="MeanParam.csv",row.names=1)
	Vparam <- MeanParam ;  Setting <- "SouthAfrica"
	source("PARAMv5.r")
	source("TIMESTEPv4.r")

### RUN BASECASE TO END 2012
 	DIAG <- 2; OutMatInit <- OutMat
	V1950 <- dget("V1950.Rdata")
	Vcurrent <- V1950/sum(V1950)*InitPop[Setting]*10^6 
	system.time( for(t in 1:(62*12)) { Out <- timestep(Vcurrent,t,ArtNdCov11,DIAG)
		Vcurrent <- Out$Vnext; OutMatInit[t,] <- Out$Vout
				if(t==732) {ArtNdCov11 <- Out$Vout["ArtNdCov"]  } } )
	V2012 <- Vcurrent 

### RUN SCENARIO 2012 TO 2042
	system.time(  for(j in 2:3) {
	DIAG <- j
	OutMati <- OutMatInit
	Vcurrent <- V2012
	for(t in (62*12+1):(92*12)) { Out <- timestep(Vcurrent,t,ArtNdCov11,DIAG);Vcurrent <- Out$Vnext; OutMati[t,] <- Out$Vout  }
	assign(paste("OutMat",j-1,sep=""),OutMati)	} )


######### GRAPHS ######### 

### READ IN VALIDATION DATA
	CheckDat	<- data.matrix(read.csv(paste("CheckDat2_",Setting,".csv",sep=""), nrows=94)[,1:29])
	colnames(CheckDat)
	colnames(OutMat)

###############111111111111111111111111111111111111###############################
###############111111111111111111111111111111111111###############################
###############111111111111111111111111111111111111###############################
	par(mfrow=c(4,2));par(mar=c(3,5,2,2))

### Demographics check
	#1 Pop Size
		plot(1:(92*12)/12+1950,OutMat1[1:(92*12),"NAll"]*10^-6, col="red", main="Adult Population (Mil)",
			ylim=c(0,50),ylab="Pop (Mil)",xlab="Year",axes=FALSE, type="l",lwd=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,45,by=5),las=1); box()
			abline(h=0:20*5,col="grey85")
			lines(1:(92*12)/12+1950,OutMat1[1:(92*12),"NAll"]*10^-6, col="red",lwd=2)
			lines(1950:2043+0.5,CheckDat[,2], col="blue",lwd=2)
			legend("topleft",c("UN Pop Division","Model"),col= c("blue","red"),lwd=3,cex=0.75,bg="white")

	#2 HIV Positive Prev
		plot((20*12):(92*12)/12+1950,OutMat1[(20*12):(92*12),"NHiv"]/OutMat1[(20*12):(92*12),"NAll"]*100, col="red", 
			main="Adult HIV Prevalence (%)",ylim=c(0,20),ylab="Prev (%)",xlab="Year",axes=FALSE, type="l",lwd=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,25,by=5),las=1); box()
			abline(h=0:5*5,col="grey85")
			lines((20*12):(92*12)/12+1950,OutMat1[(20*12):(92*12),"NHiv"]/OutMat1[(20*12):(92*12),"NAll"]*100, col="red",lwd=2)
			lines((20*12):(92*12)/12+1950,OutMat1[(20*12):(92*12),"NHiv350"]/OutMat1[(20*12):(92*12),"NAll"]*100, col="purple",lwd=2)
			lines(1950:2043+0.5,CheckDat[,3]/CheckDat[,2]*100, col="blue",lwd=2)
			lines(1950:2043+0.5,CheckDat[,5]/CheckDat[,2]*100, col="darkgreen",lwd=2)
			legend("topleft",c("UNAIDS: HIV +ve","UNAIDS: CD4<350","Model: HIV +ve","Model : CD4<500"),
			col= c("blue","darkgreen","red","purple"),lwd=3,cex=0.75,bg="white")
	#3 HAART Total
		plot((40*12):(92*12)/12+1950,OutMat1[(40*12):(92*12),"NArt"]*10^-6, col="red", main="Adult HAART Population (Mil)",
			ylim=c(0,3),ylab="Pop (Mil)",xlab="Year",axes=FALSE, type="l",lwd=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,20,by=0.5),las=1); box()
			abline(h=0:20*0.5,col="grey85")
			lines((40*12):(92*12)/12+1950,OutMat1[(40*12):(92*12),"NArt"]*10^-6, col="red",lwd=2)
			lines(1950:2043+0.5,CheckDat[,6], col="blue",lwd=2)
			legend("topleft",c("UNAIDS","Model"),col= c("blue","red"),lwd=3,cex=0.75,bg="white")

	#4 TB, Active Disease per 100,000, excluding treatment 
		plot((1*12):(92*12)/12+1950,OutMat1[(1*12):(92*12),"NUnTx"]/OutMat1[(1*12):(92*12),"NAll"]*100000,
			col="red", main="TB Prevalence (excl. TB treatment)",cex.main=1,
			ylim=c(0,2000),ylab="Active TB per 100K",xlab="Year",axes=FALSE, type="l",lwd=2)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,4000,by=1000),las=1); box()
			abline(h=0:10*500,col="grey85")
			lines((1*12):(92*12)/12+1950,OutMat1[(1*12):(92*12),"NUnTx"]/OutMat1[(1*12):(92*12),"NAll"]*100000, col="red",lwd=2)
			lines(1950:2043+0.5,CheckDat[,7], col="blue",lwd=2)
			lines(1950:2043+0.5,CheckDat[,8], col="blue",lwd=2,lty=2)
			lines(1950:2043+0.5,CheckDat[,9], col="blue",lwd=2,lty=2)
			legend("topright",c("WHO Estimates","Model Projection"),
			col= c("blue","red"),lwd=3,cex=0.75,bg="white")

	#5 TB, Active Disease per 100,000, by HIV / Non-HIV
		TBPrevInHIV		<- OutMat1[(20*12):(92*12),"NTbH"]/OutMat1[(20*12):(92*12),"NHiv"]
		TBPrevInNonHIV	<- (OutMat1[(20*12):(92*12),"NActDis"]-OutMat1[(20*12):(92*12),"NTbH"])/
					(OutMat1[(20*12):(92*12),"NAll"]-OutMat1[(20*12):(92*12),"NHiv"])

		plot((20*12):(92*12)/12+1950,OutMat1[(20*12):(92*12),"NActDis"]/OutMat1[(20*12):(92*12),"NAll"]*100000, col="black", 
			main="TB Prevalence, By HIV Status",ylim=c(0,6500),ylab="Active TB per 100K",
			xlab="Year",axes=FALSE, type="l",lwd=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,20000,by=1000),las=1); box()
			abline(h=0:20*1000,col="grey85")
			lines((20*12):(92*12)/12+1950,OutMat1[(20*12):(92*12),"NActDis"]/OutMat1[(20*12):(92*12),"NAll"]*100000, col="black",lwd=2)
			lines((20*12):(92*12)/12+1950,TBPrevInHIV*100000, col="red",lwd=2)
			lines((20*12):(92*12)/12+1950,TBPrevInNonHIV*100000, col="blue",lwd=2)
			legend("topleft",c("General","HIV","Non-HIV"),
			col= c("black","red","blue"),lwd=3,cex=0.75,bg="white")

	#6 TB Distribution across subgroups (on tx x smear status)
		TbTxSmP <- OutMat1[(20*12):(92*12),"NTxSp"]/OutMat1[(20*12):(92*12),"NAll"]*100000
		TbTxSmN <- (OutMat1[(20*12):(92*12),"NTxD"]+OutMat1[(20*12):(92*12),"NTxND"]-
				OutMat1[(20*12):(92*12),"NTxSp"])/OutMat1[(20*12):(92*12),"NAll"]*100000
		TbNoTxSmP <- (OutMat1[(20*12):(92*12),"NSmP"]-TbTxSmP)/OutMat1[(20*12):(92*12),"NAll"]*100000
		TbNoTxSmN <- (OutMat1[(20*12):(92*12),"NActDis"]-TbTxSmP-TbTxSmN-TbNoTxSmP)/OutMat1[(20*12):(92*12),"NAll"]*100000

		plot((20*12):(92*12)/12+1950,TbTxSmP, col="darkgreen", 
			main="TB Prevalence, By Smear and Tx Status",ylim=c(0,2200),ylab="No. per 100K",
			xlab="Year",axes=FALSE, type="l",lwd=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,8000,by=500),las=1); box()
			abline(h=0:10*500,col="grey85")
			lines((20*12):(92*12)/12+1950,TbTxSmP, col="darkgreen",lwd=2)
			lines((20*12):(92*12)/12+1950,TbTxSmN , col="red",lwd=2)
			lines((20*12):(92*12)/12+1950,TbNoTxSmP , col="blue",lwd=2)
			lines((20*12):(92*12)/12+1950,TbNoTxSmN , col="purple",lwd=2)
			legend("topright",c("Smear Positive, on Treatment","Smear Negative, on Treatment","Smear Positive, No Treatment",
				"Smear Negative, No Treatment"),col= c("darkgreen","red","blue","purple"),lwd=3,cex=0.75,bg="white")

	#7 TB, Incidence per 100,000 
		plot((20*12):(92*12)/12+1950,OutMat1[(20*12):(92*12),"NCase"]/OutMat1[(20*12):(92*12),"NAll"]*100000*12, col="red", 
			main="TB Incidence (per 100,000)",ylim=c(0,1400),ylab="TB Incidence per 100K",
			xlab="Year",axes=FALSE, type="l",lwd=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,3000,by=500,las=1)); box()
			abline(h=0:10*500,col="grey85")
			lines((20*12):(92*12)/12+1950,OutMat1[(20*12):(92*12),"NCase"]/OutMat1[(20*12):(92*12),"NAll"]*100000*12, col="red",lwd=2)
			lines(1950:2043+0.5,CheckDat[,10], col="blue",lwd=2)
			lines(1950:2043+0.5,CheckDat[,11], col="blue",lwd=2,lty=2)
			lines(1950:2043+0.5,CheckDat[,12], col="blue",lwd=2,lty=2)
			legend("topleft",c("WHO Projection","Model Projection"),
			col= c("blue","red"),lwd=3,cex=0.75,bg="white")

	#8 TB, Incidence, recent period, by route and HIV status
		plot((20*12):(92*12)/12+1950,OutMat1[(20*12):(92*12),"NCaseNF"]*12/1000, col="forestgreen", ylim=c(0,150),
			main="TB Incidence by Fast/Slow Breakdown and HIV Status",cex.main=1,
			ylab="Annual TB Incidence (1,000s)",xlab="Year",axes=FALSE, type="l",lwd=2)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,10000,by=20),las=1); box()
			abline(h=0:100*20,col="grey85")
			lines((20*12):(92*12)/12+1950,OutMat1[(20*12):(92*12),"NCaseNF"]*12/1000, col="forestgreen",lwd=2)
			lines((20*12):(92*12)/12+1950,OutMat1[(20*12):(92*12),"NCaseNS"]*12/1000, col="blue",lwd=2)
			lines((20*12):(92*12)/12+1950,OutMat1[(20*12):(92*12),"NCaseHF"]*12/1000, col="red",lwd=2)
			lines((20*12):(92*12)/12+1950,OutMat1[(20*12):(92*12),"NCaseHS"]*12/1000, col="purple",lwd=2)
			legend("topleft",c("HIV Neg, Fast Breakdown","HIV Neg, Slow Breakdown","HIV Pos, Fast Breakdown","HIV Pos, Slow Breakdown"),
			col= c("forestgreen","blue","red","purple"),lwd=3,cex=0.75,bg="white")

###############222222222222222222222222222222222222###############################
###############222222222222222222222222222222222222###############################
###############222222222222222222222222222222222222###############################
	par(mfrow=c(4,2));par(mar=c(3,5,2,2))

	#1 TB, Untreated TB exit rate, recent period, by route and HIV status
		plot(1:(92*12)/12+1950,OutMat1[1:(92*12),"ExTbC"]*12/(OutMat1[1:(92*12),"NUnTx"]-OutMat1[1:(92*12),"NUnTxH"]), 
			col="red", ylim=c(0,1), main="Exit Rates from Untreated TB, by Cause and HIV Status",cex.main=1,
			ylab="P(Exit)",xlab="Year",axes=FALSE, type="l",lwd=2)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,2,by=0.2),las=1); box()
			abline(h=0:20*0.1,col="grey85")
			lines(1:(92*12)/12+1950,OutMat1[1:(92*12),"ExTbC"]*12/(OutMat1[1:(92*12),"NUnTx"]-OutMat1[1:(92*12),"NUnTxH"]), col="forestgreen",lwd=2)
			lines(1:(92*12)/12+1950,OutMat1[1:(92*12),"ExTbT"]*12/(OutMat1[1:(92*12),"NUnTx"]-OutMat1[1:(92*12),"NUnTxH"]), col="blue",lwd=2)
			lines(1:(92*12)/12+1950,OutMat1[1:(92*12),"ExTbD"]*12/(OutMat1[1:(92*12),"NUnTx"]-OutMat1[1:(92*12),"NUnTxH"]), col="red",lwd=2)
			lines(1:(92*12)/12+1950,OutMat1[1:(92*12),"ExTbCH"]*12/OutMat1[1:(92*12),"NUnTxH"], col="purple",lwd=2)
			lines(1:(92*12)/12+1950,OutMat1[1:(92*12),"ExTbTH"]*12/OutMat1[1:(92*12),"NUnTxH"], col="goldenrod",lwd=2)
			lines(1:(92*12)/12+1950,OutMat1[1:(92*12),"ExTbDH"]*12/OutMat1[1:(92*12),"NUnTxH"], col="grey40",lwd=2)
			legend("topleft",c("HIV Neg, Self-Cure","HIV Neg, Init on Treatment","HIV Neg, Death",
			"HIV Pos, Self-Cure","HIV Pos, Init on Treatment","HIV Pos, Death"),
			col= c("forestgreen","blue","red","purple","goldenrod","grey40"),lwd=3,cex=0.75,bg="white")

	#2 TB Notifications per 100,000 
		TxNotifD   <- OutMat1[(20*12):(92*12),"NotifD"]/OutMat1[(20*12):(92*12),"NAll"]*100000*12
		TxInitAll <- (OutMat1[(20*12):(92*12),"NotifD"]+OutMat1[(20*12):(92*12),"NotifND"])/OutMat1[(20*12):(92*12),"NAll"]*100000*12 

		plot((20*12):(92*12)/12+1950,TxNotifD,col="red", main="DOTS Notifications (per 100,000)",ylim=c(0,1000),ylab="Tx Initiations per 100K",
			xlab="Year",axes=FALSE, type="l",lwd=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,3000,by=250),las=1); box()
			abline(h=0:20*250,col="grey85")
			lines((20*12):(92*12)/12+1950,TxNotifD,col="red",lwd=2,lty=2)
			lines((20*12):(92*12)/12+1950,TxInitAll, col="black",lwd=2,lty=2)
			points(1950:2043+0.5,CheckDat[,13]/CheckDat[,2]/10, col="blue",lwd=2,pch=4,cex=1.4)
			legend("topleft",c("WHO Reported Data","Model Projection (DOTS)","Model Projection (All)"),
			col= c("blue","red","black"),lwd=3,bg="white",cex=0.75)

	#3 Rate of Testing in DOTS/NON-DOTS Programs (Parameter Input) 
		plot(1:(92*12)/12+1950,DTestt[1:(92*12)],col="darkgreen", main="Rate of Testing, for Active TB (Model Input)",
			ylim=c(0,3),ylab="Rate",xlab="Year",axes=FALSE, type="l",lwd=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,20,by=1),las=1); box()
			abline(h=0:20*1,col="grey85")
			lines(1:(92*12)/12+1950,DTestt[1:(92*12)],col="darkgreen",lwd=2,lty=2)
			lines(1:(92*12)/12+1950,NDTestt[1:(92*12)], col="purple",lwd=2,lty=2)
			legend("topleft",c("DOTS Testing Rate","Non-DOTS Testing Rate"),
			col= c("darkgreen","purple"),lwd=3,cex=0.75,bg="white")

	#4 Annual No. New Suspects, per 100,000
 		plot((20*12):(92*12)/12+1950,OutMat1[(20*12):(92*12),"SuspctD"]/OutMat1[(20*12):(92*12),"NAll"]*100000*12,col="blue", 
			main="Annual No. TB Suspects per 100,000",cex.main=1,
			ylim=c(0,9000),ylab="per 100,000",xlab="Year",axes=FALSE, type="l",lwd=2)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,300000,by=2000),las=1); box()
			abline(h=0:20*2000,col="grey85")
			lines((20*12):(92*12)/12+1950,OutMat1[(20*12):(92*12),"SuspctD"]/OutMat1[(20*12):(92*12),"NAll"]*100000*12,col="blue",lwd=2)
			lines((20*12):(92*12)/12+1950,(OutMat1[(20*12):(92*12),"SuspctND"]+OutMat1[(20*12):(92*12),"SuspctD"])/
				OutMat1[1:(92*12),"NAll"]*100000*12, col="purple",lwd=2)
			points(1950:2043+0.5,CheckDat[,15]/CheckDat[,2]/10, col="blue",lwd=2,pch=4,cex=1.4)
			legend("topleft",c("DOTS Program","All Programs"),cex=0.75,
			col= c("blue","purple"),lwd=3,cex=1,bg="white")
			legend("topright",c("WHO Reported Data","Model Result"),lwd=3,cex=0.75,pch=c(4,NA),pt.cex=1.5,
				lty=c(0,1),bg="white")

	#5 Ratio of notifications to suspects
		plot((40*12):(92*12)/12+1950,OutMat1[(40*12):(92*12),"NotifD"]/OutMat1[(40*12):(92*12),"SuspctD"],axes=FALSE, xlab="Year", type="l",lwd=2,cex.main=1,
			col="navy", main="Ratio of Tb Notifications to Suspects (DOTS)",ylab="Notifications/Suspects", ylim=c(0,0.5))
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,200,by=0.1),las=1); box()
			abline(h=0:30*0.1,col="grey85")
			lines((40*12):(92*12)/12+1950,OutMat1[(40*12):(92*12),"NotifD"]/OutMat1[(40*12):(92*12),"SuspctD"], col="navy",lwd=2)
			points(1950:2043+0.5,CheckDat[,13]/CheckDat[,15], col="navy",lwd=2,pch=4,cex=1.5)
			legend("topleft",c("WHO Reported Data","Model Result"),lwd=3,cex=0.75,pch=c(4,NA),pt.cex=1.5,
				lty=c(0,1),bg="white", col="navy")

	#6 Probabilities of failure, default, death, DOTS programs 
		plot((40*12+6):(92*12)/12+1950,OutMat1[(40*12+6):(92*12),"PfailDtx"]*100,axes=FALSE, xlab="Year", type="l",lwd=2,
			col="red", main="Treatment Outcomes in DOTS Programs (Cures not shown)",ylab="Percent",ylim=c(0,25),cex.main=1)
			axis(1, at=seq(1990,2050,10),las=1); axis(2, at=seq(0,100,by=5),las=1); box()
			abline(h=0:30*5,col="grey85")
			lines((40*12+6):(92*12)/12+1950,OutMat1[(40*12+6):(92*12),"PfailDtx"]*100, col="red",lwd=2)
			lines((40*12+6):(92*12)/12+1950,OutMat1[(40*12+6):(92*12),"PdfltDtx"]*100, col="purple",lwd=2)
			lines((40*12+6):(92*12)/12+1950,OutMat1[(40*12+6):(92*12),"PmortDtx"]*100, col="forestgreen",lwd=2)
			points(1990:2043+0.5,CheckDat[-(1:40),26]*100, col="red",lwd=2,pch=4,cex=1.5)
			points(1990:2043+0.5,CheckDat[-(1:40),27]*100, col="purple",lwd=2,pch=4,cex=1.5)
			points(1990:2043+0.5,CheckDat[-(1:40),28]*100, col="forestgreen",lwd=2,pch=4,cex=1.5)
			legend("topright",c("Tx Outcome: Failed","Tx Outcome: Dead", "Tx Outcome: Defaulted"),
			col= c("red","forestgreen","purple"),lwd=3,cex=0.75,pch=4,pt.cex=1.5,bg="white")
			legend("topleft",c("WHO Reported Rata","Model Result"),lwd=3,cex=0.75,pch=c(4,NA),pt.cex=1.5,
				lty=c(0,1),bg="white")

	#7 Annual Mortality, in 1,000s, by group
		xMortNN <- (OutMat1[1:(92*12),"NMort"]-OutMat1[1:(92*12),"NHivMort"]
				-OutMat1[1:(92*12),"NTbMort"]+OutMat1[1:(92*12),"NTbHMort"])*12/1000
		xMortTN <- (OutMat1[1:(92*12),"NTbMort"]-OutMat1[1:(92*12),"NTbHMort"])*12/1000
		xMortNH <- (OutMat1[1:(92*12),"NHivMort"]-OutMat1[1:(92*12),"NTbHMort"])*12/1000
		xMortTH <- (OutMat1[1:(92*12),"NTbHMort"])*12/1000
		
		plot(1:(92*12)/12+1950,xMortNN,axes=FALSE, xlab="Year", type="l",lwd=2,
			col="red", main="Mortality (1,000s)",ylab="Mortality (1000s)", ylim=c(0,800),cex.main=1)
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,1200,by=100),las=1); box()
			abline(h=0:20*100,col="grey85")
			polygon(c(1950,1:(92*12)/12+1950,2042),c(0,xMortNN+xMortTN+xMortNH+xMortTH,0), col="red")  
			polygon(c(1950,1:(92*12)/12+1950,2042),c(0,xMortNN+xMortTN+xMortNH,0), col="darkgreen")  
			polygon(c(1950,1:(92*12)/12+1950,2042),c(0,xMortNN+xMortTN,0), col="blue")  
			polygon(c(1950,1:(92*12)/12+1950,2042),c(0,xMortNN,0), col="purple")  
			legend("topleft",c("Uninfected","TB Only", "HIV Only", "TB-HIV"),
			col= c("purple","blue", "darkgreen", "red"),lwd=3,cex=0.75,bg="white")

	#8 ART Coverage
		plot((40*12+1):(92*12)/12+1950,OutMat1[(40*12+1):(92*12),"ArtCov"]*100,col="blue", main="ART Coverage",
			ylim=c(0,100),ylab="Percent",xlab="Year",axes=FALSE, type="l",lwd=2,cex.main=1)  
			axis(1, at=seq(1990,2050,10),las=1); axis(2, at=seq(0,100,by=10),las=1); box()
			abline(h=0:10*10,col="grey85")
			lines((20*12+1):(92*12)/12+1950,OutMat1[(20*12+1):(92*12),"ArtCov"]*100,col="blue",lwd=2)  
			lines((20*12+1):(92*12)/12+1950,OutMat1[(20*12+1):(92*12),"ArtNdCov"]*100, col="red",lwd=2)  
			lines((20*12+1):(92*12)/12+1950,OutMat1[(20*12+1):(92*12),"Art200Cov"]*100, col="darkgreen",lwd=2)  
			legend("bottomright",c("Of All HIV Infected","Of All CD4 <350","Of All CD4<200" ),
			col= c("blue","red","darkgreen"),lwd=3,cex=0.75,bg="white")


###############333333333333333333333333333333333333###############################
###############333333333333333333333333333333333333###############################
###############333333333333333333333333333333333333###############################
	par(mfrow=c(4,2));par(mar=c(3,5,2,2))

	#1 Mortality Rate, by group
		xMortNN <- (OutMat1[1:(92*12),"NMort"]-OutMat1[1:(92*12),"NHivMort"]
				-OutMat1[1:(92*12),"NTbMort"]+OutMat1[1:(92*12),"NTbHMort"])*12/
				(OutMat1[1:(92*12),"NAll"]-OutMat1[1:(92*12),"NHiv"]-
				OutMat1[1:(92*12),"NActDis"]+OutMat1[1:(92*12),"NTbH"])		
		xMortTN <- (OutMat1[1:(92*12),"NTbMort"]-OutMat1[1:(92*12),"NTbHMort"])*12/
				(OutMat1[1:(92*12),"NActDis"]-OutMat1[1:(92*12),"NTbH"])
		xMortNH <- (OutMat1[1:(92*12),"NHivMort"]-OutMat1[1:(92*12),"NTbHMort"])*12/
				(OutMat1[1:(92*12),"NHiv"]-OutMat1[1:(92*12),"NTbH"])
		xMortTH <- (OutMat1[1:(92*12),"NTbHMort"])*12/OutMat1[1:(92*12),"NTbH"]
		
		plot(1:(92*12)/12+1950,xMortNN,axes=FALSE, xlab="Year", type="l",lwd=2,cex.main=1,
			col="purple", main="Annual Mortality Rate, Overall and by Sub-Group",ylab="Annual Mortality Rate", ylim=c(0,0.65))
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,1,by=0.1),las=1); box()
			abline(h=0:30/10,col="grey85")
			lines(1:(92*12)/12+1950,xMortNN, col="purple",lwd=2)  
			lines(1:(92*12)/12+1950,xMortNH, col="darkgreen",lwd=2)  
			lines(1:(92*12)/12+1950,xMortTN, col="blue",lwd=2)  
			lines(1:(92*12)/12+1950,xMortTH, col="red",lwd=2)  
			lines(1:(92*12)/12+1950,OutMat1[1:(92*12),"NMort"]/OutMat1[1:(92*12),"NAll"]*12, col="black",lwd=2,lty=2)
			legend("topleft",c("Uninfected","TB Only", "HIV Only", "TB-HIV","Overall"),
			col= c("purple","blue", "darkgreen", "red","black"),lwd=3,cex=0.75,bg="white")

	#2 Resistance in New Patients, % of total
		Nnaive	<- OutMat1[(20*12+1):(92*12),"NStr1n"]+OutMat1[(20*12+1):(92*12),"NStr2n"]+OutMat1[(20*12+1):(92*12),"NStr3n"]+
					OutMat1[(20*12+1):(92*12),"NStr4n"]+OutMat1[(20*12+1):(92*12),"NStr5n"]
							
		AnyInhPctn  <- (OutMat1[(20*12+1):(92*12),"NStr2n"]+OutMat1[(20*12+1):(92*12),"NStr4n"]+OutMat1[(20*12+1):(92*12),"NStr5n"])
		AnyRifPctn 	<- (OutMat1[(20*12+1):(92*12),"NStr3n"]+OutMat1[(20*12+1):(92*12),"NStr4n"]+OutMat1[(20*12+1):(92*12),"NStr5n"])
		MDRPctn 	<- (OutMat1[(20*12+1):(92*12),"NStr4n"]+OutMat1[(20*12+1):(92*12),"NStr5n"])
		XDRPctn 	<- (OutMat1[(20*12+1):(92*12),"NStr5n"])

		plot((20*12+1):(92*12)/12+1950,AnyInhPctn/Nnaive*100,axes=FALSE, xlab="Year", type="l",lwd=2,cex.main=1,
			col="red", main="Prevalence of TB Drug Resistance, Treatment Naive",ylab="Percent of All Cases", ylim=c(0,15))
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,50,by=2),las=1); box()
			abline(h=0:20*2,col="grey85")
			lines((20*12+1):(92*12)/12+1950,AnyInhPctn/Nnaive*100, col="red",lwd=2)  
			lines((20*12+1):(92*12)/12+1950,XDRPctn/Nnaive*100, col="purple",lwd=2)  
			lines((20*12+1):(92*12)/12+1950,AnyRifPctn/Nnaive*100, col="darkgreen",lwd=2)  
			lines((20*12+1):(92*12)/12+1950,MDRPctn/Nnaive*100, col="blue",lwd=2)  
			points(1970:2043+0.5,CheckDat[-(1:20),16]*100, col="red",lwd=2,pch=4,cex=1.5)
			points(1970:2043+0.5,CheckDat[-(1:20),17]*100, col="darkgreen",lwd=2,pch=4,cex=1.5)
			points(1970:2043+0.5,CheckDat[-(1:20),18]*100, col="blue",lwd=2,pch=4,cex=1.5)
			legend("topleft",c("Any INH Res","Any RIF Res", "MDR-TB","XDR-TB"),
			col= c("red","darkgreen", "blue","purple"),lwd=3,cex=0.75,pch=4,pt.cex=1.5,bg="white")
			legend("topright",c("2002 Resistance Survey","Model Result"),lwd=3,cex=0.75,pch=c(4,NA),pt.cex=1.5,
				lty=c(0,1),bg="white")
			
	#3 Resistance in Tx Experienced Patients, % of total
		Nexpd	<- OutMat1[(20*12+1):(92*12),"NStr1e"]+OutMat1[(20*12+1):(92*12),"NStr2e"]+OutMat1[(20*12+1):(92*12),"NStr3e"]+
					OutMat1[(20*12+1):(92*12),"NStr4e"]+OutMat1[(20*12+1):(92*12),"NStr5e"]
							
		AnyInhPcte  <- (OutMat1[(20*12+1):(92*12),"NStr2e"]+OutMat1[(20*12+1):(92*12),"NStr4e"]+OutMat1[(20*12+1):(92*12),"NStr5e"])
		AnyRifPcte 	<- (OutMat1[(20*12+1):(92*12),"NStr3e"]+OutMat1[(20*12+1):(92*12),"NStr4e"]+OutMat1[(20*12+1):(92*12),"NStr5e"])
		MDRPcte 	<- (OutMat1[(20*12+1):(92*12),"NStr4e"]+OutMat1[(20*12+1):(92*12),"NStr5e"])
		XDRPcte 	<- (OutMat1[(20*12+1):(92*12),"NStr5e"])

		plot((20*12+1):(92*12)/12+1950,AnyInhPcte/Nexpd*100,axes=FALSE, xlab="Year", type="l",lwd=2,cex.main=1,
			col="red", main="Prevalence of TB Drug Resistance, Treatment Experienced",ylab="Percent of All Cases",ylim=c(0,40))
			axis(1, at=seq(1970,2050,10),las=1); axis(2, at=seq(0,60,by=5),las=1); box()
			abline(h=0:10*5,col="grey85")
			lines((20*12+1):(92*12)/12+1950,AnyInhPcte/Nexpd*100, col="red",lwd=2) 
			lines((20*12+1):(92*12)/12+1950,XDRPcte/Nexpd*100, col="purple",lwd=2) 
			lines((20*12+1):(92*12)/12+1950,AnyRifPcte/Nexpd*100, col='darkgreen',lwd=2)  
			lines((20*12+1):(92*12)/12+1950,MDRPcte/Nexpd*100, col="blue",lwd=2)  
			points(1970:2043+0.5,CheckDat[-(1:20),19]*100, col="red",lwd=2,pch=4,cex=1.5)
			points(1970:2043+0.5,CheckDat[-(1:20),20]*100, col="darkgreen",lwd=2,pch=4,cex=1.5)
			points(1970:2043+0.5,CheckDat[-(1:20),21]*100, col="blue",lwd=2,pch=4,cex=1.5) 
			legend("topleft",c("Any INH Res","Any RIF Res", "MDR-TB","XDR-TB"),
			col= c("red","darkgreen", "blue","purple"),lwd=3,cex=0.75,pch=4,pt.cex=1.5,bg="white")
			legend("topright",c("2002 Resistance Survey","Model Result"),lwd=3,cex=0.75,pch=c(4,NA),pt.cex=1.5,
				lty=c(0,1),bg="white")

	#4 Resistance in All TB Cases, % of total
			plot((20*12+1):(92*12)/12+1950,(AnyInhPctn+AnyInhPcte)/(Nnaive+Nexpd)*100,axes=FALSE, xlab="Year", type="l",lwd=2,
			col="red", main="Prevalence of TB Drug Resistance, All TB Cases",ylab="%", ylim=c(0,15),cex.main=1)
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,60,by=2),las=1); box()
			abline(h=0:10*2,col="grey85")
			lines((20*12+1):(92*12)/12+1950,(AnyInhPctn+AnyInhPcte)/(Nnaive+Nexpd)*100, col="red",lwd=2) 
			lines((20*12+1):(92*12)/12+1950,(XDRPctn+XDRPcte)/(Nnaive+Nexpd)*100, col="purple",lwd=2) 
			lines((20*12+1):(92*12)/12+1950,(AnyRifPctn+AnyRifPcte)/(Nnaive+Nexpd)*100, col="darkgreen",lwd=2)  
			lines((20*12+1):(92*12)/12+1950,(MDRPctn+MDRPcte)/(Nnaive+Nexpd)*100, col="blue",lwd=2)  
			points(1970:2043+0.5,(CheckDat[-(1:20),16]*Nnaive[6+0:73*12]+CheckDat[-(1:20),19]*Nexpd[6+0:73*12])/(Nnaive[6+0:73*12]+Nexpd[6+0:73*12])*100, col="red",lwd=2,pch=4,cex=1.5)
			points(1970:2043+0.5,(CheckDat[-(1:20),17]*Nnaive[6+0:73*12]+CheckDat[-(1:20),20]*Nexpd[6+0:73*12])/(Nnaive[6+0:73*12]+Nexpd[6+0:73*12])*100, col="darkgreen",lwd=2,pch=4,cex=1.5)
			points(1970:2043+0.5,(CheckDat[-(1:20),18]*Nnaive[6+0:73*12]+CheckDat[-(1:20),21]*Nexpd[6+0:73*12])/(Nnaive[6+0:73*12]+Nexpd[6+0:73*12])*100, col="blue",lwd=2,pch=4,cex=1.5) 
			legend("topleft",c("Any INH Res","Any RIF Res", "MDR-TB","XDR-TB"),
			col= c("red","darkgreen", "blue","purple"),lwd=3,cex=0.75,pch=4,pt.cex=1.5,bg="white")
			legend("topright",c("2002 Resistance Survey","Model Result"),lwd=3,cex=0.75,pch=c(4,NA),pt.cex=1.5,
				lty=c(0,1),bg="white")

	#5 Duration of Infectiousness, by Group
			plot((20*12+1):(92*12)/12+1950,OutMat1[(20*12+1):(92*12),"DurInfAll"],axes=FALSE, xlab="Year", type="l",lwd=2,
			col="black", main="Duration Of Infectiousness, By Smear Status (years)",ylab="Years",ylim=c(0,3),cex.main=1)
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,100,by=0.5),las=1); box()
			abline(h=0:5*0.5,col="grey85")
			lines((20*12+1):(92*12)/12+1950,OutMat1[(20*12+1):(92*12),"DurInfAll"] , col="black",lwd=2)  
			lines((20*12+1):(92*12)/12+1950,OutMat1[(20*12+1):(92*12),"DurInfSn"] , col="blue",lwd=2)  
			lines((20*12+1):(92*12)/12+1950,OutMat1[(20*12+1):(92*12),"DurInfSp"], col="red",lwd=2)  
			legend("topright",c("General","Smear-negative","Smear-positive"),
			col= c("black","blue","red"),lwd=3,cex=0.75,bg="white")

	# Case Detection Rate, All Active Cases
		CdrDots 	<- (OutMat1[(20*12+1):(92*12),"NCdIpD"]+OutMat1[(20*12+1):(92*12),"NCdInD"])/OutMat1[(20*12+1):(92*12),"NCase"]*100
		CdrNonDots 	<- (OutMat1[(20*12+1):(92*12),"NCdIpND"]+OutMat1[(20*12+1):(92*12),"NCdInND"])/OutMat1[(20*12+1):(92*12),"NCase"]*100
		CdrAll 	<- CdrDots + CdrNonDots

	#6 Case Detection Rate, Smear Positive
		CdrDotsSp 		<- OutMat1[(20*12+1):(92*12),"NCdIpD"]/OutMat1[(20*12+1):(92*12),"NCaseIp"]*100
		CdrNonDotsSp 	<- OutMat1[(20*12+1):(92*12),"NCdIpND"]/OutMat1[(20*12+1):(92*12),"NCaseIp"]*100
		CdrAllSp 		<- CdrDotsSp + CdrNonDotsSp

		plot((20*12+1):(92*12)/12+1950,CdrNonDotsSp,	col="purple", main="Case Detection Rate, Smear Positives",ylim=c(0,160),ylab="Rate",
			xlab="Year",axes=FALSE, type="l",lwd=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,400,by=25),las=1); box()
			abline(h=0:10*20,col="grey85")
			lines((20*12+1):(92*12)/12+1950,CdrNonDotsSp, col="purple",lwd=2)  
			lines((20*12+1):(92*12)/12+1950,CdrDotsSp , col="darkgreen",lwd=2)  
			lines((20*12+1):(92*12)/12+1950,CdrAllSp, col="black",lwd=2)  
			points(1970:2043+0.5,CheckDat[-(1:20),25], col="darkgreen",lwd=2,pch=4,cex=1.4)

			legend("topleft",c("DOTS","Non-DOTS", "All"),
			col= c("darkgreen","purple","black"),lwd=3,cex=0.75,bg="white")
			legend("topright",c("WHO Estimate","Model Result"),lwd=3,cex=0.75,pch=c(4,NA),pt.cex=1.5,
				lty=c(0,1),bg="white")

	#7 Effective Contact Rate
		plot(1:(92*12)/12+1950,OutMat1[1:(92*12),"EffContRate"],col="navy", main="Effective Contact Rate (untreated cases)",
			ylim=c(0,10),ylab="Rate",xlab="Year",axes=FALSE, type="l",lwd=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,40,by=2),las=1); box()
			abline(h=0:10*2,col="grey85")
			lines(1:(92*12)/12+1950,OutMat1[1:(92*12),"EffContRate"],col="navy",lwd=2)  

	#8 Positive & Negative Predictive Value, DOTS TB Diagnosis
		plot((40*12+6):(92*12)/12+1950,OutMat1[(40*12+6):(92*12),"PPVTb"]*100,col="blue", main="PPV/NPV of TB Diagnosis (DOTS)",
			ylim=c(30,100),ylab="Percent",xlab="Year",axes=FALSE, type="l",lwd=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,100,by=10),las=1); box()
			abline(h=0:10*10,col="grey85")
			lines((40*12+6):(92*12)/12+1950,OutMat1[(40*12+6):(92*12),"PPVTb"]*100,col="blue",lwd=2)  
			lines((40*12+6):(92*12)/12+1950,OutMat1[(40*12+6):(92*12),"NPVTb"]*100, col="red",lwd=2)  
			legend("bottomleft",c("Positive Predictive Value","Negative Predictive Value"),
			col= c("blue","red"),lwd=3,cex=0.75,bg="white")

###############444444444444444444444444444444444444###############################
###############444444444444444444444444444444444444###############################
###############444444444444444444444444444444444444###############################
	par(mfrow=c(4,2));par(mar=c(3,5,2,2))

	#1 Positive& Negative Predictive Value, DOTS TB Diagnosis
		plot((62*12+1):(92*12)/12+1950,OutMat2[(62*12+1):(92*12),"PPVTb"]*100,col="blue", main="Positive/Negative Predictive Value of TB Diagnosis",
			ylim=c(60,100),ylab="Percent",xlab="Year",axes=FALSE, type="l",lwd=2,lty=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,100,by=10),las=1); box()
			abline(h=0:10*10,col="grey85")
			lines((62*12+1):(92*12)/12+1950,OutMat2[(62*12+1):(92*12),"PPVTb"]*100,col="blue",lwd=2,lty=2)  
			lines((62*12+1):(92*12)/12+1950,OutMat1[(62*12+1):(92*12),"PPVTb"]*100,col="blue",lwd=2)  
			lines((62*12+1):(92*12)/12+1950,OutMat2[(62*12+1):(92*12),"NPVTb"]*100, col="red",lwd=2,lty=2)  
			lines((62*12+1):(92*12)/12+1950,OutMat1[(62*12+1):(92*12),"NPVTb"]*100, col="red",lwd=2)  
			legend("bottomright",c("Positive Predictive Value","Negative Predictive Value"),
			col= c("blue","red"),lwd=3,cex=0.75,bg="white")
			legend("bottomleft",c("BaseCase","Xpert Algorithm 1"),lwd=2,cex=0.75,
				lty=c(1,2),bg="white")

	#2 Positive & Negative Predictive Value, Xpert RIF Result
		plot((62*12+1):(92*12)/12+1950,OutMat2[(62*12+1):(92*12),"PPVRif"]*100,col="blue", main="Positive/Negative Predictive Value of Xpert RIF Result",
			ylim=c(50,100),ylab="Percent",xlab="Year",axes=FALSE, type="l",lwd=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,100,by=10),las=1); box()
			abline(h=0:100*5,col="grey85")
			lines((62*12+1):(92*12)/12+1950,OutMat2[(62*12+1):(92*12),"PPVRif"]*100,col="blue",lwd=2)  
			lines((62*12+1):(92*12)/12+1950,OutMat2[(62*12+1):(92*12),"NPVRif"]*100, col="red",lwd=2)  
			legend("bottomright",c("Positive Predictive Value","Negative Predictive Value"),
			col= c("blue","red"),lwd=3,,cex=0.75,bg="white")

	#3 Annual Notifications
		plot((62*12+1):(92*12)/12+1950,OutMat2[(62*12+1):(92*12),"NotifD"]*12/1000,col="blue", main="No. Annual Notifications",
			ylim=c(0,250),ylab="Notifications (1,000s)",xlab="Year",axes=FALSE, type="l",lwd=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,1000,by=100),las=1); box()
			abline(h=0:100*100,col="grey85")
			lines((62*12+1):(92*12)/12+1950,OutMat2[(62*12+1):(92*12),"NotifD"]*12/1000,col="blue",lwd=2)  
			lines((62*12+1):(92*12)/12+1950,OutMat1[(62*12+1):(92*12),"NotifD"]*12/1000, col="red",lwd=2)  
			legend("topright",c("Base Case","Xpert Algorithm 1"),
			col= c("red","blue"),lwd=3,cex=0.75,bg="white",lty=c(1,1))

	#4 Annual True Positive Diagnoses
		plot((62*12+1):(92*12)/12+1950,OutMat2[(62*12+1):(92*12),"NotifTBD"]*12/1000,col="blue", main="No. Annual True-Positive Notifications",
			ylim=c(0,250),ylab="Notifications (1,000s)",xlab="Year",axes=FALSE, type="l",lwd=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,1000,by=100),las=1); box()
			abline(h=0:100*100,col="grey85")
			lines((62*12+1):(92*12)/12+1950,OutMat2[(62*12+1):(92*12),"NotifTBD"]*12/1000,col="blue",lwd=2)  
			lines((62*12+1):(92*12)/12+1950,OutMat1[(62*12+1):(92*12),"NotifTBD"]*12/1000, col="red",lwd=2)  
			legend("topright",c("Base Case","Xpert Algorithm 1"),
			col= c("red","blue"),lwd=3,cex=0.75,bg="white",lty=c(1,1))

	#5 Annual True Positive Tx Initiations
		plot((62*12+1):(92*12)/12+1950,rowSums(OutMat2[(62*12+1):(92*12),c("NCdInD","NCdIpD")])*12/1000,col="blue", main="No. Annual True-Positive Tx Initiations",
			ylim=c(0,250),ylab="Tx Initiations (1,000s)",xlab="Year",axes=FALSE, type="l",lwd=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,1000,by=100),las=1); box()
			abline(h=0:100*100,col="grey85")
			lines((62*12+1):(92*12)/12+1950,rowSums(OutMat2[(62*12+1):(92*12),c("NCdInD","NCdIpD")])*12/1000,col="blue",lwd=2)  
			lines((62*12+1):(92*12)/12+1950,rowSums(OutMat1[(62*12+1):(92*12),c("NCdInD","NCdIpD")])*12/1000, col="red",lwd=2)  
			legend("topright",c("Base Case","Xpert Algorithm 1"),
			col= c("red","blue"),lwd=3,cex=0.75,bg="white",lty=c(1,1))

	#6 TB, Active Disease per 100,000, recent period 
		plot((62*12+1):(92*12)/12+1950,(OutMat2[(62*12+1):(92*12),"NActDis"]-OutMat2[(62*12+1):(92*12),"NTxD"]-
			OutMat2[(62*12+1):(92*12),"NTxND"])/OutMat2[(62*12+1):(92*12),"NAll"]*100000, col="red", 
			main="TB Prevalence (Untreated Individuals, per 100,000)",ylim=c(0,1000),ylab="Active TB per 100K",
			xlab="Year",axes=FALSE, type="l",lwd=2,lty=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,4000,by=500),las=1); box()
			abline(h=0:10*250,col="grey85")
			lines((62*12+1):(92*12)/12+1950,(OutMat2[(62*12+1):(92*12),"NActDis"]-OutMat2[(62*12+1):(92*12),"NTxD"]-
				OutMat2[(62*12+1):(92*12),"NTxND"])/OutMat2[(62*12+1):(92*12),"NAll"]*100000, col="red",lwd=2)
			lines((62*12+1):(92*12)/12+1950,(OutMat1[(62*12+1):(92*12),"NActDis"]-OutMat1[(62*12+1):(92*12),"NTxD"]-
				OutMat1[(62*12+1):(92*12),"NTxND"])/OutMat1[(62*12+1):(92*12),"NAll"]*100000, col="blue",lwd=2)
			legend("topright",c("Base Case","Xpert Algorithm 1"),
			col= c("blue","red"),lwd=3,cex=0.75,bg="white",lty=c(1,1))

	#7 TB, Incidence per 100,000, recent period 
		plot((62*12+1):(92*12)/12+1950,OutMat2[(62*12+1):(92*12),"NCase"]/OutMat2[(62*12+1):(92*12),"NAll"]*100000*12, col="red", 
			main="TB Incidence (per 100,000)",ylim=c(0,750),ylab="Incidence per 100K",
			xlab="Year",axes=FALSE, type="l",lwd=2,lty=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,4000,by=250),las=1); box()
			abline(h=0:10*250,col="grey85")
			lines((62*12+1):(92*12)/12+1950,OutMat2[(62*12+1):(92*12),"NCase"]/OutMat2[(62*12+1):(92*12),"NAll"]*100000*12,
				 col="red",lwd=2)
			lines((62*12+1):(92*12)/12+1950,OutMat1[(62*12+1):(92*12),"NCase"]/OutMat1[(62*12+1):(92*12),"NAll"]*100000*12,
				 col="blue",lwd=2)
			legend("topright",c("Base Case","Xpert Algorithm 1"),
			col= c("blue","red"),lwd=3,cex=0.75,bg="white",lty=c(1,1))

	#8 Annual risk of infection (ARI), recent period 
		plot((62*12+1):(92*12)/12+1950,(1-(1-OutMat2[(62*12+1):(92*12),"NInf"]/(OutMat2[(62*12+1):(92*12),"NAll"]-
			OutMat2[(62*12+1):(92*12),"NAnyTb"]))^12)*100, col="red", main="Annual Risk of Infection"
			,ylim=c(0,4),ylab="ARI (%)",xlab="Year",axes=FALSE, type="l",lwd=2,lty=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,40,by=2),las=1); box()
			abline(h=0:10*2,col="grey85")
			lines((62*12+1):(92*12)/12+1950,(1-(1-OutMat2[(62*12+1):(92*12),"NInf"]/(OutMat2[(62*12+1):(92*12),"NAll"]-
			OutMat2[(62*12+1):(92*12),"NAnyTb"]))^12)*100, col="red",lwd=2)
			lines((62*12+1):(92*12)/12+1950,(1-(1-OutMat1[(62*12+1):(92*12),"NInf"]/(OutMat1[(62*12+1):(92*12),"NAll"]-
			OutMat1[(62*12+1):(92*12),"NAnyTb"]))^12)*100, col="blue",lwd=2)
			legend("topright",c("Base Case","Xpert Algorithm 1"),
			col= c("blue","red"),lwd=3,cex=0.75,bg="white",lty=c(1,1))


###############555555555555555555555555555555555555###############################
###############555555555555555555555555555555555555###############################
###############555555555555555555555555555555555555###############################
	par(mfrow=c(4,2));par(mar=c(3,5,2,2))

	# Resistance in Tx Naive Patients, % of total
		NnaiveAlg1	<- OutMat1[(62*12+1):(92*12),"NStr1n"]+OutMat1[(62*12+1):(92*12),"NStr2n"]+OutMat1[(62*12+1):(92*12),"NStr3n"]+
					OutMat1[(62*12+1):(92*12),"NStr4n"]+OutMat1[(62*12+1):(92*12),"NStr5n"]			
		AnyInhPctnAlg1  <- (OutMat1[(62*12+1):(92*12),"NStr2n"]+OutMat1[(62*12+1):(92*12),"NStr4n"]+
					OutMat1[(62*12+1):(92*12),"NStr5n"])
		MDRPctnAlg1 	<- (OutMat1[(62*12+1):(92*12),"NStr4n"]+OutMat1[(62*12+1):(92*12),"NStr5n"])

		NnaiveAlg2	<- OutMat2[(62*12+1):(92*12),"NStr1n"]+OutMat2[(62*12+1):(92*12),"NStr2n"]+OutMat2[(62*12+1):(92*12),"NStr3n"]+
					OutMat2[(62*12+1):(92*12),"NStr4n"]+OutMat2[(62*12+1):(92*12),"NStr5n"]			
		AnyInhPctnAlg2  <- (OutMat2[(62*12+1):(92*12),"NStr2n"]+OutMat2[(62*12+1):(92*12),"NStr4n"]+
					OutMat2[(62*12+1):(92*12),"NStr5n"])
		MDRPctnAlg2 	<- (OutMat2[(62*12+1):(92*12),"NStr4n"]+OutMat2[(62*12+1):(92*12),"NStr5n"])
		NexpcdAlg1	<- OutMat1[(62*12+1):(92*12),"NStr1e"]+OutMat1[(62*12+1):(92*12),"NStr2e"]+OutMat1[(62*12+1):(92*12),"NStr3e"]+
					OutMat1[(62*12+1):(92*12),"NStr4e"]+OutMat1[(62*12+1):(92*12),"NStr5e"]			
		AnyInhPcteAlg1  <- (OutMat1[(62*12+1):(92*12),"NStr2e"]+OutMat1[(62*12+1):(92*12),"NStr4e"]+
					OutMat1[(62*12+1):(92*12),"NStr5e"])
		MDRPcteAlg1 	<- (OutMat1[(62*12+1):(92*12),"NStr4e"]+OutMat1[(62*12+1):(92*12),"NStr5e"])

		NexpcdAlg2	<- OutMat2[(62*12+1):(92*12),"NStr1e"]+OutMat2[(62*12+1):(92*12),"NStr2e"]+OutMat2[(62*12+1):(92*12),"NStr3e"]+
					OutMat2[(62*12+1):(92*12),"NStr4e"]+OutMat2[(62*12+1):(92*12),"NStr5e"]			
		AnyInhPcteAlg2  <- (OutMat2[(62*12+1):(92*12),"NStr2e"]+OutMat2[(62*12+1):(92*12),"NStr4e"]+
					OutMat2[(62*12+1):(92*12),"NStr5e"])
		MDRPcteAlg2 	<- (OutMat2[(62*12+1):(92*12),"NStr4e"]+OutMat2[(62*12+1):(92*12),"NStr5e"])

	#1 Resistance in All Tb Cases, % of total

		plot((62*12+1):(92*12)/12+1950,(AnyInhPctnAlg1+AnyInhPcteAlg1)/(NnaiveAlg1+NexpcdAlg1)*100,axes=FALSE, xlab="Year", type="l",lwd=2,
			col="forestgreen", main="Prevalence of TB Drug Resistance, All TB Cases (%)",ylab="Percent of All Cases", ylim=c(0,15),cex.main=1)
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,50,by=2),las=1); box()
			abline(h=0:10*2,col="grey85")
			lines((62*12+1):(92*12)/12+1950,(AnyInhPctnAlg1+AnyInhPcteAlg1)/(NnaiveAlg1+NexpcdAlg1)*100, col="forestgreen",lwd=2)  
			lines((62*12+1):(92*12)/12+1950,(AnyInhPctnAlg2+AnyInhPcteAlg2)/(NnaiveAlg2+NexpcdAlg2)*100, col="forestgreen",lwd=2,lty=2)  
			lines((62*12+1):(92*12)/12+1950,(MDRPctnAlg1+MDRPcteAlg1)/(NnaiveAlg1+NexpcdAlg1)*100, col="blue",lwd=2)  
			lines((62*12+1):(92*12)/12+1950,(MDRPctnAlg2+MDRPcteAlg2)/(NnaiveAlg2+NexpcdAlg2)*100, col="blue",lwd=2,lty=2)  
			legend("topleft",c("Any INH Resistance","MDR-TB"),
			col= c("forestgreen", "blue"),lwd=3,cex=0.75,bg="white")
			legend("topright",c("BaseCase","Xpert Algorithm 1"),lwd=2,cex=0.75,
				lty=c(1,2),bg="white")

	#2 Resistance in All Patients, Absolute No.

		plot((62*12+1):(92*12)/12+1950,AnyInhPctnAlg1+AnyInhPcteAlg1,axes=FALSE, xlab="Year", type="l",lwd=2,
			col="forestgreen", main="Drug Resistant TB (Absolute No. Cases)",ylab="No. Cases", 
			ylim=c(0,30000),cex.main=1)
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,500000,by=10000),las=1); box()
			abline(h=0:10*10000,col="grey85")
			lines((62*12+1):(92*12)/12+1950,AnyInhPctnAlg1+AnyInhPcteAlg1, col="forestgreen",lwd=2)  
			lines((62*12+1):(92*12)/12+1950,AnyInhPctnAlg2+AnyInhPcteAlg2, col="forestgreen",lwd=2,lty=2)  
			lines((62*12+1):(92*12)/12+1950,MDRPctnAlg1+MDRPcteAlg1, col="blue",lwd=2)  
			lines((62*12+1):(92*12)/12+1950,MDRPctnAlg2+MDRPcteAlg2, col="blue",lwd=2,lty=2)  
			legend("topleft",c("Any INH Resistance","MDR-TB"),
			col= c("forestgreen", "blue"),lwd=3,cex=0.75,bg="white")
			legend("topright",c("BaseCase","Xpert Algorithm 1"),lwd=2,cex=0.75,
				lty=c(1,2),bg="white")
	#3 TB Mortality per 100,000, by HIV + / -
		TbMort1 <- (OutMat1[(62*12+1):(92*12),"NTbMort"]-OutMat1[(62*12+1):(92*12),"NTbHMort"])*12
		TbMort2 <- (OutMat2[(62*12+1):(92*12),"NTbMort"]-OutMat2[(62*12+1):(92*12),"NTbHMort"])*12
		TbHMort1 <- OutMat1[(62*12+1):(92*12),"NTbHMort"]*12
		TbHMort2 <- OutMat2[(62*12+1):(92*12),"NTbHMort"]*12

				
	plot((62*12+1):(92*12)/12+1950,TbMort1 , col="blue", main="Annual TB Mortality, by HIV Status",
			ylim=c(0,80000),ylab="Annual Deaths",xlab="Year",axes=FALSE, type="l",lwd=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,200000,by=80000),las=1); box()
			abline(h=0:10*20000,col="grey85")
			lines((62*12+1):(92*12)/12+1950,TbMort1 , col="blue",lwd=2)
			lines((62*12+1):(92*12)/12+1950,TbMort2 , col="blue",lwd=2,lty=2)
			lines((62*12+1):(92*12)/12+1950,TbHMort1 , col="red",lwd=2)
			lines((62*12+1):(92*12)/12+1950,TbHMort2 , col="red",lwd=2,lty=2)
			legend("topright",c("HIV negative deaths","HIV positive deaths"),
			col= c("blue","red"),lwd=3,cex=0.75,bg="white")
			legend("topleft",c("BaseCase","Xpert Algorithm 1"),lwd=2,cex=0.75,
				lty=c(1,2),bg="white")

	#4 Annual Costs by Category, DOTS
		CostTestD1	<- OutMat1[(62*12+1):(92*12),"CostTestD"]*12/1000000+OutMat1[(62*12+1):(92*12),"GetXpt"]*12/1000000*30
		CostTxD1	<- (OutMat1[(62*12+1):(92*12),"CostTxD"]+OutMat1[(62*12+1):(92*12),"CostRegD"]
				+OutMat1[(62*12+1):(92*12),"CostFalsRegD"]+OutMat1[(62*12+1):(92*12),"CostFalsTxD"])*12/1000000
		CostTestD2	<- OutMat2[(62*12+1):(92*12),"CostTestD"]*12/1000000+OutMat2[(62*12+1):(92*12),"GetXpt"]*12/1000000*30
		CostTxD2	<- (OutMat2[(62*12+1):(92*12),"CostTxD"]+OutMat2[(62*12+1):(92*12),"CostRegD"]
				+OutMat2[(62*12+1):(92*12),"CostFalsRegD"]+OutMat2[(62*12+1):(92*12),"CostFalsTxD"])*12/1000000

		plot((62*12+1):(92*12)/12+1950,CostTestD1,col="forestgreen", main="Total Annual TB Costs by Category (USD Million))",
			ylim=c(0,125),ylab="Cost (USD, Millions)",xlab="Year",axes=FALSE, type="l",lwd=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,500,by=25),las=1); box()
			abline(h=0:10*25,col="grey85")
			lines((62*12+1):(92*12)/12+1950,CostTestD1, col="forestgreen",lwd=2)  
			lines((62*12+1):(92*12)/12+1950,CostTestD2, col="forestgreen",lwd=2,lty=2)  
			lines((62*12+1):(92*12)/12+1950,CostTxD1, col="purple",lwd=2)  
			lines((62*12+1):(92*12)/12+1950,CostTxD2, col="purple",lwd=2,lty=2)  
			legend("topleft",c("Diagnosis Costs","Treatment Costs"),
			col= c("forestgreen", "purple"),lwd=3,cex=0.75,bg="white")
			legend("topright",c("BaseCase","Xpert Algorithm 1"),lwd=2,cex=0.75,
				lty=c(1,2),bg="white")

	#5 Annual Incremental Costs by Category, DOTS
		IncCostTestD	<- CostTestD2 - CostTestD1 
		IncCostTxD		<- CostTxD2 - CostTxD1 
		IncCostHAART	<- (OutMat2[(62*12+1):(92*12),"CostART"]-OutMat1[(62*12+1):(92*12),"CostART"])*12/1000000

		plot((62*12+1):(92*12)/12+1950,IncCostTestD,col="forestgreen", main="Total Annual Incremental Costs by Category (USD Million)",
			ylim=c(0,100),ylab="Cost (USD, Millions)",xlab="Year",axes=FALSE, type="l",lwd=2,cex.main=1)  
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,500,by=25),las=1); box()
			abline(h=0:10*25,col="grey85")
			lines((62*12+1):(92*12)/12+1950,IncCostTxD, col="navy",lwd=2)  
			lines((62*12+1):(92*12)/12+1950,IncCostHAART, col="brown",lwd=2)  
			legend("topleft",c("DOTS Diagnosis Costs","DOTS Treatment Costs","HAART Costs"),
			col= c("forestgreen", "navy","brown"),lwd=3,cex=0.75,bg="white")

	#6 Duration of Infectiousness, by Group
			plot((62*12+1):(92*12)/12+1950,OutMat1[(62*12+1):(92*12),"DurInfAll"],axes=FALSE, xlab="Year", type="l",lwd=2,
			col="black", main="Duration Of Infectiousness, By Smear Status (years)",ylab="Years", ylim=c(0,3),cex.main=1)
			axis(1, at=seq(1950,2050,10),las=1); axis(2, at=seq(0,100,by=1),las=1); box()
			abline(h=0:10*0.5,col="grey85")
			lines((62*12+1):(92*12)/12+1950,OutMat1[(62*12+1):(92*12),"DurInfAll"] , col="black",lwd=2)  
			lines((62*12+1):(92*12)/12+1950,OutMat1[(62*12+1):(92*12),"DurInfSn"] , col="blue",lwd=2)  
			lines((62*12+1):(92*12)/12+1950,OutMat1[(62*12+1):(92*12),"DurInfSp"], col="red",lwd=2)  
			lines((62*12+1):(92*12)/12+1950,OutMat2[(62*12+1):(92*12),"DurInfAll"] , col="black",lwd=2,lty=2)  
			lines((62*12+1):(92*12)/12+1950,OutMat2[(62*12+1):(92*12),"DurInfSn"] , col="blue",lwd=2,lty=2)  
			lines((62*12+1):(92*12)/12+1950,OutMat2[(62*12+1):(92*12),"DurInfSp"], col="red",lwd=2,lty=2)  
			legend("topright",c("General","Smear-negative","Smear-positive"),
			col= c("black","blue","red"),lwd=3,cex=0.75,bg="white")
			legend("topleft",c("BaseCase","Xpert Algorithm 1"),lwd=2,cex=0.75,
				lty=c(1,2),bg="white")

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################


	## Total Annual Costs by Category, DOTS, since 2011
		CostsDx1	<- OutMat1[(62*12+1):(92*12),c("CostTestD")]*12+OutMat1[(62*12+1):(92*12),"GetXpt"]*12*20
		CostsDx2	<- OutMat2[(62*12+1):(92*12),c("CostTestD")]*12+OutMat2[(62*12+1):(92*12),"GetXpt"]*12*20
		CostsDOTS1	<- apply(OutMat1[(62*12+1):(92*12),c("CostTestD","CostTxD","CostFalsRegD","CostFalsTxD","CostRegD")],1,sum)*12+
					OutMat1[(62*12+1):(92*12),"GetXpt"]*12*30
		CostsDOTS2	<- apply(OutMat2[(62*12+1):(92*12),c("CostTestD","CostTxD","CostFalsRegD","CostFalsTxD","CostRegD")],1,sum)*12+
					OutMat2[(62*12+1):(92*12),"GetXpt"]*12*30
		CostsTb1	<- CostsDOTS1 + apply(OutMat1[(62*12+1):(92*12),c("CostTestND","CostTxND","CostFalsRegND","CostFalsTxND","CostRegND")],1,sum)*12
		CostsTb2	<- CostsDOTS2 + apply(OutMat2[(62*12+1):(92*12),c("CostTestND","CostTxND","CostFalsRegND","CostFalsTxND","CostRegND")],1,sum)*12
		CostsAll1	<- CostsTb1 + OutMat1[(62*12+1):(92*12),"CostART"]*12
		CostsAll2	<- CostsTb2 + OutMat2[(62*12+1):(92*12),"CostART"]*12

	## Health outcomes
		Notif1	<- apply(OutMat1[(62*12+1):(92*12),c("NCdIpD","NCdInD")],1,sum)*12
		Notif2	<- apply(OutMat2[(62*12+1):(92*12),c("NCdIpD","NCdInD")],1,sum)*12
		LY1		<- OutMat1[(62*12+1):(92*12),c("NAll")]
		LY2		<- OutMat2[(62*12+1):(92*12),c("NAll")]
		Daly1		<- OutMat1[(62*12+1):(92*12),c("Ndaly")]
		Daly2		<- OutMat2[(62*12+1):(92*12),c("Ndaly")]

	ResultsMat1  <- matrix(NA,ncol=2,nrow=8)
	colnames(ResultsMat1) <- c("10 Years","30 Years")
	rownames(ResultsMat1) <- c("$/DALY,Prog Persp","$/DALY,HS Persp","$/LY,Prog Persp","$/LY,HS Persp",
						"$/DALY,Prog Persp,Disc.","$/DALY,HS Persp,Disc.","$/LY,Prog Persp,Disc.","$/LY,HS Persp,Disc.")

# ICER for DALYs Averted at 5,10,20, 30 years (DOTS Program Perspective)

	ResultsMat1[1,1] <- sum((CostsDOTS2-CostsDOTS1)[1:(12*10)])/sum((Daly2-Daly1)[1:(12*10)])
	ResultsMat1[1,2] <- sum((CostsDOTS2-CostsDOTS1)[1:(12*30)])/sum((Daly2-Daly1)[1:(12*30)])

# ICER for DALYs Averted at 5,10,20, 30 years (Health System Perspective)

	ResultsMat1[2,1] <- sum((CostsAll2-CostsAll1)[1:(12*10)])/sum((Daly2-Daly1)[1:(12*10)])
	ResultsMat1[2,2] <- sum((CostsAll2-CostsAll1)[1:(12*30)])/sum((Daly2-Daly1)[1:(12*30)])

# ICER for LY Saved at 5,10,20, 30 years (DOTS Program Perspective)

	ResultsMat1[3,1] <- sum((CostsDOTS2-CostsDOTS1)[1:(12*10)])/sum((LY2-LY1)[1:(12*10)])
	ResultsMat1[3,2] <- sum((CostsDOTS2-CostsDOTS1)[1:(12*30)])/sum((LY2-LY1)[1:(12*30)])

# ICER for LY Saved at 5,10,20, 30 years (Health System Perspective)

	ResultsMat1[4,1] <- sum((CostsAll2-CostsAll1)[1:(12*10)])/sum((LY2-LY1)[1:(12*10)])
	ResultsMat1[4,2] <- sum((CostsAll2-CostsAll1)[1:(12*30)])/sum((LY2-LY1)[1:(12*30)])

# ICER for DALYs Averted at 5,10,20, 30 years (DOTS Program Perspective)

	Vdiscount <- rep(1,30*12); for(i in 1:length(Vdiscount)) { Vdiscount[i] <- 1.03^(-(i-1)/12) }

	ResultsMat1[5,1] <- sum(((CostsDOTS2-CostsDOTS1)*Vdiscount)[1:(12*10)])/sum(((Daly2-Daly1)*Vdiscount)[1:(12*10)])
	ResultsMat1[5,2] <- sum(((CostsDOTS2-CostsDOTS1)*Vdiscount)[1:(12*30)])/sum(((Daly2-Daly1)*Vdiscount)[1:(12*30)])

# ICER for DALYs Averted at 5,10,20, 30 years (Health System Perspective)

	ResultsMat1[6,1] <- sum(((CostsAll2-CostsAll1)*Vdiscount)[1:(12*10)])/sum(((Daly2-Daly1)*Vdiscount)[1:(12*10)])
	ResultsMat1[6,2] <- sum(((CostsAll2-CostsAll1)*Vdiscount)[1:(12*30)])/sum(((Daly2-Daly1)*Vdiscount)[1:(12*30)])

# ICER for LY Saved at 5,10,20, 30 years (DOTS Program Perspective)

	ResultsMat1[7,1] <- sum(((CostsDOTS2-CostsDOTS1)*Vdiscount)[1:(12*10)])/sum(((LY2-LY1)*Vdiscount)[1:(12*10)])
	ResultsMat1[7,2] <- sum(((CostsDOTS2-CostsDOTS1)*Vdiscount)[1:(12*30)])/sum(((LY2-LY1)*Vdiscount)[1:(12*30)])

# ICER for LY Saved at 5,10,20, 30 years (Health System Perspective)

	ResultsMat1[8,1] <- sum(((CostsAll2-CostsAll1)*Vdiscount)[1:(12*10)])/sum(((LY2-LY1)*Vdiscount)[1:(12*10)])
	ResultsMat1[8,2] <- sum(((CostsAll2-CostsAll1)*Vdiscount)[1:(12*30)])/sum(((LY2-LY1)*Vdiscount)[1:(12*30)])

	ResultsMat1 <- round(ResultsMat1)
	write.table(ResultsMat1, file=paste("ICERs1",Setting,"NEW.csv",sep=""), sep = ",")







