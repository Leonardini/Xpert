##########   TB XPERT DIAGNOSTIC MODEL 2012  ##########   
##########            BURNIN FILE            ##########  

### Burn-in for 100 years
	rm(list = ls()) 
	# setwd("C:/Users/nick/Documents/Harvard/TB Diagnostics/ANALYSIS")
	MeanParam <- read.csv(file="MeanParam.csv",row.names=1)
	Vparam <- MeanParam ;  Setting <- "SouthAfrica"
	source("PARAMv5.r")
	source("TIMESTEPv5.r")
	DIAG <- 2
	OutBurnIn <- OutMat
	Vcurrent <- c(990,0,10,rep(0,502))
	for(t in 1:(100*12)) { 	Out <- timestep(Vcurrent,1,ArtNdCov11,DIAG); Vcurrent <- Out$Vnext; OutBurnIn[t,] <- Out$Vout }
	V1950 <- Vcurrent

	prev <- as.integer(OutBurnIn[1200,"NUnTx"]/OutBurnIn[1200,"NAll"]*100000); prev 
			# TB prevalence = 2000
	incid <- as.integer(OutBurnIn[1200,"NCase"]/OutBurnIn[1200,"NAll"]*100000*12); incid
		  	# TB Incidence = 965
	mort <- as.integer(OutBurnIn[1200,"NTbMort"]/OutBurnIn[1200,"NAll"]*100000*12); mort
		  	# TB Mortality = 568
	anytb <- as.integer(OutBurnIn[1200,"NAnyTb"]/OutBurnIn[1200,"NAll"]*100); anytb 
		  	# MTB Carriage = 85% 
	infrisk <- as.integer((1-(1-OutBurnIn[1200,"NInf"]/(OutBurnIn[1200,"NAll"]-OutBurnIn[1200,"NAnyTb"]))^12)*100); infrisk 
		  	# ARI = 16% 
	par(mfrow=c(3,2));par(mar=c(5,5,2,2))
	plot(1:1200/12,OutBurnIn[1:1200,"NUnTx"]/OutBurnIn[1:1200,"NAll"]*100000,
		main="TB Prevalence",type="l",col="darkgreen",lwd=2,ylab="Cases per 100,000",xlab="Years")
		legend("bottomright",paste("Equilibrium Prev =",prev,"per 100,000"),col="darkgreen",cex=0.9,lwd=2)

	plot(1:1200/12,OutBurnIn[1:1200,"NCase"]/OutBurnIn[1:1200,"NAll"]*100000*12,
		main="TB Incidence",type="l",col="blue",lwd=2,ylab="Cases per 100,000",xlab="Years")
		legend("bottomright",paste("Equilibrium Incid =",incid,"per 100,000"),col="blue",cex=0.9,lwd=2)

	plot(1:1200/12,OutBurnIn[1:1200,"NAnyTb"]/OutBurnIn[1:1200,"NAll"]*100,
		main="Percent Any TB (incl. latent)",type="l",col="purple",lwd=2,ylab="Percent",xlab="Years")
		legend("bottomright",paste("Equilibrium Latent TB =",anytb ,"%"),col="purple",cex=0.9,lwd=2)

	plot(1:1200/12,(1-(1-OutBurnIn[1:1200,"NInf"]/(OutBurnIn[1:1200,"NAll"]-OutBurnIn[1:1200,"NAnyTb"]))^12)*100,
		main="Annual Infection Risk",type="l",col="purple",lwd=2,ylab="Percent",xlab="Years")
		legend("bottomright",paste("Equilibrium Inf risk  =",infrisk ,"%"),col="purple",cex=0.9,lwd=2)

	plot(1:1200/12,OutBurnIn[1:1200,"NTbMort"]/OutBurnIn[1:1200,"NAll"]*100000*12,
		main="TB Mortality",type="l",col="blue",lwd=2,ylab="Cases per 100,000",xlab="Years")
		legend("bottomright",paste("Equilibrium Mortality =",mort,"per 100,000"),col="blue",cex=0.9,lwd=2)


	# dput(V1950,paste("V1950_",Setting,".Rdata",sep=""))
	dput(V1950,"V1950.Rdata")



