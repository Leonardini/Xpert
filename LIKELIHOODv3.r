##########   TB XPERT DIAGNOSTIC MODEL 2012  ##########   
##########             CALIB FILE            ########## 

library(mvtnorm)

############
####  Bring in data #########
############

### WHO Prevalence and Incidence Estimates, 5 year increments
	PrevMean 	<- matrix(NA,5,5)
	rownames(PrevMean) <- c("1990","1995","2000","2005","2010"); colnames(PrevMean) <- c("Bot","Les","Nam","Sou","Swa")
	PrevSD <- IncidMean <- IncidSD <- PrevMean
	PrevMean[,1] <- c(357,352,457,561,531); PrevSD[,1] <- c(141,119,153,181,148); IncidMean[,1] <- c(307,444,640,770,694); IncidSD[,1] <- c(58,42,58,79,69)
	PrevMean[,2] <- c(209,253,364,416,405); PrevSD[,2] <- c(85 ,92 ,134,155,141); IncidMean[,2] <- c(184,323,553,639,634); IncidSD[,2] <- c(28,22,38,58,45)
	PrevMean[,3] <- c(490,486,525,616,588); PrevSD[,3] <- c(205,172,190,221,204); IncidMean[,3] <- c(322,465,671,808,727); IncidSD[,3] <- c(71,47,55,58,69)
	PrevMean[,4] <- c(428,359,535,801,808); PrevSD[,4] <- c(176,119,167,252,236); IncidMean[,4] <- c(301,317,576,925,971); IncidSD[,4] <- c(55,32,59,94,96)
	PrevMean[,5] <- c(326,308,506,643,673); PrevSD[,5] <- c(141,110,185,236,206); IncidMean[,5] <- c(267,337,801,1141,1257); IncidSD[,5] <- c(61,34,82,116,125)

### MDR in new and existing cases
	MdrNpos	<- matrix(NA,3,5)
	rownames(MdrNpos) <- c("1995","2002","2009"); colnames(MdrNpos) <- c("Bot","Les","Nam","Sou","Swa")
	MdrRsamp <- MdrRpos <- MdrNsamp <- MdrNpos
	MdrNpos[,1] <- c(NA,10,NA); MdrNsamp[,1] <- c(NA,1182,NA); MdrRpos[,1] <- c(NA,11,NA); MdrRsamp[,1] <- c(NA,106,NA)
	MdrNpos[,2] <- c(3,NA,NA); MdrNsamp[,2] <- c(330,NA,NA); MdrRpos[,2] <- c(3,NA,NA); MdrRsamp[,2] <- c(53,NA,NA)
	MdrNpos[,3] <- c(NA,0.01*50,NA); MdrNsamp[,3] <- c(NA,50,NA); MdrRpos[,3] <- c(NA,0.08*50,NA); MdrRsamp[,3] <- c(NA,50,NA)
	MdrNpos[,4] <- c(NA,77,NA); MdrNsamp[,4] <- c(NA,4243,NA); MdrRpos[,4] <- c(NA,98,NA); MdrRsamp[,4] <- c(NA,1465,NA)
	MdrNpos[,5] <- c(3,NA,NA); MdrNsamp[,5] <- c(334,NA,NA); MdrRpos[,5] <- c(4,NA,NA); MdrRsamp[,5] <- c(44,NA,NA)
 
### Notifications data
	NotifLikDat <- dget("NotifLikDat.rData")

############
####  Log-Likelihoods #########
############

	PrevLnLik <- function(Dat,Setting,n)  { # Dat = full vector of prevalance projections, n= no time points
		c <- substr(Setting,1,3); z <- rep(NA,5)
		z[1]	<- dnorm(Dat[40],PrevMean[1,c],PrevSD[1,c],log=T)
		z[2]	<- dnorm(Dat[45],PrevMean[2,c],PrevSD[2,c],log=T)
		z[3]	<- dnorm(Dat[50],PrevMean[3,c],PrevSD[3,c],log=T)
		z[4]	<- dnorm(Dat[55],PrevMean[4,c],PrevSD[4,c],log=T)
		z[5]	<- dnorm(Dat[59],PrevMean[5,c],PrevSD[5,c],log=T)
	
		if(n==2) {  return(sum(z[c(1,5)],na.rm=T))  }
		if(n==3) {  return(sum(z[c(1,3,5)],na.rm=T))  }
		if(n==5) {  return(sum(z,na.rm=T))  }  }

	IncidLnLik <- function(Dat,Setting,n)  { # Dat = full vector of incidence projections, n= no time points
		c <- substr(Setting,1,3); z <- rep(NA,5)
		z[1]	<- dnorm(Dat[40],IncidMean[1,c],IncidSD[1,c],log=T)
		z[2]	<- dnorm(Dat[45],IncidMean[2,c],IncidSD[2,c],log=T)
		z[3]	<- dnorm(Dat[50],IncidMean[3,c],IncidSD[3,c],log=T)
		z[4]	<- dnorm(Dat[55],IncidMean[4,c],IncidSD[4,c],log=T)
		z[5]	<- dnorm(Dat[59],IncidMean[5,c],IncidSD[5,c],log=T)

		if(n==2) {  return(sum(z[c(1,5)],na.rm=T))  }
		if(n==3) {  return(sum(z[c(1,3,5)],na.rm=T))  }
		if(n==5) {  return(sum(z,na.rm=T))  }  }

	MdrLnLik <- function(Dat1,Dat2,Setting,d)  {  # Dat1 = MDR proportion in New TB, Dat2 = MDR proportio in Retreatment TB, d = design effect
		c <- substr(Setting,1,3); z <- rep(NA,6)
		z[1]		<- dbeta(Dat1[45],MdrNpos[1,c]/d,(MdrNsamp[1,c]-MdrNpos[1,c])/d,log=T)	
		z[2]		<- dbeta(Dat1[52],MdrNpos[2,c]/d,(MdrNsamp[2,c]-MdrNpos[2,c])/d,log=T)
		z[3]		<- dbeta(Dat1[59],MdrNpos[3,c]/d,(MdrNsamp[3,c]-MdrNpos[3,c])/d,log=T)

		z[4]		<- dbeta(Dat2[45],MdrRpos[1,c]/d,(MdrRsamp[1,c]-MdrRpos[1,c])/d,log=T)
		z[5]		<- dbeta(Dat2[52],MdrRpos[2,c]/d,(MdrRsamp[2,c]-MdrRpos[2,c])/d,log=T)
		z[6]		<- dbeta(Dat2[59],MdrRpos[3,c]/d,(MdrRsamp[3,c]-MdrRpos[3,c])/d,log=T)

		return(sum(z,na.rm=T))  }  

	NotifLnLik <- function(Dat,Setting,extra)  { # Dat = vector of Notifications per 100,000 pop, extra = additional nromal noise, with sd=extra*mean
		NotifFit 	<- NotifLikDat[-94,substr(Setting,1,3),]
		VNlik 	<- rep(NA,93)
		for(f in 1:93) { if(is.na(NotifFit[f,1])!=TRUE) { 
			VNlik[f] <- dnorm(Dat[f],NotifFit[f,1],NotifFit[f,3]*extra,log=T)  }  }
		
		return(sum(VNlik,na.rm=T))    }	



