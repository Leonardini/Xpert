########################################################
  ##                                                ##
  ##         TB XPERT DIAGNOSTIC MODEL 2011         ##
  ##                                                ##
########################################################

	library(logitnorm)
	library(foreign)
	library(base)
	library(IMIS)
	library(mvtnorm)
	library(lhs)

	rm(list=ls())
	setwd("C:/Users/nick/Documents/Harvard/TB Diagnostics/ANALYSIS")
	options(digits=4)

################## 
### FUNCTIONS FOR GETTING PSA PARAMETER VALUES  ########################
################## 

# Function for calculating lognormal parameters 

	lnormpar <- function(tgt) {
		tgt <- as.numeric(tgt)
		xu <- tgt[1]; pctr <- (tgt[3]-tgt[2])
		if(pctr>0) {
		xopt <- function(xsd,xu,pctr) {
			xq <- qlnorm(c(0.025,0.975),meanlog=log(xu)-log(1+xsd^2/exp((2*log(xu))))/2,sdlog=log(1+xsd^2/exp((2*log(xu))))^0.5)
			cir <- (xq[2]-xq[1]); return((cir-pctr)^2)  }
		xsd <- as.numeric(optimise(xopt,interval=c(0,xu),pctr=pctr,xu=xu)[1])
		mu <- log(xu)-log(1+xsd^2/exp((2*log(xu))))/2
		sd <- log(1+xsd^2/exp((2*log(xu))))^0.5
		return(c(mu,sd))  }
		else { return(c(log(xu),0)) }   } 

# Function for calculating logitnormal parameters 

	lgtnormpar <- function(tgt) {
		if((tgt[3]-tgt[2])>0) {
		#  qfunct calculates the mean and CI range of a logitnorm with given parameters
		qfunct <- function(mean,sd) {    
		zz <- qnorm(randomLHS(10000,1),mean,sd)
		zz <- exp(zz)/(1+exp(zz))
		c(mean(zz),quantile(zz,0.975)-quantile(zz,0.025)) }

		xopt <- function(x,tgt) {
			tgt <- as.numeric(tgt)
			t1 <- tgt[1]; t2 <- (tgt[3]-tgt[2])
			z <- qfunct(x[1],x[2])
			return(c(z[1]-t1,z[2]-t2)%*%c(z[1]-t1,z[2]-t2))  }
		jp <- optim(c(0.5,0.1),xopt,method="L-BFGS-B", control=list(maxit=1000),
			lower = c(-100,0.0001), upper = c(100,100),tgt=tgt)
		return(jp$par[1:2])  } 
		else { return(c(log(tgt[1]/(1-tgt[1])),0)) }}

################## 
### CREATING TABLE OF PARAMETERS FOR POINT ESTIMATE RUNS  ########################
################## 

	ParamInit <- as.data.frame(read.csv("ParamInit4.csv")[,2:6]); rownames(ParamInit) <- read.csv("ParamInit4.csv")[,1]
	MeanParam <- 1:nrow(ParamInit); names(MeanParam) <- rownames(ParamInit)
	for (i in 1:nrow(ParamInit)) {
		if(ParamInit[i,5]==1) 	{ MeanParam[i] <- log(ParamInit[i,1]/(1-ParamInit[i,1])) }
			else 			{ MeanParam[i] <- log(ParamInit[i,1])  }
		if(MeanParam[i]==Inf) 	{ MeanParam[i] <- 100   }
		if(MeanParam[i]==-Inf) { MeanParam[i] <- -100  }
		}

	write.csv(MeanParam, file="MeanParam.csv")

################## 
### CREATING TABLE OF PARAMETERS FOR PSA  ########################
################## 

	ParamInit <- as.data.frame(read.csv("ParamInit4.csv")[,2:6]); rownames(ParamInit) <- read.csv("ParamInit4.csv")[,1]
	PSAparam <- cbind(ParamInit,rep(NA,nrow(ParamInit)),rep(NA,nrow(ParamInit)))
	colnames(PSAparam) <- c(colnames(ParamInit),"Par1","Par2")
	rownames(PSAparam) <- read.csv("ParamInit4.csv")[,1]

	for (i in 1:nrow(PSAparam)) {
		if(PSAparam[i,5]==1) 	{ PSAparam[i,6:7] <- lgtnormpar(as.matrix(PSAparam[i,2:4])) }
			else 			{ PSAparam[i,6:7] <- lnormpar(as.matrix(PSAparam[i,2:4]))  } }
	PSAparam[PSAparam[,6]==-Inf,6] <- -100;	PSAparam[PSAparam[,6]==Inf,6] <- 100
	
	write.table(PSAparam, file="PSAparam.csv")

## Output and save SA ranges ################## 

	logitfunct <- function(mean,sd) {    
		zz <- qnorm(randomLHS(10000,1),mean,sd)
		zz <- exp(zz)/(1+exp(zz))
		as.numeric(c(mean(zz),quantile(zz,0.025),quantile(zz,0.975))) }

	PSArange <- matrix(NA,nrow(PSAparam),3)
		colnames(PSArange) <- c("mean","CIlow","CIhigh")
		rownames(PSArange) <- rownames(PSAparam) 

	for (i in 1:nrow(PSArange)) {
	if(PSAparam[i,5]==1) { PSArange[i,1:3] <- logitfunct(PSAparam[i,6],PSAparam[i,7]) }
	else { PSArange[i,1:3] <- c(exp(PSAparam[i,6]+0.5*PSAparam[i,7]^2),qlnorm(c(0.025,0.975),PSAparam[i,6],PSAparam[i,7])) }  }

	write.table(PSArange, file="PSArange.csv")

################## 
## CREATING sample.prior(n) WHICH DRAWS n SAMPLES FROM PRIOR ########################### 
################## 
## Latin hypercube sample design:
	sample.prior1 <- function(n)  {
		if(n>1) { normdraw <- unifdraw <- randomLHS(n,nrow(PSAparam)) }
		else normdraw <- unifdraw <- as.matrix(t(runif(nrow(PSAparam))))
		for (i in 1:ncol(unifdraw)) {
			normdraw[,i] <- qnorm(unifdraw[,i],PSAparam[i,6],PSAparam[i,7])  }
		colnames(normdraw) <- rownames(PSAparam); normdraw		
		}

# Random sample design:
	sample.prior2 <- function(n)  { rmvnorm(n,PSAparam[,6],diag(PSAparam[,7]))  }

# Using latin hypercube for the moment...	
	sample.prior <- sample.prior1

################## 
## GETTING MATRIX OF 20000 PRIOR DRAWS  ########################### 
################## 
	
	set.seed(123)
	PriorDraws <- sample.prior(20000)
	save(PriorDraws, file="PriorDraws3-19.rData")

################## 
## AND A MATRIX OF THE UNTRANSFORMED DRAWS  ########################### 
################## 
	PriorDrawsU <- PriorDraws
	for (i in 1:ncol(PriorDrawsU)) {
	if(PSAparam[i,5]==1) { 
		PriorDrawsU[,i] <- exp(PriorDraws[,i])/(1+exp(PriorDraws[,i]))
		}
	else { 
		PriorDrawsU[,i] <- exp(PriorDraws[,i])
	 	}  }

	save(PriorDrawsU, file="PriorDrawsU3-19.rData")

########################### DONE ########################### 

