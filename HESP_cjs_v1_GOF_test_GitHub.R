#######################################
#Goodness-of-fit tests for HESP   #
#Survival analysis			  #
#Simulating data for posterior	  #
#predictive checks		        #
#E. Hunter, May 2021  		  #
#######################################

rm(list=ls())


##################################
##########  LOAD PACKAGES

library(nimble)


###################################
##########   READ IN DATA

load("nimble_data_2022-01-11.RData")

#MCMC for all 4 models
load("nimble_models_2021-06-18.RData")

x <- mcmc.out.p.phi$samples
mcmc.out.samp <- as.data.frame(rbind(x$chain1, x$chain2, x$chain3))
mcmc.out.samp <- mcmc.out.samp[,c(1:33, 10726:10795)]


###################################
##########   SIMULATE CAPTURE HISTORY
#
# (Following Ergon and Gardner 2013):
# The function can be used for posterior predictive checks by iteratively
# sampling 'par' from the posterior distribution, generating new data (keeping
# the other function arguments (design variables) corresponding to the original 
# data), and comparing aspects of the simulated data with the original data.

# ARGUMENTS:
# - par: List of parameters with the same definition as in the Nimble model
# - others: Design variables and data with same definition as in Nimble model

logit <- function(p) {
	log(p/(1-p))
	}

sim.cap <- function (par.sample, nsites, ninds, firsts, nperiods, nsessions, mass, interval, effort) {
	
	#Set up "data holding" arrays
	alive <- array(NA, dim=c(nsites, max(ninds), nperiods))
	mu.phi <- array(NA, dim=c(nsites, max(ninds), nperiods))
	phi <- array(NA, dim=c(nsites, max(ninds), nperiods))
	transition <- array(NA, dim=c(nsites, max(ninds), nperiods))
	mualive <- array(NA, dim=c(nsites, max(ninds), nperiods))
	thisp0 <- array(NA, dim=c(nsites, max(ninds), nperiods, max(nsessions)))
	muy <- array(NA, dim=c(nsites, max(ninds), nperiods, max(nsessions)))
	y.new <- array(NA, dim=c(nsites, max(ninds), nperiods, max(nsessions)))
	
	#Fill in missing data probabalistically:
	miss.mass <- which(is.na(mass), arr.ind=T)
	mass.f <- mass
	for(m in 1:nrow(miss.mass)){
		mass.f[miss.mass[m,][1], miss.mass[m,][2], miss.mass[m,][3]] <- rnorm(1, 0, 1/par.sample$mass.prec)
	}


	#Survival
    for(site in 1:nsites){
    for(ind in 1:ninds[site]){  # loop through individuals
		alive[site,ind,firsts[site,ind]] <- 1  #At first capture, individual is alive 
		
		if(firsts[site,ind] < 11){  #Only for individuals caught before last period
		for(period in firsts[site,ind]:(nperiods-1)){
			mu.phi[site,ind,period] <- logit(par.sample$phi0)  + par.sample$phi.mass.eff*mass.f[site,ind,period] #+ unlist(par.sample$phi.site.prd.eff[site,period]) 
			phi[site,ind,period] <- 1 / (1 + exp(-1 * mu.phi[site,ind,period]))
		}
		
		
    		for(period in (firsts[site,ind] + 1):nperiods){
    			transition[site,ind,period] <- phi[site,ind,(period-1)] ^ interval[period-1]
    			mualive[site,ind,period] <- alive[site,ind,(period-1)] * transition[site,ind,period]
			alive[site,ind,period] <- rbinom(1, 1, mualive[site,ind,period])
		}
		}
		
	#Capture
    		for(period in firsts[site,ind]:nperiods){
    			for(session in 1:nsessions[period]){
      			thisp0[site,ind,period,session] <- 1 / (1 + exp(-1 * (logit(par.sample$p0) + par.sample$p.effort.eff*effort[site,period,session] #+ unlist(par.sample$p.site.prd.eff[site,period])
      			)))
				muy[site,ind,period,session] <- alive[site,ind,period] * thisp0[site,ind,period,session]
				y.new[site,ind,period,session] <- rbinom(1, 1, muy[site,ind,period,session])
    			}
    		}
	}
	}

	y.new
	
}

#############################
#Generate new capture history from draws from posterior distribution

dat.props <- NULL
sim.props <- NULL
bayes.p <- NULL
sims <- 100

for (i in 1:sims){
	
#Parameter samples from MCMC
psamp <- mcmc.out.samp[sample(1:nrow(mcmc.out.samp),1),]

#for site-year effects, need to put in the correct indexing
phi.site.prd <- psamp[ , grepl( "phi.site.prd.eff" , names( psamp ) ) ]
phi.site.prd.eff <- matrix(phi.site.prd, nrow=nsites, ncol=nperiods-1)
p.site.prd <- psamp[ , grepl( "p.site.prd.eff" , names( psamp ) ) ]
p.site.prd.eff <- matrix(p.site.prd, nrow=nsites, ncol=nperiods)

par.sample <- list(
	mass.prec = psamp$mass.prec,
	phi0 = psamp$phi0,
	phi.mass.eff = psamp$phi.mass.eff,
	phi.site.prd.eff = phi.site.prd.eff,
	p0 = psamp$p0,
	p.effort.eff = psamp$p.effort.eff,
	p.site.prd.eff = p.site.prd.eff
	)
	
dat.new <- sim.cap(par.sample = par.sample, nsites, ninds, firsts, nperiods, nsessions, mass, interval, effort)


#Test 2 (from Ergon and Gardner). Survival: For each primary session, count number of individuals known to be alive (captured at least once)
#and proportion of those individuals that were seen at any later primary session
dat.prop <- array(0, dim=c(nsites, nperiods-2))
sim.prop <- array(0, dim=c(nsites, nperiods-2))
	
for (site in 1:nsites){
for (period in 1:(nperiods-2)){
	temp <- y[site,,period,]
	temp.sum <- apply(temp, 1, function(x) sum(x, na.rm=TRUE))
	temp.count <- length(temp.sum[temp.sum>0])
	temp.id <- which(temp.sum>0)
	temp2 <- y[site,temp.id,(period+1):nperiods,]
	temp2.sum <- apply(temp2, 1, function(x) sum(x, na.rm=TRUE))
	temp2.count <- length(temp2.sum[temp2.sum>0])
	dat.prop[site,period] <- temp2.count / temp.count


	temp <- dat.new[site,,period,]
	temp.sum <- apply(temp, 1, function(x) sum(x, na.rm=TRUE))
	temp.count <- length(temp.sum[temp.sum>0])
	temp.id <- which(temp.sum>0)
	temp2 <- dat.new[site,temp.id,(period+1):nperiods,]
	temp2.sum <- apply(temp2, 1, function(x) sum(x, na.rm=TRUE))
	temp2.count <- length(temp2.sum[temp2.sum>0])
	sim.prop[site,period] <- ifelse(temp2.count>0, temp2.count / temp.count, 0)
}
}

	
#Bayesian p-value (proportion of data points > 1:1 ratio (this will be averaged across many draws from posteriors)
dat.props <- c(dat.props, dat.prop)
sim.props <- c(sim.props, sim.prop)
bayes.p <- c(bayes.p, mean(sim.prop > dat.prop))

}


mean(bayes.p)
plot(sim.props ~ dat.props)
abline(0, 1)

#Bayesian p-value = 0.46





#