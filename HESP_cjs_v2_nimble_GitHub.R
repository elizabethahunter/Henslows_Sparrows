########################################
#Henslow's Sparrow			   #
#CJS mark-recapture analysis       #
#E. Hunter, Apr 2021       	       #
########################################

rm(list=ls())


##################################
##########  LOAD PACKAGES

library(nimble)

###################################
##########   READ IN DATA

load("nimble_data_2022-01-11.RData")


###############
# NIMBLE MODEL - baseline (no time-varying parameters for p or phi)
###############

hesp.code <- nimbleCode({
    
    #############
    # PRIORS
    #############
    
    ###########
    # PROBABILITY OF CAPTURE IN SURVEYED AREA
    ###########
    
    p0 ~ dunif(0,1)           # mean/intercept detection prob within the surveyed area
    logit.p0 <- log(p0/(1-p0))   # convert to logit
    p.effort.eff ~ dunif(-5,5)
    
    #### mean capture probability
    for(site in 1:nsites){
      #p.site.eff[site] ~ dnorm(0, p.site.prec)
    for(ind in 1:ninds[site]){  # loop through (data augmented) individuals
    for(period in firsts[site,ind]:nperiods){
    for(session in 1:nsessions[period]){
      logit(thisp0[site,ind,period,session]) <- logit.p0 + p.effort.eff*effort[site,period,session] 
	muy[site,ind,period,session] <- alive[site,ind,period] * thisp0[site,ind,period,session]
	y[site,ind,period,session] ~ dbern(muy[site,ind,period,session])
    }
    }
    }
    }
    
    ###########
    # PROBABILITY OF SURVIVING/AGING 
    ###########
    
	#Priors
    phi0 ~ dunif(0.05,1)                 
    phi0.logit <- log(phi0/(1-phi0)) 
    phi.mass.eff ~ dunif(-3,3)

    #############
    # DEAL WITH MISSING DATA
    #############

	#Mass  
	mass.prec ~ dgamma(0.01,0.01)
	mass.sd <- pow((1/mass.prec),0.5)
	for(site in 1:nsites){
	for(ind in 1:ninds[site]){
	for(period in 1:nperiods){
		mass[site,ind,period] ~ dnorm(0, mass.prec)
	}}}

    
    #############
    # LIKELIHOOD
    #############

    for(site in 1:nsites){
    for(ind in 1:ninds[site]){  # loop through individuals
    	for(period in 1:firsts[site,ind]){
			alive[site,ind,period] ~ dbern(1)  #At first capture, individual is alive 
		}
    	for(period in (firsts[site,ind] + 1):nperiods){
			alive[site,ind,period] ~ dbern(mualive[site,ind,period])
			mualive[site,ind,period] <- alive[site,ind,(period-1)] * transition[site,ind,period]
			transition[site,ind,period] <- pow(phi[site,ind,(period-1)], interval[period-1])
		}

		for(period in firsts[site,ind]:(nperiods-1)){
			logit(phi[site,ind,period]) <- mu.phi[site,ind,period]
			mu.phi[site,ind,period] <- phi0.logit  + phi.mass.eff*mass[site,ind,period]
		}

    }
}   
    
    #################
    ####  CALCULATE ABUNDANCE USING HORVITZ-THOMPSON ESTIMATOR....

      for(site in 1:nsites){
	for(ind in 1:ninds[site]){
      for(period in 1:nperiods){
        mu.p2[site,ind,period,1] <- logit.p0 + p.effort.eff*effort[site,period,1] 
        p2[site,ind,period,1] <- 1/(1+exp(-1*mu.p2[site,ind,period,1]))   # back to prob. scale
        pncap[site,ind,period,1] <- 1-p2[site,ind,period,1]        # pncap refers to the probability of not capturing for session
        for(session in 2:nsessions[period]){
          mu.p2[site,ind,period,session] <- logit.p0 + p.effort.eff*effort[site,period,session] 
          p2[site,ind,period,session] <- 1/(1+exp(-1*mu.p2[site,ind,period,session]))   # back to prob. scale
          pncap[site,ind,period,session] <- pncap[site,ind,period,(session-1)]*(1-p2[site,ind,period,session])
        }
        pcap[site,ind,period] <- 1-pncap[site,ind,period,nsessions[period]]     # pcap refers to the prob of being captured at least once in any session
        invpcap[site,ind,period] <- 1/pcap[site,ind,period]    # for H-T estimator
      } 
    }


    for(period in 1:nperiods){
      Ntot[site,period] <- inprod(invpcap[site,1:ninds[site],period],caphist3d[site,1:ninds[site],period])    # estimate of total abundance within sampled region... 
    }
  }
    
    })   ## end Nimble model
    

######################
# PREPARE DATA FOR NIMBLE
######################

data.for.bugs <- list(
  y = y,
  caphist3d = caphist3d,
  effort = effort,
  mass = mass,
  interval = interval
)

constants <- list(
  nsites = nsites,
  ninds = ninds,
  nperiods = nperiods,
  nsessions = nsessions,
  firsts = firsts
)

#init for "alive"
Z <- array(NA, dim=c(nsites, max(ninds), nperiods))  
for(s in 1:nsites){
	ch = apply(caphist3d[s,,], c(1,2), function(i) any(i!=0))
	first.last = apply(ch, 1, function(i) range(which(i)))  #First and last periods trapped
	for (i in 1:ninds[s]){    # nindG
		Z[s,i,first.last[1,i]:first.last[2,i]] = 1    #1 when known to be alive, 0 otherwise
	}
}

initz.bugs<-function(){
  list(
    alive=Z,    
    p0=runif(1,0.1,0.15),
    p.effort.eff=rnorm(1,0,1),
    phi0=runif(1,0.1,0.5),  
    phi.mass.eff = runif(1, -1, 1),
    mass.prec = runif(1, 0.01, 0.1)
  )
}


mcmc.out <- nimbleMCMC(code = hesp.code, constants = constants,
                       data = data.for.bugs, inits = initz.bugs,
                       nchains = 3, nburnin = 5000, niter = 10000, thin=5,
                       summary = TRUE, WAIC = TRUE,
                       monitors = c("p0", "p.effort.eff", "phi0", "phi.mass.eff", "mass.prec", "Ntot"))


names(mcmc.out)

mcmc.out$WAIC

rownames(mcmc.out$summary$all.chains)
mcmc.out$summary$all.chains




###############
# NIMBLE MODEL - model 2 (time-varying parameters for p but not phi)
###############

hesp.code.p <- nimbleCode({
    
    #############
    # PRIORS
    #############
    
    ###########
    # PROBABILITY OF CAPTURE IN SURVEYED AREA
    ###########
    
    p0 ~ dunif(0,1)           # mean/intercept detection prob within the surveyed area
    logit.p0 <- log(p0/(1-p0))   # convert to logit
    p.effort.eff ~ dunif(-5,5)
    p.site.prd.prec ~ dgamma(0.01, 0.01)
    p.site.prd.sd <- pow((1/p.site.prd.prec), 0.5)

    for(site in 1:nsites){
    for(period in 1:nperiods){
	p.site.prd.eff[site,period] ~ dnorm(0, p.site.prd.prec)
    }
    }
    
    #### mean capture probability
    for(site in 1:nsites){
    for(ind in 1:ninds[site]){  # loop through (data augmented) individuals
    for(period in firsts[site,ind]:nperiods){
    for(session in 1:nsessions[period]){
      logit(thisp0[site,ind,period,session]) <- logit.p0 + p.effort.eff*effort[site,period,session] + p.site.prd.eff[site,period] 
	muy[site,ind,period,session] <- alive[site,ind,period] * thisp0[site,ind,period,session]
	y[site,ind,period,session] ~ dbern(muy[site,ind,period,session])
    }
    }
    }
    }
   
    
    ###########
    # PROBABILITY OF SURVIVING/AGING 
    ###########
    
	#Priors
    phi0 ~ dunif(0.05,1)                 
    phi0.logit <- log(phi0/(1-phi0)) 
    phi.mass.eff ~ dunif(-3,3)
    
    
    #############
    # DEAL WITH MISSING DATA
    #############

	#Mass  
	mass.prec ~ dgamma(0.01,0.01)
	mass.sd <- pow((1/mass.prec),0.5)
	for(site in 1:nsites){
	for(ind in 1:ninds[site]){
	for(period in 1:nperiods){
		mass[site,ind,period] ~ dnorm(0, mass.prec)
	}}}

    
    #############
    # LIKELIHOOD
    #############

    for(site in 1:nsites){
    for(ind in 1:ninds[site]){  # loop through individuals
    	for(period in 1:firsts[site,ind]){
			alive[site,ind,period] ~ dbern(1)  #At first capture, individual is alive 
		}
    	for(period in (firsts[site,ind] + 1):nperiods){
			alive[site,ind,period] ~ dbern(mualive[site,ind,period])
			mualive[site,ind,period] <- alive[site,ind,(period-1)] * transition[site,ind,period]
			transition[site,ind,period] <- pow(phi[site,ind,(period-1)], interval[period-1])
		}

		for(period in firsts[site,ind]:(nperiods-1)){
			logit(phi[site,ind,period]) <- mu.phi[site,ind,period]
			mu.phi[site,ind,period] <- phi0.logit  + phi.mass.eff*mass[site,ind,period] 
		}

    }
}
    
    
    #################
    ####  CALCULATE ABUNDANCE USING HORVITZ-THOMPSON ESTIMATOR....

      for(site in 1:nsites){
	for(ind in 1:ninds[site]){
      for(period in 1:nperiods){
        mu.p2[site,ind,period,1] <- logit.p0 + p.effort.eff*effort[site,period,1] + p.site.prd.eff[site,period] 
        p2[site,ind,period,1] <- 1/(1+exp(-1*mu.p2[site,ind,period,1]))   # back to prob. scale
        pncap[site,ind,period,1] <- 1-p2[site,ind,period,1]        # pncap refers to the probability of not capturing for session
        for(session in 2:nsessions[period]){
          mu.p2[site,ind,period,session] <- logit.p0 + p.effort.eff*effort[site,period,session] + p.site.prd.eff[site,period] 
          p2[site,ind,period,session] <- 1/(1+exp(-1*mu.p2[site,ind,period,session]))   # back to prob. scale
          pncap[site,ind,period,session] <- pncap[site,ind,period,(session-1)]*(1-p2[site,ind,period,session])
        }
        pcap[site,ind,period] <- 1-pncap[site,ind,period,nsessions[period]]     # pcap refers to the prob of being captured at least once in any session
        invpcap[site,ind,period] <- 1/pcap[site,ind,period]    # for H-T estimator
      } 
    }


    for(period in 1:nperiods){
      Ntot[site,period] <- inprod(invpcap[site,1:ninds[site],period],caphist3d[site,1:ninds[site],period])    # estimate of total abundance within sampled region... 
    }
  }
    
    })   ## end Nimble model
    



######################
# PREPARE DATA FOR NIMBLE
######################
#Initial values are different...

initz.bugs<-function(){
  list(
    alive=Z,    
    p0=runif(1,0.1,0.15),
    p.effort.eff=rnorm(1,0,1),
    p.site.prd.prec=runif(1, 0.01, 0.1),
    phi0=runif(1,0.1,0.5),  
    phi.mass.eff = runif(1, -1, 1),
    mass.prec = runif(1, 0.01, 0.1)
  )
}


mcmc.out.p <- nimbleMCMC(code = hesp.code.p, constants = constants,
                       data = data.for.bugs, inits = initz.bugs,
                       nchains = 3, nburnin = 5000, niter = 10000, thin=5,
                       summary = TRUE, WAIC = TRUE,
                       monitors = c("p0", "p.effort.eff", "p.site.prd.sd", "p.site.prd.eff", 
               "phi0", 
			"phi.mass.eff", "mass.prec", "alive", "Ntot"))


mcmc.out.p$WAIC



###############
# NIMBLE MODEL - model 3 (time-varying parameters for phi but not p)
###############

hesp.code.phi <- nimbleCode({
    
    #############
    # PRIORS
    #############
    
    ###########
    # PROBABILITY OF CAPTURE IN SURVEYED AREA
    ###########
    
    p0 ~ dunif(0,1)           # mean/intercept detection prob within the surveyed area
    logit.p0 <- log(p0/(1-p0))   # convert to logit
    p.effort.eff ~ dunif(-5,5)
    
    #### mean capture probability
    for(site in 1:nsites){
      #p.site.eff[site] ~ dnorm(0, p.site.prec)
    for(ind in 1:ninds[site]){  # loop through (data augmented) individuals
    for(period in firsts[site,ind]:nperiods){
    for(session in 1:nsessions[period]){
      logit(thisp0[site,ind,period,session]) <- logit.p0 + p.effort.eff*effort[site,period,session] 
	muy[site,ind,period,session] <- alive[site,ind,period] * thisp0[site,ind,period,session]
	y[site,ind,period,session] ~ dbern(muy[site,ind,period,session])
    }
    }
    }
    }
    
    ###########
    # PROBABILITY OF SURVIVING/AGING 
    ###########
    
	#Priors
    phi0 ~ dunif(0.05,1)                 
    phi0.logit <- log(phi0/(1-phi0)) 
    phi.mass.eff ~ dunif(-3,3)

   	#Survival varies by period 
	phi.site.prd.prec ~ dgamma(0.01,0.01)
    phi.site.prd.sd <- pow((1/phi.site.prd.prec),0.5)
    for(site in 1:nsites){
   	for(period in 1:(nperiods-1)){
     		 phi.site.prd.eff[site,period] ~ dnorm(0,phi.site.prd.prec)
   	}
	}
    
    #############
    # DEAL WITH MISSING DATA
    #############

	#Mass  
	mass.prec ~ dgamma(0.01,0.01)
	mass.sd <- pow((1/mass.prec),0.5)
	for(site in 1:nsites){
	for(ind in 1:ninds[site]){
	for(period in 1:nperiods){
		mass[site,ind,period] ~ dnorm(0, mass.prec)
	}}}


    
    #############
    # LIKELIHOOD
    #############

    for(site in 1:nsites){
    for(ind in 1:ninds[site]){  # loop through individuals
    	for(period in 1:firsts[site,ind]){
			alive[site,ind,period] ~ dbern(1)  #At first capture, individual is alive 
		}
    	for(period in (firsts[site,ind] + 1):nperiods){
			alive[site,ind,period] ~ dbern(mualive[site,ind,period])
			mualive[site,ind,period] <- alive[site,ind,(period-1)] * transition[site,ind,period]
			transition[site,ind,period] <- pow(phi[site,ind,(period-1)], interval[period-1])
		}

		for(period in firsts[site,ind]:(nperiods-1)){
			logit(phi[site,ind,period]) <- mu.phi[site,ind,period]
			mu.phi[site,ind,period] <- phi0.logit  + phi.mass.eff*mass[site,ind,period] + phi.site.prd.eff[site,period] 
		}

    }
}
    
    
    #################
    ####  CALCULATE ABUNDANCE USING HORVITZ-THOMPSON ESTIMATOR....

      for(site in 1:nsites){
	for(ind in 1:ninds[site]){
      for(period in 1:nperiods){
        mu.p2[site,ind,period,1] <- logit.p0 + p.effort.eff*effort[site,period,1] 
        p2[site,ind,period,1] <- 1/(1+exp(-1*mu.p2[site,ind,period,1]))   # back to prob. scale
        pncap[site,ind,period,1] <- 1-p2[site,ind,period,1]        # pncap refers to the probability of not capturing for session
        for(session in 2:nsessions[period]){
          mu.p2[site,ind,period,session] <- logit.p0 + p.effort.eff*effort[site,period,session] 
          p2[site,ind,period,session] <- 1/(1+exp(-1*mu.p2[site,ind,period,session]))   # back to prob. scale
          pncap[site,ind,period,session] <- pncap[site,ind,period,(session-1)]*(1-p2[site,ind,period,session])
        }
        pcap[site,ind,period] <- 1-pncap[site,ind,period,nsessions[period]]     # pcap refers to the prob of being captured at least once in any session
        invpcap[site,ind,period] <- 1/pcap[site,ind,period]    # for H-T estimator
      } 
    }


    for(period in 1:nperiods){
      Ntot[site,period] <- inprod(invpcap[site,1:ninds[site],period],caphist3d[site,1:ninds[site],period])    # estimate of total abundance within sampled region... 
    }
  }
    
    })   ## end Nimble model
    



######################
# PREPARE DATA FOR NIMBLE
######################
#Initial values are different...

initz.bugs<-function(){
  list(
    alive=Z,    
    p0=runif(1,0.1,0.15),
    p.effort.eff=rnorm(1,0,1),
    phi0=runif(1,0.1,0.5),  
    phi.site.prd.prec=runif(1,0.01, 0.1),
    phi.mass.eff = runif(1, -1, 1),
    mass.prec = runif(1, 0.01, 0.1)
  )
}


mcmc.out.phi <- nimbleMCMC(code = hesp.code.phi, constants = constants,
                       data = data.for.bugs, inits = initz.bugs,
                       nchains = 3, nburnin = 5000, niter = 10000, thin=5,
                       summary = TRUE, WAIC = TRUE,
                       monitors = c("p0", "p.effort.eff", 
                                      "phi0", "phi.site.prd.sd", "phi.site.prd.eff", 
                                      			"phi.mass.eff", "mass.prec", "alive", "Ntot"))


names(mcmc.out.phi)

mcmc.out.phi$WAIC




###############
# NIMBLE MODEL - model 4 (time-varying parameters for p and phi)
###############

hesp.code.p.phi <- nimbleCode({
    
    #############
    # PRIORS
    #############
    
    ###########
    # PROBABILITY OF CAPTURE IN SURVEYED AREA
    ###########
    
    p0 ~ dunif(0,1)           # mean/intercept detection prob within the surveyed area
    logit.p0 <- log(p0/(1-p0))   # convert to logit
    p.effort.eff ~ dunif(-5,5)
    p.site.prd.prec ~ dgamma(0.01, 0.01)
    p.site.prd.sd <- pow((1/p.site.prd.prec), 0.5)

    for(site in 1:nsites){
    for(period in 1:nperiods){
	p.site.prd.eff[site,period] ~ dnorm(0, p.site.prd.prec)
    }
    }
    
    #### mean capture probability
    for(site in 1:nsites){
    for(ind in 1:ninds[site]){  # loop through (data augmented) individuals
    for(period in firsts[site,ind]:nperiods){
    for(session in 1:nsessions[period]){
      logit(thisp0[site,ind,period,session]) <- logit.p0 + p.effort.eff*effort[site,period,session] + p.site.prd.eff[site,period] 
	muy[site,ind,period,session] <- alive[site,ind,period] * thisp0[site,ind,period,session]
	y[site,ind,period,session] ~ dbern(muy[site,ind,period,session])
    }
    }
    }
    }
    
    ###########
    # PROBABILITY OF SURVIVING/AGING 
    ###########
    
	#Priors
    phi0 ~ dunif(0.05,1)                 
    phi0.logit <- log(phi0/(1-phi0)) 
    phi.mass.eff ~ dunif(-3,3)

   	#Survival varies by period 
	phi.site.prd.prec ~ dgamma(0.01,0.01)
    phi.site.prd.sd <- pow((1/phi.site.prd.prec),0.5)
    for(site in 1:nsites){
   	for(period in 1:(nperiods-1)){
     		 phi.site.prd.eff[site,period] ~ dnorm(0,phi.site.prd.prec)
   	}
	}
   
    #############
    # DEAL WITH MISSING DATA
    #############

	#Mass  
	mass.prec ~ dgamma(0.01,0.01)
	mass.sd <- pow((1/mass.prec),0.5)
	for(site in 1:nsites){
	for(ind in 1:ninds[site]){
	for(period in 1:nperiods){
		mass[site,ind,period] ~ dnorm(0, mass.prec)
	}}}

    #############
    # LIKELIHOOD
    #############

    for(site in 1:nsites){
    for(ind in 1:ninds[site]){  # loop through individuals
    	for(period in 1:firsts[site,ind]){
			alive[site,ind,period] ~ dbern(1)  #At first capture, individual is alive 
		}
    	for(period in (firsts[site,ind] + 1):nperiods){
			alive[site,ind,period] ~ dbern(mualive[site,ind,period])
			mualive[site,ind,period] <- alive[site,ind,(period-1)] * transition[site,ind,period]
			transition[site,ind,period] <- pow(phi[site,ind,(period-1)], interval[period-1])
		}

		for(period in firsts[site,ind]:(nperiods-1)){
			logit(phi[site,ind,period]) <- mu.phi[site,ind,period]
			mu.phi[site,ind,period] <- phi0.logit  + phi.mass.eff*mass[site,ind,period] + phi.site.prd.eff[site,period] 
		}

    }
}
    
    
    #################
    ####  CALCULATE ABUNDANCE USING HORVITZ-THOMPSON ESTIMATOR....

      for(site in 1:nsites){
	for(ind in 1:ninds[site]){
      for(period in 1:nperiods){
        mu.p2[site,ind,period,1] <- logit.p0 + p.effort.eff*effort[site,period,1] + p.site.prd.eff[site,period] 
        p2[site,ind,period,1] <- 1/(1+exp(-1*mu.p2[site,ind,period,1]))   # back to prob. scale
        pncap[site,ind,period,1] <- 1-p2[site,ind,period,1]        # pncap refers to the probability of not capturing for session
        for(session in 2:nsessions[period]){
          mu.p2[site,ind,period,session] <- logit.p0 + p.effort.eff*effort[site,period,session] + p.site.prd.eff[site,period] 
          p2[site,ind,period,session] <- 1/(1+exp(-1*mu.p2[site,ind,period,session]))   # back to prob. scale
          pncap[site,ind,period,session] <- pncap[site,ind,period,(session-1)]*(1-p2[site,ind,period,session])
        }
        pcap[site,ind,period] <- 1-pncap[site,ind,period,nsessions[period]]     # pcap refers to the prob of being captured at least once in any session
        invpcap[site,ind,period] <- 1/pcap[site,ind,period]    # for H-T estimator
      } 
    }


    for(period in 1:nperiods){
      Ntot[site,period] <- inprod(invpcap[site,1:ninds[site],period],caphist3d[site,1:ninds[site],period])    # estimate of total abundance within sampled region... 
    }
  }
    
    })   ## end Nimble model
    



######################
# PREPARE DATA FOR NIMBLE
######################
#Initial values are different...

initz.bugs<-function(){
  list(
    alive=Z,    
    p0=runif(1,0.1,0.15),
    p.effort.eff=rnorm(1,0,1),
    p.site.prd.prec=runif(1, 0.01, 0.1),
    phi0=runif(1,0.1,0.5),  
    phi.site.prd.prec=runif(1,0.01, 0.1),
    phi.mass.eff = runif(1, -1, 1),
    mass.prec = runif(1, 0.01, 0.1)
  )
}


mcmc.out.p.phi <- nimbleMCMC(code = hesp.code.p.phi, constants = constants,
                       data = data.for.bugs, inits = initz.bugs,
                       nchains = 3, nburnin = 5000, niter = 10000, thin=5,
                       summary = TRUE, WAIC = TRUE,
                       monitors = c("p0", "p.effort.eff", "p.site.prd.sd", "p.site.prd.eff", 
                                      "phi0", "phi.site.prd.sd", "phi.site.prd.eff", 
                                      			"phi.mass.eff", "mass.prec", "alive", "Ntot"))


names(mcmc.out.p.phi)

mcmc.out.p.phi$WAIC
mcmc.out.p$WAIC
mcmc.out.phi$WAIC
mcmc.out$WAIC

#