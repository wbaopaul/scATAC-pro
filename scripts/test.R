


#######################################

### raw code to build SIBER package

#######################################



### (1) Explicitely deal with 0-inflation in NB and GP models

###     Rather than using BIC, we detect 0-inflated genes by empirical approach based on observed percentage of 0s.



### (2) actually 3 models are fitted: E, V and 0-inflated. When zeroPercentThr is not achieved, 0-inflated model is inactive. However,

###     when this is achieved, E and V are deactivated while 0-inflated model is chosen.



### (3) For output, all results return c('mu1', 'mu2', 'sigma1', 'sigma2', 'pi1', 'logLik', 'BIC') or 

###     c('mu1', 'mu2', 'phi1', 'phi2', 'pi1', 'logLik', 'BIC')

###     For 0-inflated models, mu1=phi1=0 or mu1=sigma1=0



### (4) note for the parameter inits:

###     This parameter is assumed to be a vector of length 5 corresponding to c('mu1', 'mu2', 'phi1', 'phi2', 'pi1'). Note LN model needs no initial values.

###     For E model, phi1=phi2

###     For 0-inflated model, only mu2, phi2, pi1 is used to initialize. The other elements can be NA or anything else which are not used.



## compute BIC

getBIC <- function(logLik, nPar, nObs) {
  
  -2*logLik+nPar*log(nObs)
  
}

# convert mean mu and variance Var to overdispersion phi by: Var=mu+phi*mu^2 (NB)

# or Var=mu phi (GP)

# modified on 07/31/12: add scenarios where dist=LN,Normal or Gaussian: just

# return sqrt(Var) since we only focus on parameterization on log scale in

# LNmix or normalMix

muVarToPhi <- function(mu, Var, dist='NB') {
  
  if(dist=='NB') {
    
    res <- (Var-mu)/mu^2
    
    res[which(res<0)] <- 0 # under-dispersion, really rare	
    
  }
  
  if(dist=='GP') {
    
    res <- Var/mu
    
  }
  
  if(dist=='LN'|dist=='Gaussian'|dist=='Normal') { # no action on normal or LN
    
    # since we've only focus on parameterization on log scale in LNmix or
    
    # normalMix: just return SD
    
    res <- sqrt(Var) # returned phi is sd itself
    
  }
  
  res
  
}



# given mu and phi, calculates variance. i.e. in NB(mu, phi), compute its variance

# as mu+phi*mu^2; for GP, var=phi*mu

muPhiToVar <- function(mu, phi, dist='NB') {
  
  if(dist=='NB') {
    
    res <- mu+phi*mu^2
    
  } 
  
  if(dist=='GP') {
    
    res <- phi*mu
    
  }
  
  res
  
}



# given mu and phi, calculates variance. i.e. in NB(mu, phi), compute its variance

# as mu+phi*mu^2; for GP, var=phi*mu

muPhiToVar <- function(mu, phi, dist='NB') {
  
  if(dist=='NB') {
    
    res <- mu+phi*mu^2
    
  } 
  
  if(dist=='GP') {
    
    res <- phi*mu
    
  }
  
  res
  
}





########################### NB mix ###############################



logLikNegBinOneCompByOverdisp <- function(mu=100, phi=10, y, d) {
  
  lgamma(1/phi+y)-lgamma(y+1)-lgamma(1/phi)-
    
    1/phi*log(phi*mu*d+1)+y*log(phi*mu*d)-y*log(phi*mu*d+1)
  
}

##### check:

## logLikNegBinOneCompByOverdisp(mu=200, phi=100, datasetsGP_H1[[1]][1, ], d=d)



logLikNegBinMixByOverdisp <- function(par, y=y, d, model='V') {
  
  if(model=='V') {
    
    mu1 <- par[1]; mu2 <- par[2]
    
    phi1 <- par[3]; phi2 <- par[4]
    
    pi1 <- par[5]
    
    res <- sum(log(pi1*exp(logLikNegBinOneCompByOverdisp(mu1, phi1, y, d)) +
                     
                     (1-pi1)*exp(logLikNegBinOneCompByOverdisp(mu2, phi2, y, d))))
    
  } else { # for the E model, only 4 parameters
    
    mu1 <- par[1]; mu2 <- par[2]
    
    phi <- par[3];
    
    pi1 <- par[4]	
    
    res <- sum(log(pi1*exp(logLikNegBinOneCompByOverdisp(mu1, phi, y, d)) +
                     
                     (1-pi1)*exp(logLikNegBinOneCompByOverdisp(mu2, phi, y, d))))
    
  }
  
  if(is.infinite(res)) {
    
    res <- -1e100 # very important
    
  }
  
  res
  
}



## logLik for 0-inflated model: 3 parameters: mu2, phi2, pi1

logLikNegBin0inflByOverdisp <- function(par, y=y, d) {
  
  mu2 <- par[1]
  
  phi2 <- par[2]
  
  pi1 <- par[3]
  
  res <- sum(log(pi1*as.numeric(y==0)+(1-pi1) *
                   
                   exp(logLikNegBinOneCompByOverdisp(mu2, phi2, y, d))))
  
  if(is.infinite(res)) {
    
    res <- -1e100 # very important
    
  }
  
  res
  
}



## add 0-inflation option zeroPercentThr=0.2

fitNB <- function(y, d=NULL, inits=NULL, model='V', zeroPercentThr=0.2) {
  
  model <- try(match.arg(model, c('E', 'V'), several.ok=FALSE), silent=TRUE)
  
  # stop if model not recognizable.
  
  if(class(model)=='try-error') stop('Only model E or V can be specified!\n')
  
  percent0 <- mean(y==0, na.rm=TRUE) 
  
  res <- rep(NA, 7)	# res=c(mu1, mu2, phi1, phi2, pi1, logLik, BIC)
  
  names(res) <- c('mu1', 'mu2', 'phi1', 'phi2', 'pi1', 'logLik', 'BIC')
  
  nPar <- ifelse(model=='V', 5, 4) # number of parameters
  
  
  
  if(is.null(d)) d <- rep(1, length(y)) #	default, no normalization
  
  
  
  # Case I: no severe 0-inflation, fit 2-comp E or V models
  
  if(percent0<=zeroPercentThr) { 
    
    if(model=='V') { # V model
      
      # initial value
      
      initials <- rep(NA, 5)
      
      if(is.null(inits)) {
        
        temp <- quantile(y, prob=c(1/8, 7/8), na.rm=T)
        
        initials[1] <- temp[1]
        
        initials[2] <- temp[2]
        
        overallPhi <- muVarToPhi(mean(y), var(y), dist='NB')
        
        initials[3] <- overallPhi*0.5
        
        initials[4] <- overallPhi*2
        
        initials[5] <- 0.5
        
      } else {
        
        initials <- inits
        
      }
      
      constrL <- c(1e-5, 1e-5, 1e-5, 1e-5, 0) # box constraint.
      
      #phi not set to exactly 0, mu not 0. otherwise, pdf not defined
      
      constrU <- c(5e7, 5e7, 1e3, 1e3, 1)
      
    } else { # E model
      
      # initial value
      
      initials <- rep(NA, 4)
      
      if(is.null(inits)) {
        
        temp <- quantile(y, prob=c(1/8, 7/8), na.rm=T)
        
        initials[1] <- temp[1]
        
        initials[2] <- temp[2]
        
        overallPhi <- muVarToPhi(mean(y), var(y), dist='NB')
        
        initials[3] <- overallPhi
        
        initials[4] <- 0.5
        
      } else {
        
        if(length(inits)==5) { # in case mu1, mu2, phi, phi, pi1 is specified,
          
          # we only need 4 elements
          
          initials <- inits[c(1:3, 5)]
          
        } else {
          
          initials <- inits # assumed 4 elements
          
        }
        
      }
      
      constrL <- c(0, 0, 1e-5, 0) # box constraint
      
      constrU <- c(5e7, 5e7, 1e3, 1)
      
    }
    
    optimRes <- try(optim(par=initials, logLikNegBinMixByOverdisp,
                          
                          y=y, d=d, model=model, lower=constrL, upper=constrU,
                          
                          control=list(fnscale=-1, pgtol=1e-16, factr=1e3,
                                       
                                       maxit=3000), method="L-BFGS-B"),
                    
                    silent=TRUE)
    
    if(class(optimRes)!='try-error') {
      
      if(model=='V') {
        
        res <- c(optimRes$par, optimRes$value)
        
      } else { # E model, formulate result such that it is like mu1, mu2,
        
        # phi1, phi2, pi1, logLik (phi1=phi2)
        
        temp <- optimRes$par
        
        res <- c(temp[1:3], temp[3], temp[4], optimRes$value)
        
      }
      
      res[7] <- getBIC(logLik=optimRes$value, nPar=nPar, nObs=length(y)) 
      
    }
    
  } else {
    
    # Case II: 0-inflation, override E or V models	
    
    # initial value
    
    initials <- rep(NA, 3)
    
    if(is.null(inits)) {
      
      initials[1] <- mean(y[y!=0], na.rm=TRUE)
      
      initials[2] <- muVarToPhi(mean(y[y!=0]), var(y[y!=0]), dist='NB')
      
      initials[3] <- 0.5
      
    } else {
      
      initials <- inits[, c(2, 4, 5)] # pick mu2, phi2, pi1 as initials
      
    }
    
    constrL <- c(1e-5, 1e-5, 0) 
    
    constrU <- c(5e7, 1e3, 1)
    
    optimRes <- try(optim(par=initials, logLikNegBin0inflByOverdisp,
                          
                          y=y, d=d, lower=constrL, 
                          
                          upper=constrU, control=list(fnscale=-1, pgtol=1e-16,
                                                      
                                                      factr=1e3, maxit=3000),
                          
                          method="L-BFGS-B"), silent=TRUE)
    
    if(class(optimRes)!='try-error') {
      
      temp <- optimRes$par
      
      res[1:6] <- c(0, temp[1], 0, temp[2], temp[3], optimRes$value)
      
      res[7] <- getBIC(logLik=optimRes$value, nPar=3, nObs=length(y)) 
      
    }
    
  }	
  
  res
  
}

##################### GP mix ###########################

## logLik for GenPoi parameterized by overdispersion and mean. this is a 

## temporary function since it is only for one component

# phi: needs phi>1

# d: normalization constant, a vector of length(y)

# returns a vector of logLik of length(y)



logLikGenPoiOneCompByOverdisp <- function(mu=100, phi=10, y, d) {
  
  log(mu*d/sqrt(phi)) + (y-1)*log(mu*d/sqrt(phi) + (1-1/sqrt(phi))*y) -
    
    mu*d/sqrt(phi) - (1-1/sqrt(phi))*y-lfactorial(y)
  
}



## logLik for GenPoi mixture distribution parameterized by overdispersion and V model

# par: in the order of (mu1, mu2, phi1, phi2, pi1)

# returns a scalar

logLikGenPoiMixByOverdisp <- function(par, y=y, d, model='V') {
  
  if(model=='V') {
    
    mu1 <- par[1]; mu2 <- par[2]
    
    phi1 <- par[3]; phi2 <- par[4]
    
    pi1 <- par[5]
    
    res <- sum(log(pi1*exp(logLikGenPoiOneCompByOverdisp(mu1, phi1, y, d))+
                     
                     (1-pi1)*exp(logLikGenPoiOneCompByOverdisp(mu2, phi2, y, d))))
    
  } else { # for the E model, only 4 parameters
    
    mu1 <- par[1]; mu2 <- par[2]
    
    phi <- par[3];
    
    pi1 <- par[4]	
    
    res <- sum(log(pi1*exp(logLikGenPoiOneCompByOverdisp(mu1, phi, y, d))+
                     
                     (1-pi1)*exp(logLikGenPoiOneCompByOverdisp(mu2, phi, y, d))))
    
  }
  
  if(is.infinite(res)) {
    
    res <- -1e100 # very important
    
  }
  
  res
  
}



logLikGenPoi0inflByOverdisp <- function(par, y=y, d) {
  
  mu2 <- par[1]
  
  phi2 <- par[2]
  
  pi1 <- par[3]
  
  res <- sum(log(pi1*as.numeric(y==0) +
                   
                   (1-pi1)*exp(logLikGenPoiOneCompByOverdisp(mu2, phi2, y, d))))
  
  if(is.infinite(res)) {
    
    res <- -1e100 # very important
    
  }
  
  res
  
}





###################

### fitting GP mixture of 2-component by direct optimization. E and V model

### are available parameterized by over-dispersion

###################

### suppressWarn: whether to suppress warning that comes from optim()

### convergence criterion

fitGP <- function(y, d=NULL, inits=NULL, model='V', zeroPercentThr=0.2) {
  
  model <- try(match.arg(model, c('E', 'V'), several.ok=FALSE), silent=TRUE)
  
  # stop if model not recognizable.
  
  if(class(model)=='try-error') stop('Only model E or V can be specified!\n')
  
  percent0 <- mean(y==0, na.rm=TRUE) 
  
  
  
  res <- rep(NA, 7)	# res=c(mu1, mu2, phi1, phi2, pi1, logLik, BIC)
  
  names(res) <- c('mu1', 'mu2', 'phi1', 'phi2', 'pi1', 'logLik', 'BIC')
  
  nPar <- ifelse(model=='V', 5, 4) # number of parameters
  
  
  
  if(is.null(d)) d <- rep(1, length(y)) #	default, no normalization
  
  
  
  # Case I: no severe 0-inflation, fit 2-comp E or V models
  
  if(percent0<=zeroPercentThr) { 
    
    if(model=='V') { # V model
      
      # initial value
      
      initials <- rep(NA, 5)
      
      if(is.null(inits)) {
        
        temp <- quantile(y, prob=c(1/8, 7/8), na.rm=T)
        
        initials[1] <- temp[1]
        
        initials[2] <- temp[2]
        
        overallPhi <- muVarToPhi(mean(y), var(y), dist='GP')
        
        initials[3] <- overallPhi*0.5
        
        initials[4] <- overallPhi*2
        
        initials[5] <- 0.5
        
      } else {
        
        initials <- inits
        
      }
      
      constrL <- c(1e-5, 1e-5, 1, 1, 0) # box constraint
      
      constrU <- c(5e7, 5e7, 1e10, 1e10, 1)
      
    } else { # E model
      
      # initial value
      
      initials <- rep(NA, 4)
      
      if(is.null(inits)) {
        
        temp <- quantile(y, prob=c(1/8, 7/8), na.rm=T)
        
        initials[1] <- temp[1]
        
        initials[2] <- temp[2]
        
        overallPhi <- muVarToPhi(mean(y), var(y), dist='GP')
        
        initials[3] <- overallPhi
        
        initials[4] <- 0.5
        
      } else {
        
        if(length(inits)==5) { # in case mu1, mu2, phi, phi, pi1 is specified, we only need 4 elements
          
          initials <- inits[c(1:3, 4)]
          
        } else {
          
          initials <- inits
          
        }
        
      }
      
      constrL <- c(1e-5, 1e-5, 1, 0) # box constraint
      
      constrU <- c(5e7, 5e7, 1e10, 1)
      
    }
    
    optimRes <- try(optim(par=initials, logLikGenPoiMixByOverdisp, y=y, d=d,
                          
                          model=model, lower=constrL, upper=constrU, 
                          
                          control=list(fnscale=-1, pgtol=1e-16, factr=1e3,
                                       
                                       maxit=3000), method="L-BFGS-B"),
                    
                    silent=TRUE)
    
    if(class(optimRes)!='try-error') {
      
      if(model=='V') {
        
        res <- c(optimRes$par, optimRes$value)
        
      } else { # E model, formulate result such that it is like mu1, mu2,
        
        # phi1, phi2, pi1, logLik (phi1=phi2)
        
        temp <- optimRes$par
        
        res <- c(temp[1:3], temp[3], temp[4], optimRes$value)
        
      }
      
      
      
      res[7] <- getBIC(logLik=optimRes$value, nPar=nPar, nObs=length(y)) 
      
      names(res) <- c('mu1', 'mu2', 'phi1', 'phi2', 'pi1', 'logLik', 'BIC')
      
    }
    
  } else {
    
    # Case II: 0-inflation, override E or V models	
    
    # initial value
    
    initials <- rep(NA, 3)
    
    if(is.null(inits)) {
      
      initials[1] <- mean(y[y!=0], na.rm=TRUE)
      
      initials[2] <- muVarToPhi(mean(y[y!=0]), var(y[y!=0]), dist='GP')
      
      initials[3] <- 0.5
      
    } else {
      
      initials <- inits[, c(2, 4, 5)] # pick mu2, phi2, pi1 as initials
      
    }
    
    constrL <- c(1e-5, 1, 0) 
    
    constrU <- c(5e7, 1e10, 1)
    
    optimRes <- try(optim(par=initials, logLikGenPoi0inflByOverdisp, y=y, d=d,
                          
                          lower=constrL, upper=constrU,
                          
                          control=list(fnscale=-1, pgtol=1e-16, factr=1e3,
                                       
                                       maxit=3000), method="L-BFGS-B"),
                    
                    silent=TRUE)
    
    if(class(optimRes)!='try-error') {
      
      temp <- optimRes$par
      
      res[1:6] <- c(0, temp[1], 0, temp[2], temp[3], optimRes$value)
      
      res[7] <- getBIC(logLik=optimRes$value, nPar=3, nObs=length(y)) 
      
    }
    
  }
  
  res
  
}







###################### LN mix ###########################

# d: normalization constant. The definition is different from TMM or RLE.

# Our d is the scale factor from true expression level to observed count.

# Hence, Cs~d*mu_{c(s)}. However, TMM is the scale factor from  observed

# count to true expression level, which is the opposite direction:

# Cs*TMM~mu_{c(s)}. Thus, d=1/TMM For LN model, we need fit normal mixture

# on log((Cs+eps)/d)=log(Cs+eps)-log(d)

##### capable dealing with 0-inflation

fitLN <- function(y, base=10, eps=10, d=NULL, model='E',
                  
                  zeroPercentThr=0.2, logLikToLN=TRUE) {
  
  if(is.null(d)) d <- rep(1, length(y)) #	default, no normalization
  
  model <- try(match.arg(model, c('E', 'V'), several.ok=FALSE), silent=TRUE)
  
  # stop if model not recognizable.
  
  if(class(model)=='try-error') stop('Only model E or V can be specified!\n')
  
  percent0 <- mean(y==0, na.rm=TRUE)
  
  res <- rep(NA, 7) # mu1, mu2, sigma1, sigma2, pi1, logLik, BIC
  
  names(res) <- c('mu1', 'mu2', 'sigma1', 'sigma2', 'pi1', 'logLik', 'BIC') # 1:7
  
  if(percent0<=zeroPercentThr) { # no severe 0-inflation, fit 2-comp
    
    Dat <- (log(y+eps)-log(d))/log(base) # transform data;
    
    # change of base formula: log_b(x)=log_e(x)/log_e(b)
    
    mc <- try(Mclust(Dat, G = 2, modelNames = model), silent=TRUE)
    
    res[1:7] <- extractMclustPar(mc, modelName=model, logLikToLN=logLikToLN, dat=Dat)
    
  } else { # severe 0-inflation. 1-comp model. need control mu and sd for RPKM
    
    nonzero <- y[y!=0]
    
    #if(mean(nonzero)<1) nonzero <- nonzero+max(c(1, eps)) ## when count<1, i.e. in RPKM, usually 0.00001, force it at least 1 to avoid FP
    
    #### median is more robust: not susceptible to single outlier expression
    
    if(median(nonzero)<1) nonzero <- nonzero+max(c(1, eps)) ## when count<1, i.e. in RPKM, usually 0.00001, force it at least 1 to avoid FP
    
    Dat <- (log(nonzero)-log(d[y!=0]))/log(base) ###### no need to add eps since no 0 now!
    
    #mc_1comp <- try(Mclust(Dat, G = 1), silent=TRUE)
    
    #Modified to MLE such that c(00000, 11111) wouldn't fail
    
    mu2 <- mean(Dat)
    
    sigma2 <- max(sd(Dat), 1e-2) ## in case sd=0.0004
    
    pi1 <- percent0 # 
    
    res[1:5] <- c(0, mu2, 0, sigma2, pi1) # mu1=0, sigma1=0, logLik=BIC=NA
    
    res[6] <- logLik0inflatedLN(y, mu=mu2, sigma=sigma2, pi1=pi1,
                                
                                logLikToLN=logLikToLN)
    
    res[7] <- getBIC(logLik=res[6], nPar=3, nObs=length(y)) 
    
  }
  
  res
  
}



# calculate log-lik for 0-inflated normal or log normal model. log normal model is achieved by setting logLikToLN=TRUE

logLik0inflatedLN <- function(y, mu, sigma, pi1, logLikToLN=TRUE) {
  
  if(logLikToLN==TRUE) {
    
    res <- sum(log(pi1*as.numeric(y==0)+(1-pi1)*dlnorm(y, meanlog=mu, sdlog=sigma)))
    
  } else {
    
    res <- sum(log(pi1*as.numeric(y==0)+(1-pi1)*dnorm(y, mean=mu, sd=sigma)))
    
  }
  
  res
  
}



# a function to extract parameters estimated by mclust with G=2, compatible to try-error of Mclust() function

## modified on 07/22/12: add option logLikToLN and dat. by default, it's disabled so that it is compatible with previous version

## add option logLikToLN on 7/22/2012 so that logLik and BIC can be defined on the exponent LN scale. 

## This is needed for pseudo-LN mixture fitting		

# dat: the data exactly used to fit the normal mixture. it is defined on the log scale and hence can be negative		

extractMclustPar <- function(mc, modelName='E', logLikToLN=FALSE, dat=NA) {
  
  res <- rep(NA, 7)
  
  nPar <- ifelse(modelName=='V', 5, 4) # number of parameters
  
  if(class(mc)!="try-error") {
    
    # extract mu1, mu2
    
    res[1:2] <- mc$parameters$mean
    
    # extract sigma1, sigma2
    
    temp <- sqrt(mc$parameters$variance$sigmasq)
    
    if(length(temp)==1) { # E model, 1 sigma
      
      res[3:4] <- rep(temp, 2)
      
    } else {
      
      res[3:4] <- temp
      
    }
    
    # extract p1
    
    res[5] <- mc$parameters$pro[1]
    
    if(logLikToLN) { # modify logLik and BIC to the exponent LN scale.
      
      # extract logLik
      
      res[6] <- logLikLN(y=exp(dat), theta=res[1:5]) # data is now transformed to the exponent scale
      
      # extract BIC
      
      res[7] <- getBIC(logLik=res[6], nPar=nPar, nObs=length(dat)) ## confirmed by BIC function in R
      
    } else {
      
      # extract logLik
      
      res[6] <- mc$loglik
      
      # extract BIC
      
      res[7] <- ifelse(modelName=="V", -bic(modelName="V",
                                            
                                            loglik=mc$loglik, n=mc$n, d=1, G=2), 
                       
                       -bic(modelName="E", loglik=mc$loglik, n=mc$n, d=1, G=2))
      
    }								 
    
  }
  
  res
  
}







## logLik for LN data. used to compute BIC of LN since we fit normalMixture whose logLik is defined on the log-transformed data

# theta: a vector. if length is 2, theta=c(meanlog , sdlog),logLik of 1-comp data is obtained. 

#		 if length=5, logLik of 2-comp data is obtained. theta=c(meanlog1, meanlog2, sdlog1, sdlog2, pi1)

## bug fixed at 07/27/12:  dlnorm(y, meanlog=theta[2], sdlog=theta[5], log=FALSE) --->  dlnorm(y, meanlog=theta[2], sdlog=theta[4], log=FALSE)

logLikLN <- function(y, theta) {
  
  if(length(theta)==2) { # 1-comp logLik
    
    res <- sum(dlnorm(y, meanlog=theta[1], sdlog=theta[2], log=TRUE))
    
  } else { # 2-comp logLik
    
    res <- sum(log(dlnorm(y, meanlog=theta[1], sdlog=theta[3], log=FALSE)*theta[5]+
                     
                     dlnorm(y, meanlog=theta[2], sdlog=theta[4],
                            
                            log=FALSE)*(1-theta[5])))
    
  }
  
  res
  
}









############################ NL mix ##########################

# for completeness, we also incorporate normal mixture for microarray data. We use NL to denote normal.

fitNL <- function(y, d=NULL, model='E') {
  
  if(is.null(d)) d <- rep(1, length(y)) #	default, no normalization
  
  model <- try(match.arg(model, c('E', 'V'), several.ok=FALSE), silent=TRUE)
  
  # stop if model not recognizable.
  
  if(class(model)=='try-error') stop('Only model E or V can be specified!\n')
  
  res <- rep(NA, 7) # mu1, mu2, sigma1, sigma2, pi1, logLik, BIC
  
  names(res) <- c('mu1', 'mu2', 'sigma1', 'sigma2', 'pi1', 'logLik', 'BIC') # 1:7
  
  Dat <- y/d # normalization
  
  mc <- try(Mclust(Dat, G = 2, modelNames = model), silent=TRUE)
  
  res[1:7] <- extractMclustPar(mc, modelName=model, logLikToLN=FALSE, dat=Dat)
  
  res
  
}



# given mu1, mu2, phi1, phi2, p1, logLik (optional), BIC (optional), transform phi to sigma

# transformNBfit() is its special case

transformPhiToSigmaInMixFit <- function(est, dist='NB') {
  
  sigma1 <- sqrt(muPhiToVar(mu=est[1], phi=est[3], dist=dist))
  
  sigma2 <- sqrt(muPhiToVar(mu=est[2], phi=est[4], dist=dist))
  
  res <- est
  
  res[3:4] <- c(sigma1, sigma2)
  
  res
  
}



# phi2Var: when mu1, mu2, phi1, phi2, pi1 is specified, we need to set this so that var can be computed

# modified on 08/03/12: deal with 0-inflation: input=c(0, mu2, 0, sigma2, pi1, ...), delta=mu2/sigma2 

# modified on 08/03/12: !any(is.na(est))----->!any(is.na(est[1:5])). necessary since 0-inflation always get logLik and BIC as NA

# modified on 08/04/12: for 0-inflation case, general BI2 formula also works. instead sqrt((1-pi))*mu2/sigma2 is wrong: no detection of such genes!

parToBI <- function(est, phi2Var=FALSE, dist='NB') {
  
  res <- c(NA, NA)
  
  names(res) <- c('delta', 'BI')
  
  if(!any(is.na(est[1:5]))) {
    
    if(phi2Var==TRUE) {
      
      est[3] <- sqrt(muPhiToVar(est[1], est[3], dist=dist))
      
      est[4] <- sqrt(muPhiToVar(est[2], est[4], dist=dist))
      
    }
    
    #est[5] <- restrictPi1(est[5]) # deal with numeric issue from optim fit (GP). i.e. pi= -2.775558e-17
    
    # BI2 automatically deal with 0-inflation 
    
    res[1] <- abs(diff(est[1:2]))/sqrt((1-est[5])*est[3]^2+est[5]*est[4]^2)
    
    res[2] <- sqrt(est[5]*(1-est[5]))*res[1]
    
  }
  
  res
  
}



########################### SIBER ###########################

# the final product presented to the user

# parameters (zeroPercentThr=0.1, base=10, eps=10) are only relevant to LN model.

SIBER <- function(y, d=NULL, model=c('LN', 'NB', 'GP', 'NL'),
                  
                  zeroPercentThr=0.2, base=exp(1), eps=10) {
  
  model <- try(match.arg(model, c('LN', 'NB', 'GP', 'NL'),
                         
                         several.ok=FALSE), silent=TRUE)
  
  # stop if model not recognizable.
  
  if(class(model)=='try-error') {
    
    stop('Only model LN, NB, GP or NL can be specified!\n')
    
  }
  
  res <- rep(NA, 7) 
  
  names(res) <- c('mu1', 'mu2', 'sigma1', 'sigma2', 'pi1', 'delta', 'BI') 
  
  if(model=='LN') {
    
    fit <- fitLN(y, d=d, model='E', zeroPercentThr=zeroPercentThr,
                 
                 base=base, eps=eps)[1:5]
    
  } else if(model=='NB') {
    
    fitPhiScale <- fitNB(y, d=d, model='E')[1:5]
    
    fit <- transformPhiToSigmaInMixFit(fitPhiScale, dist='NB')
    
  } else if(model=='GP') {
    
    fitPhiScale <- fitGP(y, d=d, model='E')[1:5]
    
    fit <- transformPhiToSigmaInMixFit(fitPhiScale, dist='GP')
    
  } else {
    
    fit <- fitNL(y, d=d, model='E')[1:5]
    
  }
  
  BIinfo <- parToBI(fit, phi2Var=FALSE)
  
  res[1:5] <- fit
  
  res[6:7] <- BIinfo
  
  res
  
}



