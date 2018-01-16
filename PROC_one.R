

Model_single_ISE_cal_only<-function()
{
  
  ### Calculate x from log x
  for (i in 1:N) {
    x[i] <- pow(10, log10x[i])
  }
  
  ### Calibration data
  # Each ISE emf response is based on x and a, b, c, and tau for the particular ISE
  for (i in 1:N) {
    emf[i] ~ dnorm(mu.emf[i], Tau[i])
    mu.emf[i] <- a + b *log(x[i] + c)/log(10)
    Tau[i] <- tau
  }
  
  ### Priors for the ISE model
  a ~ dnorm(0, 0.000001)
  
  ###########################
  # ISE-specific priors     #
  b ~ dnorm(mu.b, 0.01)     #
  cstar ~ dunif(0.1, 0.7)   #
  sigma ~ dunif(0,50)       #
  ###########################
  
  # Logical nodes
  c <- pow(cstar, 10)	
  logsigma <- log(sigma)
  tau <- 1/(sigma*sigma)
  
}

Model_single_ISE_second_step<-function()
{
  
  ### Calculate x from log x
  for (i in 1:N) {
    x[i] <- pow(10, log10x[i])
  }
  
  ### Calibration data
  # Each ISE emf response is based on x and a, b, c, and tau for the particular ISE
  for (i in 1:N) {
    emf[i] ~ dnorm(mu.emf[i], Tau[i])
    mu.emf[i] <- a + b *log(x[i] + c)/log(10)
    Tau[i] <- tau
  }
  
  ### Priors for the ISE model
  a ~ dnorm(0, 0.000001)
  
  ###########################
  # ISE-specific priors     #
  b ~ dnorm(mu.b, 0.01)     #
  cstar ~ dunif(0.1, 0.7)   #
  sigma ~ dunif(0,50)       #
  ###########################
  
  # Logical nodes
  c <- pow(cstar, 10)	
  logsigma <- log(sigma)
  tau <- 1/(sigma*sigma)
  
  # Separate ISE parameter estimation from sample calibration
  a.cut <- cut(a)
  b.cut <- cut(b)
  c.cut <- cut(c)
  tau.cut <- cut(tau)
  
  ###  Experimental Samples
  for (i in 1:M) {
    ##################################
    # Prior on log x                 #
    log10x.exp[i] ~ dunif(-9, -1)    #
    ################################## 
    x.exp[i] <- pow(10, log10x.exp[i])
  }
  
  for (i in 1:M) {	
    emf.exp[i] ~ dnorm(mu.emf.exp[i], Tau.exp[i])
    mu.emf.exp[i] <- a.cut + b.cut *log(x.exp[i] + c.cut)/log(10)
    Tau.exp[i] <- tau.cut
    SD.exp[i] <- 1/sqrt(Tau.exp[i])				
  }
}



Model_multi_ISEs<-function()
{
  ### Calculate x from log x
  for (i in 1:N) {
    x[i] <- pow(10, log10x[i])
  }
  
  ### Calibration data
  # Each ISE emf response is based on x and a, b, c, and tau for the particular ISE
  for (i in 1:N) {
    emf[i] ~ dnorm(mu.emf[i], Tau[i])
    mu.emf[i] <- a[ISEID[i]] + b[ISEID[i]] *log(x[i] + c[ISEID[i]])/log(10)
    Tau[i] <- tau[ISEID[i]]
  }
  
  ### Priors for each ISE
  for (j in 1:R) {
    a[j] ~ dnorm(0, 0.000001)
    #############################
    # ISE-specific priors       #
    b[j] ~ dnorm(mu.b, 0.01)    #
    cstar[j] ~ dunif(0.1, 0.8)  #
    sigma[j] ~ dunif(0,50)      #
    #############################
    
    # Logical nodes
    c[j] <- pow(cstar[j], 10)
    logsigma[j] <- log(sigma[j])
    tau[j] <- 1/(sigma[j]*sigma[j])
    
    # Separate ISE parameter estimation from sample calibration
    a.cut[j] <- cut(a[j])
    b.cut[j] <- cut(b[j])
    c.cut[j] <- cut(c[j])
    tau.cut[j] <- cut(tau[j])
  }
  
  ### Experimental samples
  for (i in 1:M) {
    ##################################
    # Prior on log x                 #
    log10x.exp[i] ~ dunif(-12, -1)  #
    ################################## 
    x.exp[i] <- pow(10, log10x.exp[i])
  }
  
  for (i in 1:M.obs) {		
    emf.exp[i] ~ dnorm(mu.emf.exp[i], Tau.exp[i])
    mu.emf.exp[i] <- a.cut[ISEID.exp[i]] + b.cut[ISEID.exp[i]] *log(x.exp[xID.exp[i]] + c.cut[ISEID.exp[i]])/log(10)
    Tau.exp[i] <- tau.cut[ISEID.exp[i]]
    SD.exp[i] <- 1/sqrt(Tau.exp[i])			
  }
}

Model_multi_ISEs_cal_only<-function()
{
  ### Calculate x from log x
  for (i in 1:N) {
    x[i] <- pow(10, log10x[i])
  }
  
  ### Calibration data
  # Each ISE emf response is based on x and a, b, c, and tau for the particular ISE
  for (i in 1:N) {
    emf[i] ~ dnorm(mu.emf[i], Tau[i])
    mu.emf[i] <- a[ISEID[i]] + b[ISEID[i]] *log(x[i] + c[ISEID[i]])/log(10)
    Tau[i] <- tau[ISEID[i]]
  }
  
  ### Priors for each ISE
  for (j in 1:R) {
    a[j] ~ dnorm(0, 0.000001)
    #############################
    # ISE-specific priors       #
    b[j] ~ dnorm(mu.b, 0.01)    #
    cstar[j] ~ dunif(0.1, 0.8)  #
    sigma[j] ~ dunif(0,50)      #
    #############################
    
    # Logical nodes
    c[j] <- pow(cstar[j], 10)
    logsigma[j] <- log(sigma[j])
    tau[j] <- 1/(sigma[j]*sigma[j])
    
  }
  
}

Model_single_ISE_second_step<-function()
{
  
  ### Calculate x from log x
  for (i in 1:N) {
    x[i] <- pow(10, log10x[i])
  }
  
  ### Calibration data
  # Each ISE emf response is based on x and a, b, c, and tau for the particular ISE
  for (i in 1:N) {
    emf[i] ~ dnorm(mu.emf[i], Tau[i])
    mu.emf[i] <- a + b *log(x[i] + c)/log(10)
    Tau[i] <- tau
  }
  
  ### Priors for the ISE model
  a ~ dnorm(0, 0.000001)
  
  ###########################
  # ISE-specific priors     #
  b ~ dnorm(mu.b, 0.01)     #
  cstar ~ dunif(0.1, 0.7)   #
  sigma ~ dunif(0,50)       #
  ###########################
  
  # Logical nodes
  c <- pow(cstar, 10)	
  logsigma <- log(sigma)
  tau <- 1/(sigma*sigma)
  
  # Separate ISE parameter estimation from sample calibration
  a.cut <- cut(a)
  b.cut <- cut(b)
  c.cut <- cut(c)
  tau.cut <- cut(tau)
  
  ###  Experimental Samples
  for (i in 1:M) {
    ##################################
    # Prior on log x                 #
    log10x.exp[i] ~ dunif(-12, -1)    #
    ################################## 
    x.exp[i] <- pow(10, log10x.exp[i])
  }
  
  for (i in 1:M) {	
    emf.exp[i] ~ dnorm(mu.emf.exp[i], Tau.exp[i])
    mu.emf.exp[i] <- a.cut + b.cut *log(x.exp[i] + c.cut)/log(10)
    Tau.exp[i] <- tau.cut
    SD.exp[i] <- 1/sqrt(Tau.exp[i])				
  }
}

gen_inits_single_PROC_one<- function(data=data.subset, a.init=NA, b.init=NA, cstar.init=NA, 
                                     logc.limits = c(-8.9, -1.1), sigma.upper=NA, stdadd = F, offset = 1, calibration.only = F)
{
  x <- 10^data$log10x
  order.x = order(x)
  x = sort(x)
  log10x = log10(x)
  emf <- data$emf[order.x]
  
  # Assignments
  n = length(x)
  b = (emf[n] - emf[n-offset])/(log10x[n] - log10x[n-offset])	
  a = emf[n] - b*log10x[n]
  C = min(max(10^logc.limits[1], 10^((emf[1]-a)/b) - x[1]), 10^logc.limits[2])
  cstar = 10^(log10(C)/10)
  
  # over-ride if initial values are provided
  if(!is.na(a.init)) { a = a.init  }
  if(!is.na(b.init)) { b = b.init  }
  if(!is.na(cstar.init)) { cstar = cstar.init  }
  
  Preds = a + b*log10(x + C)
  Resids = emf - Preds
  logsigma.temp = log(sqrt( sum(Resids^2)/(length(Resids) - 3) ))
  logsigma <- logsigma.temp
  sigma = exp(logsigma)
  # If sigma is too large, it may be out of bounds.
  # Use sigma.upper to place an upper limit on sigma
  if(!is.na(sigma.upper)) {
    sigma = min(sigma.upper, exp(logsigma))
  }
  
  sigma.sq = sigma^2
  tau = 1/sigma.sq
  
  # Calculate LOD using s/n = 3
  mu0 <- a + b*log10(cstar^10)
  mu0.3sig <- mu0 + 3*sigma
  LOD <- 10^( (mu0.3sig - a)/b ) - cstar^10
  
  
  
  init.stretch = 0.01
  inits <- list(a=a*runif(1, 1 - init.stretch, 1 + init.stretch),
                b=b*runif(1, 1 - init.stretch, 1 + init.stretch),
                cstar=cstar*runif(1, 1 - init.stretch, 1 + init.stretch),
                sigma=sigma*runif(1, 1 - init.stretch, 1 + init.stretch))
  
  
  return(inits)
  
}


gen_inits_multiple_PROC_one <- function(data=data.subset, a.init= NA, b.init= NA, cstar.init= NA, offset = 1,
                                        logc.limits = c(-8.9, -1.9), sigma.upper = 14, stdadd = F, calibration.only=F) {
  ###
  # Similar to gen.inits.single, but for multiple ISEs
  # If initial values are specified, it should be vectors with length equal to the number of ISEs
  ###
  # Initialise vectors
  a <- rep(NA, data$R)
  b <- rep(NA, data$R)
  cstar <- rep(NA, data$R)
  logsigma <- rep(NA, data$R)
  sigma <- rep(NA, data$R)
  sigma.sq <- rep(NA, data$R)
  tau <- rep(NA, data$R)
  
  # Generate initial values for calibration data
  for (i in 1:data$R) {
    
    x <- 10^data$log10x[data$ISEID==i]
    order.x = order(x) 			#new
    x = sort(x)  #new
    log10x = log10(x)
    emf <- data$emf[data$ISEID==i][order.x] #new
    
    # Assignments
    n = length(x)
    b[i] = (emf[n] - emf[n-offset])/(log10x[n] - log10x[n-offset])	
    a[i] = emf[n] - b[i]*log10x[n]
    C = min(max(10^logc.limits[1], 10^((emf[1]-a[i])/b[i]) - x[1]), 10^logc.limits[2])
    cstar[i] = 10^(log10(C)/10)
    
    # over-ride if initial values are provided
    if(!is.na(a.init[1])) { a[i] = a.init[i]  }
    if(!is.na(b.init[1])) { b[i] = b.init[i]  }
    if(!is.na(cstar.init[1])) { cstar[i] = cstar.init[i]  }
    
    Preds = a[i] + b[i]*log10(x + C)
    Resids = emf - Preds
    logsigma.temp = log(sqrt( sum(Resids^2)/(length(Resids) - 3) ))
    
    logsigma[i] <- logsigma.temp
    sigma[i] = exp(logsigma[i])
    
    # If sigma is too large, it may be out of bounds.
    # Use sigma.upper to place an upper limit on sigma
    if(!is.na(sigma.upper)) {
      sigma[i] = min(sigma.upper, exp(logsigma[i]))
    }
    sigma.sq[i] = sigma[i]^2
    tau[i] = 1/sigma.sq[i]
  }
  
  # Calculate LOD using s/n = 3
  mu0 <- a + b*log10(cstar^10)
  mu0.3sig <- mu0 + 3*sigma*sign(b)
  LOD <- 10^( (mu0.3sig - a)/b ) - cstar^10
  
  
  
  #INITS <- list(a=a, b=b, cstar=cstar, logsigma=logsigma)
  init.stretch = 0.01
  inits <- list(a=a*runif(data$R, 1 - init.stretch, 1 + init.stretch),
                b=b*runif(data$R, 1 - init.stretch, 1 + init.stretch),
                cstar=cstar*runif(data$R, 1 - init.stretch, 1 + init.stretch),
                sigma=sigma*runif(data$R, 1 - init.stretch, 1 + init.stretch))
  
  return(inits)
}



library(R2OpenBUGS)
PROC_one = function(filename.calibration, params,Z) {
  ###
  
  # This function takes Calibration data, the parameters of interest, and Z (valance)
  # The output will be an R list with 4 elements. Eachelement of the list is an M \times R matrix containing elements for the
  #	specified model parameters. Columns of each matrix correspond to different ISEs and each row of the matrix correspond to 
  #an iteration of the Bayesian MCMC output for model parameters.
  #
  #Input:		filename.calibration, The location of calibration data
  #			params, a vector of model parameters that needed to be estimated (Here it's \beta_0,\beta_1,\beta_2 and \sigma)
  #			Z, the ionic valence ,e.g.  Z = 2
  ###
  
  # Load calibration data
  data.calib = read.delim(filename.calibration,
                          header = TRUE, sep = "\t", quote="\"", dec=".", fill = TRUE, comment.char="")
  
  # Format data from the calibration file
  N = nrow(data.calib)       # Number of calibration data
  R = max(data.calib$ISEID)  # Number of ISEs in the calibration data
  ISEID = data.calib$ISEID   # A vector of the ISEs
  log10x = data.calib$log10x # A vector of concentration values
  emf = data.calib$emf 		 # A vector of signal values
  mu.b = 1000*log(10)*8.31433*(21 + 273.15)/(96487*Z) # The theoritical slope
  if (R == 1) { single.ISE = T} # If true, use single ISE 
  if (R > 1) { single.ISE = F } # Otherwise use multiple ISEs
  # Format the data as a list for initial values and to be used in OpenBUGS 
  data.subset = list(N=N, R=R, ISEID=ISEID, log10x=log10x, emf=emf, mu.b=mu.b)	
  # If model.name is not specified, revert to appropriate defaults
  
  if ( R == 1) {
    model.name = Model_single_ISE_cal_only
  }
  if ( R > 1 ){
    model.name = Model_multi_ISEs_cal_only
  }
  
  
  # Initial values
  if (single.ISE==T) {
    init <- rep(list(NA), 1)
    for (i in 1:1) {
      # Generate initial values
      init[[i]]<-gen_inits_single_PROC_one(data=data.subset, a.init= NA, b.init= NA, cstar.init= NA, offset = 1,
                                           logc.limits = c(-8.9, -1.9), sigma.upper = 14, stdadd = F, calibration.only=F)
    }}
  
  if (single.ISE==F) {
    # Generate initial values
    init <- rep(list(NA), 1)
    for (i in 1:1) {
      init[[i]]<-gen_inits_multiple_PROC_one(data=data.subset, a.init= NA, b.init= NA, cstar.init= NA, offset = 1,
                                             logc.limits = c(-8.9, -1.9), sigma.upper = 14, stdadd = F, calibration.only=F)
    }
  }
  
  #return(init)}
  
  # Call BUGS
  require(R2OpenBUGS)
  my.BUGS<-bugs(data=data.subset, inits=init, parameters.to.save=params,
                n.iter=100000, n.chains=1, n.burnin=98000, n.thin=1,
                model.file=model.name, debug=FALSE, codaPkg=FALSE)
  
  
  if(single.ISE==T) {
    # Other vars
    ahat.coda = as.vector(my.BUGS$sims.array[,,1])
    coda.length = length(ahat.coda)
    bhat.coda = as.vector(my.BUGS$sims.array[,,2])
    cstarhat.coda = as.vector(my.BUGS$sims.array[,,4])
    chat.coda = cstarhat.coda^10
    sigmahat.coda = as.vector(my.BUGS$sims.array[,,5])
    
    #result=matrix(c(ahat.coda=ahat.coda,bhat.coda=bhat.coda,cstarhat.coda=cstarhat.coda,sigmahat.coda=sigmahat.coda),nrow = 2000,byrow = FALSE)
    # result=list(ahat.coda=ahat.coda,bhat.coda=bhat.coda,cstarhat.coda=cstarhat.coda,sigmahat.coda=sigmahat.coda)
    result=list(matrix(c(ahat.coda=ahat.coda)),matrix(c(bhat.coda=bhat.coda)),matrix(c(cstarhat.coda=cstarhat.coda)),matrix(c(sigmahat.coda=sigmahat.coda)))
    
  }
  
  
  if (single.ISE == F) {
    # Extract paramaters information needed
    coda.length = length(as.vector(my.BUGS$sims.array[,,1]))
    ahat.coda = matrix(NA, nrow = coda.length, ncol = data.subset$R)
    bhat.coda <- cstarhat.coda <- chat.coda <- sigmahat.coda <- LOD.SN3.coda <- ahat.coda
    
    ahat = rep(NA, data.subset$R)
    ahat.lcl = rep(NA, data.subset$R)
    ahat.ucl = rep(NA, data.subset$R)
    bhat = rep(NA, data.subset$R)
    bhat.lcl = rep(NA, data.subset$R)
    bhat.ucl = rep(NA, data.subset$R)
    chat = rep(NA, data.subset$R)
    chat.lcl = rep(NA, data.subset$R)
    chat.ucl = rep(NA, data.subset$R)
    cstarhat = rep(NA, data.subset$R)
    cstarhat.lcl = rep(NA, data.subset$R)
    cstarhat.ucl = rep(NA, data.subset$R)
    sigmahat = rep(NA, data.subset$R)
    sigmahat.lcl = rep(NA, data.subset$R)
    sigmahat.ucl = rep(NA, data.subset$R)
    
    
    for (j in 1:data.subset$R) {
      ahat.coda[,j] = as.vector(my.BUGS$sims.array[,,j])
      bhat.coda[,j] = as.vector(my.BUGS$sims.array[,,data.subset$R + j])
      cstarhat.coda[,j] = as.vector(my.BUGS$sims.array[,,3*data.subset$R + j])
      sigmahat.coda[,j] = as.vector(my.BUGS$sims.array[,,4*data.subset$R + j])
      chat.coda[,j] = cstarhat.coda[,j]^10
      
      
    }
    
    result=list(ahat.coda=ahat.coda,bhat.coda=bhat.coda,cstarhat.coda=cstarhat.coda,sigmahat.coda=sigmahat.coda)
    # result=list(finalresults=finalresults)
    
  }
  
  return(result)
  
}

parameters_matrix_ISE16<-PROC_one("/home/balsaedi/new_method/6obs_simulation/6cal_sim_4ISEs.txt",params=c("a", "b", "c", "cstar", "sigma"),Z=2)

parameters_matrix_ISE16[[1]][,1]
