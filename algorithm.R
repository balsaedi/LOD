
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

PROC_two=function(E,results,filename.calibration,Z)
{
  ###
  # This function takes Calibration data, and the output of PROC_one
  # The output will be a vector of detection threshold in the log concentration. Eachelement of the vector correspond to a combination 
  # of parameter value for a sample. 
  #
  #Input:	E, indictor correspond to a row of parameter values	
  # 		filename.calibration, The location of calibration data
  #			results, PROC_one output   (model parameters are \beta_0,\beta_1,\beta_2 and \sigma)
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
  
  rejection_region=rep(list(NA),100)
  # Create an empty matrix for storage
  emf.2<-replicate(100, matrix(NA, nrow = 200, ncol = R), simplify = F)
  
  emf.1<-matrix(NA, nrow = 200, ncol = R)
  Sample_1<-matrix(NA, nrow = 200, ncol = R)
  ISE_1<-matrix(NA, nrow = 200, ncol = R)
  
  for(j in 1:R)
  {
    # Create new experimental response that correspond to x=0 
    emf1=results[[1]][E,j]+results[[2]][E,j]*log10(results[[3]][E,j]^10)+rnorm(200,0,results[[4]][E,j])
    emf.2[[E]][,j]<-emf1
    SampleID_2<-1:length(emf1)
    Sample_1[,j]<-SampleID_2
    ISEEX_2<- rep(j:j, each=200)
    ISE_1[,j]<- ISEEX_2
  }
  
  #Reshape the data into the appropriate format
  emf.exp<-as.vector(emf.2[[E]])
  xID.exp <- as.vector(Sample_1)
  ISEID.exp <- as.vector(ISE_1)
  M = max(xID.exp)
  M.obs = length(emf.exp)
  
  # Format the data as a list for initial values and to be used in OpenBUGS
  data.subset2= list(N=N, R=R, ISEID=ISEID, log10x=log10x, emf=emf, mu.b=mu.b,
                     ISEID.exp=ISEID.exp, xID.exp=xID.exp, emf.exp=emf.exp,M=M,M.obs=M.obs)
  #if single ISE, use this model}
  if ( R == 1) {
    model.name = Model_single_ISE_second_step
  }
  # If multiple ISEs use this model
  if ( R > 1 ){
    model.name = Model_multi_ISEs
  }
  
  
  # Initial values for the single ISE case
  if (single.ISE==T) {
    init <- rep(list(NA), 1)
    for (i in 1:1) {
      # Generate initial values
      init[[i]]<-gen.inits.single1(data=data.subset2, a.init= NA, b.init= NA, cstar.init= NA, offset = 1,
                                   logc.limits = c(-8.9, -1.9), sigma.upper = 14, stdadd = F, calibration.only=F)
    }}
  # Initial values for the multiple ISEs case
  if (single.ISE==F) {
    # Generate initial values
    init <- rep(list(NA), 1)
    for (i in 1:1) {
      init[[i]]<-gen.inits.multiple(data=data.subset2, a.init= NA, b.init= NA, cstar.init= NA, offset = 1,
                                    logc.limits = c(-8.9, -1.9), sigma.upper = 14, stdadd = F, calibration.only=F)
    }
  }
  
  # call OpenBUGS
  my.BUGS = bugs(data=data.subset2, inits=init, parameters.to.save=c("a", "b", "c", "cstar", "sigma","log10x.exp"),OpenBUGS.pgm="/usr/bin/OpenBUGS",
                 n.iter=2000, n.chains=1, n.burnin=1000, n.thin=1,
                 model.file=model.name, debug=FALSE, codaPkg=FALSE)
  # format the results for single ISE
  if(single.ISE==T) {
    # Calculate the median for xhat
    log10x.exp <- matrix(NA, nrow=data.subset2$M, ncol=3)
    for (i in 1:data.subset2$M) {
      tmp <- as.vector(my.BUGS$sims.array[,,i+5])
      quantile.tmp <- quantile(tmp, c(0.025, 0.5, 0.975, 0.25, 0.75))
      log10x.exp[i,] <- quantile.tmp[c(2,1,3)]
    }
    # This calculate detection threshold for each E 
    rejection_region[[E]]<-quantile(log10x.exp[,1],0.95)
    
  }
  # format the results for multiple ISEs 
  if (single.ISE == F) {
    # Calculate the median for xhat
    log10x.exp <- matrix(NA, nrow=data.subset2$M, ncol=3)
    for (i in 1:data.subset2$M) {
      tmp = as.vector(my.BUGS$sims.array[,,i+5*data.subset2$R])
      quantile.tmp <- quantile(tmp, c(0.025, 0.5, 0.975, 0.25, 0.75))
      log10x.exp[i,] <- quantile.tmp[c(2,1,3)]
      
    }
    # This calculate detection threshold for each E 
    rejection_region[[E]]<-quantile(log10x.exp[,1],0.95)
    
    
  }
  return(rejection_region[[E]])
}



PROC_three = function(filename.calibration,Z,fun.ahat,fun.bhat,fun.cstarhat,fun.sigmahat,newx,rejection) {
  ###
  
  # This function takes Calibration data, the output of PROC_one, the output fro PROC_two (DT) and 
  # an x value provided by PROC_four.
  # The output will be a vector of detection limit (LOD) in the log concentration. Eachelement of the 
  # result is calculated conditionally on parameter values (PROC_one) and detection threshold (PROC_two)
  #
  #Input:	 filename.calibration, The location of calibration data
  #			   fun.ahat, element one of the list from PROC_one (\beta_0)
  #        fun.bhat, element two of the list from PROC_one (\beta_1)
  #        fun.cstarhat, element three of the list from PROC_one (\beta_2)
  #        fun.sigmahat, element four of the list from PROC_one (\sigma)
  #			   Z, the ionic valence ,e.g.  Z = 2
  #        newx, is an x value provided by PROC_four
  #        rejection, is a vector of detection threshold values (PROC_two)
  ###
  
  # Load calibration data
  data.calib = read.delim(filename.calibration,
                          header = TRUE, sep = "\t", quote="\"", dec=".", fill = TRUE, comment.char="")
  # Format data from the calibration file
  N = nrow(data.calib)
  R = max(data.calib$ISEID)
  ISEID = data.calib$ISEID
  log10x = data.calib$log10x
  emf = data.calib$emf
  mu.b = 1000*log(10)*8.31433*(21 + 273.15)/(96487*Z)
  if (R == 1) { single.ISE = T}
  if (R > 1) { single.ISE = F }
  
  #if single ISE, use this model}
  if ( R == 1) {
    model.name = Model_single_ISE_second_step
  }
  
  # If multiple ISEs use this model
  if ( R > 1 ){
    model.name = Model_multi_ISEs
  }
  
  finalresults <- rep(NA, length(newx))
  emf.new=matrix(NA,nrow = 200,ncol = R)
  samp=matrix(NA,nrow = 200,ncol = R)
  ISE_1=matrix(NA,nrow = 200,ncol = R)
  
  
  for (i in 1:length(newx)) {
    for (m in 1:R) {
      # Create new experimental response that correspond to an x value (newx)
      newemf=fun.ahat[m]+fun.bhat[m]*log10(newx[i]+fun.cstarhat[m]^10)+ rnorm(200,0,fun.sigmahat[m])
      emf.new[,m]<-newemf
      SampleID_3<-1:length(newemf)
      samp[,m]<-SampleID_3
      ISEIDexp_3<-rep(m:m, each=200)
      ISE_1[,m]<-ISEIDexp_3
    }
    # Reshape the data into the appropriate format
    emf.exp=as.vector(emf.new)
    xID.exp = as.vector(samp)
    M = max(xID.exp)
    ISEID.exp = as.vector(ISE_1)
    M.obs = length(emf.exp)
    
    
    data.subset3=list(N=N, R=R, ISEID=ISEID, log10x=log10x, emf=emf, mu.b=mu.b,
                      ISEID.exp=ISEID.exp, xID.exp=xID.exp, emf.exp=emf.exp,M=M,M.obs=M.obs)
    
    # Initial values for the single ISE case
    if (single.ISE==T) {
      init <- rep(list(NA), 1)
      for (i in 1:1) {
        # Generate initial values
        init[[i]]<-gen.inits.single1(data=data.subset3, a.init= NA, b.init= NA, cstar.init= NA, offset = 1,
                                     logc.limits = c(-8.9, -1.9), sigma.upper = 14, stdadd = F, calibration.only=F)
      }}
    
    # Initial values for the multiple ISEs case
    if (single.ISE==F) {
      # Generate initial values
      init <- rep(list(NA), 1)
      for (i in 1:1) {
        init[[i]]<-gen.inits.multiple(data=data.subset3, a.init= NA, b.init= NA, cstar.init= NA, offset = 1,
                                      logc.limits = c(-8.9, -1.9), sigma.upper = 14, stdadd = F, calibration.only=F)
      }
    }
    # call OpenBUGS
    my.BUGS = bugs(data=data.subset3, inits=init, parameters.to.save=c("a", "b", "c", "cstar", "sigma","log10x.exp"),OpenBUGS.pgm="/usr/bin/OpenBUGS",
                   n.iter=2000, n.chains=1, n.burnin=1000, n.thin=1,
                   model.file=Model_multi_ISEs, debug=FALSE, codaPkg=FALSE)
    
    # format the results for single ISE
    if(single.ISE==T) {
      # Calculate quantiles for x
      log10x.exp <- matrix(NA, nrow=data.subset3$M, ncol=1)
      for (i in 1:data.subset3$M) {
        tmp <- as.vector(my.BUGS$sims.array[,,i+5])
        quantile.tmp <- quantile(tmp, c(0.025, 0.5, 0.975, 0.25, 0.75))
        log10x.exp[i,] <- quantile.tmp[c(2)]
      }
    }
    # format the results for multiple ISEs 
    if (single.ISE == F) {
      # Calculate quantiles for x
      log10x.exp <- matrix(NA, nrow=data.subset3$M, ncol=1)
      for (t in 1:data.subset3$M) {
        tmp = as.vector(my.BUGS$sims.array[,,t+5*data.subset3$R])
        quantile.tmp <- quantile(tmp, c(0.025, 0.5, 0.975, 0.25, 0.75))
        log10x.exp[t,] <- quantile.tmp[c(2)]
      }
    }
    # Calculate the probability of x being above DT
    dg<-(sum(log10x.exp > rejection))/200
    #finalresults[i]<-dg  
    
  }
  result=(dg=dg)
  return(result)
}

PRC_four=function(E,filename.calibration,results,rejection)
{
  ###
  # This bisection function takes Calibration data, and the output of PROC_one and the output fro PROC_two (DT) 
  # The output will the detection limit in the log concentration.
  #
  #Input:	E, indictor correspond to a row of parameter values	
  # 		filename.calibration, The location of calibration data
  #			results, PROC_one output   (model parameters are \beta_0,\beta_1,\beta_2 and \sigma)
  #        rejection, is a vector of detection threshold values (PROC_two)
  ###
  # A is the desire beta level, here beta=0.087
  A=0.913
  # logxl is a starting value where we know LOD is above it
  logxl=-12
  # evaluate the value of logxl using PROC_three
  flogxl<-PROC_three_multiple_ISEs(filename.calibration,Z=2,fun.ahat=results[[1]][E,],fun.bhat=results[[2]][E,],fun.cstarhat=results[[3]][E,],fun.sigmahat=results[[4]][E,],newx=10^logxl, rejection=DT[[E]])
  fun.probability_fl= flogxl-A
  
  # logxr is a starting value where we know LOD is below it
  logxr=-3
  # evaluate the value of logxr using PROC_three
  flogxr<-PROC_three_multiple_ISEs(filename.calibration,Z=2,fun.ahat=results[[1]][E,],fun.bhat=results[[2]][E,],fun.cstarhat=results[[3]][E,],fun.sigmahat=results[[4]][E,],newx=10^logxr, rejection=DT[[E]])
  
  fun.probability_fr= flogxr-A
  logxnew=rep(NA,1)
  fun.probability_fnew=rep(NA,1)
  flr=rep(NA,1)
  loopmax=100 # maximum number of bisection iterations
  # t
  for (i in 1:loopmax){
    # pick the point in the middle between logxl and logxr as a new x to be evaluated
    logxnew[i]=(logxl+logxr)/2
    # evaluate the value of logxnew using PROC_three
    flogxnew<- PROC_three_multiple_ISEs(filename.calibration,Z=2, fun.ahat=results[[1]][E,],fun.bhat=results[[2]][E,],fun.cstarhat=results[[3]][E,],fun.sigmahat=results[[4]][E,] ,newx=10^logxnew[i], rejection=DT[[E]])
    
    fun.probability_fnew=flogxnew-A
    # replace logxl by logxnew if flr is greater than zero
    flr=(fun.probability_fnew)*(fun.probability_fl)
    if(flr>0){
      logxl=logxnew[i]
      fun.probability_fl=fun.probability_fnew
    }
    else{
      # replace logxr by logxnew if flr is less than zero
      
      logxr=logxnew[i]
      fun.probability_fr=fun.probability_fnew
      
    }
    if(abs(logxl-logxr)<5*10^-7) break # The bisection loop will stop on the 24th iteration
  }
  # the last value in the bisection loop is LOD
  res=logxnew[[i]]
  return(res)
  
  
  
  
}

# This is the model is for PROC_one (single ISE)
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

# This is the model is for PROC_one (multipe ISEs)
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

# Generating initial values for PROC_one (single ISE)
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


# Generating initial values for PROC_one (multipe ISEs)
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

# Generating initial values for PROC_two and PROC_three (single ISE)
gen.inits.single1<- function(data, a.init=NA, b.init=NA, cstar.init=NA, 
                             logc.limits = c(-8.9, -1.1), sigma.upper=14, stdadd = F, offset = 1, calibration.only = F)
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
  
  # Generate inital values for experimental data
  # Values with ISEs below S/N=3 are set to E(ISE|x=0) + 3*sigma
  if (!calibration.only)	{
    x.exp = rep(NA, data$M)
    if (stdadd != T) {
      emf.tmp <- data$emf.exp
      emf.tmp2 <- pmax(emf.tmp, mu0)
      x.tmp <- pmax((10^( (emf.tmp2 - a)/b ) - cstar^10), 10^logc.limits[1])
    }
    if (stdadd == T) {
      delta.emf.tmp = data$delta.emf
      V.add.tmp = data$V.add
      conc.add.tmp = data$conc.add
      V.s.tmp = data$V.s
      #	x.tmp = (V.add.tmp*conc.add.tmp/V.s.tmp)/(10^(delta.emf.tmp/b))
      x.tmp = pmax(10^logc.limits[1],
                   (V.add.tmp*conc.add.tmp/V.s.tmp)/(10^(delta.emf.tmp/b)))
    }
    
    # Initial values are between 90% and 100% of the average estimate
    x.exp = x.tmp*runif(1, 0.9, 1)
    log10x.exp = log10(x.exp)
    INITS <- list(a=a, b=b, cstar=cstar, logsigma=logsigma)
    init.stretch = 0.01
    
    inits <- list(a=a*runif(1, 1 - init.stretch, 1 + init.stretch),
                  b=b*runif(1, 1 - init.stretch, 1 + init.stretch),
                  cstar=cstar*runif(1, 1 - init.stretch, 1 + init.stretch),
                  sigma=sigma*runif(1, 1 - init.stretch, 1 + init.stretch),
                  log10x.exp=log10x.exp)
    
    
    
    
  }
  if (calibration.only) {
    INITS <- list(a=a, b=b, cstar=cstar, logsigma=logsigma)
    init.stretch = 0.01
    inits <- list(a=a*runif(1, 1 - init.stretch, 1 + init.stretch),
                  b=b*runif(1, 1 - init.stretch, 1 + init.stretch),
                  cstar=cstar*runif(1, 1 - init.stretch, 1 + init.stretch),
                  sigma=sigma*runif(1, 1 - init.stretch, 1 + init.stretch))
  }
  
  return(inits)
  
}


# Generating initial values for PROC_two and PROC_three (multiple ISEs)
gen.inits.multiple <- function(data, a.init= NA, b.init= NA, cstar.init= NA, offset = 1,
                               logc.limits = c(-8.9, -1.9), sigma.upper = 14, stdadd = F, calibration.only=F) {
  ###
  # Similar to gen.inits.single, but for multiple ISEs
  # If initial values are specified, the should be vectors with length equal to the number of ISEs
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
  
  if (!calibration.only) {
    # Generate inital values for experimental data
    # Values with ISEs below S/N=3 are set to E(ISE|x=0) + 3*sigma
    x.exp = rep(NA, data$M)
    
    if (stdadd != T) {
      for (i in 1:data$M) {
        emf.tmp <- data$emf.exp[data$xID.exp==i]
        ISEID.tmp <- data$ISEID.exp[data$xID.exp==i]
        emf.tmp2 <- pmax(emf.tmp, mu0[ISEID.tmp])
        x.tmp <- 10^( (emf.tmp2 - 
                         a[ISEID.tmp])/b[ISEID.tmp] ) - 
          cstar[ISEID.tmp]^10
        # Initial values are between 90% and 100% of the average estimate
        # but at least equal to the lower limit for log c
        x.exp[i] = max(10^logc.limits[1], mean(x.tmp))*runif(1, 0.9, 1)
      }
    }
    
    if (stdadd == T) {
      for (i in 1:data$M) {
        delta.emf.tmp = data$delta.emf[data$xID.exp==i]
        ISEID.tmp <- data$ISEID.exp[data$xID.exp==i]
        V.add.tmp = data$V.add[data$xID.exp==i]
        conc.add.tmp = data$conc.add[data$xID.exp==i]
        V.s.tmp = data$V.s[data$xID.exp==i]
        x.tmp = (V.add.tmp*conc.add.tmp/V.s.tmp)/(10^(delta.emf.tmp/b[ISEID.tmp]))
        # Initial values are between 90% and 100% of the average estimate
        x.exp[i] = mean(x.tmp)*runif(1, 0.9, 1)
      }
    }
    
    log10x.exp = log10(x.exp)
    
    INITS <- list(a=a, b=b, cstar=cstar, logsigma=logsigma)
    init.stretch = 0.01
    inits <- list(a=a*runif(data$R, 1 - init.stretch, 1 + init.stretch),
                  b=b*runif(data$R, 1 - init.stretch, 1 + init.stretch),
                  cstar=cstar*runif(data$R, 1 - init.stretch, 1 + init.stretch),
                  sigma=sigma*runif(data$R, 1 - init.stretch, 1 + init.stretch),
                  log10x.exp=log10x.exp)
  }
  if (calibration.only) {
    INITS <- list(a=a, b=b, cstar=cstar, logsigma=logsigma)
    init.stretch = 0.01
    inits <- list(a=a*runif(data$R, 1 - init.stretch, 1 + init.stretch),
                  b=b*runif(data$R, 1 - init.stretch, 1 + init.stretch),
                  cstar=cstar*runif(data$R, 1 - init.stretch, 1 + init.stretch),
                  sigma=sigma*runif(data$R, 1 - init.stretch, 1 + init.stretch))
  }
  return(inits)
}

# This is the model is for PROC_two and PROC_three (single ISE)
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

# This is the model is for PROC_two and PROC_three (multiple ISEs)
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