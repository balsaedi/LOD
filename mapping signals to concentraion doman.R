
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
                n.iter=100000, n.chains=1, n.burnin=90000, n.thin=1,
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


IterationCheck<-function(filename.experimental,results,SampleNumber)
{
  ###
  
  # This function takes Experimental data,  the output of PROC_one, and sample number (the sample needed to be estimated)
  # The output will be a vector of iteration numbers that can give an correct bisection procedure. Any combination of parameters
  # values (from PROC ONE) that gives the same signal in the insial evaluation of the bisection method, will be excluded.
  #Input:	filename.experimental, The location of calibration data
  #			results, PROC_one output   (model parameters are \beta_0,\beta_1,\beta_2 and \sigma)
  #			SampleNumber, is the number of the sample if interest, ,e.g.  SampleNumber = 1
  ###

  # Load experimental data
  data.exp = read.delim(filename.experimental,
                        header = TRUE, sep = "\t", quote="\"", dec=".", fill = TRUE, comment.char="")
  # Format data from the experimental file
  emf.exp=data.exp$emf
  xID.exp = data.exp$SampleID
  y1=emf.exp[xID.exp==SampleNumber][1]
  y2=emf.exp[xID.exp==SampleNumber][2]
  # Create an empty vector to store, the iteration number that passes the test
  E1=NA
  ## Here i choose 10000 becasue i made the output of PROC ONE to be equal to that
  for (E in 1:10000)
  {
    # The bisection method
    # IF logxl=-12 and logxr=-1 gives the same signals for the derivative, then the bisection method will not work
    # very well
    logxl=-12
    # the evaluation function, which is the derivavtive of the product of the distributions of two signals
    # with respect to x (concentration value)
    flogxl<-(-1/(results[[4]][E,1]^2))*(y1-(results[[1]][E,1]+(results[[2]][E,1]*log(10^logxl+(results[[3]][E,1])^10)/log(10))))*(-(results[[2]][E,1]/(log(10)*(10^logxl+(results[[3]][E,1])^10))))+(-1/(results[[4]][E,2]^2))*(y2-(results[[1]][E,2]+(results[[2]][E,2]*log(10^logxl+(results[[3]][E,2])^10)/log(10))))*(-(results[[2]][E,2]/(log(10)*(10^logxl+(results[[3]][E,2])^10))))
    fun.probability_fl= flogxl
    logxr=-1
    flogxr<-(-1/(results[[4]][E,1]^2))*(y1-(results[[1]][E,1]+(results[[2]][E,1]*log(10^logxr+(results[[3]][E,1])^10)/log(10))))*(-(results[[2]][E,1]/(log(10)*(10^logxr+(results[[3]][E,1])^10))))-(1/(results[[4]][E,2]^2))*(y2-(results[[1]][E,2]+(results[[2]][E,2]*log(10^logxr+(results[[3]][E,2])^10)/log(10))))*(-(results[[2]][E,2]/(log(10)*(10^logxr+(results[[3]][E,2])^10))))
    fun.probability_fr= flogxr
    # The condition, if both (flogxl and flogxr) have the same signals, ignore that iteration
    if (fun.probability_fl*fun.probability_fr>0) next
    E1[E]=E
  }
  return(E1[!is.na(E1)])
}

ss1<-TR("C:/Users/balsaedi/Desktop/ch3/DATA/generate/exp_data_7obs2.txt",results=ddd,SampleNumber=1)


PRC_four=function(E,filename.experimental,results,SampleNumber)
{
  ###
  # This bisection function takes experimental data, and the output of PROC_one and the output from IterationCheck 
  # The output will be the concentration value in the log concentration that correspond to the specified SampleNumber.
  #
  #Input:	E, indictor correspond to a row of parameter values	that is supplied by IterationCheck
  # 		filename.experimental, The location of experimental data
  #			results, PROC_one output   (model parameters are \beta_0,\beta_1,\beta_2 and \sigma)
  #       SampleNumber, is the number of the sample if interest, ,e.g.  SampleNumber = 1
  ###
  
  # Load experimental data
    data.exp = read.delim(filename.experimental,
                        header = TRUE, sep = "\t", quote="\"", dec=".", fill = TRUE, comment.char="")

     # Format data from the experimental file
  emf.exp=data.exp$emf
  xID.exp = data.exp$SampleID
  y1=emf.exp[xID.exp==SampleNumber][1]
  y2=emf.exp[xID.exp==SampleNumber][2]

  # The bisection method
  #  logxl=-12 and logxr=-1 will gives different signals for the derivative, 
  logxl=-12
  flogxl<-(-1/(results[[4]][E,1]^2))*(y1-(results[[1]][E,1]+(results[[2]][E,1]*log(10^logxl+(results[[3]][E,1])^10)/log(10))))*(-(results[[2]][E,1]/(log(10)*(10^logxl+(results[[3]][E,1])^10))))+(-1/(results[[4]][E,2]^2))*(y2-(results[[1]][E,2]+(results[[2]][E,2]*log(10^logxl+(results[[3]][E,2])^10)/log(10))))*(-(results[[2]][E,2]/(log(10)*(10^logxl+(results[[3]][E,2])^10))))
  
  fun.probability_fl= flogxl
  logxr=-1
  flogxr<-(-1/(results[[4]][E,1]^2))*(y1-(results[[1]][E,1]+(results[[2]][E,1]*log(10^logxr+(results[[3]][E,1])^10)/log(10))))*(-(results[[2]][E,1]/(log(10)*(10^logxr+(results[[3]][E,1])^10))))-(1/(results[[4]][E,2]^2))*(y2-(results[[1]][E,2]+(results[[2]][E,2]*log(10^logxr+(results[[3]][E,2])^10)/log(10))))*(-(results[[2]][E,2]/(log(10)*(10^logxr+(results[[3]][E,2])^10))))
  fun.probability_fr= flogxr
  logxnew=rep(NA,1)
  fun.probability_fnew=rep(NA,1)
  flr=rep(NA,1)
  loopmax=10
  # the bisection algorithm
  for (Y in 1:loopmax){
    
    logxnew[Y]=(logxl+logxr)/2
    flogxnew<-(-1/(results[[4]][E,1]^2))*(y1-(results[[1]][E,1]+(results[[2]][E,1]*log(10^logxnew[Y]+(results[[3]][E,1])^10)/log(10))))*(-(results[[2]][E,1]/(log(10)*(10^logxnew[Y]+(results[[3]][E,1])^10))))-(1/(results[[4]][E,2]^2))*(y2-(results[[1]][E,2]+(results[[2]][E,2]*log(10^logxnew[Y]+(results[[3]][E,2])^10)/log(10))))*(-(results[[2]][E,2]/(log(10)*(10^logxnew[Y]+(results[[3]][E,2])^10))))
    fun.probability_fnew=flogxnew
    flr=(fun.probability_fnew)*(fun.probability_fl)
    if(flr>0){
      logxl=logxnew[Y]
      fun.probability_fl=fun.probability_fnew
    }
    else{
      
      logxr=logxnew[Y]
      fun.probability_fr=fun.probability_fnew
      
    }
    
  }
  # the last value in the bisection loop is our choosen value
  res=logxnew[[Y]]
  return(res)
  
}

ss<-sapply(ss1,PRC_four,"C:/Users/balsaedi/Desktop/ch3/DATA/generate/exp_data_7obs2.txt",results=ddd,SampleNumber=1)
