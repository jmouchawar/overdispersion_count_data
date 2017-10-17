###############
# Jason Mouchawar
# Overdispersion Simulation
###############
library(MASS)

dispersion = function(model){
  #Calculate the dispersion statistic
  return(sum(residuals(model, type = "pearson")^2)/model$df.residual)
}

###############
#Global Vars
nobs = 500
b0 = 1
b1 = .2
conf.level = .95
iterations = 250
alphavec = seq(0,5, by=.1)
###############

data.sim = function(nobs,b0,b1,alpha){
  # Simulates count data:
  # if alpha = 0, will simulate data from Poisson distribution ; if alpha > 0, will simulate overdispersed poisson data
  #  
  # Args:   nobs (numeric) ------------ number of observations
  #         b0 (numeric) -------------- value of intercept parameter for Poisson link function
  #         b1 (numeric) -------------- value of slope parameter for Poisson link function
  #         alpha (numeric/vector) ---- value(s) of overdispersion parameter
  #         
  # Returns: dataframe of linked response values (y), and predictor values (x)
  
  if(alpha == 0){
    
    x = runif(nobs)
    xb = b0 + b1*x
    y = rpois(nobs,exp(xb))
    df = as.data.frame(cbind(y,x))
    
  }else{
    
    x = runif(nobs)
    xb = b0 + b1*x
    y = rnbinom(nobs, size = 1/alpha, mu = (exp(xb))) #see ?rnibom var parameterization
    df = as.data.frame(cbind(y,x))
    
  }
  
  return(df)
}

poismodel.sim = function(conf.level,alpha){
  # Fits Poisson model to data:
  # Args:   conf.level (numeric) ---- confidence level for computing confidence intervals
  #         alpha (numeric/vector) -- value(s) of overdispersion parameter
  #         
  # Returns: (vector) of estimated values of:
  #   estimate b1, SE b1, pval, CI LL, CI UL (length=5)
  
  pois = glm(y ~ x, data = data.sim(nobs,b0,b1,alpha), family = poisson)
  
  output = c(summary(pois)$coefficients[2,1], #b1
             summary(pois)$coefficients[2,2], #se b1
             summary(pois)$coefficients[2,4], #pval
             confint(pois, 'x', level = conf.level))
  
  names(output)[1] <- "b1"
  names(output)[2] <- "se b1"
  names(output)[3] <- "pval"
  
  return(output)
  
}

simulation = function(iterations,alpha){
  # Runs simulation
  # Examine effect of overdispersion parameter on coverage probability of confidence intervals for the slope parameter
  # Args:   iterations (numeric) ---- number of iterations (denominator for coverage probability calculation)
  #                                       **NOTE** large values may take several minutes
  #         alpha (numeric/vector) -- value(s) of overdispersion parameter
  #
  # Returns: (matrix)   columns are values returned from poismodel.sim and contains. 
  #                     contains is binary indicator if b1 is contained in the calculated confidence interval.
  #                     each row = a new iteration. 
  #                     
  
  contains=NULL
  data.frame=NULL
  
  for(i in 1:iterations){
    
    temp = poismodel.sim(conf.level,alpha)
    
    if(temp[4] < b1 & temp[5] > b1){
      contains = 1
    }else{
      contains = 0
    }
    
    temp = c(temp,contains)
    
    data.frame = rbind(temp, data.frame)
    
  }
  
  rownames(data.frame) <- 1:iterations
  
  return(data.frame)
  
}


#Calculation of coverage probability
cov.prob = NULL
for(i in 1:length(alphavec)){
  
  tempsim = simulation(iterations, alpha=alphavec[i])
  cov.prob[i] = sum(tempsim[,6])/dim(tempsim)[1]
  
}

#Plot Coverage probability
plot(alphavec,cov.prob, 
     main = '95% CI Coverage Probability by Overdispersion Parameter',
     ylab = 'Coverage Probability', ylim = c(0,1),
     xlab = "Overdispersion Parameter")

