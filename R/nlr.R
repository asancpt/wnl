nlr = function(Fx, Data, pNames, IE, LB, UB, Error="A", ObjFx=ObjDef, SecNames, SecForms, Method="L-BFGS-B")
{
#  e = new.env(parent=globalenv()) # environment should exist
  e$Fx = Fx # function of structural model. Fx should return a vector of the same length to Y
  e$DATA = Data # Fx use this
  e$Y  = Data[,"DV"] # Observation values, Data should have "DV" column.
  e$nRec = length(e$Y)
  e$pNames = pNames # parameter names in the order of Fx arguments
  e$IE = IE # initial estimate of Fx arguments
  e$nTheta = length(IE)
  e$Error = UT(Error) # BasicUtil fx toupper and Trim
  e$AddErrVar = (min(e$Y[e$Y > 0])/4)^2 # initial estimate of addtive error varaince
  e$PoisErrVar = 1    # initial estimate of proportional error varaince
  e$PropErrVar = 0.1  # initial estimate of proportional error varaince
  e$Obj = ObjFx

  if (e$Error == "A") {
    e$nEps = 1
    e$IE = c(e$IE, e$AddErrVar)
    e$pNames = c(e$pNames, "AddErrVar")
  } else if (e$Error == "POIS") {
    e$nEps = 1
    e$IE = c(e$IE, e$PoisErrVar)
    e$pNames = c(e$pNames, "PoisErrVar")
  } else if (e$Error == "P") {
    e$nEps = 1
    e$IE = c(e$IE, e$PropErrVar)
    e$pNames = c(e$pNames, "PropErrVar")
  } else if (e$Error == "C") {
    e$nEps = 2
    e$IE = c(e$IE, e$AddErrVar, e$PropErrVar)
    e$pNames = c(e$pNames, "AddErrVar", "PropErrVar")
  } else {
    print("Using custom objective function!\nSet e$nEps, e$IE, e$pNames yourself!")
  }

  e$nPara = e$nTheta + e$nEps # number of parameters
  if (missing(LB)) {
    e$LB = rep(0, e$nPara) # lower bound
  } else {
    e$LB = c(LB, rep(0, e$nEps))
  }
  if (missing(UB)) {
    e$UB = rep(1e+6, e$nPara) # upper bound
  } else {
    e$UB = c(UB, rep(1e+6, e$nEps))
  }

  e$alpha  = 0.1 - log((e$IE - e$LB)/(e$UB - e$LB)/(1 - (e$IE - e$LB)/(e$UB - e$LB))) # scaling constant
  e$SGindex = (e$nTheta + 1):(e$nTheta + e$nEps) # index of error variance(s)

  e$r = optim(rep(0.1, e$nPara), ObjEst, method=Method)
  e$PE = exp(e$r$par - e$alpha)/(exp(e$r$par - e$alpha) + 1)*(e$UB - e$LB) + e$LB

  e$InvCov = hessian(e$Obj, e$PE)/2  # FinalEst from EstStep()
  e$Cov = solve(e$InvCov)
  colnames(e$Cov) = e$pNames
  rownames(e$Cov) = e$pNames

  e$SE = sqrt(diag(e$Cov))
  e$Correl = cov2cor(e$Cov)
  e$EigenVal = eigen(e$Correl)$values
  
  e$PE = c(e$PE, sqrt(e$PE[e$SGindex]))
  e$SE = c(e$SE, e$SE[e$SGindex]/2/sqrt(e$PE[e$SGindex])) # Delta method, See Wackerly p484, gr = (x^0.5)' = 0.5x^(-0.5), gr^2 = 1/(4*x)
  e$RSE = e$SE / e$PE * 100
  e$Est = rbind(e$PE, e$SE, e$RSE)

#  colnames(e$Est) = c(e$pNames)
  if(e$Error == "A") {
    colnames(e$Est) = c(e$pNames, "AddErrSD")
  } else if (e$Error == "POIS") {
    colnames(e$Est) = c(e$pNames, "PoisErrSD")
  } else if (e$Error == "P") {
    colnames(e$Est) = c(e$pNames, "PropErrSD")
  } else if (e$Error == "C") {
    colnames(e$Est) = c(e$pNames, "AddErrSD", "PropErrSD")
  }   
  rownames(e$Est) = c("PE", "SE", "RSE")

  if (!missing(SecNames)) {
    nSec = length(SecNames)
    tRes = matrix(nrow=3, ncol=nSec)
    colnames(tRes) = SecNames
    rownames(tRes) = c("PE", "SE", "RSE")
    for (i in 1:nSec) {
      tRes[,i] = t(Secondary(SecForms[[i]], e$Est["PE",1:e$nPara], e$Cov))
    }
    e$Est = cbind(e$Est, tRes)
  }

  e$Residual = e$Y - e$Fx(e$PE[1:e$nTheta])
  e$run = run.test(e$Residual)
  e$AIC = e$nRec*log(2*pi) + e$r$value + 2*e$nPara
  e$AICc = e$AIC + 2*e$nPara*(e$nPara + 1)/(e$nRec - e$nPara - 1)
  
  Result = list(e$Est, e$Cov, e$run, e$AIC, e$AICc)
  names(Result) = c("Est", "Cov", "run", "AIC", "AICc")
  return(Result)
}

