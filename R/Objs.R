ObjEst = function(vPara)
{
  b0 = exp(vPara - e$alpha)
  vPara = b0/(b0 + 1)*(e$UB - e$LB) + e$LB
  return(e$Obj(vPara))
}

ObjDef = function(vPara) # Default Obj for Covariance Step (original obj)
{
  Fi = e$Fx(vPara[1:e$nTheta])
  Ri = e$Y - Fi

  if (e$Error == "A") {
    Ci = rep(vPara[e$SGindex], e$nRec)
  } else if (e$Error == "POIS") {
    Ci = vPara[e$SGindex]*Fi     # Caution on Fi is zero
    Ci[Fi == 0] = 1 # Weight cannot be calculated with zero values
    Ri[Fi == 0] = 0
  } else if (e$Error == "P") {
    Ci = vPara[e$SGindex]*Fi*Fi  # Caution on Fi is zero
    Ci[Fi == 0] = 1 # Weight cannot be calculated with zero values
    Ri[Fi == 0] = 0
  } else if (e$Error == "C") {
    Ci = rep(vPara[e$SGindex[1]], e$nRec) + vPara[e$SGindex[2]]*Fi*Fi
  }

  return(sum(log(Ci) + Ri*Ri/Ci))
}

ObjLS = function(vPara) # Default Obj for Covariance Step (original obj)
{
  Fi = e$Fx(vPara)
  Ri = e$Y - Fi
  if (e$Error == "POIS") {
    Ri[Fi != 0] = Ri[Fi != 0] / sqrt(Fi[Fi != 0])   # Fi should contain zero.
    Ri[Fi == 0] = 0
  } else if (e$Error == "PROP") {
    Ri[Fi != 0] = Ri[Fi != 0] / Fi[Fi != 0]         # Fi should contain zero.
    Ri[Fi == 0] = 0
  }

  return(sum(Ri*Ri))
}
