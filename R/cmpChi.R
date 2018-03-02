cmpChi = function(r1, r2)
{
  npar1 = dim(r1$Cov)[1]
  npar2 = dim(r2$Cov)[1]
  ofv1 = r1$'Objective Function Value'
  ofv2 = r2$'Objective Function Value'
  if ((npar2 - npar1)*(ofv2 - ofv1) < 0) {
    p.val = 1 - pchisq(abs(ofv2 - ofv1), abs(npar2 - npar1))
  } else {
    p.val = 0 
  }
  return(p.val)
}