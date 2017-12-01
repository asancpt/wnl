\name{nlr}
\alias{nlr}
\title{Nonlinear Regression in R}
\description{
It performs nonlinear regression usually for pharmacokinetic and pharmacodynamic models.
}
\usage{
nlr(Fx, Data, pNames, IE, LB, UB, Error="A", ObjFx=ObjDef, SecNames, SecForms, 
    Method="L-BFGS-B")
}
\arguments{
  \item{Fx}{Function for structural model. It should return a vector of same length to observations.}
  \item{Data}{Data table which will be used in Fx. Fx should access this with \code{e$DATA}.}
  \item{pNames}{Parameter names in the order of Fx arguments}
  \item{IE}{Initial estimates of parameters}
  \item{LB}{Lower bound for \code{optim} function. Default value is 0.}
  \item{UB}{Upper bound for \code{optim} function. Default value is 1e+06.}
  \item{Error}{Error model. One of \code{"A"} for additive error, \code{"POIS"} for poisson error, \code{"P"} for proportional error, and \code{"C"} for combined error model.}
  \item{ObjFx}{Objective function to be minimized. Default is maximum likelihood estimation function.}
  \item{SecNames}{Names of secondary parameter estimates}
  \item{SecForms}{Formula to calculate the secondary parameter estimates}
  \item{Method}{\code{"L-BFGS-B"} is default. See \code{optim} for more detail.}
}
\details{
This uses scaled transformed parameters and environment \code{e} internally.
}
\value{
This returns point estimate, standard error, relative standard error, covariance matrix of estimates, run test results on residuals, AIC, and AICc. Details are stored in the environment \code{e}.
}
\author{Kyun-Seop Bae <k@acr.kr>}
\examples{
tData = Theoph
colnames(tData) = c("ID", "BWT", "DOSE", "TIME", "DV")

fPK = function(THETA) # Prediction function
{
  DOSE = 320000 # in microgram
  TIME = e$DATA[,"TIME"] # use data in e$DATA

  K    = THETA[1]
  Ka   = THETA[2]
  V    = THETA[3]

  P  = DOSE/V*Ka/(Ka - K) * (exp(-K*TIME) - exp(-Ka*TIME))
  return(P)
}

IDs = unique(tData[,"ID"])
nID = length(IDs)
for (i in 1:nID) {
  Data = tData[tData$ID == IDs[i],]
  Res = nlr(fPK, Data, pNames=c("k", "ka", "V"), IE=c(0.1, 3, 500), 
            SecNames=c("CL", "Thalf", "MRT"), SecForms=c(~V*k, ~log(2)/k, ~1/k))
  print(paste("## ID =", i, "##"))
  print(Res)
}

}
\keyword{nlr}
