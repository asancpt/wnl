\name{wnl5}
\alias{wnl5}
\title{Old type WinNonlin - Least Square not MLE}
\description{
It performs old type Winnonlin regression.
}
\usage{
wnl5(Fx, Data, pNames, IE, LB, UB, Error="A", ObjFx=ObjLS)
}
\arguments{
  \item{Fx}{Function for structural model. It should return a vector of same length to observations.}
  \item{Data}{Data table which will be used in Fx. Fx should access this with \code{e$DATA}.}
  \item{pNames}{Parameter names in the order of Fx arguments}
  \item{IE}{Initial estimates of parameters}
  \item{LB}{Lower bound for \code{optim} function. Default value is 0.}
  \item{UB}{Upper bound for \code{optim} function. Default value is 1e+06.}
  \item{Error}{Error model. One of \code{"POIS"} for poisson error, \code{"PROP"} for proportional error, and others for additive error model.}
  \item{ObjFx}{Objective function to be minimized. Default is least square function.}
}
\details{
This uses scaled transformed parameters and environment \code{e} internally. Here we do not provide standard error. If you want standard error, use \code{nlr}.
}
\value{
This returns point estimates, run test results on residuals, AIC, and SBC. Details are stored in the environment \code{e}.
}
\author{Kyun-Seop Bae <k@acr.kr>}
\examples{
tData = Theoph
colnames(tData) = c("ID", "BWT", "DOSE", "TIME", "DV")

fPK = function(THETA) # Prediction function
{
  DOSE = 320000 # in microgram
  TIME = e$DATA[,"TIME"] # use data in e$DATA

  K  = THETA[1]
  Ka = THETA[2]
  V  = THETA[3]
  Cp  = DOSE/V*Ka/(Ka - K)*(exp(-K*TIME) - exp(-Ka*TIME))
  return(Cp)
}

IDs = unique(tData[,"ID"])
nID = length(IDs)
for (i in 1:nID) {
  Data = tData[tData$ID == IDs[i],]
  Res = wnl5(fPK, Data, pNames=c("k", "ka", "V"), IE=c(0.1, 3, 500))
  print(paste("## ID =", i, "##"))
  print(Res)
}

}
\keyword{wnl5}
