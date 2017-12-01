\name{wnl-package}
\alias{wnl-package}
\alias{wnl}
\docType{package}
\title{Minimization Tool for Pharmacokinetic-Pharmacodynamic Data Analysis}
\description{This is a minimization tool to solve the examples in the book Gabrielsson J, Weiner D. 'Pharmacokinetic and Pharmacodynamic Data Analysis - Concepts and Applications' 5th ed. 2016. (ISBN:9198299107).}
\details{
Description: This is a minimization tool to solve all the examples(PK1-PK53, PD1-PD52) in the book 'Pharmacokinetic and Pharmacodynamic Data Analysis - Concepts and Applications'.
}
\author{Kyun-Seop Bae <k@acr.kr>}
\references{Gabrielsson J, Weiner D. Pharmacokinetic and Pharmacodynamic Data Analysis - Concepts and Applications. 5th ed. 2016.}
\keyword{Pharmacokinetics}
\keyword{Pharmacodynamics}
\examples{
tData = Theoph
colnames(tData) = c("ID", "BWT", "DOSE", "TIME", "DV")

fPK = function(THETA) # Prediction function
{
  DOSE = 320000 # in microgram
  TIME = e$DATA[,"TIME"]  # use data in e$DATA

  K    = THETA[1]
  Ka   = THETA[2]
  V    = THETA[3]

  Cp  = DOSE/V*Ka/(Ka - K)*(exp(-K*TIME) - exp(-Ka*TIME))
  return(Cp)
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
