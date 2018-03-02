\name{cmpChi}
\alias{cmpChi}
\title{Compare model with Chi-square test}
\description{
It performs chi-square test for two models comparison.
}
\usage{
cmpChi(r1, r2)
}
\arguments{
  \item{r1}{A result from \code{nlr}}
  \item{r2}{Another result from \code{nlr}}
}
\details{
One model should include the other model.
}
\value{
Returns a p-value from \code{pchisq}
}
\author{Kyun-Seop Bae <k@acr.kr>}

\keyword{Model Comparison}
