\name{sim.subject.data}
\alias{sim.subject.data}
\title{Simulates the individual effect related to heterogeneity in baseline disease risk}
\usage{
sim.subject.data(num.obs = NULL, baseline.OR = NULL)
}
\arguments{
  \item{num.obs}{number of observations to simulate.}

  \item{baseline.OR}{baseline odds ratio for subject on 95
  percent population centile versus 5 percentile. This
  parameter reflects the heterogeneity in disease risk
  arising from determinantes that have not been measured or
  have not been included in the model.}
}
\value{
a numerical vector.
}
\description{
The variation in baseline disease risk is assumed to be
normally distributed on a logistic scale. If this parameter
is set to 10, the implication is that a 'high risk' subject
(someone at the upper 95 percent entile of population risk)
is, all else being equal, at 10 times the offs of
developing disease compared to someone else who is at 'low
risk' (at the lower 5 percent centile of population risk).
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

