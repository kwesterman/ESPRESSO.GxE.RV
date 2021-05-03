\name{power.calc}
\alias{power.calc}
\title{Calculates the empirical and theoretical power}
\usage{
power.calc(pval = NULL, z.values = NULL, mean.model.z = NULL)
}
\arguments{
  \item{pval}{cut-off p-value defining statistical
  significance.}

  \item{z.values}{z-statistic of the determinant.}

  \item{mean.model.z}{mean z-statistic of the environmental
  determinant.}
}
\value{
a list that contains the computed empirical power and
theoretical power.
}
\description{
The function determines the empirical and theoretical
power. The empirical power is the proportion of simulations
in which the z-statistic for the parameter of interest
exceeds the z-statistic for the desured level if
statistical significance. The theoretical power is the
power of the study.
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

