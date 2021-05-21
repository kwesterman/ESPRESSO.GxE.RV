\name{power.calc}
\alias{power.calc}
\title{Calculates the empirical and theoretical power}
\usage{
power.calc(pval = NULL, p.values = NULL)
}
\arguments{
  \item{pval}{cut-off p-value defining statistical
  significance.}

  \item{p.values}{Vector of p-values (one per genetic variant).}

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

