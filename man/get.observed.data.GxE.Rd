\name{get.observed.data.GxE}
\alias{get.observed.data.GxE}
\title{Generates exposure data with some error}
\usage{
get.observed.data.GxE(data = NULL, g.error = NULL, g.model = NULL,
  freq = NULL, e.error = NULL, e.model = NULL, e.prev = NULL,
  e.sd = NULL, e.reliability = NULL)
}
\arguments{
  \item{data}{Input table of simulated data considered as
  true data}

  \item{g.error}{Misclassification rates in genetic
  assessment: 1-sensitivity and 1-specificity}

  \item{g.model}{Genetic model; 0 for binary and 1 for
  additive}

  \item{freq}{Minor allele frequency}

  \item{e.error}{Misclassification rates in environmental
  exposures assessment: 1-sensitivity and 1-specificity}

  \item{e.model}{Model of the exposure: binary=0,
  quantitative-normal=1 or quantitative-uniform=2}

  \item{e.prev}{Prevalence of environmental exposure}

  \item{e.sd}{Standard Deviation}

  \item{e.reliability}{Reliability of the assessment of
  quantitative exposure}
}
\value{
A matrix
}
\description{
Uses functions make.obs.geno and make.obs.env to generate
effect data with a set level of error
}
\details{

}
\author{
Gaye A
}
\keyword{internal}

