\name{get.obs.env}
\alias{get.obs.env}
\title{Generates the observed environmental exposure data}
\usage{
get.obs.env(env.data = NULL, env.model = NULL, env.sd = NULL,
  env.prev = NULL, env.error = NULL, env.reliability = NULL)
}
\arguments{
  \item{env.data}{a vector of environmental measures that
  represents the true data.}

  \item{env.model}{distribution of the exposure: binary=0 ,
  normal=1 or uniform=2.}

  \item{env.prev}{prevalence of the environmental
  exposure.}

  \item{env.error}{misclassification rates: 1-sensitivity
  and 1-specificity.}

  \item{env.reliability}{reliability of the assessment of
  quantitative exposure.}
}
\value{
a dataframe with two coloumns/variables:
\code{true.environment} the error free exposure data (true
data). \code{observed.environment} the true esposure data
with some added error (observed data).
}
\description{
Adds a set level of error to binary or quantitative data
(the true data) to obtain data with a larger variance (the
observed data). The level of error is determined by the
misclassification rates in binary exposure and by the set
level of variance in the quantitative exposure.
}
\details{

}
\author{
Amadou Gaye
}
\keyword{internal}

