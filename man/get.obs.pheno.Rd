\name{get.obs.pheno}
\alias{get.obs.pheno}
\title{Generates observed outcome data}
\usage{
get.obs.pheno(phenotype = NULL, pheno.model = NULL, pheno.sd = NULL,
  pheno.error = NULL, pheno.reliability = NULL)
}
\arguments{
  \item{seed}{}

  \item{phenotype}{outcome status.}

  \item{pheno.model}{distribution of the outcome variable:
  binary=0, normal=1}

  \item{pheno.sd}{standard deviation of the outcome in the
  study population}

  \item{pheno.model}{distribution of the outcome variable:
  binary=0, normal=1 or uniform=2.}

  \item{pheno.sd}{standard deviation of the outcome in the
  study the population}

  \item{pheno.error}{misclassification rates: 1-sensitivity
  and 1-specificity}

  \item{pheno.reliability}{reliability of the assessment
  for a quantitative outcome.}
}
\value{
A dataframe containing: \code{true.phenotype} the error
free outcome data (true data). \code{observed.phenotype}
the true outcome data with some added error (observed
data).
}
\description{
Adds a set level of error to error free binary or
quantitative data (the true phenotype data) to obtain data
with a larger variance (the observed phenotype data).
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

