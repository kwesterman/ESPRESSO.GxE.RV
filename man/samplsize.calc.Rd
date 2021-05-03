\name{samplsize.calc}
\alias{samplsize.calc}
\title{Calculates the sample size required to achieve the desired power}
\usage{
samplsize.calc(numcases = NULL, numcontrols = NULL, num.subjects = NULL,
  pheno.model = NULL, pval = NULL, power = NULL, mean.model.z = NULL)
}
\arguments{
  \item{numcases}{number of cases when outcome is binary.}

  \item{numcontrols}{number of controls when outcome is
  binary.}

  \item{num.subjects}{number of subjects when outcome is
  quantitative.}

  \item{pheno.model}{outcome model, binary=0 and
  quantitafive=1.}

  \item{pval}{cut-off p-value defining statistical
  significance.}

  \item{power}{desired power}

  \item{mean.model.z}{ratio of mean beta estimate over mean
  se estimate}
}
\value{
A table containing: \code{numcases.required} number of
cases required to achieve the desired power under binary
outcome model \code{numcontrols.required} number of
controls required to achieve the desired power under binary
outcome model \code{numsubjects.required} number of
subjects required to achieve the desired power under a
quantitative outcome model
}
\description{
Estimates by how much the simulated study size needs to be
inflated or shrank in order to obtain the specified level
of power. The ratio of z-statistic required for desired
power to mean model z-statistic obtained indicates the
relative changes in standard error required. This
corresponds to the relative change on scale and square root
of sample size.
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

