\name{sim.QTL.data.GxE}
\alias{sim.QTL.data.GxE}
\title{Simulates subjects for continuous outcome}
\usage{
sim.QTL.data.GxE(n = NULL, ph.mean = NULL, ph.sd = NULL, freq = NULL,
  g.model = NULL, g.efkt = NULL, e.model = NULL, e.efkt = NULL,
  e.prev = NULL, e.mean = NULL, e.sd = NULL, e.low.lim = NULL,
  e.up.lim = NULL, i.efkt = NULL, pheno.rel = NULL)
}
\arguments{
  \item{n}{Number of subjects to simulate}

  \item{ph.mean}{statistical mean}

  \item{ph.sd}{standard deviation}

  \item{freq}{Minor allele frequencies of the genetic
  variant}

  \item{g.model}{Genetic model; 0 for binary and 1 for
  continuous}

  \item{g.efkt}{Effects of the genetic variant}

  \item{e.model}{Model of the environmental exposure}

  \item{e.efkt}{Effects of the environment determinats}

  \item{e.prev}{Prevalences of the environmental
  determinants}

  \item{e.mean}{Mean under quantitative-normal model}

  \item{e.sd}{Standard deviation under quantitative-normal
  model}

  \item{e.low.lim}{lower limit under quantitative-uniform
  model}

  \item{e.up.lim}{upper limit under quantitative-uniform
  model}

  \item{i.efkt}{Interaction effect}

  \item{pheno.rel}{reliability of the assessment for a
  quantitative outcome.}
}
\value{
A matrix
}
\description{
Generates the specified number of subjects
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

