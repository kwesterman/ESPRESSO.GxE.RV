\name{sim.CC.data.GxE}
\alias{sim.CC.data.GxE}
\title{Simulates case and controls}
\usage{
sim.CC.data.GxE(n = NULL, ncases = NULL, ncontrols = NULL,
  max.sample.size = NULL, pheno.prev = NULL, freq = NULL,
  g.model = NULL, g.OR = NULL, e.model = NULL, e.prev = NULL,
  e.mean = NULL, e.sd = NULL, e.low.lim = NULL, e.up.lim = NULL,
  e.OR = NULL, i.OR = NULL, b.OR = NULL, ph.error = NULL)
}
\arguments{
  \item{n}{Number of observations to generate per
  iteration}

  \item{ncases}{Number of cases to simulate}

  \item{ncontrols}{Number of controls to simulate}

  \item{max.sample.size}{Maximum number of observations
  allowed}

  \item{pheno.prev}{Prevalence of the binary outcome}

  \item{freq}{Minor allele frequency}

  \item{g.model}{Genetic model; 0 for binary and 1 for
  additive}

  \item{g.OR}{Odds ratios of the genetic determinant}

  \item{e.model}{Model of the environmental exposure}

  \item{e.prev}{Prevelance of the environmental
  determinates}

  \item{e.mean}{Mean under quantitative-normal model}

  \item{e.sd}{Standard deviation under quantitative-normal
  model}

  \item{e.low.lim}{Lower limit under quantitative-uniform
  model}

  \item{e.up.lim}{Upper limit under quantitative-uniform
  model}

  \item{e.OR}{Odds ratios of the environmental
  determinants}

  \item{i.OR}{Odds ration of the interaction}

  \item{b.OR}{Baseline odds ratio for subject on 95 percent
  population centile versus 5 percentile. This parameter
  reflects the heterogeneity in disease risk arising from
  determinates that have not been measured or have not been
  included in the model.}

  \item{ph.error}{misclassification rates: 1-sensitivity
  and 1-specificity}
}
\value{
A matrix
}
\description{
Generates affected and non-affected subjects until the set
sample size is achieved.
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

