\name{sim.pheno.qtl.GxE}
\alias{sim.pheno.qtl.GxE}
\title{Simulates continuous outcome data}
\usage{
sim.pheno.qtl.GxE(numsubjects = NULL, pheno.mean = NULL, pheno.sd = NULL,
  genotype = NULL, geno.efkt = NULL, environment = NULL,
  env.efkt = NULL, interaction = NULL, int.efkt = NULL)
}
\arguments{
  \item{num.subjects}{Number of subjects to simulate}

  \item{pheno.mean}{statistical mean}

  \item{pheno.sd}{standard deviation}

  \item{genotype}{Genotype}

  \item{geno.efkt}{Effects of the genetic variant}

  \item{environment}{Exposure data for environment}

  \item{env.efkt}{Effects of the environmental determiants}

  \item{interaction}{data}

  \item{int.efkt}{Interaction effect}
}
\value{
A dataframe of phenotype
}
\description{
Uses the effects data of the determinants to construct a
linear predictor(LP). The outcome is normally distributed
variable generated with a mean equal to LP and a standard
deviation of 1. Some error is then added to the simulated
outcome to obtained the observed outcome.
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

