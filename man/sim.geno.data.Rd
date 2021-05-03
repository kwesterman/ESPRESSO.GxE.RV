\name{sim.geno.data}
\alias{sim.geno.data}
\title{Generates genotypes for a genetic variant}
\usage{
sim.geno.data(num.obs = NULL, geno.model = NULL, MAF = NULL)
}
\arguments{
  \item{num.obs}{number of observations to generate.}

  \item{geno.model}{genetic model; binary=0 and
  additive=1.}

  \item{MAF}{minor allele frequency of the SNP (in ESPRESSO
  this is the frequency of the 'at risk' allele)}
}
\value{
a dataframe that contains the following variables:
\code{allele.A} major allele \code{allele.B} minor allele
\code{genotype} genotype
}
\description{
Generates two alleles and combines them to form the
genotype of a SNP under a binary or additive genetic model.
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

