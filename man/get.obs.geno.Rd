\name{get.obs.geno}
\alias{get.obs.geno}
\title{Adds some error to genotype data}
\usage{
get.obs.geno(allele.A = NULL, allele.B = NULL, geno.model = NULL,
  MAF = NULL, geno.error = NULL)
}
\arguments{
  \item{allele.A}{Allele A}

  \item{allele.B}{Allele B}

  \item{geno.model}{genetic model; binary=0 and additive=1}

  \item{MAF}{minor allele frequency of the SNP (in ESPRESSO
  this is the frequency of the 'at risk' allele)}

  \item{geno.error}{a vector with two values, the
  misclassififcation rates related to the sensitivity and
  specificity of the assessment of the alleles, i.e.
  1-sensitivity and 1-specificity}
}
\value{
a dataframe that contains the below data:
\code{observed.genotype} observed genotypes
\code{observed.allele.A} observed A alleles
\code{observed.allele.B} observed B alleles
}
\description{
Simulates errors and adds it to the true data to obtain
observed data.  The alleles simulated by the function
sim.geno.data are randomly misclassified and used to form
new genotypes that represent the observed genotypes.
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

