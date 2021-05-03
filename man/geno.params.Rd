\name{geno.params}
\alias{geno.params}
\docType{data}
\title{Parameters to generate genetic exposures data}
\description{A table that contains parameters to set the genetic exposure}
\usage{data(geno.params)}
\format{
  A data frame with 4 observations on the following 7 variables.
  \describe{
    \item{\code{scenario.id}}{scenario ID}
    \item{\code{geno.model}}{genetic model: binary (0) or additive (1)}
    \item{\code{MAF}}{minor allele frequency of the SNP}
    \item{\code{geno.OR}}{odds ratio if outcome is binary}
    \item{\code{geno.efkt}}{effect size if outcome is quantitative}
    \item{\code{geno.sensitivity}}{sensitivity of the assessment of the SNP alleles}
    \item{\code{geno.specificity}}{specificity of the assessment of the SNP alleles}
  }
}
\examples{data(geno.params)}

