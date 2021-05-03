\name{pheno.params}
\alias{pheno.params}
\docType{data}
\title{Parameters for the outcome variable}
\description{A table with 6 columns: scenario.id, pheno.model, disease.prev, pheno.sensitivity, pheno.specificity, pheno.reliability}
\usage{data(pheno.params)}
\format{
  A table with 8 columns and 1 row for each scenario
  \describe{
    \item{\code{scenario.id}}{scenario ID}
    \item{\code{pheno.model}}{type of outcome; 0 for binary and 1 for continuous}
    \item{\code{disease.prev}}{prevalence of the binary outcome}
    \item{\code{pheno.mean}}{statistical mean in the population the study is based on}
    \item{\code{pheno.sd}}{standard deviation in the population the study is based on}
    \item{\code{pheno.sensitivity}}{sensitivity of the assessment of binary outcome}
    \item{\code{pheno.specificity}}{specificity of the assessment of binary outcome}
    \item{\code{pheno.reliability}}{reliability of the assessment of continuous outcome}
  }
}

\examples{data(pheno.params)}
