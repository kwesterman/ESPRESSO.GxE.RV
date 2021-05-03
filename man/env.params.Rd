\name{env.params}
\alias{env.params}
\docType{data}
\title{Parameters to generate environmental exposure data}
\description{A table that contains parameters to set the environmental exposure}
\usage{data(env.params)}
\format{
  A data frame with 4 observations on the following 12 variables.
  \describe{
    \item{\code{scenario.id}}{scenario ID}
    \item{\code{env.model}}{model of the environmental exposure; 0 for binary and 1 for continuous}
    \item{\code{env.prevalence}}{prevelance of the 'at risk' environmental exposure}
    \item{\code{env.OR}}{odds ratio related to the environmental exposure (when exposure is binary)}
    \item{\code{env.efkt}}{effect size related to the environmental exposure (when exposure is quantitative)}
    \item{\code{env.mean}}{statistical mean of the normally distributed environment exposure}
    \item{\code{env.sd}}{standard deviation of the normally distributed environment exposure}
    \item{\code{env.low.lim}}{Lower limit of the uniformly distributed environment exposure}
    \item{\code{env.up.lim}}{Upper limit of the uniformly distributed environment exposure}
    \item{\code{env.sensitivity}}{sensitivity of the assessment of the binary exposure}
    \item{\code{env.specificity}}{specificity of the assessment of the binary exposure}
    \item{\code{env.reliability}}{reliability of the assessment of the quantative exposure}
  }
}

\examples{data(env.params)}
