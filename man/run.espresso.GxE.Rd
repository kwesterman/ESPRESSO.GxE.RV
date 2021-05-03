\name{run.espresso.GxE}
\alias{run.espresso.GxE}
\title{Runs a full ESPRESSO analysis}
\usage{
run.espresso.GxE(simulation.params = NULL, pheno.params = NULL,
  geno.params = NULL, env.params = NULL, scenarios2run = 1)
}
\arguments{
  \item{simulation.params}{general parameters for the
  scenario(s) to analyse}

  \item{pheno.params}{paramaters for the outcome variables}

  \item{geno.params}{parameters for the genetic
  determinant}

  \item{env.params}{parameters for the environmental
  determinant}

  \item{scenarios2run}{the indices of the scenarios one
  wish to analyse}
}
\value{
a summary table that contains both the input parameters and
the results of the analysis
}
\description{
This function calls the functions required to run a full
ESPRESSO analysis where the model consists of an outcome
(binary or continuous) determined by two interacting
covariates (a SNP and an environmental exposure)
}
\details{

}
\examples{
{

# load the table that hold the input parameters; each of the table
# hold parameters for 4 scenarios:
# scenario 1: a binary outcome determined by a binary SNP and binary exposure
# and an interaction between the the genetic variant and the environmental exposure
# scenario 2: a binary outcome determined by an additive SNP and continuous
# exposure and an interaction between the the genetic variant and the environmental exposure
# scenario 3: a quantitative outcome determined by a binary SNP and binary exposure,
# and an interaction between the the genetic variant and the environmental exposure
# scenario 4: a quantitative outcome determined by an additive SNP and continuous
# exposure and an interaction between the the genetic variant and the environmental exposure
data(simulation.params)
data(pheno.params)
data(geno.params)
data(env.params)

# run the function for the first two scenarios, two binomial models
run.espresso.GxE(simulation.params, pheno.params, geno.params, env.params, scenarios2run=c(1,2))

# run the function for the last two scenarios, two gaussian models
run.espresso.GxE(simulation.params, pheno.params, geno.params, env.params, scenarios2run=c(3,4))

}
}
\author{
Gaye A.
}

