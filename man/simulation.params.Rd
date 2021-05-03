\name{simulation.params}
\alias{simulation.params}
\docType{data}
\title{Main parameters of the simulations}
\description{A table with 11 columns: scenario.id, seed.val, numsims, numcases, numcontrols, numsubjects, interaction.OR, interaction.efkt, RR.5.95, p.val, power}
\usage{data(simulation.params)}
\format{
  A table with 11 columns and 1 row for each scenario
  \describe{
    \item{\code{scenario.id}}{The ID of the scenario}
    \item{\code{seed.val}}{Seed Value}
    \item{\code{numsims}}{Number of runs for each simulation}
    \item{\code{numcases}}{Number of cases under binary outcome}
    \item{\code{numcontrols}}{Number of controls under binary outcome}
    \item{\code{numsubjects}}{Number of subjects under continuous outcome}
    \item{\code{interaction.OR}}{Odds ration of the interaction}
    \item{\code{interaction.efkt}}{Interaction effect}
    \item{\code{RR.5.95}}{The baseline odds ratio for subjects on 95 percent population centile versus 5 percentile. This parameter reflects the heterogeneity in disease risk arising from determinants that have not been measured or have not been included in the model. If this parameter is set to 10, the implication is that a high risk subject (someone at the upper 95 percent centile of population risk) is, all else being equal, at 10 times the odds of developing disease compared to someone else who is at low risk (individual at the lower 5 percent centile of population risk)}
    \item{\code{p.val}}{Cut-off p-value defining statistical significance}
    \item{\code{power}}{Desired power}
  }
}


\examples{data(simulation.params)}
