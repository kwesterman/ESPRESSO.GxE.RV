#' 
#' @title Simulates continuous outcome data
#' @description Uses the effects data of the determinants to construct a linear predictor(LP). The outcome is normally distributed variable generated with a mean equal to LP and a standard deviation of 1. Some error is then added to the simulated outcome to obtained the observed outcome.
#' @param num.subjects Number of subjects to simulate
#' @param pheno.mean statistical mean
#' @param pheno.sd standard deviation
#' @param genotype Genotype matrix (N x g.nvar)
#' @param geno.efkt Effect of the genetic variant
#' @param environment Exposure data for environment
#' @param env.efkt Effect of the environmental determinants
#' @param interaction Interaction matrix (N x g.nvar)
#' @param int.efkt Interaction effect
#' @return A dataframe of phenotype
#' @keywords internal
#' @author Gaye A.; Westerman K.
#'
sim.pheno.qtl.GxE <- function(
   numsubjects=NULL,pheno.mean=NULL,pheno.sd=NULL,
   genotype=NULL,geno.efkt=NULL,geno.frac_causal=NULL,
   environment=NULL,env.efkt=NULL,interaction=NULL,int.efkt=NULL
) {  

  # ALPHA IS EQUAL TO THE MEAN OF THE TRAIT, WHICH IS 0
   num.obs <- numsubjects
   alpha <- pheno.mean
   
   # SET FINAL VARIANT EFFECTS BASED ON THE CAUSAL FRACTION
   num_causal <- ceiling(geno.frac_causal * ncol(genotype))  # Round up
   geno.efkts <- c(rep(geno.efkt, num_causal), 
                   rep(0, ncol(genotype) - num_causal))
   int.efkts <- c(rep(int.efkt, num_causal), 
                  rep(0, ncol(genotype) - num_causal))

   # GENERATE THE LINEAR PREDICTOR
   lp <- alpha + 
      (genotype %*% geno.efkts) + 
      (environment * env.efkt) + 
      (interaction %*% int.efkts)

   # GENERATE THE TRUE PHENOTYPE DATA
   phenotype <- rnorm(num.obs, mean=lp, sd=sqrt(pheno.sd^2 - var(lp)))
   
   return(phenotype)
}

