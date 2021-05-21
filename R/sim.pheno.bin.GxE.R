#'
#' @title Generates phenotype status
#' @description Generates affected and non-affected subjects
#' @param num.obs Number of observations to generate per iteration
#' @param disease.prev Prevalence of the binary outcome
#' @param genotype Exposure data for genetic determinant
#' @param environment Exposure data for environment
#' @param interaction data
#' @param subject.effect.data Subject effect data, reflects the heterogenity 
#' in baseline disease risk
#' @param geno.OR Odds ratios of the genetic determinant
#' @param geno.frac_causal Fraction of causal genetic variants
#' @param env.OR Odds ratios of the two environments
#' @param int.OR Odds ration of the interaction
#' @return A dataframe of phenotype
#' @keywords internal
#' @author Gaye A.; Westerman K.
#'
sim.pheno.bin.GxE <- function(
   num.obs=NULL, disease.prev=NULL, genotype=NULL, environment=NULL, interaction=NULL, 
   subject.effect.data=NULL, geno.OR=NULL, geno.frac_causal=NULL, env.OR=NULL, int.OR=NULL
) { 
   # GET THE ALPHA AND BETA VALUES
   alpha <- log(disease.prev/(1-disease.prev))
   geno.beta <- log(geno.OR)
   env.beta <-	log(env.OR)
   int.beta <- log(int.OR)
   
   # SET FINAL VARIANT EFFECTS BASED ON THE CAUSAL FRACTION
   num_causal <- ceiling(geno.frac_causal * ncol(genotype))  # Round up
   geno.betas <- c(rep(geno.beta, num_causal), 
                   rep(0, ncol(genotype) - num_causal))
   int.betas <- c(rep(int.beta, num_causal), 
                  rep(0, ncol(genotype) - num_causal))
   
   # GENERATE THE LINEAR PREDICTOR
   # lp <- alpha + (geno.beta*genotype) + (env.beta*environment) + (int.beta*interaction) + subject.effect.data
   
   lp <- alpha + 
      (genotype %*% geno.betas) + 
      (environment * env.beta) + 
      (interaction %*% int.betas) +
      subject.effect.data
   
   # GET THE 'mu' THE PROBABILITY OF DISEASE THROUGH LOGISTIC TRANSFORMATION
   mu <- exp(lp)/(1 + exp(lp))
   
   # GENERATE THE PHENOTYPE DATA AND RETURN IT AS A DATAFRAME
   phenotype <- rbinom(num.obs,1,mu)
   
   return(phenotype)
}

