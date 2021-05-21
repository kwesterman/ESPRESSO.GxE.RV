#'
#' @title Simulates subjects for continuous outcome
#' @description Generates the specified number of subjects
#' @param n Number of subjects to simulate
#' @param ph.mean statistical mean
#' @param ph.sd standard deviation
#' @param g.nvar Integer number of genetic variants to simulate
#' @param g.frac_causal Fraction of causal genetic variants
#' @param freq Minor allele frequencies of the genetic variants
#' @param g.model Genetic model; 0 for binary and 1 for continuous
#' @param g.efkt Effects of the genetic variant
#' @param e.model Model of the environmental exposure
#' @param e.efkt Effects of the environment determinats
#' @param e.prev Prevalences of the environmental determinants
#' @param e.mean Mean under quantitative-normal model
#' @param e.sd Standard deviation under quantitative-normal model
#' @param e.low.lim lower limit under quantitative-uniform model
#' @param e.up.lim upper limit under quantitative-uniform model
#' @param i.efkt Interaction effect
#' @param pheno.rel reliability of the assessment for a quantitative outcome.
#' @return A matrix
#' @keywords internal
#' @author Gaye A.; Westerman K.

sim.QTL.data.GxE <- function(
   n=NULL,ph.mean=NULL,ph.sd=NULL,
   g.nvar=NULL,g.frac_causal=NULL,freq=NULL,g.model=NULL,g.efkt=NULL,
   e.model=NULL,e.efkt=NULL,e.prev=NULL,e.mean=NULL,e.sd=NULL,e.low.lim=NULL,
   e.up.lim=NULL,i.efkt=NULL,pheno.rel=NULL
) {
   
   # GENERATE THE TRUE GENOTYPE DATA
   geno <- matrix(nrow=n, ncol=g.nvar)
   for (j in 1:g.nvar) {
      geno.data <- sim.geno.data(num.obs=n, geno.model=g.model, MAF=freq)
      allele.A <- geno.data$allele.A  # Is this info ultimately needed? Can it be the same across all variants?
      allele.B <- geno.data$allele.B
      geno[, j] <- geno.data$genotype
   }
   
   # GENERATE THE TRUE ENVIRONMENTAL EXPOSURE DATA			
   env <- sim.env.data(num.obs=n,env.model=e.model,env.prev=e.prev,env.mean=e.mean,
                       env.sd=e.sd,env.low.lim=e.low.lim,env.up.lim=e.up.lim)
   
   # GENERATE THE TRUE INTERACTION DATA           
   interaction <- geno * env  # N x g.nvar matrix
   
   # GENERATE THE TRUE OUTCOME DATA
   pheno.data <- sim.pheno.qtl.GxE(numsubjects=n,pheno.mean=ph.mean,pheno.sd=ph.sd,
                                   genotype=geno,geno.efkt=g.efkt,environment=env,env.efkt=e.efkt,
                                   interaction=interaction,int.efkt=i.efkt,
                                   geno.frac_causal=g.frac_causal)
   true.phenotype <- pheno.data
   
   # GENERATE THE OBSERVED OUTCOME DATA 
   obs.phenotype <- get.obs.pheno(phenotype=true.phenotype, pheno.model=1, 
                                  pheno.sd=ph.sd, pheno.reliability=pheno.rel)
   pheno <- obs.phenotype
   
   # STORE THE GENERATED TRUE DATA INTO AN OUTPUT MATRIX 
   sim.matrix <- cbind(pheno,env,geno,interaction,allele.A,allele.B)
   
   # ADD IDs (JUST A ROW COUNT)
   sim.matrix <- cbind(1:nrow(sim.matrix), sim.matrix)
   
   # ADD COLUMN NAMES AND RETURN A DATAFRAME
   colnames(sim.matrix) <- c("id","phenotype","environment",
                             paste0("genotype", 1:ncol(geno)),
                             paste0("interaction", 1:ncol(interaction)),
                             "allele.A","allele.B")
   mm <- data.frame(sim.matrix)
}

