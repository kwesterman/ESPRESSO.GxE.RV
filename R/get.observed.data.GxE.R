#'
#' @title Generates exposure data with some error
#' @description Uses functions make.obs.geno and make.obs.env to generate effect data with a set level of error
#' @param data Input table of simulated data considered as true data
#' @param g.error Misclassification rates in genetic assessment: 1-sensitivity and 1-specificity
#' @param g.model Genetic model; 0 for binary and 1 for additive
#' @param freq Minor allele frequency
#' @param e.error Misclassification rates in environmental exposures assessment: 1-sensitivity and 1-specificity
#' @param e.model Model of the exposure: binary=0, quantitative-normal=1 or quantitative-uniform=2
#' @param e.prev Prevalence of environmental exposure
#' @param e.sd Standard Deviation
#' @param e.reliability Reliability of the assessment of quantitative exposure
#' @return A matrix
#' @keywords internal
#' @author Gaye A.; Westerman K.
#'
get.observed.data.GxE <- function(
    data=NULL,g.error=NULL,g.model=NULL,freq=NULL,
    e.error=NULL,e.model=NULL,e.prev=NULL,e.sd=NULL,
    e.reliability=NULL
) {
    
    sim.df <- data
    
    # # GET THE OBSERVED GENOTYPES
    # # KW comment: ESPRESSO.GxE.RV won't add any genotyping error for now
    # for (genotype in grep("genotype", names(sim.df), value=TRUE)) {
    #     true.genotype <- sim.df[[genotype]]
    #     obs.genotype <- get.obs.geno(allele.A=sim.df$allele.A,allele.B=sim.df$allele.B,
    #                                  geno.model=g.model,MAF=freq,geno.error=g.error)
    # }
    obs.genotype <- sim.df[, grepl("genotype", names(sim.df)), drop=FALSE]
    
    # GET THE OBSERVED ENVIRONMENTAL EXPOSURE DATA
    true.environment <- sim.df$environment
		
    obs.environment <- get.obs.env(env.data=true.environment,env.model=e.model,env.prev=e.prev,
                                   env.sd=e.sd,env.error=e.error,env.reliability=e.reliability)
    
    # GET THE OBSERVED INTERACTION DATA
    obs.interaction <- obs.genotype * obs.environment
    colnames(obs.interaction) <- gsub("genotype", "interaction",
                                      colnames(obs.interaction))

    # REPLACE THE TRUE DATA BY THE NOW GENERATED OBSERVED GENOTYPES
    # IN THE INITIAL MATRIX THAT HELD THE TRUE DATA
    # sim.df$genotype <- obs.genotype$observed.genotype
    # sim.df$allele.A <- obs.genotype$observed.allele.A
    # sim.df$allele.B <- obs.genotype$observed.allele.B
    sim.df$environment <- obs.environment
    for (field in names(obs.interaction)) {
        sim.df[[field]] <- obs.interaction[[field]]
    }
    # sim.df$interaction <- obs.interaction
    
    # RETURN THE MATRIX WHICH NOW CONTAINS ONLY THE OBSERVED DATA TO ANALYSE BY GLM
    # colnames(sim.df) <- c("id", "phenotype", "genotype", "allele.A", "allele.B", "environment", "interaction")
    return(sim.df)
}

