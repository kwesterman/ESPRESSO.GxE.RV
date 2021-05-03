#' 
#' @title Generates observed outcome data
#' @description Adds a set level of error to error free binary or quantitative data (the true phenotype data) 
#' to obtain data with a larger variance (the observed phenotype data).
#' @param seed 
#' @param phenotype outcome status.
#' @param pheno.model distribution of the outcome variable: binary=0, normal=1
#' @param pheno.sd standard deviation of the outcome in the study population
#' @param pheno.model distribution of the outcome variable: binary=0, normal=1 or uniform=2.
#' @param pheno.sd standard deviation of the outcome in the study the population
#' @param pheno.error misclassification rates: 1-sensitivity and 1-specificity
#' @param pheno.reliability reliability of the assessment for a quantitative outcome.
#' @return A dataframe containing:
#' \code{true.phenotype} the error free outcome data (true data).
#' \code{observed.phenotype} the true outcome data with some added error (observed data).
#' @keywords internal
#' @author Gaye A.
#'
get.obs.pheno <- function (phenotype=NULL, pheno.model=NULL, pheno.sd=NULL, pheno.error=NULL, pheno.reliability=NULL){ 
  
  # GET THE OBSERVED OUTCOME DATA
  if(pheno.model==0){ # IF THE OUTCOME IS BINARY
    
    pheno.original<-phenotype
    disease.missed<-rbinom(length(phenotype),1,pheno.error[1])
    non.disease.missed<-rbinom(length(phenotype),1,pheno.error[2])
    pheno.U<-
      ((pheno.original==0) * (disease.missed==0) * (non.disease.missed==0))*pheno.original+
      ((pheno.original==1) * (disease.missed==0) * (non.disease.missed==0))*pheno.original+
      ((pheno.original==0) * (disease.missed==0) * (non.disease.missed==1))*(1-pheno.original)+
      ((pheno.original==1) * (disease.missed==0) * (non.disease.missed==1))*pheno.original+
      ((pheno.original==0) * (disease.missed==1) * (non.disease.missed==0))*pheno.original+
      ((pheno.original==1) * (disease.missed==1) * (non.disease.missed==0))*(1-pheno.original)+
      ((pheno.original==0) * (disease.missed==1) * (non.disease.missed==1))*(1-pheno.original)+
      ((pheno.original==1) * (disease.missed==1) * (non.disease.missed==1))*(1-pheno.original)
    observed.phenotype <- pheno.U
    
  }else{ # IF THE OUTCOME IS CONTINUOUS NORMAL
    # USE THE RELIABITLITY OF PHENOTYPE ASSESSMENT TO COMPUTE THE VARIANCE OF MEASURED PHENOTYPES.
    # RELIABITY = (VAR.true.data/VAR.obs.data) AND VAR.obs.data = (VAR.true.data + VAR.measurement)
    # IT FOLLOWS THAT VAR.true.data + VAR.of.estimate = VAR.true.data/RELIABILITY AND THEN:
    # VAR.measurement = (VAR.true.data/RELIABILITY) - VAR.true.data
    var.m <- (pheno.sd^2/pheno.reliability)-(pheno.sd^2)
    
    # ADD THE ERROR TO ORIGINAL PHENOTYPES TO GENERATE THE OBSERVED PHENOTYPE DATA
    n <- length(phenotype)
    observed.phenotype <- rnorm(n, phenotype, sqrt(var.m))
  } 
  
  # RETURN THE TRUE AND OBSERVED PHENOTYPE DATA AS A DATAFRAME
  df <- observed.phenotype
  
}
