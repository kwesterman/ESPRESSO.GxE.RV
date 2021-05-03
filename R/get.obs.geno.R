#' 
#' @title Adds some error to genotype data
#' @description Simulates errors and adds it to the true data to obtain observed data. 
#'  The alleles simulated by the function sim.geno.data are randomly misclassified and used to 
#'  form new genotypes that represent the observed genotypes.
#' @param allele.A Allele A
#' @param allele.B Allele B
#' @param geno.model genetic model; binary=0 and additive=1
#' @param MAF minor allele frequency of the SNP (in ESPRESSO this is the frequency of the 'at risk' allele)
#' @param geno.error a vector with two values, the misclassification rates related to the sensitivity and 
#' specificity of the assessment of the alleles, i.e. 1-sensitivity and 1-specificity
#' @return a dataframe that contains the below data:
#' \code{observed.genotype} observed genotypes
#' \code{observed.allele.A} observed A alleles
#' \code{observed.allele.B} observed B alleles
#' @keywords internal
#' @author Gaye A.; Westerman K.
#' 
get.obs.geno <- function (allele.A=NULL, allele.B=NULL, geno.model=NULL, MAF=NULL, geno.error=NULL){
  
  mean.add <- (2*MAF*(1-MAF)+2*(MAF^2))
  mean.bin <- (2*MAF*(1-MAF)+(MAF^2))
  
  numsubs <- length(allele.A)
  allele.A.1.missed <- rbinom(numsubs, 1, geno.error[1])
  allele.A.0.missed <- rbinom(numsubs, 1, geno.error[2])
  allele.B.1.missed <- rbinom(numsubs, 1, geno.error[1])
  allele.B.0.missed <- rbinom(numsubs, 1, geno.error[2])
  
  observed.allele.A <-
    ((allele.A==0) * (allele.A.1.missed==0) * (allele.A.0.missed==0) * allele.A)+
    ((allele.A==1) * (allele.A.1.missed==0) * (allele.A.0.missed==0) * allele.A)+
    ((allele.A==0) * (allele.A.1.missed==0) * (allele.A.0.missed==1) * (1-allele.A))+
    ((allele.A==1) * (allele.A.1.missed==0) * (allele.A.0.missed==1) * allele.A)+
    ((allele.A==0) * (allele.A.1.missed==1) * (allele.A.0.missed==0) * allele.A)+
    ((allele.A==1) * (allele.A.1.missed==1) * (allele.A.0.missed==0) * (1-allele.A))+
    ((allele.A==0) * (allele.A.1.missed==1) * (allele.A.0.missed==1) * (1-allele.A))+
    ((allele.A==1) * (allele.A.1.missed==1) * (allele.A.0.missed==1) * (1-allele.A))
  
  observed.allele.B <-
    ((allele.B==0) * (allele.B.1.missed==0) * (allele.B.0.missed==0) * allele.B)+
    ((allele.B==1) * (allele.B.1.missed==0) * (allele.B.0.missed==0) * allele.B)+
    ((allele.B==0) * (allele.B.1.missed==0) * (allele.B.0.missed==1) * (1-allele.B))+
    ((allele.B==1) * (allele.B.1.missed==0) * (allele.B.0.missed==1) * allele.B)+
    ((allele.B==0) * (allele.B.1.missed==1) * (allele.B.0.missed==0) * allele.B)+
    ((allele.B==1) * (allele.B.1.missed==1) * (allele.B.0.missed==0) * (1-allele.B))+
    ((allele.B==0) * (allele.B.1.missed==1) * (allele.B.0.missed==1) * (1-allele.B))+
    ((allele.B==1) * (allele.B.1.missed==1) * (allele.B.0.missed==1) * (1-allele.B))
  
  observed.genotype <- observed.allele.A + observed.allele.B
  
  if(geno.model==0){
    genotyp.U <- observed.genotype > 0
    observed.genotype <- genotyp.U - mean.bin
  }
  if(geno.model==1){
    genotyp.U <- observed.genotype
    observed.genotype <- genotyp.U - mean.add
  }
  
  df <- data.frame(observed.genotype, observed.allele.A, observed.allele.B)
}

