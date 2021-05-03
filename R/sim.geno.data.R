#' 
#' @title Generates genotypes for a genetic variant
#' @description Generates two alleles and combines them to form the genotype of 
#' a SNP under a binary or additive genetic model.
#' @param num.obs number of observations to generate.
#' @param geno.model genetic model; binary=0 and additive=1.
#' @param MAF minor allele frequency of the SNP (in ESPRESSO this is the frequency of the 'at risk' allele)
#' @return a dataframe that contains the following variables:
#' \code{allele.A} major allele
#' \code{allele.B} minor allele
#' \code{genotype} genotype
#' @keywords internal
#' @author Gaye A.
#' 
sim.geno.data <- function(num.obs=NULL, geno.model=NULL, MAF=NULL){
  

  # CORRECTION TERM FOR MEAN CENTERING FOR ADDITIVE 
  mean.geno.additive <- (2*MAF*(1-MAF)+2*(MAF^2))

  # CORRECTION TERM FOR MEAN CENTERING FOR BINARY GENE
  mean.geno.binary <- (2*MAF*(1-MAF)+(MAF^2))

  # CREATE, CENTRE AND ROUND AN ADDITIVE GENETIC GENOTYPE COVARIATE WITH
  # APPROPRIATE MAF 
  allele.A <- rbinom(num.obs,1,MAF)
  allele.B <- rbinom(num.obs,1,MAF)

  # ACTUAL GENOTYPE IS SUM OF ALLELES IN ADDITIVE GENETIC MODEL
  # AND IS 1 IF SUM OF ALLELES IS 1 OR GREATER IN THE BINARY MODEL 
  genotype <- allele.A+allele.B

		if(geno.model==0){
				geno.U <- genotype > 0
				genotype <- geno.U - mean.geno.binary
		}
		if(geno.model==1){
				geno.U <- genotype
				genotype <- geno.U - mean.geno.additive
		}

  # return the data as a dataframe
  dtframe <- data.frame(allele.A, allele.B, genotype)
}

