#'
#' @title Summarizes the main results
#' @description Gets the number of cases and controls or subjects and the empirical and theoretical power under each model and prints a summary on the screen
#' @param scenario Scenario number
#' @param pheno.model Type of the outcome; 0 for binary and 2 for continuous
#' @param geno.model Genetic model; 0 for binary and 1 for additive
#' @param env.model Model of the enviromental explosure
#' @param empirical.power Estimated empirical power
##' @param sample.sizes.required Number of cases and controls or number of subjects required to achieve the desired power
#' @return A table containing the following items:
#' \code{genetic.model} Model of the genetic determinant
#' \code{environment.model} Model of the environmental determinant
##' \code{number.of.subjects.required} Number of subjects required to achieve the desired power under a quantitative outcome model
#' \code{empirical.power} Estimated empirical power under each model
#' @keywords internal
#' @author Gaye A.; Westerman K.
#'
get.critical.results.GxE <- function(
  scenario=NULL, pheno.model=NULL,geno.model=NULL,env.model=NULL,
  empirical.power=NULL#,sample.sizes.required=NULL
  # modelled.power=NULL,mean.beta=NULL)
) {
		 
 	 if(geno.model==0){
 	   g.model <- "binary"
 	 }else{
 	   g.model <- "additive"
 	 }
  
   if(env.model==0){
 	   e.model <- "binary"
 	 }else{
 	   if(env.model==1){
 	     e.model <- "normally distributed"
 	   }else{
 	     e.model <- "uniformly distributed"
 	   }
 	 }
    #if(pheno.model == 0){
    #    				numcases <- sample.sizes.required[[1]]
    #    			 numcontrols <- sample.sizes.required[[2]]
    #}else{
    #    				numsubjects <- sample.sizes.required[[1]]
    #}

  if(pheno.model==0) {
     # estimated ORs
     # estimated.OR <- exp(mean.beta)
     # estimated.effect <- 'NA'

					cat("\n---- SUMMARY OF SCENARIO",scenario,"----\n")
					cat("\nModels\n")
					cat("------\n")
					cat("  Outcome: binary \n")
					cat("  Genetic determinant:",g.model,"\n")
					cat("  Environmental determinant:",e.model)

					#cat("\n\nNumber of cases required\n")
					#cat("------------------------\n")
					#cat(" ", numcases)
					#
					#cat("\n\nNumber of controls required\n")
					#cat("---------------------------\n")
					#cat(" ", numcontrols)

					cat("\n\nEmpirical power\n")
					cat("---------------\n")
					cat(" ",empirical.power)

					# cat("\n\nModel power\n")
					# cat("-----------\n")
					# cat(" ",round(modelled.power,2))

					# cat("\n\nEstimated interaction OR\n")
					# cat("-----------\n")
					# cat(" ",round(estimated.OR,2))

					cat("\n\n---- END OF SUMMARY ----\n")

		   # crit.res <- c(g.model,e.model,numcases,numcontrols,round(empirical.power,2),round(modelled.power,2),
		   #               round(estimated.OR,2), estimated.effect)
		   # return(list(genetic.model=crit.res[1], environment.model=crit.res[2],number.of.cases.required=crit.res[3],
		   #             number.of.controls.required=crit.res[4],empirical.power=crit.res[5], 
		   #             modelled.power=crit.res[6],estimated.OR=crit.res[7], estimated.effect=crit.res[8]))
		   return(list(genetic.model=g.model, environment.model=e.model, #numer.of.subjects.required=numcases+numcontrols,
		               empirical.power=round(empirical.power,2)))

  } else {
     # estimated ORs
     # estimated.effect <- mean.beta
     # estimated.OR <- 'NA'

					cat("\n---- SUMMARY OF SCENARIO",scenario,"----\n")
					cat("\nModels\n")
					cat("------\n")
					cat(" Outcome: quantitative; ")
					cat(" Genetic determinant: ",g.model,"\n")
					cat(" Environmental determinant: ",e.model)

					#cat("\n\nNumber of subjects required\n")
					#cat("------------------------\n")
					#cat(" ",numsubjects)

					cat("\n\nEmpirical power\n")
					cat("---------------\n")
					cat(" ",empirical.power)

					# cat("\n\nModel power\n")
					# cat("-----------\n")
					# cat(" ",round(modelled.power,2))
					# 
					# cat("\n\nEstimated interaction effect\n")
					# cat("-----------\n")
					# cat(" ",estimated.effect)

					cat("\n\n---- END OF SUMMARY ----\n")

		   # crit.res <- c(g.model,e.model,numsubjects,round(empirical.power,2),round(modelled.power,2),estimated.OR,round(estimated.effect,2))
		   # return(list(genetic.model=crit.res[1], environment.model=crit.res[2],number.of.subjects.required=crit.res[3],
		   #            empirical.power=crit.res[4], modelled.power=crit.res[5], estimated.OR=crit.res[6],estimated.effect=crit.res[7]))
    return(list(genetic.model=g.model, environment.model=e.model, #number.of.subjects.required=numsubjects,
                empirical.power=round(empirical.power,2)))
  }
}

