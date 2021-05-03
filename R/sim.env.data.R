#' 
#' @title Simulates cases and controls
#' @description Generates data for a binary, quantitative-normal, or quantitative-uniform environmental determinant.
#' @param num.obs number of observations to simulate.
#' @param env.model model of the exposure: binary=0, quantitative-normal=1, quantitative-uniform=2.
#' @param env.prev prevalence of the environmental exposure.
#' @param env.mean statisitical man under quantitative-normal model.
#' @param env.sd standard deviation under quantitative-normal model.
#' @param env.low.lim lower limit under quantitative-uniform model.
#' @param env.up.lim upper limit under quantitative-uniform model.
#' @return a vector of continuous or binary values.
#' @keywords internal
#' @author Gaye A.
#' 
sim.env.data <- function(num.obs=NULL, env.model=NULL, env.prev=NULL, env.mean=NULL,
                         env.sd=NULL,env.low.lim=NULL,env.up.lim=NULL){

  # CREATE THE FIRST ENVIRONMENTAL COVARIATE 
  if(env.model==0){   # BINARY DISTRIBUTION
     env.U <- rbinom(num.obs, 1, env.prev)
     e.mean <- env.prev
     env.U <- env.U-e.mean
  }
  if(env.model==1){   # NORMAL DISTRIBUTION
     env.U <- rnorm(num.obs, env.mean, env.sd)
     env.U <- env.U-mean(env.U)     # mean centering
  }
  if(env.model==2){  # UNIFORM DISTRIBUTION
     if(env.low.lim >= env.up.lim){
       stop("\n\nALERT!\n Uniform Distribution: The upper limit must be greater than the lower limit\n\n")
     }else{
       env.U <- runif(num.obs, env.low.lim, env.up.lim)
       env.U <- env.U-mean(env.U)     # mean centering
     }
  }

  # return a vector 
  return(env.U)
}

