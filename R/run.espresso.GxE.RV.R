#' 
#' @title Runs a full ESPRESSO.GxE.RV analysis
#' @description This function calls the functions required to run a full ESPRESSO.GxE.RV analysis 
#'  where the model consists of an outcome (binary or continuous) determined by two interacting
#'  covariates (a SNP  and an environmental exposure)
#' @param simulation.params general parameters for the scenario(s) to analyse
#' @param pheno.params paramaters for the outcome variables
#' @param geno.params parameters for the genetic determinant(s)
#' @param env.params parameters for the environmental determinant
#' @param scenarios2run the indices of the scenarios one wish to analyse
#' @return a summary table that contains both the input parameters and 
#' the results of the analysis
#' @export
#' @author Gaye A.; Westerman K.
#' @examples {
#'   
#' # load the table that hold the input parameters; each of the table
#' # hold parameters for 4 scenarios:
#' # scenario 1: a binary outcome determined by a binary SNP and binary exposure 
#' # and an interaction between the the genetic variant and the environmental exposure
#' # scenario 2: a binary outcome determined by an additive SNP and continuous 
#' # exposure and an interaction between the the genetic variant and the environmental exposure
#' # scenario 3: a quantitative outcome determined by a binary SNP and binary exposure, 
#' # and an interaction between the the genetic variant and the environmental exposure
#' # scenario 4: a quantitative outcome determined by an additive SNP and continuous 
#' # exposure and an interaction between the the genetic variant and the environmental exposure
#' data(simulation.params) 
#' data(pheno.params)
#' data(geno.params)
#' data(env.params)
#' 
#' # run the function for the first two scenarios, two binomial models
#' run.espresso.GxE.RV(simulation.params, pheno.params, geno.params, env.params, scenarios2run=c(1,2))
#'
#' # run the function for the last two scenarios, two gaussian models
#' run.espresso.GxE.RV(simulation.params, pheno.params, geno.params, env.params, scenarios2run=c(3,4))
#' 
#' }
#'
run.espresso.GxE.RV <- function(simulation.params=NULL, pheno.params=NULL, geno.params=NULL, 
                                env.params=NULL, scenarios2run=1) {
   
   # IF AN INPUT FILE IS NOT SUPPLIED LOAD THE DEFAULT TABLES WARNING
   if (is.null(simulation.params)) {
      cat("\n WARNING!\n")
      cat(" No simulation parameters supplied\n")
      cat(" The default simulation parameters will be used\n")
      simulation.params <- data(simulation.params)
   }
   
   if (is.null(pheno.params)) {
      cat("\n WARNING!\n")
      cat(" No outcome parameters supplied\n")
      cat(" The default outcome parameters will be used\n")
      pheno.params <- data(pheno.params)
   }
   
   if (is.null(geno.params)) {
      cat("\n WARNING!\n")
      cat(" No genotype parameters supplied\n")
      cat(" The default genotype parameters will be used\n")
      geno.params <- data(geno.params)
   }
   
   if (is.null(env.params)) {
      cat("\n WARNING!\n")
      cat(" No environmental parameters supplied\n")
      cat(" The default environmental parameters will be used\n")
      env.params <- data(env.params)
   }
   
   # MERGE INPUT FILES TO MAKE ONE TABLE OF PARAMETERS
   s.temp1 <- merge(simulation.params, pheno.params)
   s.temp2 <- merge(s.temp1, geno.params)
   s.parameters <- merge(s.temp2, env.params)
   
   
   #----------LOAD SET UP UP INITIAL PARAMETERS------------#
   
   # PRINT TRACER CODE EVERY Nth ITERATION
   # THIS ENSURES THAT YOU CAN SEE IF THE PROGRAM GRINDS TO A HALT FOR SOME REASON (IT SHOULDN'T)
   trace.interval <- 10
   
   
   # CREATE UP TO 20M SUBJECTS IN BLOCKS OF 20K UNTIL REQUIRED NUMBER OF
   # CASES AND CONTROLS IS ACHIEVED. IN GENERAL THE ONLY PROBLEM IN ACHIEVING THE
   # REQUIRED NUMBER OF CASES WILL OCCUR IF THE DISEASE PREVALENCE IS VERY LOW
   allowed.sample.size <- 20000000
   block.size <- 20000
   
   
   # DECLARE MATRIX THAT STORE THE RESULTS FOR EACH SCENARIO (ONE PER SCENARIO PER ROW)
   output.file <- "output.RV.csv"
   output.matrix <- matrix(numeric(0), ncol=ncol(s.parameters) + 2)
   column.names <- c(colnames(s.parameters), "exceeded.sample.size?", "empirical.power")#, "numsubjects.required")
   write(t(column.names), output.file, dim(output.matrix)[2], sep=";")
   
   #-----------LOOP THROUGH THE SCENARIOS - DEALS WITH ONE SCENARIO AT A TIME-------------
   
   for (j in scenarios2run) {
      
      # RANDOM NUMBER GENERATOR STARTS WITH SEED SET AS SPECIFIED 
      set.seed(s.parameters$seed.val[j])
      
      # SIMULATION PARAMETERS
      scenario.id <- s.parameters$scenario.id[j]         
      seed.val <- s.parameters$seed.val[j]               
      numsims <- s.parameters$numsims[j]                 
      numcases <- s.parameters$numcases[j]               
      numcontrols <- s.parameters$numcontrols[j]  
      numsubjects <- s.parameters$numsubjects[j]
      int.OR <- s.parameters$interaction.OR[j]
      int.efkt <- s.parameters$interaction.efkt[j]            
      baseline.OR <- s.parameters$RR.5.95[j]                
      pval <- s.parameters$p.val[j]                      
      power <- s.parameters$power[j]
      
      # OUTCOME PARAMETERS
      pheno.model <- s.parameters$pheno.model[j]
      pheno.mean <- s.parameters$pheno.mean[j]
      pheno.sd <- s.parameters$pheno.sd[j]
      disease.prev <- s.parameters$disease.prev[j]
      pheno.error <- c(1-s.parameters$pheno.sensitivity[j], 1-s.parameters$pheno.specificity[j])
      pheno.reliability <- s.parameters$pheno.reliability[j]    
      
      # GENETIC DETERMINANT PARAMETERS
      geno.model <- s.parameters$geno.model[j]
      MAF <-  s.parameters$MAF[j]    
      geno.OR <- s.parameters$geno.OR[j]
      geno.efkt <- s.parameters$geno.efkt[j]
      geno.error <- c(1-s.parameters$geno.sensitivity[j], 1-s.parameters$geno.specificity[j])
      geno.nvariants <- s.parameters$geno.nvariants[j]
      geno.frac_causal <- s.parameters$geno.frac_causal[j]
      
      # ENVIRONMENTAL DETERMINANT PARAMETERS
      env.model <- s.parameters$env.model[j]
      env.prev <-  s.parameters$env.prevalence[j]    
      env.OR <- s.parameters$env.OR[j]
      env.efkt <- s.parameters$env.efkt[j]
      env.mean <- s.parameters$env.mean[j]
      env.sd <- s.parameters$env.sd[j]
      env.low.lim <- s.parameters$env.low.lim[j]
      env.up.lim <- s.parameters$env.up.lim[j]
      env.error <- c(1-s.parameters$env.sensitivity[j], 1-s.parameters$env.specificity[j])
      env.reliability <- s.parameters$env.reliability[j]
      
      # VECTORS TO HOLD BETA, SE AND Z VALUES AFTER EACH RUN OF THE SIMULATION
      ##beta.values <- rep(NA, numsims)
      ##se.values <- rep(NA, numsims)
      ##z.values<-rep(NA, numsims)
      p.values <- rep(NA, numsims)
      
      # TRACER TO DETECT EXCEEDING MAX ALLOWABLE SAMPLE SIZE
      sample.size.excess <- 0
      
      # GENERATE AND ANALYSE DATASETS ONE AT A TIME 
      for(s in 1:numsims) {  # s from 1 to total number of simulations
         
         #----------------------------------GENERATE "TRUE" DATA-----------------------------#
         
         if (pheno.model == 0) { # UNDER BINARY OUTCOME MODEL
            # GENERATE CASES AND CONTROLS UNTILL THE REQUIRED NUMBER OF CASES, CONTROLS IS ACHIEVED 
            sim.data <- sim.CC.data.GxE(n=block.size, ncases=numcases, ncontrols=numcontrols, 
                                        max.sample.size=allowed.sample.size, pheno.prev=disease.prev,
                                        g.nvar=geno.nvariants, g.frac_causal=geno.frac_causal,
                                        freq=MAF, g.model=geno.model, g.OR=geno.OR, e.model=env.model, 
                                        e.prev=env.prev, e.mean=env.mean, e.sd=env.sd, e.low.lim=env.low.lim, 
                                        e.up.lim=env.up.lim, e.OR=env.OR, i.OR=int.OR, b.OR=baseline.OR, 
                                        ph.error=pheno.error)
            
            true.data <- sim.data$data
            
         } else { # UNDER QUANTITATIVE OUTCOME MODEL
            # GENERATE THE SPECIFIED NUMBER OF SUBJECTS
            true.data <- sim.QTL.data.GxE(n=numsubjects, ph.mean=pheno.mean, ph.sd=pheno.sd, 
                                          g.nvar=geno.nvariants, g.frac_causal=geno.frac_causal, freq=MAF,
                                          g.model=geno.model, g.efkt=geno.efkt, e.model=env.model,
                                          e.efkt=env.efkt, e.prev=env.prev, e.mean=env.mean, e.sd=env.sd,
                                          e.low.lim=env.low.lim, e.up.lim=env.up.lim, i.efkt=int.efkt, 
                                          pheno.rel=pheno.reliability)
         }
         
         #------------SIMULATE ERRORS AND ADD THEM TO THE TRUE COVARIATES DATA TO OBTAIN OBSERVED COVARIATES DATA-----------#
         
         # ADD APPROPRIATE ERRORS TO PRODUCE OBSERVED GENOTYPES 
         observed.data <- get.observed.data.GxE(data=true.data, g.error=geno.error, g.model=geno.model, freq=MAF, 
                                                e.error=env.error, e.model=env.model, e.prev=env.prev, e.sd=env.sd, 
                                                e.reliability=env.reliability)
         
         #--------------------------DATA ANALYSIS ----------------------------#
         
         variant_results_list <- list()
         for (g.idx in 1:geno.nvariants) {
            glm.estimates.variant <- glm.analysis.GxE(pheno.model, observed.data, g.idx)
            variant_results_list[[g.idx]] <- glm.estimates.variant
         }
         p.combined <- combine.variants(variant_results_list)
         
         # beta.values[s] <- 100 #glm.estimates[[1]]
         # se.values[s] <- 100 #glm.estimates[[2]]
         # z.values[s] <- 100 #glm.estimates[[3]]
         p.values[s] <- p.combined
         
         # PRINT TRACER AFTER EVERY Nth DATASET CREATED
         if (s %% trace.interval == 0)cat("\n", s, "of", numsims, "runs completed in scenario", scenario.id)
         
      }
      cat("\n\n")
      
      #------------------------ SUMMARISE RESULTS ACROSS ALL SIMULATIONS---------------------------#
      
      # # SUMMARISE PRIMARY PARAMETER ESTIMATES
      # # COEFFICIENTS ON LOG-ODDS SCALE
      # mean.beta <- round(mean(beta.values, na.rm=T), 3)
      # mean.se <- round(sqrt(mean(se.values^2, na.rm=T)), 3)
      # mean.model.z <- mean.beta/mean.se
      
      
      #---------------------------POWER AND SAMPLE SIZE CALCULATIONS----------------------#
      
      # CALCULATE EMPIRICAL POWER AND THE MODELLED POWER 
      # THE EMPIRICAL POWER IS SIMPLY THE PROPORTION OF SIMULATIONS IN WHICH
      # THE Z STATISTIC FOR THE PARAMETER OF INTEREST EXCEEDS THE Z STATISTIC
      # FOR THE DESIRED LEVEL OF STATISTICAL SIGNIFICANCE
      empirical.power <- power.calc(pval, p.values)

      ## CALCULATE THE SAMPLE SIZE REQUIRED UNDER EACH MODEL
      #mean.model.z <- mean(qnorm(1 - p.values / 2))  # This is a rough approximation of the mean z (precisely: mean of back-transformed z-values after combining multi-genotype p-values)
      #sample.sizes.required <- samplsize.calc(numcases, numcontrols, numsubjects, pheno.model, pval, power, mean.model.z)
      
      #------------------MAKE FINAL A TABLE THAT HOLDS BOTH INPUT PARAMETERS AND OUTPUT RESULTS---------------#
      
      critical.res <- get.critical.results.GxE(
         j, pheno.model, geno.model, env.model, empirical.power#, sample.sizes.required
      )
      # power$modelled, mean.beta)
      
      #  WHEN OUTCOME IS BINARY INFORM IF RECORD EXCEEDED MAXIMUM SAMPLE SIZE
      if (pheno.model == 0) {
         sample.size.excess <- sim.data$allowed.sample.size.exceeded
         if (sample.size.excess == 1) {
            excess <- "yes"
            cat("\nTO GENERATE THE NUMBER OF CASES SPECIFIED AT OUTSET\n")
            cat("THE SIMULATION EXCEEDED THE MAXIMUM POPULATION SIZE OF ", allowed.sample.size, "\n")
         } else {
            excess <- "no"
         }
      }
      inparams <- s.parameters[j, ]
      cat(colnames(inparams))
      # all_cols <- c(
      #    "scenario.id", "seed.val", "numsims", "numcases", "numcontrols", 
      #    "numsubjects", "interaction.OR", "interaction.efkt", "p.val", "power", 
      #    "pheno.model", "disease.prev", "pheno.mean", "pheno.sd", 
      #    "pheno.sensitivity", "pheno.specificity", "pheno.reliability", 
      #    "geno.model", "MAF", "geno.efkt", "geno.sensitivity", 
      #    "geno.specificity", "env.model", "env.prevalence", "env.efkt", 
      #    "env.mean", "env.sd", "env.low.lim", "env.up.lim", "env.reliability", 
      #    "exceeded.sample.size?", "numcases.required", "numcontrols.required", 
      #    "numsubjects.required", "empirical.power"
      #    # "modelled.power", "estimated.OR", "estimated.effect"
      # )
      if (pheno.model == 0) {
         # mod <- "binary"
         if (env.model == 0) {
            irrelevant_params <- c(
               "geno.efkt", "env.efkt", "interaction.efkt", 
               "env.mean", "env.sd",
               "env.low.lim", "env.up.lim", "env.reliability"
            )
            inparams[irrelevant_params] <- "NA"
            inputs <- inparams
         } else if (env.model == 1) {
            irrelevant_params <- c(
               "geno.efkt", "env.efkt", "interaction.efkt", 
               "env.prevalence", "env.OR", "env.sensitivity", "env.specificity"
            )
            inputs <- inparams
         } 
         # else {
         #    inparams [c(6,8,18,22,27,29:31,34,35)] <- "NA"
         #    inputs <- inparams          
         # }
         outputs <- c(sample.size.excess, critical.res$empirical.power)#, critical.res$number.of.subjects.required)
      } else {
         mod <- "quantitative"
         if (env.model == 0) {
            irrelevant_params <- c(
               "geno.OR", "env.OR", "interaction.OR",
               "env.mean", "env.sd",
               "env.low.lim", "env.up.lim", "env.reliability"
            )
            inputs <- inparams
         } else if (env.model == 1) {
            irrelevant_params <- c(
               "geno.OR", "env.OR", "interaction.OR",
               "env.prevalence", "env.OR", "env.sensitivity", "env.specificity"
            )
            inputs <- inparams
            # else {
            #    inparams [c(4,5,7,15:17,21,26,27,29,30,33,35)] <- "NA"
            #    inputs <- inparams          
            # }
         }
         outputs <- c(sample.size.excess, critical.res$empirical.power)#, critical.res$number.of.subjects.required)
      }
      jth.row <- as.character(c(inputs, outputs))
      write(t(jth.row), output.file, dim(output.matrix)[2], append=TRUE, sep=";")
   }
}
