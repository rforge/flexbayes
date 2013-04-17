{cat("----- No extraneous info printed while fitting AE model --\n");T}
{
# Functions: runAEmodel
# Description: Make sure that no extraneous info is printed 
# while fitting an adverse event model
  
  # save printed output to a file
  sink( "temp.txt" )
  
  # run the Berry and Berry (2004) model on the drug 
  # trial adverse event data. Create convergence 
  # diagnostics, then summarize the posterior 
  # inferences.

  prepareAEdata <- function(bodySystem, controlCounts, 
    treatmentCounts, totalCtl, totalTrt){
        
    nBodySystems <- max (bodySystem) 
    if (length(unique(bodySystem)) != nBodySystems)
      stop("The maximum body system index must be 
        equal to the number of body systems")

    nAEs <- length(bodySystem)
    if ((length(controlCounts) != nAEs) || 
      (length(treatmentCounts) != nAEs) ||
      (length(totalCtl) != nAEs) || 
      (length(totalTrt) != nAEs))
      stop("The body system index vector, control 
        counts vector, and treatment counts vector must 
        all be the same length")

    # Count the number of AEs for each body system
    # Create a matrix with each row equal to the vector 
    # 1:nBodySystems
    bodySystemVals <- t(array (1:nBodySystems, 
      dim = c(nBodySystems, nAEs)))
    # sum the indicators that each AE is for a 
    # particular body system
    numAEsPerBodySystem <- 
      colSums (bodySystem == bodySystemVals)              
    maxNumAEs <- max(numAEsPerBodySystem)

    # calculate the number of exposures for each AE in 
    # the control and treatment groups. There should be 
    # the same number of exposures in the control group 
    # for all AEs, and the same number of exposures 
    # in the treatment group for all AEs.  
    # First remove NAs, because sometimes the user will 
    # create these counts by dividing the AE count by 
    # the AE percentage.
    totalCtl <- totalCtl [!is.na(totalCtl)]
    totalTrt <- totalTrt [!is.na(totalTrt)]
    nExposedControl <- unique(totalCtl)
    nExposedTreatment <- unique(totalTrt)
    if (length (nExposedControl) != 1)
      stop("The number of exposures in the control 
        group is not the same for all AEs")
    if (length (nExposedTreatment) != 1)
      stop("The number of exposures in the treatment 
        group is not the same for all AEs")
        
    # Create the AE counts and exposure counts arrays; 
    # one row for each body system.  The number of 
    # columns is the maximum number of different 
    # adverse events for a single body system
    controlCountsMat <- array (dim = 
      c(nBodySystems, maxNumAEs))
    treatmentCountsMat <- array (dim = 
      c(nBodySystems, maxNumAEs))

    for (i in 1:nBodySystems){
      nAEsThisBodySystem <- numAEsPerBodySystem[i]
      controlCountsMat [i, 1:nAEsThisBodySystem] <-
        controlCounts [bodySystem == i]
      treatmentCountsMat [i, 1:nAEsThisBodySystem] <- 
        treatmentCounts [bodySystem == i]
    }
         
    aeData <- list("B"=nBodySystems, 
      "Nt"=nExposedTreatment, "Nc"=nExposedControl, 
      "Nobs"=numAEsPerBodySystem, 
      "Y"= treatmentCountsMat,    "X" = controlCountsMat)
                
    aeData
  }

  # Load the adverse event data set 
  data(AE.PTwSOC) 

  # Sort the data set in alphabetical order on the body 
  # system name
  AE.PTwSOC <- AE.PTwSOC[order(AE.PTwSOC$AESOC),]
  controlData <- AE.PTwSOC[AE.PTwSOC$TRT == "TRT B",]
  treatmentData <- AE.PTwSOC[AE.PTwSOC$TRT == "TRT A",]

  nAEs <- numRows(treatmentData)
  b <- if (nAEs != numRows(controlData))
    stop("The number of AEs is different in the 
      treatment and control groups")

  # Preprocess the data 
  aeCounts <- prepareAEdata(bodySystem = 
    as.numeric(controlData$AESOC), 
    controlCounts = controlData$Count, 
    treatmentCounts = treatmentData$Count, 
    totalCtl = controlData$N, 
    totalTrt = treatmentData$N)    

  ##### Run the Berry and Berry model.  THIS TAKES A 
  ##### FEW MINUTES.  

  aeSim <- runAEModel(aeCounts, engine = "WinBUGS")

  ##### Perform convergence diagnostics

  # Choose a representative subset of parameters
  diagnosticParameters <- aeSim[,c("theta.latent[3,1]", 
    "mu.gamma[3]", "mu.theta[3]", "tau.theta", "tau.gamma", 
    "sigma.gamma", "sigma.theta[1]", "sigma.theta[3]", "t[3,1]", 
    "c[3,1]", "mu.theta.0", "mu.gamma.0", "pi[1]", 
    "pi[3]", "alpha.pi", "beta.pi")]
  # The effective sample sizes are large
  e <- effectiveSize(diagnosticParameters)
  # The Geweke and Heidelberger-Welch diagnostics
  # pass
  g <- geweke.diag(diagnosticParameters)
  h <- heidel.diag(diagnosticParameters)

  ##### The convergence diagnostics are acceptable, so 
  ##### summarize the AE results.

  # Only the first 30 parameters are displayed by 
  # the summary function by default.
  s <- summary(aeSim)

  # Calculate the probability that the relative risk is 
  # less than or equal to one for treatment vs. 
  # control, for the first AE.  This is the "Bayesian
  # p-value" for this AE.
  t.samples <- getSamples(aeSim, "t[1,1]")
  c.samples <- getSamples(aeSim, "c[1,1]")
  rr.samples <- t.samples / c.samples
  m <- mean(rr.samples <= 1)
  
  sink()
  ( file.lengths("temp.txt")[1] == 0 )
}
