{cat("----- No extraneous info printed while fitting dose-finding model --\n");T}
{
# Functions: chooseNextDoseTR
# Description: Make sure that no extraneous info is printed 
# while fitting a dosage model
  
  # save printed output to a file
  sink( "temp.txt" )
  
  cohSize <- 3  # the cohort size  
  nCohorts <- 10  # the number of cohorts

  ## Define the acceptable tox / eff probabilities to 
  ## be used in the Thall and Russell (1998) 
  ## dose-finding method
  lowestEff <- 0.5
  highestTox <- 0.1
  pEff <- 0.1
  pTox <- 0.1

  # Pick a true model and a model to fit. 
  # "CR": continuation ratio or 
  # "PO": proportional odds
  fitModel <- "CR"
  trueModel <- "PO"

  ## Define and transform the dose levels.
  doseValues <- c(0, 2.5, 5, 7.5, 10)
  nDose <- length(doseValues)
  # The dose levels here need to be ordered
  # lowest to highest
  doseValues <- doseValues[order(doseValues)]
  # log-transform and center the dose values
  # for the regression
  ifVal <- if (doseValues[1] == 0){
    doseTrans <- doseValues + doseValues[2]
  }
  doseTrans <- log(doseTrans)
  doseTrans <- doseTrans - mean(doseTrans)

  ifVal <- if (fitModel == "CR"){
    # For the continuation ratio model,
    # define the priors for the coefficients in the
    # toxicity and efficacy-given-no-toxicity models.
    coefPriors <- list( toxCoefMean = c( -1.966, 1.059 ),
      toxCoefSD = c( 1.791, 1.791 ), 
      effCoefMean = c( 0.464, 0.968 ), 
      effCoefSD = c( 0.332, 0.333 ) )
  } else {
    # For the proportional odds model,
    # define the priors for the coefficients.
    coefPriors <- list( betaLower = 0,
      betaUpper = 2, muLower = -6,
      muUpper = 0, alphaLower = 0,
      alphaUpper = 4 )
  }

  # Define the dose efficacy and toxicity rates for the
  # data to be simulated.  
  ifVal <- if (trueModel == "CR"){
    plogitTox <- -2.5 + doseTrans*2
    doseProbTox <- exp(plogitTox) / (1+exp(plogitTox))
    plogitEff <- 0.15 + doseTrans*1
    doseProbEff <- exp(plogitEff) / (1+exp(plogitEff))
    doseProbEff <- doseProbEff * (1-doseProbTox)
  } else {
    plogitTox <- -2.5 + doseTrans*2
    doseProbTox <- exp(plogitTox) / (1+exp(plogitTox))
    plogitEff <- -2.5 + 2.5 + doseTrans*2
    doseProbEff <- exp(plogitEff) / (1+exp(plogitEff))
    doseProbEff <- doseProbEff * (1-doseProbTox)
  }

  # initialize the vectors of patient assigned doses
  # and outcomes
  patientDose <- numeric(0)
  outcome <- numeric(0)
  # start at the minimum dose
  currentDoseInd <- 1  

  # initialize the cohort index
  coh <- 1
  # Run the simulation, assigning dose to patients 
  # using the adaptive method.
  # After each cohort, plot point and interval estimates
  # for the probability of efficacy and toxicity at 
  # each dose.
  while ((coh <= nCohorts) && (currentDoseInd > 0)) {
          
    patientDose <- c(patientDose, 
      rep(currentDoseInd,cohSize))
    trueProbEff <- doseProbEff[currentDoseInd]
    trueProbTox <- doseProbTox[currentDoseInd]

    # simulate patient response for the new cohort
    for (i in (1:cohSize)){
      tmp <- runif(n=1)
      if (tmp < trueProbEff)  
        outcome <- c(outcome, 1)  # efficacy event
      else if (tmp < trueProbEff + trueProbTox)  
        outcome <- c(outcome, 2)  # toxicity event
      else outcome <- c(outcome, 0)  # no event
    }
        
    currentDoseInd <- chooseNextDoseTR (
      patientDose = patientDose, 
      outcome = outcome, doseTrans = doseTrans,
      doseValues = doseValues,
      currentDoseInd = currentDoseInd,
      lowestEff = lowestEff, highestTox = highestTox,
      pEff = pEff, pTox = pTox,
      coefPriors = coefPriors, model = fitModel ) 

    coh <- coh + 1
  }

  sink()
  ( file.lengths("temp.txt")[1] == 0 )
}
