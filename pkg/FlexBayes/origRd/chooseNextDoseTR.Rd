\name{chooseNextDoseTR}
\alias{chooseNextDoseTR}
\title{Choose the Next Drug Dose Assignment}

\usage{
chooseNextDoseTR (patientDose, outcome, 
  doseValues, doseTrans, currentDoseInd,
  lowestEff, highestTox, pEff, pTox, coefPriors,
  model = "CR", plot = T)
}
\arguments{
\item{patientDose}{A vector of the patient dose assignments so far in the trial.
It should have length equal to the number of patients who have been assigned doses,
and whose outcome is known.  Each element should be an index from 1 to length(doseLevels),
indicating which dose level was used.}
\item{outcome}{A vector of the patient outcomes so far in the trial.  It should have length 
equal to the number of patients who have been assigned doses,
and whose outcome is known.  Each element should be 1 if the patient had an efficacy outcome,
2 if a toxicity outcome, or 0 if neither.}
\item{doseValues}{A vector of the dose levels on the original scale, ordered from smallest to largest.
These are used only for labelling the plots.}
\item{doseTrans}{A vector of the dose levels, centered and transformed as described in 
Thall and Cook (2004), p.686.  Should be ordered from smallest to largest.}
\item{currentDoseInd}{The index of the most recent dose assignment.  The next recommended dose
assignment will not exceed (1+currentDoseInd).}
\item{lowestEff}{The quantity theta*_1 as defined on p. 255 of Thall and Russell (1998).}
\item{highestTox}{The quantity theta*_2 as defined on p. 255 of Thall and Russell (1998).}
\item{pEff}{The quantity pi_1 as defined on p. 255 of Thall and Russell (1998).}
\item{pTox}{The quantity pi_2 as defined on p. 255 of Thall and Russell (1998).}
\item{coefPriors}{A list defining the priors for the fixed effects in the continuation
  ratio models.  Should include elements: 
  \bold{toxCoefMean:} a vector of length 2, giving the prior means for the intercept and 
  slope parameters in the model for toxicity.
  \bold{toxCoefSD:} a vector of length 2, giving the prior standard deviations for the intercept and 
  slope parameters in the model for toxicity.
  \bold{effCoefMean:} a vector of length 2, giving the prior means for the intercept and 
  slope parameters in the model for efficacy given no toxicity.
  \bold{effCoefSD:} a vector of length 2, giving the prior standard deviations for the intercept and 
  slope parameters in the model for efficacy given no toxicity.
  }
\item{model}{Which model to use for the joint probability of efficacy and toxicity.  Either "CR",
  meaning the continuation ratio model of Thall and Cook (2004), or "PO", meaning the proportional odds model
  of Thall and Russell (1998).}
\item{plot}{Whether to plot the posterior intervals of the probability of efficacy and of toxicity at
  each dose level.}
}
\description{
	Use the adaptive dose-finding method of Thall and Russell (1998) to choose the
	dose assignment for the next cohort in the trial.
  For the joint probability of efficacy and toxicity use either the proportional odds model
  of Thall and Russell (1998), or the
  continuation ratio model described in Thall and Cook (2004), p.686.}
\value{
The dose assignment indicator for the next cohort (an index from 1 to length(doseLevels)).
}
  
\source{
O'Quigley, J, Hughes, MD and Fenton, T (2001). Dose-finding designs for HIV studies. 
Biometrics 57, 1018-1029.

TThall, PF and Cook, JD (2004). Dose-finding based on efficacy-toxicity trade-offs. 
Biometrics 60, 684-693.

hall, PF and Russell, KE (1998). A strategy for dose-finding and safety monitoring based
on efficacy and adverse outcomes in phase I/II clinical trials.
Biometrics 54, 251-264.
}

\details{
For an example, see the \code{effToxStroke} help file. 
}
