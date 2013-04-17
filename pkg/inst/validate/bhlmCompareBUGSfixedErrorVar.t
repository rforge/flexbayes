{cat("----- Compare bhlm, BUGS for non-informative priors ----------\n");T}
{
# Functions: posteriorSamples(engine='WinBUGS'), bhlm
# Description: Compare the results for the bhlm function
#  to those from the BUGS engine, for a linear 
#  regression model fit to the orthodont data as 
#  described in the S+Bayes manual.

  bhlmCompareBUGSfixedErrorVar();
}
