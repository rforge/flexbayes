{cat("----- bhlm: Orthodont data, Wishart prior, second-level effects -------\n");T}
{
# Functions: posteriorSamples(engine='WinBUGS'), bhlm
# Description: Compare the results for the bhlm function
#  to those from the BUGS engine, for a linear 
#  regression model fit to the orthodont data.
#  Use an inverse Wishart prior for the random effect 
#  covariance matrix, and include gender as a second-
#  level effect.

  bhlmCompareBUGS10();
}
