{cat("----- bhlm: Orthodont data, one subject, Wishart prior -------\n");T}
{
# Functions: posteriorSamples(engine='WinBUGS'), bhlm
# Description: Compare the results for the bhlm function
#  to those from the BUGS engine, for a linear 
#  regression model fit to just one subject in the orthodont data.
#  Use an inverse Wishart prior for the random effect variance.

  bhlmCompareBUGS6();
}
