{cat("----- Compare bhlm, WinBUGS on the orthodont data, with only a random effect intercept -------\n");T}
{
# Functions: posteriorSamples(engine='WinBUGS'), bhlm
# Description: Compare the results for the bhlm function
#  to those from the BUGS engine, for a linear 
#  regression model fit to the orthodont data with
#  just a random effect intercept at the first level
#  and an intercept and sex effect at the second level.

  bhlmCompareBUGS5();
}
