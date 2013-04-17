.First.lib <- function(lib.loc, section)
{
 if( interactive() ) {
   cat("Welcome to FlexBayes\n  Examples are available by calling help(FlexBayes) from the command line,\n  or by opening the Help->Available Help->FlexBayes menu item,\n  then navigating to the FlexBayes help file.\n")
 } 
 # detach libraries that could cause a conflict
 if ("coda" %in% search())
   detach("coda")
 if ("R2WinBUGS" %in% search())
   detach("R2WinBUGS")
 if ("BRugs" %in% search())
   detach("BRugs")
 invisible()
}
