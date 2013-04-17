loadCoda <- function(){
  # if coda is before FlexBayes in the search path, detach coda
  if( ( "coda" %in% search() ) &&
      ( codaLoc() < FlexBayesLoc() ) )
    detach("coda")
  # attempt to attach coda
  loaded <- require("coda")
  if( !loaded )
    stop("The S-PLUS open-source package 'coda' is required for use of this function.  Obtain coda from http://csan.insightful.com")
  invisible( loaded )
}

loadR2WinBUGS <- function(){
	# if R2WinBUGS is before FlexBayes in the search path, detach it
  if( ( "R2WinBUGS" %in% search() ) &&
      ( R2WinBUGSLoc() < FlexBayesLoc() ) )
    detach("R2WinBUGS")
  # attempt to attach R2WinBUGS
  loaded <- require("R2WinBUGS")
  if( !loaded )
    stop("The S-PLUS open-source package 'R2WinBUGS' is required for use of this function.  Obtain R2WinBUGS from http://csan.insightful.com")
  
  #check that WinBUGS is installed
  
  invisible( loaded )
}

loadBRugs <- function(){
  # if BRugs is before FlexBayes in the search path, detach BRugs
  if( ( "BRugs" %in% search() ) &&
      ( BRugsLoc() < FlexBayesLoc() ) )
    detach("BRugs")
  # attempt to attach BRugs
  loaded <- require("BRugs")
  if( !loaded )
    stop("The S-PLUS open-source package 'BRugs' is required for use of this function.  Obtain BRugs from http://csan.insightful.com")
  invisible( loaded )
}

FlexBayesLoc <- function(){
  ##  Finds the location of FlexBayes on the search path in a case-insensitive fashion
  getSearchPathLocation( "posteriorSamples" )
}

codaLoc <- function(){
  ##  Finds the location of coda on the search path in a case-insensitive fashion
  getSearchPathLocation( "read.coda" )
}

R2WinBUGSLoc <- function(){
  ##  Finds the location of R2WinBUGS on the search path in a case-insensitive fashion
  getSearchPathLocation( "read.bugs" )
}

BRugsLoc <- function(){
  ##  Finds the location of BRugs on the search path in a case-insensitive fashion
  getSearchPathLocation( "writeModel" )
}

getSearchPathLocation <- function( lib.fn.name ){
  ##  Finds the location on the search path of the library that contains the function lib.fn.name.  If there is more than
  ##  one such library, then the location of the first is returned.
  
  lib.locations <- names( objects( where = "*", pattern = lib.fn.name ) ) [ 1 ]
  which( search() == lib.locations )
}
