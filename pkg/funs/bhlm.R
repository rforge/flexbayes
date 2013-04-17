##################################################
# bhlm class in the Bayes Module                 #
# Insightful: NIH Bayes II                       #
# Alejandro Murua 1/03                           #
##################################################




#####################################################################################
##
##  LIKELIHOOD
##
#####################################################################################

bhlm.likelihood <- function( type = "normal",
                             df = 3 )
{
  errorCov <- NULL
  type <- match.arg( type, choices = c( "normal", "t" ) )
  list( type = type, errorCov = errorCov, df = df )
}#end





#####################################################################################
##
##  PRIOR
##
#####################################################################################

bhlm.prior <- function( error.var = bayes.nonInfoPower( -1.0 ),
                        fixed.coef = "non-informative",
                        level2.coef = "non-informative",
                        random.var = bayes.nonInfoPower( -1.0 ),
                        common.error.var = 2 )
{

  # create a list of indicators of whether the prior has been specified 
  # by the user for each parameter.  This is used in bhlm to check that the user 
  # has not specified a prior for a parameter that is not in the model.
  priorSpec <- list( error.var = !missing( error.var ),
    fixed.coef = !missing( fixed.coef ),
    level2.coef = !missing( level2.coef ),
    random.var = !missing( random.var ) )

  random.coef <- bayes.normal(mean=zero, cov=identity)

  # for common sigma = error.var
  # 0 ==> independent sigmas; different hyper prior parameters
  # 1 ==> independent sigmas; same hyper prior parameters
  # 2 ==> common sigma

  if ( !is.null( error.var ) && class( error.var ) == "bayes.distribution" )
  {
    if ( !is.element( common.error.var, c( 1, 2 ) ) )
    {
      stop( "bhlm.prior: need a list of prior parameters for error variance in this model." )
    }


    if ( !is.element( error.var@name, c( "invChisq", "uniform shrinkage", "du Mouchel", "non-informative power", "mass point", "invWishart" ) ) )
    {
      stop( "bhlm.prior: Invalid prior distribution for error variance." )
    }
  }
  else if ( !is.null( error.var ) && is.list( error.var ) )
  {
    if ( common.error.var != 0 )
    {
      stop( "bhlm.prior: need a single prior for error variance in this model." )
    }

    for ( i in seq( 1, length( error.var ) ) )
    {
      if ( class( error.var[[i]] ) == "bayes.distribution" )
      {
        if ( !is.element( error.var[[i]]@name, c( "invChisq", "uniform shrinkage", "du Mouchel", "non-informative power", "mass point", "invWishart" ) ) )
        {
          stop( "bhlm.prior: Invalid prior distribution for error variance." )
        }
      }
      else
      {
        stop( "bhlm.prior: Invalid prior distribution for error variance." )
      } 
    }#end for loop
  }
  else
  {
    stop( "bhlm.prior: Invalid prior distribution for error variance." )
  }  

  ##
  ## Now check the variance for the random coefficients
  ##

  if ( !is.null( random.var ) && class( random.var ) == "bayes.distribution" )
  {
    if ( !is.element( random.var@name, c( "invChisq", "uniform shrinkage", "du Mouchel", "invWishart", "non-informative power" ) ) )
    {
      stop( "bhlm.prior: Invalid prior distribution for random effects coefficients variance." )
    }
  }
  else if ( !is.null( random.var ) )
  {
    stop( "bhlm.prior: Invalid prior distribution for random effects coefficients variance." )
  }



  ##
  ## check the prior for random coefficients
  ## if a list then check each element of the list
  ##

  if ( !is.null( random.coef ) && class( random.coef ) == "bayes.distribution" )
  {
    if ( !is.element( random.coef@name, c( "normal", "t", "non-informative" ) ) )
    {
      stop( "bhlm.prior: Invalid prior distribution for random effects coefficients." )
    }
  }
  else if ( !is.null( random.coef ) && is.list( random.coef ) )
  {
    for ( i in seq( 1, length( random.coef ) ) )
    {
      if ( class( random.coef[[i]] ) == "bayes.distribution" )
      {
        if ( !is.element( random.coef[[i]]@name, c( "normal", "t", "non-informative" ) ) )
        {
          stop( "bhlm.prior: Invalid prior distribution for random effects coefficients." )
        }        
      }
      else if ( random.coef[[i]] != "non-informative" )
      {
        stop( "bhlm.prior: Invalid prior distribution for random effects coefficients." )
      }
      else if ( random.coef[[i]] == "non-informative" )  
      {
        random.coef[[i]] <- bayes.nonInformative()
      }
    }#end for loop
  }
  else if ( !is.null( random.coef )  && random.coef != "non-informative" )
  {
    stop( "bhlm.prior: Invalid prior distribution for random effects coefficients." )
  }
  else if ( !is.null( random.coef )  && random.coef == "non-informative" )  
  {
    random.coef <- bayes.nonInformative()
  }


  ##
  ## check prior for fixed effects
  ##


  if ( !is.null( fixed.coef ) && class( fixed.coef ) == "bayes.distribution" )
  {
    if ( !is.element( fixed.coef@name, c( "normal", "t", "non-informative" ) ) )
    {
      stop( "bhlm.prior: Invalid prior distribution for fixed effects coefficients." )
    }
  }
  else if ( !is.null( fixed.coef ) && fixed.coef != "non-informative" )
  {  
    stop( "bhlm.prior: Invalid prior distribution for fixed effects coefficients." )
  }
  else if ( fixed.coef == "non-informative" )  
  {
    fixed.coef <- bayes.nonInformative()
  }



  ##
  ## check prior for sencond level regression parameters
  ##


  if ( !is.null( level2.coef ) && class( level2.coef ) == "bayes.distribution" )
  {
    if ( !is.element( level2.coef@name, c( "normal", "t", "non-informative" ) ) )
    {
      stop( "bhlm.prior: Invalid prior distribution for second level coefficients." )
    }
  }
  else if ( !is.null( level2.coef ) && level2.coef != "non-informative" )
  {  
    stop( "bhlm.prior: Invalid prior distribution for second level coefficients." )
  }
  else if ( !is.null( random.coef )  && level2.coef == "non-informative" )  
  {
    level2.coef <- bayes.nonInformative()
  }



  ##
  ## create prior object for hierarchical model
  ##

  list( error.var = error.var, 
        random.coef = random.coef,
        fixed.coef = fixed.coef,
        level2.coef = level2.coef,
        random.var = random.var,
        common.error.var = common.error.var,
        priorSpec = priorSpec )

}#end



#####################################################################################
##
##  SAMPLER
##
#####################################################################################

bhlm.sampler <- function( nBurnin = 1000, nSamples = 1000, nThin = 1, 
                          nChains = 1,
                          init.point = "prior",
                          params.init = NULL)
{
  number.chains <- nChains
  random.effects <- T
  fixed.effects <- T   
  
  #assume params.init is a list of list with initial points for each chain
  if ( number.chains <= 0 )
  {
    warning( "bhlm.sampler: zero or negative number of chains to simulate. Assuming 1 chain." )
    number.chains <- 1
  }

  #just in case this is not an integer as it should be  
  if ( ceiling( number.chains ) - number.chains > 0 )
  {
    new.number.chains <- ceiling( number.chains )
    warning( "bhlm.sampler: Non-integer number of chains specified [", number.chains, "]. Set number of chains to [", new.number.chains, "]." )
    number.chains <- new.number.chains
  }

  init.point <- match.arg( init.point, choices = c( "prior", "prior + likelihood", "likelihood", "user's choice" ) )
  if ( init.point == "user's choice" )
  {
    if ( is.null( params.init ) )
    {
      stop( "bhlm.sampler: initial points must be specified." )
    }
    else
    { 
      if ( number.chains == 1 && !is.list( params.init[[1]] ) )
      {
        params.init <- list( params.init )
      }

      if ( length( params.init ) < number.chains )
      {
        stop( "bhlm.sampler: initial values for all parameters in all chains must be specified." )
      }
      else
      {
        for ( i in seq( 1, number.chains ) )
        {
          #cat("init values are\n")
          #print( params.init[[i]] )

          error.var <- params.init[[ i ]]$error.var
          random.coef <- params.init[[ i ]]$random.coef
          fixed.coef <- params.init[[ i ]]$fixed.coef
          level2.coef <-  params.init[[ i ]]$level2.coef
          random.var <- params.init[[ i ]]$random.var

          #check error variance
          if ( is.null( error.var ) )
          {
            stop( "bhlm.sampler: initial value for error variance must be provided." )
          }    
          else if ( !is.vector( error.var ) && !is.matrix( error.var ) && !is.list( error.var ) )
          {
            stop( "bhlm.sampler: Initial value for error variance must be provided. " )
          }
          else if ( is.vector( error.var ) )
          {
            if ( any( error.var ) <= 0 )
            {
              stop( "bhlm.sampler: initial value for error variance must be positive." )
            }
          }    
          else if ( is.list( error.var ) )
          {
            if ( any ( diag( error.var[[1]] < 0 ) ) )
            {
              stop( "bhlm.sampler: initial value for error variances must be positive." )
            }

            if ( error.var[[1]] - t( error.var[[1]] ) != 0 )
            {
              stop( "bhlm.sampler: initial value for error variance must be a symmetric matrix." )
            }

            this.error.var <- error.var[[1]]
            if ( length( error.var ) > 1 )
            {
              for (  j in seq( 2, length( error.var ) ) )
              {
                if ( any ( diag( error.var[[j]] < 0 ) ) )
                {
                  stop( "bhlm.sampler: initial value for error variances must be positive." )
                }

                if ( error.var[[j]] - t( error.var[[j]] ) != 0 )
                {
                  stop( "bhlm.sampler: initial value for error variance must be a symmetric matrix." )
                }

                this.error.var <- cbind( this.error.var, error.var[[j]] )
              }
            }
            params.init[[ i ]]$error.var <- this.error.var
          }
          
          if ( random.effects && is.null( level2.coef ) )
          {     
            if ( is.null( random.coef ) )
            {
              stop( "bhlm.sampler: initial value for random effects coefficients must be provided." )
            }
          }

          if ( fixed.effects && is.null( fixed.coef ) )
          {
            stop( "bhlm.sampler: initial values for the fixed effects coefficients must be provided." )
          }

          #check parameters variance
          if ( is.null( random.var ) )
          {
            stop( "bhlm.sampler: initial value for random effects coefficients scale must be provided." )
          }    
          else if ( is.vector( random.var ) )
          {
            if ( any( random.var ) <= 0 )
            {
              stop( "bhlm.sampler: Initial value for random effects coefficients scale must be positive. " )
            }
          }
          else if ( is.matrix( random.var ) )
          {
            if ( any( diag( random.var ) ) <= 0 )
            {
              stop( "bhlm.sampler: Initial value for random effects coefficients variances must be positive. " )
            }


            if ( any( random.var - t( random.var ) != 0 ) )
            {
              stop( "bhlm.sampler: Initial value for random effects coefficients variances must be a symmetric matrix." )
            }
          }  
          else
          {
            stop( "bhlm.sampler: Initial value for random effects coefficients scale are not valid." )
          }
        }#end for i
      }
    }
  }#end user's choice

  type <- "Gibbs"
  
  list( control = sampler.control(nSamples=nSamples, nThin=nThin, nBurnin=nBurnin), 
        number.chains = number.chains,
        init.point = list( values = params.init, type = init.point ), 
        sampler = type )
}#end



#####################################################################################
##
##  VALIDATION FOR MEAN AND COVARIANCE
##
#####################################################################################


valid.mean.specification <- function( x, dim )
{

  if ( is.function( x ) )
  {
    x <- x( dim )
  }

  if ( is.name( x ) )
  {
    x <- eval( call( x, p = dim ) )
  }

  if ( is.vector( x ) && length( x ) != dim )
  {
    stop( "bhlm: mean specification has wrong dimension." )
  }

  x
}#end


valid.Cov.specification <- function( x, dim )
{

  if ( is.function( x ) )
  {
    x <- x( dim )
  }

  if ( is.name( x ) )
  {
    x <- eval( call( x, p = dim ) )
  }


  if ( is.vector( x ) )
  {
    if ( any( x <= 0 ) )
    {
      stop( "bhlm: covariance main diagonal must be positive." )
    }
    
    if ( length( x ) == 1 && dim > 1 )
    {
      x <- x * diag( dim )
    }
    else if ( length( x ) != dim )
    {
      stop( "bhlm: covariance specification has wrong dimension." )
    }
    else if ( length( x ) == dim && dim > 1 )
    {
      x <- diag( x )
    }
    else if ( length( x ) == dim && dim == 1 )
    {
      x <- matrix( x, nrow = 1 )
    }  
  }

  if ( is.matrix( x ) )
  {
    if ( ncol( x ) != dim || nrow( x ) != dim )
    {
      stop( "bhlm: covariance specification has wrong dimensions." )
    }
    else if ( any( x - t(x) != 0 ) )
    {
      stop( "bhlm: covariance specification is not symmetric." )
    }
    else if ( any( diag( x ) <= 0 ) )
    {
      stop( "bhlm: covariance specification contains negative or zero variances." )
    }    
  }
  else
  {
    stop( "bhlm: covariance specification is not a matrix!." )
  }

  x
}#end



##############################################################################################
##
##VALIDATION OF INITIAL POINTS
##
##############################################################################################

validate.initial.points.bhlm <- function( n.groups, dimX, dimM, dimZ, init.point, random.effects, 
fixed.effects, second.effects, error.var.common, error.var.type, random.coef.type, random.var.type )
{
  s.error.var <- NULL
  s.random.coef <- NULL
  s.fixed.coef <- NULL
  s.random.var <- NULL
  s.level2.coef <- NULL

  #check that is not invWishart nor known
  if ( !is.element( error.var.type, c( 4, 5 ) ) )
  {
    if ( error.var.common != 0 )
    {
      if ( length( init.point$error.var ) > 1 )
      {
        stop( "bhlm: initial value for error variance should not be an array." )
      }
      else
      {
        s.error.var <- matrix( init.point$error.var, nrow = 1 )
      }
    }
    else 
    {
      if ( length( init.point$error.var ) != n.groups )
      {
        stop( "bhlm: initial values for error variance should be an array of number of groups elements." )
      }
      else
      {
        s.error.var <- matrix( init.point$error.var, nrow = n.groups )  
      }
    }
  }

  if ( random.effects )
  {
    if ( second.effects )
    {
      if ( is.null( init.point$level2.coef ) )
      {
        stop( "bhlm: initial values for second stage regression coefficients must be provided." )
      }
      else if ( !is.vector( init.point$level2.coef ) )
      {
        stop( "bhlm: invalid initial values for second stage regression coefficients." )
      }
      else if ( length( init.point$level2.coef ) != dimZ )
      {
        stop( "bhlm: initial values for second stage regression coefficients have wrong dimension." )
      }       
      else
      {
        s.level2.coef <- matrix( init.point$level2.coef, nrow = dimZ )
      }
    }

    #if beta has an informative prior
    if ( random.coef.type != 3 )
    {
      if ( random.var.type == 4 )
      {
        #inv wishart case
        if ( is.null( init.point$random.var ) )
        {
          stop( "bhlm: initial value for the covariance matrix associated to the random coefficients must be provided." )
        }

        s.random.var <- valid.Cov.specification( init.point$random.var, dimX )
      }
      else
      {
        if ( is.null( init.point$random.var ) )
        {
          stop( "bhlm: initial value for the variance associated to the random coefficients must be provided." )
        }

        if ( is.vector( init.point$random.var ) )
        {
          if ( any( init.point$random.var <= 0 ) )
          {
            stop( "bhlm: initial value for the variance associated to the random coefficients must be positive." )
          }

          s.random.var <- matrix( init.point$random.var, nrow = 1 )
        }
        else
        {
          stop( "bhlm: wrong specification for initial value for the variance associated to the random coefficients." ) 
        }

      }#end a scalar variance
    }


    if ( random.coef.type != 0 )
    {
      if ( is.null( init.point$random.coef ) )
      {
        stop( "bhlm: initial values for random coefficients must be specified." )
      }

      if ( random.coef.type == 1 )
      {
        if ( is.vector( init.point$random.coef ) && length( init.point$random.coef ) == dimX )
        {
          s.random.coef <- matrix( init.point$random.coef, nrow = dimX )
        }
        else
        {
          stop( "bhlm: wrong specification for initial values for random coefficients." )
        }
      }
      else if ( random.coef.type == 2 )
      {
        if ( is.list( init.point$random.coef ) && length( init.point$random.coef ) == n.groups )
        {
          s.random.coef <- apply( matrix( seq( 1, n.groups ), nrow = 1 ), 
                                   MARGIN = 2, 
                                   FUN = function( x, y, dim )
                                   {
                                     r.mean <- y[[ x ]]
                                     r.mean <- valid.mean.specification( r.mean, dim )
                                     r.mean
                                   },
                                   y = init.point$random.coef,
                                   dim = dimX )            

        }
        else if ( is.matrix( init.point$random.coef ) && ncol( init.point$random.coef ) == n.groups && nrow( init.point$random.coef ) == dimX )
        {
          s.random.coef <- init.point$random.coef
        }
        else
        {
          stop( "bhlm: wrong specification for initial values of random effects coefficients." )
        }
      }
      else #non-informative case
      {
        if ( is.vector( init.point$random.coef ) && length( init.point$random.coef ) == dimX )
        {
          s.random.coef <- matrix( init.point$random.coef, nrow = dimX )
        }
        else
        {
          stop( "bhlm: wrong specification for initial values of random effects coefficients." )
        }
      }
    }#end no second level regression

  }#end if random effects

  if ( fixed.effects )
  {
    if ( is.null( init.point$fixed.coef ) )
    {
      stop( "bhlm: initial values of fixed effects coefficients must be provided." )
    }

    if ( is.vector( init.point$fixed.coef ) && length( init.point$fixed.coef ) == dimM )
    {
      s.fixed.coef <- matrix( init.point$fixed.coef, nrow = dimM )
    }
    else
    {
      stop( "bhlm: wrong specification for initial values of fixed effects coefficients." )
    }
  }#end if fixed effects

  list( error.var = s.error.var, random.coef = s.random.coef, fixed.coef = s.fixed.coef, random.var = s.random.var, level2.coef = s.level2.coef )

}#end validate initial point



##############################################################################################
##
##VALIDATION OF INITIAL POINTS (MultiVariate and Imputing Means case
##
##############################################################################################

validate.initial.points.bhlmMV <- function( n.groups, dimR, dimX, dimM, dimZ, init.point, random.effects, fixed.effects, second.effects, error.var.common, random.coef.type, error.var.type, random.var.type )
{
  s.error.var <- NULL
  s.random.coef <- NULL
  s.fixed.coef <- NULL
  s.random.var <- NULL
  s.level2.coef <- NULL


  if ( error.var.type == 4 )
  {
    if ( error.var.common != 0 )
    {
      #inv wishart case
      if ( is.null( init.point$error.var ) )
      {
        stop( "bhlm: initial value for the covariance matrix associated to the response variables must be provided." )
      }

      s.error.var <- valid.Cov.specification( init.point$error.var, dimR )
    }
    else if ( is.list( init.point$error.var ) )
    {
      if ( length( init.point$error.var ) != n.groups )
      {
        stop( "bhlm: initial values for error variance should be a list of number of groups matrices." )
      }
      else
      {
        for ( i in seq( 1, n.groups ) )
        {
          if ( any( diag( init.point$error.var[[i]] ) < 0 ) )
          {
            stop( "bhlm: Measurement error variances must be positive.\n" )
          }
          else if ( t( init.point$error.var[[i]] ) != init.point$error.var[[i]] )
          {
            stop( "bhlm: Measurement error covariance matrix must be symmetric.\n" )
          }
          else
          {
            if ( i == 1)
            {
              s.error.var <- init.point$error.var[[i]]
            }
            else
            {
              s.error.var <- cbind( s.error.var, init.point$error.var[[i]] )
            }
          }
        }#end for i
      }
    }
    else if ( is.matrix( init.point$error.var ) )
    {
      if ( ncol( init.point$error.var ) != dim.response * n.groups || nrow( init.point$error.var ) != dim.response )
      {
        stop( "bhlm: initial values for error variance should be a matirx of number of groups times number of response-variable columns and number of response-variable rows." )
      }
      else
      {
        s.error.var <- init.point$error.var
      }
    }
  }
  else if ( error.var.common != 0 )
  {
    if ( length( init.point$error.var ) > 1 )
    {
      stop( "bhlm: initial value for error variance should not be an array." )
    }
    else
    {
      s.error.var <- matrix( init.point$error.var, nrow = 1 )
    }
  }
  else 
  {
    if ( length( init.point$error.var ) != n.groups )
    {
      stop( "bhlm: initial values for error variance should consists of number of groups elements." )
    }
    else
    {
      s.error.var <- matrix( init.point$error.var, nrow = n.groups )  
    }
  }

  if ( random.effects )
  {
    if ( second.effects )
    {
      if ( is.null( init.point$level2.coef ) )
      {
        stop( "bhlm: initial values for second stage regression coefficients must be provided." )
      }
      else if ( !is.vector( init.point$level2.coef ) )
      {
        stop( "bhlm: invalid initial values for second stage regression coefficients." )
      }
      else if ( length( init.point$level2.coef ) != dimZ )
      {
        stop( "bhlm: initial values for second stage regression coefficients have wrong dimension." )
      }       
      else
      {
        s.level2.coef <- matrix( init.point$level2.coef, nrow = dimZ )
      }
    }

    #if beta has an informative prior
    if ( random.coef.type != 3 )
    {
      if ( random.var.type == 4 )
      {
        #inv wishart case
        if ( is.null( init.point$random.var ) )
        {
          stop( "bhlm: initial value for the covariance matrix associated to the random coefficients must be provided." )
        }

        s.random.var <- valid.Cov.specification( init.point$random.var, dimX )
      }
      else
      {
        if ( is.null( init.point$random.var ) )
        {
          stop( "bhlm: initial value for the variance associated to the random coefficients must be provided." )
        }

        if ( is.vector( init.point$random.var ) && length( init.point$random.var ) == 1 )
        {
          if ( init.point$random.var <= 0 )
          {
            stop( "bhlm: initial value for the variance associated to the random coefficients must be positive." )
          }

          s.random.var <- matrix( init.point$random.var, nrow = 1 )
        }
        else
        {
          stop( "bhlm: wrong specification for initial value for the variance associated to the random coefficients." ) 
        }

      }#end a scalar variance
    }


    if ( random.coef.type != 0 )
    {
      if ( is.null( init.point$random.coef ) )
      {
        stop( "bhlm: initial values for random coefficients must be specified." )
      }

      if ( random.coef.type == 1 )
      {
        if ( is.vector( init.point$random.coef ) && length( init.point$random.coef ) == dimX )
        {
          s.random.coef <- matrix( init.point$random.coef, nrow = dimX )
        }
        else
        {
          stop( "bhlm: wrong specification for initial values for random coefficients." )
        }
      }
      else if ( random.coef.type == 2 )
      {
        if ( is.list( init.point$random.coef ) && length( init.point$random.coef ) == n.groups )
        {
          s.random.coef <- apply( matrix( seq( 1, n.groups ), nrow = 1 ), 
                                   MARGIN = 2, 
                                   FUN = function( x, y, dim )
                                   {
                                     r.mean <- y[[ x ]]
                                     r.mean <- valid.mean.specification( r.mean, dim )
                                     r.mean
                                   },
                                   y = init.point$random.coef,
                                   dim = dimX )            

        }
        else if ( is.matrix( init.point$random.coef ) && ncol( init.point$random.coef ) == n.groups && nrow( init.point$random.coef ) == dimX )
        {
          s.random.coef <- init.point$random.coef
        }
        else
        {
          stop( "bhlm: wrong specification for initial values of random effects coefficients." )
        }
      }
      else #non-informative case
      {
        if ( is.vector( init.point$random.coef ) && length( init.point$random.coef ) == dimX )
        {
          s.random.coef <- matrix( init.point$random.coef, nrow = dimX )
        }
        else
        {
          stop( "bhlm: wrong specification for initial values for random coefficients." )
        }
      }
    }#end no second level regression

  }#end if random effects

  if ( fixed.effects )
  {
    if ( is.null( init.point$fixed.coef ) )
    {
      stop( "bhlm: initial values of fixed effects coefficients must be provided." )
    }

    if ( is.vector( init.point$fixed.coef ) && length( init.point$fixed.coef ) == dimM )
    {
      s.fixed.coef <- matrix( init.point$fixed.coef, nrow = dimM )
    }
    else
    {
      stop( "bhlm: wrong specification for initial values of fixed effects coefficients." )
    }
  }#end if fixed effects

  list( error.var = s.error.var, random.coef = s.random.coef, fixed.coef = s.fixed.coef, random.var = s.random.var, level2.coef = s.level2.coef )

}#end validate initial point






#####################################################################################
##
##  BAYESIAN HIERARCHICAL LINEAR MODEL
##
#####################################################################################

#response.formula, random.formula, fixed.formula, level2.formula and group.formula are of type "formula"
bhlm <- function( fixed.formula = NULL, 
                  random.formula = NULL, 
                  level2.formula = NULL, 
                  group.formula = NULL, 
                  data,
                  prior = bhlm.prior(), 
                  likelihood  = bhlm.likelihood(),
	                sampler = bhlm.sampler(), 
                  random.seed = .Random.seed, na.action = NULL, contrasts = NULL,
                  debug = F )
{
  response.formula <- NULL

  ####-------------------------------------------------------
  ## Get Data 
  ##

  if ( is.data( data ) )
  {
    stop( "bhlm: data must be provided." )
  }


  found.response <- FALSE
  multi.variate.response <- FALSE

  if ( !is.null( response.formula ) && length( response.formula ) > 0 && !is.all.white( response.formula )[1] )
  {
    if ( length( terms( response.formula )@term.labels ) > 1 )
    {
      multi.variate.response <- T
    }
    else if ( length( terms( response.formula )@term.labels ) == 1 )
    {
      if ( terms( response.formula )@term.labels == "." )
      {
        multi.variate.response <- T
      }
      else
      {
        multi.variate.response <- F
      }
    }
    else
    {
      stop( "bhlm: Formula specifying response variable(s) is not valid.\n" )
    } 

    found.response <- T 
  }

  if ( !found.response )
  {
    if ( !is.null( random.formula ) && length( random.formula ) > 0 && !is.all.white( random.formula )[1] )
    {
      if ( terms( random.formula )@response == 1 )
      {
        found.response <- T
        if ( charmatch( "cbind(", dimnames( terms( random.formula )@factors )[[1]][1], nomatch = - 1 ) == 1 )
        {
          if ( grep( ",", dimnames( terms( random.formula )@factors )[[1]][1] ) >= 1 )
          {
            multi.variate.response <- T
          }
          else
          {
            multi.variate.response <- F
          }
        }
        else
        {
          multi.variate.response <- F
        }
      }
    }
  }

  if ( !found.response )
  {
    if ( !is.null( fixed.formula ) && length( fixed.formula ) > 0 && !is.all.white( fixed.formula )[1] )
    {
      if ( terms( fixed.formula )@response == 1 )
      {
        found.response <- T
        if ( charmatch( "cbind(", dimnames( terms( fixed.formula )@factors )[[1]][1], nomatch = - 1 ) == 1 )
        {
          if ( grep( ",", dimnames( terms( fixed.formula )@factors )[[1]][1] ) >= 1 )
          {
            multi.variate.response <- T
          }
          else
          {
            multi.variate.response <- F
          }
        }
        else
        {
          multi.variate.response <- F
        }
      }
    }
  }

  if ( !found.response )
  {
    stop( "bhlm: response variable(s) must be provided.\n" )
  }


  #output
  out.bhlm <- NULL

  if ( !multi.variate.response )
  {
    #call bhlmUV
    out.bhlm <- bhlmUV( response.formula = response.formula, 
                random.formula = random.formula, fixed = fixed.formula, 
                level2 = level2.formula, group = group.formula, data = data, 
                prior = prior, likelihood = likelihood, sampler = sampler, 
                random.seed = random.seed, na.action = na.action, contrasts = contrasts,
                debug = debug )
  }
  else
  {
    #need to check if random, fixed, level2 formulas are trivial

    means.model <- T
    if ( !is.null( random.formula ) && length( random.formula ) > 0 && !is.all.white( random.formula )[1] )
    {
      if ( length( terms( random.formula )@term.labels ) > 0 )
      {
        means.model <- F
      }
    }

    if ( means.model && !is.null( fixed.formula ) && length( fixed.formula ) > 0 && !is.all.white( fixed.formula )[1] )
    {
      if ( length( terms( fixed.formula )@term.labels ) > 0 )
      {
        means.model <- F
      }
      else if ( terms( fixed.formula )@intercept == 1 )
      {
        means.model <- F
      }
    }

    if ( !is.null( level2.formula ) && length( level2.formula ) > 0 && !is.all.white( level2.formula )[1] )
    {
      if ( length( terms( level2.formula )@term.labels ) > 0 )
      {
        means.model <- F
      }
    }


    if ( means.model )
    {
      random.effects <- F
      fixed.effects <- F
      second.effects <- F

      #call bhlmMD
      if ( !is.null( random.formula ) && length( random.formula ) > 0 && !is.all.white( random.formula )[1] )
      {
        random.effects <- T
      }

      if ( !is.null( fixed.formula ) && length( fixed.formula ) > 0 && !is.all.white( fixed.formula )[1] )
      {
        fixed.effects <- T
      }

      if ( !is.null( level2.formula ) && length( level2.formula ) > 0 && !is.all.white( level2.formula )[1] )
      {
        second.effects <- T
      }

      ##"Missing Data\n\n"

      out.bhlm <- bhlmMD( response.formula = response.formula, group = group.formula, 
                  data = data, random.effects = random.effects, fixed.effects = fixed.effects,
                  second.effects = second.effects,
                  prior = prior, likelihood = likelihood, sampler = sampler, 
                  random.seed = random.seed, na.action = na.action, contrasts = contrasts,
                  debug = debug )

    }
    else
    {
      out.bhlm <- bhlmMV( response.formula = response.formula, random.formula = random.formula, 
                  fixed.formula = fixed.formula, level2.formula = level2.formula, 
                  group = group.formula, data = data, 
                  prior = prior, likelihood = likelihood, sampler = sampler, 
                  random.seed = random.seed, na.action = na.action, contrasts = contrasts,
                  debug = debug )
    }

  }

  out.bhlm$call = match.call()
  out.bhlm

}#end bhlm


