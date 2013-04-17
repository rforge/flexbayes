##################################################
# bhlm class in the Bayes Module                 #
# Insightful: NIH Bayes II                       #
# Alejandro Murua 1/03                           #
##################################################




#####################################################################################
##
##  BAYESIAN Missing Data (Response) HIERARCHICAL LINEAR MODEL
##
#####################################################################################

#random.formula, fixed, level2 and group are of type "formula"
bhlmMD <- function( response.formula = NULL, group = NULL, data = sys.parent(),
                    random.effects = T, fixed.effects = F, second.effects = T,
                    prior = bhlm.prior(), 
                    likelihood  = bhlm.likelihood(),
	            sampler = bhlm.sampler(), 
                    random.seed = .Random.seed, na.action = NULL, contrasts = NULL,
                    debug = F )
{

  #get all stats from Gibbs sampling
  print.stats <- 1

  #formula refers to random effects

####-------------------------------------------------------
  ## Get Data 
  ##

  if ( is.null( data ) )
  {
    stop( "bhlm: data must be provided." )
  }

  data <- as.data.frame( data )

  if ( !is.null( na.action ) )
  {
    #get rid of rows with missing data if na.omit, otherwise stop and fail (na.fail)

    data <- na.action( data )
    if ( nrow( data ) == 0 )
    {
      stop( "bhlm: every observation (row) of the data set contains at least a missing value.\n" )
    }
  }
  else
  {
    na.action <- na.fail
  }


  if ( is.null( response.formula ) )
  {
    stop( "You must provide a formula with the response variables\n" )
  }
    
  data.copy <- data
  missing.response <- is.na( data )
  sum.missing.response <- sum( missing.response )
  if ( sum.missing.response > 0 )
  {
    data[ is.na( data ) ] <- 0 
  }

  model.response <- call( "model.frame", formula = response.formula, data = data )
  model.response <- eval( model.response, sys.parent() )
  Terms <- attr( model.response, "terms" )
  Terms@intercept <- 0
  Y <- model.matrix( Terms, model.response, contrasts )
  used.contrasts <- contrasts( Y )
  response.names <- dimnames( Y )[[2]]

  #response.names <- attr( Terms, "term.labels" )

  #sort data by group

  if ( is.null( Y ) )
  {
    stop( "bhlm: response variable must be provided.\n" )
  }


  idx <- seq( 1, nrow( data ) )
  if ( !is.null( group ) )
  {
    group.var <- attr( terms( group ), "term.labels" )

    idx <- order( data[[ group.var ]] )

    Y <- Y[ idx, ]

    #now sort missing response
    if ( sum.missing.response > 0 )
    {
      missing.response <- as.matrix( missing.response[ idx, ] )
    }

    #counts by groups
    n.responses <- table( data[ , group.var ] )
    n.groups <- length( n.responses )

    group.names <- data[[ group.var ]][ idx ] 
    group.names <- group.names[ cumsum( n.responses ) ]

  }#end if group
  else
  {
    n.responses <- length( Y ) / length( response.names )
    n.groups <- 1
    group.names <- 1
  }

  number.data <- length( Y )

  #transpose Y
  Y <- t(Y)
  if ( sum.missing.response > 0 )
  {
    missing.response <- t( missing.response[ , dimnames( Y )[[1]] ] )
  }

  if ( is.vector( Y ) )
  {
    Y <- matrix( Y, nrow = 1 )

    if ( sum.missing.response > 0 )
    {
      missing.response <- matrix( missing.response, nrow = 1 )
    }
  }

  dim.response <- nrow( Y )

  missing.response.components <- 0 
  missing.response.counts <- 0
  total.missing.response <- 0
  if ( sum.missing.response > 0 )
  {
    missing.response.counts <- apply( missing.response, 2, sum )
    total.missing.response <- sum( missing.response.counts )
    if ( total.missing.response > 0 )
    {
      missing.response.components <- apply( missing.response, 2, 
                                            FUN = function( x, p ) 
                                                    { seq( 1, p )[x == T] }, 
                                            p = nrow( missing.response ) )
    
      missing.response.components <- unlist( missing.response.components )
    }
  }





####end Data-------------------------------------------------------


####----------------------------------------------------------------------
  ## Now validate/interpret the arguments for the likelihood 
  ##

  unique.lklhd.Cov <- 1
  dim.error.Cov <- 0
  if ( likelihood$df <= 0.0 ) 
  {
    if ( likelihood$type == "t" )
    {
      warning( paste(  "bhlm: the degrees of freedom for the multivariate-t likelihood must be",
	               "positive.  Using the default value of 3.", sep = "\n" ) )
    }

    likelihood$df <- 3
  }

  if ( likelihood$type == "normal" )
  {
    likelihood$type <- 0
  }
  else if ( likelihood$type == "t" )
  {
    likelihood$type <- 1
  }  
  else if ( likelihood$type == "t groups" )
  {
    likelihood$type <- 2
  }  


  
  if ( is.null( likelihood$errorCov ) ) 
  {
    likelihood$errorCov <- 1
  }
  else if ( is.name( likelihood$errorCov ) )
  {
    errorCov <- as.vector( eval( call( likelihood$errorCov, p = dim.response ) ) )
    if ( n.groups > 1 )
    {
      for ( i in seq( 2, n.groups ) )
      {
        this.errorCov <- as.vector( eval( call( likelihood$errorCov, p = dim.response ) ) )
        errorCov <- c( errorCov, this.errorCov )
      }
    }
    likelihood$errorCov <- errorCov
  }#end name
  else if ( is.function( likelihood$errorCov ) )
  {
    errorCov <- as.vector( likelihood$errorCov( dim.response ) )
    if ( n.groups > 1 )
    {
      for ( i in seq( 2, n.groups ) )
      {
        this.errorCov <- as.vector( likelihood$errorCov( dim.response ) )
        errorCov <- c( errorCov, this.errorCov )
      }
    }
    likelihood$errorCov <- errorCov
  }#end function
  else if ( is.vector( likelihood$errorCov ) ) 
  {

    if ( length( likelihood$errorCov ) == 1 ) 
    {
      if ( likelihood$errorCov <= 0 )
      { 
        stop( "bhlm: The error covariance matrix must be positive." )
      }
    }
    else 
    {
      unique.lklhd.Cov <- 0
      if ( any( likelihood$errorCov ) <= 0 )
      { 
        stop ( "bhlm: The error covariance matrix must be positive." )
      }

      if ( length( likelihood$errorCov ) == n.groups )
      {
        dim.error.Cov <- 0
      }
      else if ( length( likelihood$errorCov ) == dim.response )
      {
        dim.error.Cov <- 1
      }
      else
      { 
	stop( "bhlm: The dimension of the error covariance matrix must equal the number of response variables." )
      }
    }
  }
  else if ( is.list( likelihood$errorCov ) ) 
  {
    unique.lklhd.Cov <- 0
    dim.error.Cov <- 2
    if ( length( likelihood$errorCov ) != n.groups )
    {
      stop( "bhlm: Need ", n.groups, " error variances matrices. Got ", length( likelihood$errorCov ), "." )
    }

    errorCov <- NULL
    for ( i in seq( 1, n.groups ) )
    {
      if ( any( dim( likelihood$errorCov[[ i ]]  ) != dim.response ) )
      {
        stop( "bhlm: The error covariance matrix must be square with dimension equal",
	      "to the number of observations." )
      }
      
      if ( any( likelihood$errorCov[[ i ]] - t( likelihood$errorCov[[ i ]] ) != 0 ) )
      {
        stop( "bhlm: The error covariance matrix is not symmetric." )
      }
      
      if ( i == 1 )
      {
        errorCov <- as.vector( likelihood$errorCov[[ i ]] )
      }
      else
      {
        errorCov <- c( errorCov, as.vector( likelihood$errorCov[[ i ]] ) )
      }
    }#end for loop

    likelihood$errorCov <- errorCov
  }#end list		


####end likelihood----------------------------------------------------------------------



####------------------------------------------------------------------------------------
  ## Validate/interpret the arguments specifying the prior
  ##
  
  valid.prior <- bhlm.prior( error.var = prior$error.var, random.coef = prior$random.coef, fixed.coef = prior$fixed.coef, level2.coef = prior$level2.coef, random.var = prior$random.var, common.error.var = prior$common.error.var )

  level2.coef.mean <- 0
  level2.coef.Cov <- 1
  level2.coef.df <- 0
  level2.coef.type <- 1 #non-informative

  fixed.coef.mean <- 0
  fixed.coef.Cov <- 1
  fixed.coef.df <- 0
  fixed.coef.type <- 1 #non-informative


  random.var.scale <- 1
  random.var.nu <- 0
  random.var.power <- -0.5 #Haar measure
  random.var.type <- 0  #invChisq


  random.coef.mean <- 0
  random.coef.Cov <- 1
  random.coef.df <- 0
  random.coef.type <- 0

  if ( random.effects )
  {
    if ( second.effects )
    {
      if ( valid.prior$level2.coef@name == "non-informative" )
      {
        level2.coef.mean <- rep( 0.0, nrow( Y ) )
        level2.coef.Cov <- diag( nrow( Y ) )
        level2.coef.type <- 1
      }
      else
      {
        level2.coef.mean <- valid.prior$level2.coef@parameters[["mean vector"]]
        level2.coef.Cov <- valid.prior$level2.coef@parameters[["covariance matrix"]]

        level2.coef.mean <- valid.mean.specification( level2.coef.mean, nrow( Y ) )
        level2.coef.Cov <- valid.Cov.specification( level2.coef.Cov, nrow( Y ) )

        if ( valid.prior$level2.coef@name == "t" )
        {
          level2.coef.df <- valid.prior$level2.coef@parameters[["t degrees of freedom"]]
        }
        level2.coef.type <- 0
      }

      ## beta for the full model
      random.coef.type <- 0

      if ( valid.prior$random.coef@name == "non-informative" )
      {
        random.coef.mean <- rep( 0.0, nrow( Y ) )
        random.coef.Cov <- diag( nrow( Y ) )
      }
      else
      {
        #ignore prior mean in this case

        random.coef.mean <- rep( 0.0, nrow( Y ) )
        random.coef.Cov <- valid.prior$random.coef@parameters[["covariance matrix"]]
        random.coef.Cov <- valid.Cov.specification( random.coef.Cov, nrow( Y ) )

        if ( valid.prior$random.coef@name == "t" )
        {
          random.coef.df <- valid.prior$random.coef@parameters[["t degrees of freedom"]]
        }
      }
    }#end second effects
    else
    {
      ##check prior for random coefficients

      ## beta for the mixed effects model
      random.coef.type <- 1
      
      if ( class( valid.prior$random.coef ) == "bayes.distribution" )
      {
        if ( valid.prior$random.coef@name == "non-informative" )
        {
          random.coef.type <- 3
          random.coef.mean <- rep( 0, nrow( Y ) )
          random.coef.Cov <- diag( nrow( Y ) )
        }
        else
        {
          random.coef.mean <- valid.prior$random.coef@parameters[["mean"]]
          random.coef.Cov <- valid.prior$random.coef@parameters[["covariance matrix"]]

          random.coef.mean <- valid.mean.specification( random.coef.mean, nrow( Y ) )
          random.coef.Cov <- valid.Cov.specification( random.coef.Cov, nrow( Y ) )

          if ( valid.prior$random.coef@name == "t" )
          {
            random.coef.df <- valid.prior$random.coef@parameters[["t degrees of freedom"]]
          }
        }
      }#end if a single prior
      else if ( is.list( valid.prior$random.coef ) )
      {
        random.coef.type <- 2

        if ( length( valid.prior$random.coef ) != n.groups )
        {
          stop( "bhlm: list of prior parameters for random coefficients is not of the appropriate size." )
        }

        ##all covariances are assume equal

        if ( valid.prior$random.coef[[1]]@name == "non-informative" )
        {
          random.coef.Cov <- diag( nrow( Y ) )
          random.coef.mean <- rep ( 0, nrow( Y ) )
          random.coef.type <- 1
        }
        else
        {
          random.coef.Cov <- valid.prior$random.coef[[1]]@parameters[["covariance matrix"]]
          random.coef.Cov <- valid.Cov.specification( random.coef.Cov, nrow( Y ) )         

          random.coef.mean <- apply( matrix( seq( 1, length( valid.prior$random.coef ) ), nrow = 1 ), 
                                     MARGIN = 2, 
                                     FUN = function( x, y, dim )
                                     {
                                       if ( y[[ x ]]@name == "non-informative" )
                                       {
                                         r.mean <- rep( 0, dim ) 
                                       }
                                       else
                                       {
                                         r.mean <- y[[ x ]]@parameters[["mean"]]
                                         r.mean <- valid.mean.specification( r.mean, dim )
                                       }
                                       r.mean
                                     },
                                     y = valid.prior$random.coef,
                                     dim = nrow( Y ) )

        }
 
        ## all t's are assume to have the same df's
        if ( valid.prior$random.coef[[1]]@name == "t" )
        {
          random.coef.df <- valid.prior$random.coef[[1]]@parameters[["t degrees of freedom"]]
        }

      }#end several prior parameters

    }#end no second effects

    #if beta has an informative prior
    if ( random.coef.type != 3 )
    {
      random.var.type <- switch( valid.prior$random.var@name,
      "invChisq" =
      {
        random.var.scale <- valid.prior$random.var@parameters[["sigma0.sq"]]
        random.var.nu <- valid.prior$random.var@parameters[["df"]]
        0
      },

      "non-informative power" =
      {
        random.var.power <- valid.prior$random.var@parameters[["power"]]
        1
      },

      "uniform shrinkage" =
      {
        random.var.scale <- valid.prior$random.var@parameters[["median"]]
        2
      },

      "du Mouchel" =
      {
        random.var.scale <- valid.prior$random.var@parameters[["dispersion"]]
        3
      },

      "invWishart" =
      {
        random.var.scale <- valid.Cov.specification( valid.prior$random.var@parameters[["scale"]], nrow( Y ) )
        random.var.nu <- valid.prior$random.var@parameters[["degrees of freedom"]]
        4
      },

      stop( "bhlm: distribution for random coefficients variance must be either \"invChisq\", \"non-informative power\", \"uniform shrinkage\", \"du Mouchel\", or \"invWishart\" \n" )
      )#end switch
    }

  }#end random effects prior



  ##fixed effects prior
  if ( fixed.effects )
  {
    if ( valid.prior$fixed.coef@name == "non-informative" )
    {
      fixed.coef.mean <- rep( 0.0, nrow( Y ) )
      fixed.coef.Cov <- diag( nrow( Y ) )
      fixed.coef.type <- 1
    }
    else
    {
      fixed.coef.mean <- valid.prior$fixed.coef@parameters[["mean vector"]]
      fixed.coef.Cov <- valid.prior$fixed.coef@parameters[["covariance matrix"]]

      fixed.coef.mean <- valid.mean.specification( fixed.coef.mean, nrow( Y ) )
      fixed.coef.Cov <- valid.Cov.specification( fixed.coef.Cov, nrow( Y ) )

      if ( valid.prior$fixed.coef@name == "t" )
      {
        fixed.coef.df <- valid.prior$fixed.coef@parameters[["t degrees of freedom"]]
      }
      fixed.coef.type <- 0
    }
  }#end fixed effects



  ##sigma prior
  error.var.scale <- 1
  error.var.nu <- 0
  error.var.power <- -0.5 #haar measure
  error.var.common <- valid.prior$common.error.var
  error.var.type <- 0


  if ( is.element( error.var.common, c( 1, 2 ) ) )
  {
    # common prior parameters for error variance

    error.var.type <- switch( valid.prior$error.var@name,
      "invChisq" =
      {
        error.var.scale <- valid.prior$error.var@parameters[["sigma0.sq"]]
        error.var.nu <- valid.prior$error.var@parameters[["df"]]
        0
      },

      "non-informative power" =
      {
        error.var.power <- valid.prior$error.var@parameters[["power"]]
        1
      },

      "uniform shrinkage" =
      {
        error.var.scale <- valid.prior$error.var@parameters[["median"]]
        2
      },

      "du Mouchel" =
      {
        error.var.scale <- valid.prior$error.var@parameters[["dispersion"]]
        3
      },

      "invWishart" =
      {
        error.var.scale <- valid.Cov.specification( valid.prior$error.var@parameters[["scale"]], nrow( Y ) )
        error.var.nu <- valid.prior$error.var@parameters[["degrees of freedom"]]
        4
      },

      "mass point" =
      {
        error.var.scale <- valid.prior$error.var@parameters[["value"]]
        error.var.scale <- valid.Cov.specification( error.var.scale, nrow( Y ) )
        error.var.nu <- nrow( Y )
        5
      },


      stop( "bhlm: distribution for error variance must be either \"invChisq\", \"non-informative power\", \"uniform shrinkage\", \"du Mouchel\", \"invWishart\", or \"mass point (known)\" \n" )
    )#end switch
  }
  else
  {
    # different prior parameters for error variance

    if ( length( valid.prior$error.var ) != n.groups )
    {
      stop( "bhlm: list of prior parameters for error variance is not of the appropriate size." )
    }

    if ( is.list( valid.prior$error.var ) )
    {
      error.var.scale <- vector( "list", length( valid.prior$error.var ) )
      error.var.nu <- rep( 0, length( valid.prior$error.var ) )
      error.var.power <- rep( 0, length( valid.prior$error.var ) )

      for ( i in seq( 1, length( valid.prior$error.var ) ) )
      {
        error.var.type <- switch( valid.prior$error.var[[i]]@name,
          "invChisq" =
          {
            error.var.scale[[i]] <- valid.prior$error.var[[i]]@parameters[["sigma0.sq"]]
            error.var.nu[i] <- valid.prior$error.var[[i]]@parameters[["df"]]
            0
          },

          "non-informative power" =
          {
            error.var.power[i] <- valid.prior$error.var[[i]]@parameters[["power"]]
            1
          },

          "uniform shrinkage" =
          {
            error.var.scale[[i]] <- valid.prior$error.var[[i]]@parameters[["median"]]
            2
          },

          "du Mouchel" =
          {
            error.var.scale[[i]] <- valid.prior$error.var[[i]]@parameters[["dispersion"]]
            3
          },

          "invWishart" =
          {
            error.var.scale[[i]] <- valid.Cov.specification( valid.prior$error.var[[i]]@parameters[["scale"]], nrow( Y ) )
            error.var.nu[i] <- valid.prior$error.var[[i]]@parameters[["degrees of freedom"]]
            4
          },

          "mass point" =
          {
            error.var.scale[[i]] <- valid.prior$error.var[[i]]@parameters[["value"]]
            error.var.scale[[i]] <- valid.Cov.specification( error.var.scale[[i]], nrow( Y ) )
            error.var.nu[i] <- nrow( Y )
            5
          },


          stop( "bhlm: distribution for error variance must be either \"invChisq\", \"non-informative power\", \"uniform shrinkage\", \"du Mouchel\", \"invWishart\", or \"mass point (known)\" \n" )
        )#end switch

        if ( i == 1 )
        {
          prev.type <- error.var.type
        }
        else if ( prev.type != error.var.type )
        {
          stop( "bhlm: all priors for error variances should be of the same type." )
        }
        
      }#end for loop

      #now convert error.var.scale to matrix
      if ( !is.null( error.var.scale ) && length( error.var.scale ) > 0 )
      {
        error.scale <- error.var.scale[[1]]

        if ( length( valid.prior$error.var ) > 1 )
        {
          for ( i in seq( 2, length( valid.prior$error.var ) ) )
          {
            error.scale <- cbind( error.scale, error.var.scale[[i]] )
          }
        }
        error.var.scale <- error.scale
      }
    }#end list
    else
    {
      stop( "bhlm: prior for error variance is not valid." )
    }
  }#end different prior parameters  





####end Priors----------------------------------------------------------------------------


####-----------------------------------------------------------------
  ## Get the seed
  ##

  set.seed( random.seed )

####end Seed-----------------------------------------------------------------



####--------------------------------------------------------------
  ##Get the control arguments
  ##

  burnInLength <- sampler$control$bSize
  simulationsToPerform <- sampler$control$simSize
  sampleFrequency <- sampler$control$freqSize

####end Control-----------------------------------------------------



####-------------------------------------------------------------------
  ##Finally get the starting point for the simulation
  ##

  if ( random.effects )
  {
    dim.random <- nrow( Y )
  }
  else
  {
    dim.random <- 0
  }

  if ( fixed.effects )
  {
    dim.fixed <- nrow( Y )
  }
  else
  {
    dim.fixed <- 0
  }

  if ( second.effects )
  {
    dim.level2 <- nrow( Y )
  }
  else
  {
    dim.level2 <- 0
  }


  ## flag for default starting points
  number.chains <- sampler$number.chains
  read.init.point <- 0
  init.point <- sampler$init.point
  if ( is.null( init.point ) || init.point$type == "prior" )
  {
    ## flag for specify/mixture starting points
    read.init.point = 1

    # draw number.chains initial values at random from "prior"
    starting.points <- generateInitPoints.bhlmMD( number.chains, n.responses, Y, dim.response,
                                                  dim.random, dim.fixed, dim.level2,
                                                  error.var.nu, error.var.scale, error.var.type, error.var.common,
                                                  random.coef.df, random.coef.mean, random.coef.Cov, random.coef.type,
                                                  fixed.coef.df, fixed.coef.mean, fixed.coef.Cov,
                                                  random.var.nu, random.var.scale, random.var.type,
                                                  level2.coef.df, level2.coef.mean, level2.coef.Cov,
                                                  mix.with.MLE = 0 )

  }#end prior init points
  else if ( !is.null( init.point ) && init.point$type == "user's choice" )
  {
    read.init.point <- 1

    starting.points <- validate.initial.points.bhlm( n.groups, dim.response, dim.random, dim.fixed, dim.level2, init.point$values[[ 1 ]], random.effects, fixed.effects, second.effects, error.var.common, random.coef.type, error.var.type, random.var.type )

    if ( !is.null( starting.points$random.var ) )    
    {
      s.random.var <- vector( "list", number.chains )
      s.random.var[[1]] <- starting.points$random.var
    }

    if ( !is.null( starting.points$error.var ) )    
    {
      s.error.var <- vector( "list", number.chains )
      s.error.var[[1]] <- starting.points$error.var
    }

    if ( number.chains > 1 )
    {
      for ( i in seq( 2, number.chains ) )
      {
        s.points <- validate.initial.points.bhlm( n.groups, dim.response, dim.random, dim.fixed, dim.level2, init.point$values[[ i ]], random.effects, fixed.effects, second.effects, error.var.common, random.coef.type, error.var.type, random.var.type )
        if ( !is.null( starting.points$error.var ) )
        {
          s.error.var[[i]] <- s.points$error.var
        }

        if ( !is.null( starting.points$random.coef ) )
        {
          starting.points$random.coef <- cbind( starting.points$random.coef, s.points$random.coef )
        }

        if ( !is.null( starting.points$fixed.coef ) )
        {
          starting.points$fixed.coef <- cbind( starting.points$fixed.coef, s.points$fixed.coef )
        }

        if ( !is.null( starting.points$level2.coef ) )
        {
          starting.points$level2.coef <- cbind( starting.points$level2.coef, s.points$level2.coef )
        }

        if ( !is.null( starting.points$random.var ) )
        {
          s.random.var[[i]] <- s.points$random.var
        }

      }#end for chains

      
      if (  !is.null( starting.points$random.var ) )
      {
        starting.points$random.var <- s.random.var
      }

      if ( !is.null( starting.points$error.var ) )
      {
        starting.points$error.var <- s.error.var
      }


    }#end if chains > 1
    else
    {
      if (  !is.null( starting.points$random.var ) )
      {
        starting.points$random.var <- s.random.var
      }

      if ( !is.null( starting.points$error.var ) )
      {
        starting.points$error.var <- s.error.var
      }     
    }
    ##make sure no element is null



    if ( is.null( starting.points$random.var ) )
    {
      starting.points$random.var <- vector( "list", number.chains )
      for ( i in seq( 1, number.chains ) )
      {
        starting.points$random.var[[i]] <- 0
      }
    }

    if ( is.null( starting.points$random.coef ) )
    {
      starting.points$random.coef <- matrix( 0, ncol = number.chains )
    }

    if ( is.null( starting.points$level2.coef ) )
    {
      starting.points$level2.coef <-  matrix( 0, ncol = number.chains )
    }

    if ( is.null( starting.points$fixed.coef ) )
    {
      starting.points$fixed.coef <-  matrix( 0, ncol = number.chains )
    }

    if ( is.null( starting.points$error.var ) )
    {
      starting.points$error.var <- vector( "list", number.chains )
      for ( i in seq( 1, number.chains ) )
      {
        starting.points$error.var[[i]] <- 0
      }
    }
  }#end user's choice
  else if ( !is.null( init.point ) )
  {
    stop( "bhlm: initialization procedure [", init.point$type, "] has not been implemented yet." )
  }


####end Initial Points------------------------------------------------------------------------------------------


####-------------------------------------------------------------------------------------------
  ##Fit the model
  ##

  if ( is.null( sampler$sampler ) )
  {
    sampler.type <- 0
  }
  else if ( sampler$sampler == "Gibbs" )
  {
    sampler.type <- 0
  }
  else
  {
    stop( "bhlm: sampler type [", sampler$sampler, "] has not been implemented." )
  }  

  #prepare formula outputs
  random.formula <- NULL
  fixed.formula <- NULL
  level2.formula <- NULL
  if ( random.effects )
  {
    random.formula <- ~ 1
  }
  if ( fixed.effects )
  {
    fixed.formula <- ~ 1
  }
  if ( second.effects )
  {
    level2.formula <- ~ 1
  }


  bhlmodel <- vector( "list", number.chains )


  for ( i in seq( 1, number.chains ) )
  {
    simulation.seed <- .Random.seed
    fit <- fit.bayeshlmMD( n.groups, idx, n.responses, dim.response, dim.random, dim.fixed, dim.level2, group.names, Y, response.names, total.missing.response, missing.response.counts, missing.response.components, unique.lklhd.Cov, dim.error.Cov, likelihood$errorCov, likelihood$df, likelihood$type, random.coef.mean, random.coef.Cov, random.coef.df, random.coef.type, fixed.coef.mean, fixed.coef.Cov, fixed.coef.df, fixed.coef.type, level2.coef.mean, level2.coef.Cov, level2.coef.df, level2.coef.type, error.var.nu, error.var.scale, error.var.power, error.var.common, error.var.type, random.var.nu, random.var.scale, random.var.power, random.var.type, read.init.point, starting.points$random.coef[,i], starting.points$fixed.coef[,i], starting.points$level2.coef[,i], starting.points$error.var[[i]], starting.points$random.var[[i]], sampler.type, burnInLength, simulationsToPerform, sampleFrequency, print.stats )
    bhlmodel[[i]] <- fit

  }#end for chains

	return(posterior(mcmc.list(bhlmodel)))
}#end



############################################################################################
##
## FIT BAYES HIERARCHICAL LINEAR MODEL
##
############################################################################################

fit.bayeshlmMD <- function( n.groups, sorted.by.group.idx, n.responses, dim.response, dim.random, dim.fixed, dim.level2, 
                          group.names, Y, response.names, 
                          total.missing.response, missing.response.counts, missing.response.components,
                          unique.lklhd.Cov, dim.error.Cov, likelihood.errorCov, likelihood.df, likelihood.type, 
                          random.coef.mean, random.coef.Cov, random.coef.df, random.coef.type, 
                          fixed.coef.mean, fixed.coef.Cov, fixed.coef.df, fixed.coef.type, 
                          level2.coef.mean, level2.coef.Cov, level2.coef.df, level2.coef.type, 
                          error.var.nu, error.var.scale, error.var.power, error.var.common, error.var.type, 
                          random.var.nu, random.var.scale, random.var.power, random.var.type, 
                          read.init.point, starting.random.coef, 
                          starting.fixed.coef, starting.level2.coef, 
                          starting.error.var, starting.random.var, 
                          sampler.type = 0, burnInLength, simulationsToPerform, sampleFrequency, 
                          print.stats )
{

  #if there are t distributions then report their augmented data tau ~ IvChisq( nu, 1 )
  #keep_tau2_errors <- 1

  number.data <- sum( n.responses )

  #now count how many variables there are
  n.vars <- 0
  output.dim <- 0

  if (  error.var.type != 5 )
  {
    if ( error.var.common == 0 || error.var.common == 1 )
    {
      #different error variance for each group
      n.vars <- n.groups
    }
    else
    {
      #only one common error var
      n.vars <- 1
    }

    if ( error.var.type == 4 )
    {
      dim.error.var <- dim.response * dim.response
    }
    else
    {
      dim.error.var <- 1
    }     
    output.dim <- n.vars * dim.error.var
  }


  #now account for other variables including t distributed variables
  if ( dim.random > 0 )
  {
    #there are random effects
    if ( random.coef.df == 0 )
    {
      n.vars <- n.vars + n.groups
      output.dim <- output.dim + n.groups * dim.random 
    }
    else
    {
      n.vars <- n.vars + 2 * n.groups
      output.dim <- output.dim + n.groups * ( dim.random + 1 )
    }


    if ( dim.level2 > 0 )
    {
      #there are second effects
      if ( level2.coef.df == 0 )
      {
        n.vars <- n.vars + 1
        output.dim <- output.dim + dim.level2
      }
      else
      {
        n.vars <- n.vars + 2
        output.dim <- output.dim + dim.level2 + 1 
      }
    }

    #if beta has an informative prior
    if ( random.coef.type != 3 )
    {
      #account for random coef variance
      n.vars <- n.vars + 1
      if ( random.var.type == 4 )
      {
        dim.random.var <- dim.random * dim.random
      }
      else
      {
        dim.random.var <- 1
      }
      output.dim <- output.dim + dim.random.var
    }
  }

  if ( dim.fixed > 0 )
  {
    #there are fixed effects
    if ( fixed.coef.df == 0 )
    {
      n.vars <- n.vars + 1
      output.dim <- output.dim + dim.fixed
    }
    else
    {
      n.vars <- n.vars + 2
      output.dim <- output.dim + dim.fixed + 1
    }

  }

  #account for t distributed errors
  if ( likelihood.type == 1 )
  {
    #t likelihood
    n.vars <- n.vars + number.data
    output.dim <- output.dim + number.data
  }
  else if ( likelihood.type == 2 )
  {
    #t group likelihood
    n.vars <- n.vars + n.groups
    output.dim <- output.dim + n.groups
  }

  #account for missing response
  if ( total.missing.response > 0 )
  {
    n.vars <- n.vars + sum( missing.response.counts > 0 ) 
    output.dim <- output.dim + total.missing.response
  }

      
  #now create output array with the right dimensions
  output.samples <- matrix( 0, nrow = simulationsToPerform, ncol = output.dim )

  if ( print.stats == 1 )
  {
    gibbs.stats <- rep( 0, n.vars )
  }
  else
  {
    gibbs.stats <- 0
  }

  #call the function

  fit <- .C( "fitBayesianMDHLM",
             as.integer( n.groups ),
             as.integer( n.responses ),
             as.integer( dim.response ), 
             as.integer( dim.random ), 
             as.integer( dim.fixed ), 
             as.integer( dim.level2 ), 
             #
             as.double( Y ), 
             #
             as.integer( total.missing.response ),
             as.integer( missing.response.counts ),
             as.integer( missing.response.components ),
             #
             as.integer( unique.lklhd.Cov ), 
             as.integer( dim.error.Cov ), 
             as.double( likelihood.errorCov ),
             #
             as.double( likelihood.df ), 
             as.integer( likelihood.type ), 
             #
             as.double( fixed.coef.mean ), 
             as.double( fixed.coef.Cov ), 
             as.double( fixed.coef.df ), 
             as.integer( fixed.coef.type ),
             #
             as.double( random.coef.mean ), 
             as.double( random.coef.Cov ), 
             as.double( random.coef.df ), 
             as.integer( random.coef.type ), 
             #
             as.double( level2.coef.mean ), 
             as.double( level2.coef.Cov ), 
             as.double( level2.coef.df ), 
             as.integer( level2.coef.type ),
             #
             as.double( error.var.nu ), 
             as.double( error.var.scale ), 
             as.double( error.var.power ),
             as.integer( error.var.common ),
             as.integer( error.var.type ),
             #
             as.double( random.var.nu ), 
             as.double( random.var.scale ),
             as.double( random.var.power ),
             as.integer( random.var.type ),
             #
             as.integer( read.init.point ), 
             as.double( starting.random.coef ),
             as.double( starting.fixed.coef ),
             as.double( starting.level2.coef ),
             as.double( starting.error.var ),
             as.double( starting.random.var ),
             #
             as.integer( burnInLength ), 
             as.integer( simulationsToPerform ), 
             as.integer( sampleFrequency ), 
             #
             as.integer( print.stats ),
             output.samples = as.double( output.samples ),
             gibbs.stats = as.double( gibbs.stats ) )


  # OUTPUT
  # ------
  # The simulationsToPerform output simulations are returned in the array output_simulations.
  # This should be an array of length 
  #   J * [p] + [r]  + q + < 1 | J > + < 1 | q*q > + [J] + [1] + [1] + [J] + [n] ) * simulationsToPerform
  # where 
  # J: number of groups
  # p: dimension of beta vector
  # r: dimension of gamma vector
  # q: dimension of alpha vector
  # < 1 | J >: sigma2 (either a common sigma2 or J independent sigma2's)
  # < 1 | q*q >: dimension of tau2 (either a random variable or a random matrix )
  # J: augmented variables for t-distributed betas
  # 1: augmented variable for t-distributed gamma
  # 1: augmented variable for t-distributed alpha
  # J: augmented variables for t-distributed group errors
  # n: number of observations. augmented variables for t-distributed errors.
  #
  # < .|. > denotes a choice.
  # [.] denotes an optional argument.
  #
  # The above dimensions specification also denotes the order in which the simulations are returned:
  #       beta (by group), gamma , alpha, sigma, tau | tau2, 
  #             t_beta (by group), t_gamma, t_alpha, t_groups, t_errors.
  #
  #
  # Note: Simulations for sigma (not sigma2) are returned.
  #       Simulations for tau (not tau2) are returned if tau2 is a scalar,
  #       otherwise, the matrices tau2 are returned.
  #

  out.samples <- matrix( fit$output.samples, nrow = simulationsToPerform, ncol = output.dim )

  Y.imputed <- Y
  group.names <- as.character( group.names )
  
  end.index <- 0
  if ( dim.random > 0 )
  {
    random.coef <- matrix( out.samples[ , seq( 1, n.groups * dim.random ) ], nrow = simulationsToPerform )

    random.names <- as.vector( outer( outer( c("Mean" ), response.names, paste, sep = ":" ), group.names, paste, sep = ":" ) )

    #random.names <- as.vector( outer( response.names, group.names, paste, sep = ":" ) )

    end.index <- n.groups * dim.random
  }

  if ( dim.fixed > 0 )
  {
    fixed.coef <- matrix( out.samples[ , end.index + seq( 1, dim.fixed ) ], nrow = simulationsToPerform )
    fixed.names <- response.names
    end.index <- end.index + dim.fixed
  }

  if ( dim.random > 0 && dim.level2 > 0 )
  {
    level2.coef <- matrix( out.samples[ , end.index + seq( 1, dim.level2 ) ], nrow = simulationsToPerform )
    level2.names <- response.names

    end.index <- end.index + dim.level2
  }


  if ( error.var.type != 5 )
  {
    if ( error.var.common == 0 ||  error.var.common == 1 )
    {
      #different error variance for each group
      if ( error.var.type != 4 )
      {
        error.var <- matrix( out.samples[ , end.index + seq( 1, n.groups ) ], nrow = simulationsToPerform )
        error.var.names <- as.vector( outer( "MEASUREMENT ERROR [ SIGMA ]:", group.names, paste, sep = "" ) )
        end.index <- end.index + n.groups
      }
      else
      {
        #first group
        error.var <- matrix( out.samples[ , end.index + seq( 1, dim.response * dim.response ) ], nrow = simulationsToPerform )
        end.index <- end.index + dim.response * dim.response
 
        error.var.names <- c( as.vector( t( outer( outer( paste( "MEASUREMENT ERROR [ SIGMA ]:", group.names, sep = "" ), seq( 1, dim.response ), paste, sep = ":" ), seq( 1, dim.response ), paste, sep = "." ) ) ) )

        #remainder groups 
        if ( n.groups > 1 )
        {
          for ( i in seq( 2, n.groups ) )
          {
            error.var <- cbind( error.var, matrix( out.samples[ , end.index + seq( 1, dim.response * dim.response ) ], nrow = simulationsToPerform ) )
            end.index <- end.index + dim.response * dim.response
          }
        }
      }
    }
    else
    {
      if ( error.var.type != 4 )
      {    
        #same error variance prior
        error.var <- matrix( out.samples[ , end.index + 1 ], nrow = simulationsToPerform )
        error.var.names <- "MEASUREMENT ERROR [ SIGMA ]"
        end.index <- end.index + 1
      }
      else
      {
        error.var <- matrix( out.samples[ , end.index + seq( 1, dim.response * dim.response ) ], nrow = simulationsToPerform )
        error.var.names <- as.vector( t( outer( paste( "MEASUREMENT ERROR [ SIGMA ]:", seq( 1, dim.response ), sep = "" ), seq( 1, dim.response ), paste, sep = "." ) ) )
        end.index <- end.index + dim.response * dim.response
      }
    }
  }
  else
  {
    error.var <- NULL
    error.var.names <- NULL
  }


  #if beta has an informative prior
  if ( dim.random > 0 && random.coef.type != 3 )
  {
    if ( random.var.type == 4 )
    {
      random.var <- matrix( out.samples[ , end.index + seq( 1, dim.random * dim.random ) ], nrow = simulationsToPerform )
      random.var.names <- as.vector( t( outer( paste( "RANDOM:TAU:", seq( 1, dim.random), sep = "" ), seq( 1, dim.random ), paste, sep = "." ) ) )
      end.index <- end.index + dim.random * dim.random

    }
    else
    {
      random.var <- matrix( out.samples[ , end.index + 1 ], nrow = simulationsToPerform )
      random.var.names <- "RANDOM:TAU"
      end.index <- end.index + 1
    }
  }
  else
  {
    random.var <- NULL
    random.var.names <- NULL
  }


  if ( dim.random > 0 && random.coef.df > 0 )
  {
    random.coef.tau <- matrix( out.samples[ , end.index + seq( 1, n.groups ) ], nrow = simulationsToPerform )
    end.index <- end.index + n.groups
  }

  if ( dim.fixed > 0 && fixed.coef.df > 0 )
  {
    fixed.coef.tau <- matrix( out.samples[ , end.index + 1 ], nrow = simulationsToPerform )
    end.index <- end.index + 1
  }

  if ( dim.random > 0 )
  {
    if ( dim.level2 > 0 && level2.coef.df > 0 )
    {
      level2.coef.tau <- matrix( out.samples[ , end.index + 1 ], nrow = simulationsToPerform )
      end.index <- end.index + 1
    }
  }


  likelihood.tau <- NULL
  likelihood.group.tau <- NULL
  if ( likelihood.type == 1 )
  {
    #t likelihood
    likelihood.tau <- matrix( out.samples[ , end.index + seq( 1, number.data ) ], nrow = simulationsToPerform )
    end.index <- end.index + number.data
    t.names <- paste( "TAU:ERROR:", sorted.by.group.idx, sep = "" )
  }
  else if ( likelihood.type == 2 )
  {
    #t group likelihood
    likelihood.group.tau <- matrix( out.samples[ , end.index + seq( 1, n.groups ) ], nrow = simulationsToPerform )
    end.index <- end.index + n.groups
    t.names <- as.vector( outer( "TAU:GROUP:ERROR:", group.names, paste, sep = "" ) )
  }

  #missing data
  if ( total.missing.response > 0 )
  {
    imputed.values <- matrix( out.samples[ , end.index + seq( 1, total.missing.response ) ], nrow = simulationsToPerform )
    imputed.names <- NULL
    obs.missing <- seq( 1, length( missing.response.counts ) )[ missing.response.counts > 0 ]
    k <- 1
    for ( i in seq( 1, length( obs.missing ) ) )
    {
      for ( j in seq( 1, missing.response.counts[ obs.missing[i] ] ) )
      {
        imputed.names <- c( imputed.names, paste( "IMPUTED", obs.missing[i], response.names[ missing.response.components[ k ] ], sep = ":" ) )

        #complete data
        Y.imputed[ response.names[ missing.response.components[ k ] ], obs.missing[i] ] <- mean( imputed.values[ , k ] )

        k <- k + 1
      }
    }
  }

  t.likelihood <- NULL
  level2.effects <- NULL
  random.effects <- NULL
  fixed.effects <- NULL
  imputed.components <- NULL

  if ( dim.random > 0 )
  {
    if ( dim.level2 > 0 )
    {
      if ( level2.coef.df > 0 )
      {
        level2.tau <- list( tau = level2.coef.tau, names = "TAU:LEVEL2" )
      }
      else
      {
        level2.tau <- NULL
      }
      level2.effects <- list( dim = dim.level2, names = level2.names, coef = level2.coef, tau = level2.tau )
    }

    if ( random.coef.df > 0 )
    {
      random.tau <- list( tau = random.coef.tau, names = paste( "TAU:RANDOM:", seq( 1, n.groups ), sep = "" ) )
    }
    else
    {
      random.tau <- NULL
    }

    random.effects <- list( dim = dim.random, names = random.names, orig.names = "Mean", coef = random.coef, scale = random.var, scale.names = random.var.names, tau = random.tau, level2 = level2.effects )
  }

  if ( dim.fixed > 0 )
  {
    if ( fixed.coef.df > 0 )
    {
      fixed.tau <- list( tau = fixed.coef.tau, names = "TAU:FIXED" )
    }
    else
    {
      fixed.tau <- NULL
    }
    fixed.effects <- list( dim = dim.fixed, names = fixed.names, coef = fixed.coef, tau = fixed.tau )
  }

  if ( is.element( likelihood.type, c( 1, 2 ) ) )
  {
    t.likelihood <-list( errors = likelihood.tau, groups = likelihood.group.tau, names = t.names )
  }


  if ( error.var.type != 5 )
  {
    scale <- list( scale = error.var, names = error.var.names, common = error.var.common, type = error.var.type )
  }
  else
  {
    scale <- NULL
  }


  if ( total.missing.response > 0 )
  {
    imputed.components <- list( values = imputed.values, names = imputed.names, data = t( Y.imputed ) )
  }

  ######## Create a named matrix of the posterior samples.

	post.samples <- matrix(nrow=simulationsToPerform, ncol=0)
  
  # First the fixed effects
  if ( dim.fixed > 0 ){
  	dimnames(fixed.coef) <- list(NULL, fixed.names)
		post.samples <- cbind(post.samples, fixed.coef)
		if ( fixed.coef.df > 0 ){
      fixed.tau.mat <- matrix( fixed.coef.tau, ncol=1, dimnames = list(NULL, c("TAU:FIXED")) )
      post.samples <- cbind(post.samples, fixed.tau.mat)
    }
	}
	
	# the random effects
	if ( dim.random > 0 ){
	  dimnames(random.coef) <- list(NULL, random.names)
		post.samples <- cbind(post.samples, random.coef)
		if (!is.null(random.var)){
		  dimnames(random.var) <- list(NULL, random.var.names)
		  post.samples <- cbind(post.samples, random.var)
		}
		if ( random.coef.df > 0 )
    {
      dimnames(random.coef.tau) <- list(NULL, paste( "TAU:RANDOM:", seq( 1, n.groups ), sep = "" ))
      post.samples <- cbind(post.samples, random.coef.tau)
    }
    # level 2 effects
    if ( dim.level2 > 0 )
    {
      dimnames(level2.coef) <- list(NULL, level2.names)
      post.samples <- cbind(post.samples, level2.coef)
      if ( level2.coef.df > 0 )
      {
        level2.tau.mat <- matrix( level2.coef.tau, ncol=1, dimnames = list(NULL, c("TAU:LEVEL2")) )
        post.samples <- cbind(post.samples, level2.tau.mat)
      }
    }
	}

  # The "scale" parameter
  if ( error.var.type != 5 )
  {
  	dimnames(error.var) <- list(NULL, error.var.names)
  	post.samples <- cbind(post.samples, error.var)
  }

  # Imputed values
  if ( total.missing.response > 0 )
  {
    dimnames( imputed.values ) <- list(NULL, imputed.names)
    post.samples <- cbind(post.samples, imputed.values)
  }
	
  if ( likelihood.type == 1 )
  {
    #t likelihood
    dimnames(likelihood.tau) <- list(NULL, t.names)
    post.samples <- cbind(post.samples, likelihood.tau)
  }
  else if ( likelihood.type == 2 )
  {
    #t group likelihood
    dimnames(likelihood.group.tau) <- list(NULL, t.names)
    post.samples <- cbind(post.samples, likelihood.group.tau)
  }

  return(mcmc(post.samples))

}#end



