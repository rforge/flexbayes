##################################################
# bhlm class in the Bayes Module                 #
# Insightful: NIH Bayes II                       #
# Alejandro Murua 1/03                           #
##################################################


#####################################################################################
##
##  BAYESIAN HIERARCHICAL LINEAR MODEL
##
#####################################################################################

#random.formula, fixed, level2 and group are of type "formula"
bhlmUV <- function( response.formula = NULL, random.formula = NULL, fixed = NULL, level2 = NULL, group = NULL, data = sys.parent(),
                  prior = bhlm.prior(), 
                  likelihood  = bhlm.likelihood(),
	                sampler = bhlm.sampler(), 
                  random.seed = .Random.seed, na.action = NULL, contrasts = NULL,
                  debug = FALSE )
{

  #get all stats from Gibbs sampling
  print.stats <- 1

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


  missing.data <- is.na( data )
  sum.missing.data <- sum( missing.data )
  if ( sum.missing.data > 0 )
  {
    data[ is.na( data ) ] <- 0 
  }


  #name of variables must be preserved for output
  response.names <- NULL
  random.names <- NULL
  fixed.names <- NULL
  level2.names <- NULL

  X <- 0
  M <- 0
  Z <- 0

  random.effects <- FALSE
  fixed.effects <- FALSE
  random.vars <- -1
  fixed.vars <- -1

  if ( !is.null( random.formula ) )
  {
    attrib <- attributes(terms(random.formula, data = data))

    #if ( length( terms( random.formula )@term.labels ) == 0 && terms( random.formula )@response == 0 )
    if ( length( attrib$term.labels ) == 0 && attrib$response == 0 )
    {
      #may contain only the intercept
      if ( attrib$intercept != 1 )
      {
        stop( "bhlm: formula for random effects must contain at least one variable or intercept." )
      }
      else
      {
        #formula contains the intercept
        random.vars <- 0
      }
    }
    else if ( length( attrib$term.labels ) > 0 && attrib$response == 1 )
    {
      #formula contains the response variable, and predictors other than intercept
      random.vars <- 1
    }
    #else if ( length( terms( random.formula )@term.labels ) == 0 )
    else if ( length( attrib$term.labels ) == 0 )
    {
      #formula contains response and perhaps intercept
      if ( attrib$intercept != 1 )
      {
        stop( "bhlm: formula for random effects must contain at least one variable or intercept." )
      }
      else
      {
        #formula contains the intercept
        random.vars <- 2
      }
    }
    else
    {
      #no response: only predictors (other than just intercept)
      random.vars <- 3
    }
  }

  if ( !is.null( fixed ) )
  {
    attrib <- attributes(terms(fixed, data = data))

    #if ( length( terms( fixed )@term.labels ) == 0 && terms( fixed )@response == 0 )
    if ( length( attrib$term.labels ) == 0 && attrib$response == 0 )
    {
      #may contain only the intercept
      #if ( terms( fixed )@intercept != 1 )
      if ( attrib$intercept != 1 )
      {
        stop( "bhlm: formula for fixed effects must contain at least one variable or intercept." )
      }
      else
      {
        #formula contains the intercept
        fixed.vars <- 0
      }
    }
    #else if ( length( terms( fixed )@term.labels ) > 0 && terms( fixed )@response == 1 )
    else if ( length( attrib$term.labels ) > 0 && attrib$response == 1 )
    {
      #formula contains the response variable (and predictors other than just intercept)
      fixed.vars <- 1
    }
    #else if ( length( terms( fixed )@term.labels ) == 0 )
    else if ( length( attrib$term.labels ) == 0 )
    {
      #formula contains response and perhaps intercept
      #if ( terms( fixed )@intercept != 1 )
      if ( attrib$intercept != 1 )
      {
        stop( "bhlm: formula for fixed effects must contain at least one variable or intercept." )
      }
      else
      {
        #formula contains the intercept
        fixed.vars <- 2
      }
    }
    else
    {
      fixed.vars <- 3
    }
  }
    

#used.contrasts <- NULL

  response.in.response <- FALSE
  if ( !is.element( random.vars, c( 1, 2 ) ) && !is.element( fixed.vars, c( 1, 2 ) ) )
  {
    #check response.formula
    if ( !is.null( response.formula ) && length( response.formula ) > 0 )
    {
      attrib <- attributes(terms(response.formula, data = data))

      #if ( length( terms( response.formula )@term.labels ) == 1 )
      if ( length( attrib$term.labels ) == 1 )
      {
        #found response variable
        response.in.response <- TRUE
      }
      else
      {
        stop( "bhlm: model specification is not valid. You need to provide a valid response variable." ) 
      }
    }
    else
    {
      stop( "bhlm: model specification is not valid. You need to provide a response variable." )  
    }
  }

  if ( random.vars >= 0 )
  {
    if ( random.vars > 0 )
    {
      model.random <- call( "model.frame", formula = random.formula, data = data, na.action = na.action )
      model.random <- eval( model.random, sys.parent() )
      Terms <- attr( model.random, "terms" )
  
      if ( random.vars != 2 )
      {
        X <- model.matrix( Terms, model.random, contrasts )
        #used.contrasts <- contrasts( X )
        random.names <- dimnames( X )[[2]]

      }
      else
      {
        #only intercept as predictor
        X <- matrix( rep( 1, nrow( data ) ), nrow = nrow( data ) )
        random.names <- "(Intercept)"
      }

      if ( !response.in.response && is.element( random.vars, c( 1, 2 ) ) )
      {
        Y <- model.extract( model.random, response )
        response.names <- dimnames( attr( Terms, "factors" ) )[[1]][1]
      }

      random.effects <- TRUE
    }
    else #random effects contain only the intercept
    {
      X <- matrix( rep( 1, nrow( data ) ), nrow = nrow( data ) )
      random.effects <- TRUE
      random.names <- "(Intercept)"
    }

    dimnames( X ) <- list( dimnames( data )[[1]], random.names )

  }#end if random formula
  #cat( "X is \n")
  #print(X)



  if ( fixed.vars >= 0 )
  {
    if ( fixed.vars > 0 )
    {
      model.fixed <- call( "model.frame", formula = fixed, data = data, na.action = na.action )
      model.fixed <- eval( model.fixed, sys.parent() )
      Terms <- attr( model.fixed, "terms" )

      if ( fixed.vars != 2 )
      {
        M <- model.matrix( Terms, model.fixed, contrasts )
        #used.contrasts <- c( used.contrasts, contrasts( M ) )

        fixed.names <- dimnames( M )[[2]]

      }
      else
      {
        M <- matrix( rep( 1, nrow( data ) ), nrow = nrow( data ) )
        fixed.names <- "(Intercept)"
      }

      if ( !response.in.response && is.element( fixed.vars, c( 1, 2 ) ) && !is.element( random.vars, c( 1, 2 ) ) )
      {
        Y <- model.extract( model.fixed, response ) 
        response.names <- dimnames( attr( Terms, "factors" ) )[[1]][1]
      }

      fixed.effects <- TRUE
    }    
    else #fixed effects contain only the intercept
    {
      M <- matrix( rep( 1, nrow( data ) ), nrow = nrow( data ) )
      fixed.effects <- TRUE
      fixed.names <- "(Intercept)"
    }

    dimnames( M ) <- list( dimnames( data )[[1]], fixed.names )

  }#end if fixed formula



  if ( response.in.response )
  {
    model.response <- call( "model.frame", formula = response.formula, data = data, na.action = na.action )
    model.response <- eval( model.response, sys.parent() )
    Terms <- attr( model.response, "terms" )
    attributes(Terms)$intercept <- 0
    Y <- model.matrix( Terms, model.response, contrasts )
    #used.contrasts <- c( used.contrasts, contrasts( Y ) )

    response.names <- dimnames( Y )[[2]]

  }

  second.effects <- FALSE
  if ( !is.null( level2 ) )
  {
    level2.vars <- -1
    attrib <- attributes(terms(level2, data = data))
    
    if ( length( attrib$term.labels ) == 0 )
    {
      #may contain only the intercept
      if ( attrib$intercept != 1 )
      {
        stop( "bhlm: formula for second level effects must contain at least one variable or intercept." )
      }
      else
      {
        #formula contains the intercept
        level2.vars <- 0
      }
    }
    else
    {
      level2.vars <- 2
    }
  
    if ( level2.vars == 2 )
    {
      model.level2 <- call( "model.frame", formula = level2, data = data, na.action = na.action )

      model.level2 <- eval( model.level2, sys.parent() )
      Terms <- attr( model.level2, "terms" )
  
      Z <- model.matrix( Terms, model.level2, contrasts )
      #used.contrasts <- c( used.contrasts, contrasts( Z ) )

      second.effects <- TRUE

      level2.names <- dimnames( Z )[[2]]

    }
    else if ( level2.vars == 0 )
    {
      #only the intercept
      Z <- matrix( rep( 1, nrow( data ) ), nrow = nrow( data ) )

      second.effects <- TRUE
      level2.names <- "(Intercept)"
    }

  }#end if level2 effects

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

    #now sort missing data
    if ( sum.missing.data > 0 )
    {
      missing.data <- as.matrix( missing.data[ idx, ] )
    }


    if ( random.effects )
    {
      X <- as.matrix( X[ idx, ] )
      Y <- Y[ idx ]
    }
    
    if ( fixed.effects )
    {
      M <- as.matrix( M[ idx, ] )
      if ( !random.effects )
      {
        Y <- Y[ idx ]
      }
    }

    #counts by groups
    n.responses <- table( data[ , group.var ] )
    n.groups <- length( n.responses )

    group.names <- data[[ group.var ]][ idx ] 
    group.names <- group.names[ cumsum( n.responses ) ]

   
    if ( second.effects )
    {
      Z <- as.matrix( Z[ idx, ] )
      #one row per group
      Z <- as.matrix( Z[ cumsum( n.responses ), ] )

      #now arrange Z 


      Id <- diag( ncol( X ) )
      Z.X <- matrix( 0, ncol = ncol( Z ) * ncol( X ), nrow = ncol( X ) * n.groups )
      i1 <- 1
      for ( i in seq( 1, n.groups ) )
      {
        i2 <- i1 + ncol( X ) - 1

        Z.X[ seq( i1, i2 ), ] <- t( matrix( Z[ i, ] %o% Id, nrow = ncol( Z ) * ncol( X ) ) )

        i1 <- i1 + ncol( X )

      }
      Z <- Z.X

    }    

  }#end if group
  else
  {
    n.responses <- length( Y )
    n.groups <- 1
    group.names <- 1

    if ( second.effects )
    {
      #one row per group
      Z <- Z[ 1, ]
      Z <- t( matrix( Z[ 1, ] %o% diag( ncol( X ) ), nrow = ncol( Z ) * ncol( X ) ) )
    }
  }

  number.data <- length( Y )


####end Data-------------------------------------------------------

  ##handle missing data
  Y <- t( Y )


  missing.response <- 0
  missing.random <- 0
  missing.fixed <- 0
  if ( sum.missing.data > 0 )
  {
    missing.response <- t( missing.data[ , response.names ] )

    if ( random.effects )
    {
      if ( length( dimnames( X )[[2]] ) > 1 || dimnames( X )[[2]][1] != "(Intercept)" )
      {
        missing.random <- t( missing.data[ , dimnames( X )[[2]][ dimnames( X )[[2]] != "(Intercept)" ] ] )
        if ( is.element( "(Intercept)", dimnames( X )[[2]] ) )
        {
          missing.random <- rbind( rep( F, nrow( X ) ), missing.random )
          dimnames( missing.random ) <- dimnames( t( X ) )
        }
      }
    }


    if ( fixed.effects )
    {
      if ( length( dimnames( M )[[2]] ) > 1 || dimnames( M )[[2]][1] != "(Intercept)" )
      {
        missing.fixed <- t( missing.data[ , dimnames( M )[[2]][ dimnames( M )[[2]] != "(Intercept)" ] ] )
        if ( is.element( "(Intercept)", dimnames( M )[[2]] ) )
        {
          missing.fixed <- rbind( rep( F, nrow( M ) ), missing.fixed )
          dimnames( missing.fixed ) <- dimnames( t( M ) )
        }
      }
    }

  }



  if ( length( Y ) == sum( n.responses ) )
  {
    #a vector
    Y <- matrix( Y, nrow = 1 )
    dimnames( Y ) <- list( response.names, NULL )

    if ( sum( missing.response ) > 0 )
    {
      missing.response <- matrix( missing.response, nrow = 1 )
    }
  }

  if ( is.vector( X ) )
  {
    if ( sum( missing.random ) > 0 )
    {
      missing.random <- matrix( missing.random, nrow = 1 )
    }
  }

  if ( is.vector( M ) )
  {
    if ( sum( missing.fixed ) > 0 )
    {
      missing.fixed <- matrix( missing.fixed, nrow = 1 )
    }
  }


  dim.response <- nrow( Y )

  missing.response.components <- 0 
  missing.response.counts <- 0
  total.missing.response <- 0
  if ( sum( missing.response ) > 0 )
  {
    missing.response.counts <- apply( missing.response, 2, sum )
    total.missing.response <- sum( missing.response.counts )
    if ( total.missing.response > 0 )
    {
      missing.response.components <- apply( missing.response, 2, 
                                            FUN = function( x, p ) 
                                                    { seq( 1, p )[x == TRUE] },
                                            p = nrow( missing.response ) )
    
      missing.response.components <- unlist( missing.response.components )
    }
  }

  missing.random.components <- 0 
  missing.random.counts <- 0
  total.missing.random <- 0
  if ( sum( missing.random ) > 0 )
  {
    missing.random.counts <- apply( missing.random, 2, sum )
    total.missing.random <- sum( missing.random.counts )
    if ( total.missing.random > 0 )
    {
      missing.random.components <- apply( missing.random, 2, 
                                            FUN = function( x, p ) 
                                                    { seq( 1, p )[x == TRUE] },
                                            p = nrow( missing.random ) )
    
      missing.random.components <- unlist( missing.random.components )
    }
  }

  missing.fixed.components <- 0 
  missing.fixed.counts <- 0
  total.missing.fixed <- 0
  if ( sum( missing.fixed ) > 0 )
  {
    missing.fixed.counts <- apply( missing.fixed, 2, sum )
    total.missing.fixed <- sum( missing.fixed.counts )
    if ( total.missing.fixed > 0 )
    {
      missing.fixed.components <- apply( missing.fixed, 2, 
                                            FUN = function( x, p ) 
                                                    { seq( 1, p )[x == TRUE] },
                                            p = nrow( missing.fixed ) )
    
      missing.fixed.components <- unlist( missing.fixed.components )
    }
  }

  if ( ( length( M ) == 1 && M == 0 ) && ( ncol( X ) == 1 && random.names[1] == "(Intercept)" ) )
  {
    random.names <- "MEAN"
    dimnames( X )[[2]] <- random.names

    if ( second.effects )
    {
      if ( length( level2.names ) == 1 && level2.names[1] == "(Intercept)" )
      {
        level2.names <- "MEAN"
      }
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


  if ( !is.null( likelihood$errorCov ) && is.character( likelihood$errorCov ) )
  {
    likelihood$errorCov <- parse( text = likelihood$errorCov )[[1]]
  }


  check.errorCov <- TRUE
  if ( is.null( likelihood$errorCov ) ) 
  {
    likelihood$errorCov <- 1
    check.errorCov <- FALSE
  }
  else if ( is.name( likelihood$errorCov ) )
  {
    if (  likelihood$errorCov != "identity" )
    {
      errorCov <- as.vector( eval( likelihood$errorCov ) )
      if ( is.vector( errorCov ) && length( errorCov ) != n.responses[1] && length( errorCov == n.groups ) )
      {
        proceed.with.groups <- FALSE
        check.errorCov <- TRUE
      }
      else
      {
        proceed.with.groups <- TRUE
        check.errorCov <- FALSE
      }   
    }
    else
    {
      errorCov <- as.vector( eval( call( likelihood$errorCov, p = n.responses[1] ) ) )
      proceed.with.groups <- TRUE
      check.errorCov <- FALSE
    }

    if ( proceed.with.groups && n.groups > 1 )
    {
      for ( i in seq( 2, n.groups ) )
      {
        if (  likelihood$errorCov != "identity" )
        {
          this.errorCov <- as.vector( eval( likelihood$errorCov ) )
        }
        else
        {
          this.errorCov <- as.vector( eval( call( likelihood$errorCov, p = n.responses[i] ) ) )
        }
        errorCov <- c( errorCov, this.errorCov )
      }
    }
    likelihood$errorCov <- errorCov
  }#end name
  else if ( is.function( likelihood$errorCov ) )
  {
    errorCov <- as.vector( likelihood$errorCov( n.responses[1] ) )
    if ( n.groups > 1 )
    {
      for ( i in seq( 2, n.groups ) )
      {
        this.errorCov <- as.vector( likelihood$errorCov( n.responses[i] ) )
        errorCov <- c( errorCov, this.errorCov )
      }
    }
    likelihood$errorCov <- errorCov
    check.errorCov <- FALSE
  }#end function

  if ( check.errorCov && is.vector( likelihood$errorCov ) ) 
  {

    #cat( "is vector lklh:\n" )
    #print( likelihood$errorCov )

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
      else if ( length( likelihood$errorCov ) == number.data )
      {
        dim.error.Cov <- 1
      }
      else
      { 
	stop( "bhlm: The dimension of the error covariance matrix must equal the number of observations." )
      }
    }
  }
  else if ( check.errorCov && is.list( likelihood$errorCov ) ) 
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
      if ( any( dim( likelihood$errorCov[[ i ]]  ) != n.responses[i] ) )
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
  
  # check that the user has not specified priors for non-existent parameters
  if( !random.effects && prior$priorSpec$random.var )
    stop( "bhlm: a prior has been specified for the random effect variance, but there are no random effects in the model." )
  if( !fixed.effects && prior$priorSpec$fixed.coef )
    stop( "bhlm: a prior has been specified for fixed effects, but there are no fixed effects in the model." )
  if( !second.effects && prior$priorSpec$level2.coef )
    stop( "bhlm: a prior has been specified for level 2 effects, but there are no level 2 effects in the model." )
  
  valid.prior <- bhlm.prior( error.var = prior$error.var, fixed.coef = prior$fixed.coef, 
    level2.coef = prior$level2.coef, random.var = prior$random.var, 
    common.error.var = prior$common.error.var )

  level2.coef.mean <- 0
  level2.coef.Cov <- 1
  level2.coef.df <- 0
  level2.coef.type <- 1 #nonInformative

  fixed.coef.mean <- 0
  fixed.coef.Cov <- 1
  fixed.coef.df <- 0
  fixed.coef.type <- 1 #nonInformative


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
      if ( valid.prior$level2.coef$name == "nonInformative" )
      {
        level2.coef.mean <- rep( 0.0, ncol( Z ) )
        level2.coef.Cov <- diag( ncol( Z ) )
        level2.coef.type <- 1
      }
      else
      {
        level2.coef.mean <- valid.prior$level2.coef$parameters[["mean"]]
        level2.coef.Cov <- valid.prior$level2.coef$parameters[["S"]]

        level2.coef.mean <- valid.mean.specification( level2.coef.mean, ncol( Z ) )
        level2.coef.Cov <- valid.Cov.specification( level2.coef.Cov, ncol( Z ) )

        if ( valid.prior$level2.coef$name == "t" )
        {
          level2.coef.df <- valid.prior$level2.coef$parameters[["df"]]
        }
        level2.coef.type <- 0
      }

      ## beta for the full model
      random.coef.type <- 0

      if ( valid.prior$random.coef$name == "nonInformative" )
      {
        random.coef.mean <- rep( 0.0, ncol( X ) )
        random.coef.Cov <- diag( ncol( X ) )
      }
      else
      {
        random.coef.mean <- rep( 0.0, ncol( X ) )
        random.coef.Cov <- valid.prior$random.coef$parameters[["S"]]
        random.coef.Cov <- valid.Cov.specification( random.coef.Cov, ncol( X ) )

        if ( valid.prior$random.coef$name == "t" )
        {
          random.coef.df <- valid.prior$random.coef$parameters[["df"]]
        }
      }
    }#end second effects
    else
    {
      ## beta for the mixed effects model
      random.coef.type <- 1
      
      if ( class( valid.prior$random.coef ) == "fbprior" )
      {
        if ( valid.prior$random.coef$name == "nonInformative" )
        {
          random.coef.type <- 3
          random.coef.mean <- rep( 0, ncol( X ) )
          random.coef.Cov <- diag( ncol( X ) )
        }
        else
        {
          random.coef.mean <- valid.prior$random.coef$parameters[["mean"]]
          random.coef.Cov <- valid.prior$random.coef$parameters[["S"]]

          random.coef.mean <- valid.mean.specification( random.coef.mean, ncol( X ) )
          random.coef.Cov <- valid.Cov.specification( random.coef.Cov, ncol( X ) )

          if ( valid.prior$random.coef$name == "t" )
          {
            random.coef.df <- valid.prior$random.coef$parameters[["df"]]
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

        if ( valid.prior$random.coef[[1]]$name == "nonInformative" )
        {
          random.coef.Cov <- diag( ncol( X ) )
          random.coef.mean <- rep ( 0, ncol( X ) )
          random.coef.type <- 1
        }
        else
        {
          random.coef.Cov <- valid.prior$random.coef[[1]]$parameters[["S"]]
          random.coef.Cov <- valid.Cov.specification( random.coef.Cov, ncol( X ) )         

          random.coef.mean <- apply( matrix( seq( 1, length( valid.prior$random.coef ) ), nrow = 1 ), 
                                     MARGIN = 2, 
                                     FUN = function( x, y, dim )
                                     {
                                       if ( y[[ x ]]$name == "nonInformative" )
                                       {
                                         r.mean <- rep( 0, dim ) 
                                       }
                                       else
                                       {
                                         r.mean <- y[[ x ]]$parameters[["mean"]]
                                         r.mean <- valid.mean.specification( r.mean, dim )
                                       }
                                       r.mean
                                     },
                                     y = valid.prior$random.coef,
                                     dim = ncol( X ) )

        }
 
        ## all t's are assume to have the same df's
        if ( valid.prior$random.coef[[1]]$name == "t" )
        {
          random.coef.df <- valid.prior$random.coef[[1]]$parameters[["df"]]
        }

      }#end several prior parameters

    }#end no second effects


    ##check prior for random coef. variance

    ##check beta has an informative prior
    if ( random.coef.type != 3 )
    {
      random.var.type <- switch( valid.prior$random.var$name,
      "invChisq" =
      {
        random.var.scale <- valid.prior$random.var$parameters[["sigma0.sq"]]
        random.var.nu <- valid.prior$random.var$parameters[["df"]]
        
        if( length( random.var.scale ) != ncol(X) ){
          if( length( random.var.scale ) == 1 )
            random.var.scale <- rep( random.var.scale, ncol(X) )
          else
            stop("bhlm: the length of the vector sigma0.sq specified in the prior for the random effect variance is neither the length of the random effect vector nor 1.")
        }
        if( length( random.var.nu ) != ncol(X) ){
          if( length( random.var.nu ) == 1 )
            random.var.nu <- rep( random.var.nu, ncol(X) )
          else
            stop("bhlm: the length of the vector df specified in the prior for the random effect variance is neither the length of the random effect vector nor 1.")
        }
        0
      },

      "nonInfoPower" =
      {
        random.var.power <- valid.prior$random.var$parameters[["power"]]
        1
      },

      "uniformShrinkage" =
      {
        random.var.scale <- valid.prior$random.var$parameters[["median"]]
        2
      },

      "duMouchel" =
      {
        random.var.scale <- valid.prior$random.var$parameters[["dispersion"]]
        3
      },

      "invWishart" =
      {
        random.var.scale <- valid.Cov.specification( valid.prior$random.var$parameters[["Sigma"]], ncol( X ) )
        random.var.nu <- valid.prior$random.var$parameters[["df"]]
        4
      },

      stop( "bhlm: distribution for random coefficients variance must be either \"invChisq\", \"nonInformative power\", \"uniform shrinkage\", \"du Mouchel\", or \"invWishart\" \n" )
      )#end switch
    }
  }#end random effects prior



  ##fixed effects prior
  if ( fixed.effects )
  {
    if ( valid.prior$fixed.coef$name == "nonInformative" )
    {
      fixed.coef.mean <- rep( 0.0, ncol( M ) )
      fixed.coef.Cov <- diag( ncol( M ) )
      fixed.coef.type <- 1
    }
    else
    {
      fixed.coef.mean <- valid.prior$fixed.coef$parameters[["mean"]]
      fixed.coef.Cov <- valid.prior$fixed.coef$parameters[["S"]]

      fixed.coef.mean <- valid.mean.specification( fixed.coef.mean, ncol( M ) )
      fixed.coef.Cov <- valid.Cov.specification( fixed.coef.Cov, ncol( M ) )

      if ( valid.prior$fixed.coef$name == "t" )
      {
        fixed.coef.df <- valid.prior$fixed.coef$parameters[["df"]]
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

    error.var.type <- switch( valid.prior$error.var$name,
      "invChisq" =
      {
        error.var.scale <- valid.prior$error.var$parameters[["sigma0.sq"]]
        error.var.nu <- valid.prior$error.var$parameters[["df"]]
        0
      },

      "nonInfoPower" =
      {
        error.var.power <- valid.prior$error.var$parameters[["power"]]
        1
      },

      "uniformShrinkage" =
      {
        error.var.scale <- valid.prior$error.var$parameters[["median"]]
        2
      },

      "duMouchel" =
      {
        error.var.scale <- valid.prior$error.var$parameters[["dispersion"]]
        3
      },

      "mass point" =
      {
        error.var.scale <- valid.prior$error.var$parameters[["value"]]
        4
      },

      stop( "bhlm: distribution for error variance must be either \"invChisq\", \"nonInformative power\", \"uniform shrinkage\", \"du Mouchel\", or \"mass point\" \n" )

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
      error.var.scale <- rep( 0, length( valid.prior$error.var ) )
      error.var.nu <- rep( 0, length( valid.prior$error.var ) )
      error.var.power <- rep( 0, length( valid.prior$error.var ) )

      for ( i in seq( 1, length( valid.prior$error.var ) ) )
      {
        error.var.type <- switch( valid.prior$error.var[[i]]$name,
          "invChisq" =
          {
            error.var.scale[i] <- valid.prior$error.var[[i]]$parameters[["sigma0.sq"]]
            error.var.nu[i] <- valid.prior$error.var[[i]]$parameters[["df"]]
            0
          },

          "nonInfoPower" =
          {
            error.var.power[i] <- valid.prior$error.var[[i]]$parameters[["power"]]
            1
          },

          "uniformShrinkage" =
          {
            error.var.scale[i] <- valid.prior$error.var[[i]]$parameters[["median"]]
            2
          },

          "duMouchel" =
          {
            error.var.scale[i] <- valid.prior$error.var[[i]]$parameters[["dispersion"]]
            3
          },

          "mass point" =
          {
            error.var.scale[i] <- valid.prior$error.var[[i]]$parameters[["value"]]
            4
          },

          stop( "bhlm: distribution for error variance must be either \"invChisq\", \"nonInformative power\", \"uniform shrinkage\", \"du Mouchel\", or \"mass point\" \n" )
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

#set.seed( random.seed )

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
    dim.random <- ncol( X )
  }
  else
  {
    dim.random <- 0
  }

  if ( fixed.effects )
  {
    dim.fixed <- ncol( M )
  }
  else
  {
    dim.fixed <- 0
  }

  if ( second.effects )
  {
    dim.level2 <- ncol( Z )
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
    starting.points <- generateInitPoints.bhlm( number.chains, n.responses, Y, X, M, Z, 
                                                dim.random, dim.fixed, dim.level2,
                                                error.var.nu, error.var.scale, error.var.type, error.var.common,
                                                random.coef.df, random.coef.mean, random.coef.Cov, random.coef.type,
                                                fixed.coef.df, fixed.coef.mean, fixed.coef.Cov,
                                                random.var.nu, random.var.scale, random.var.type,
                                                level2.coef.df, level2.coef.mean, level2.coef.Cov,
                                                mix.with.MLE = 0 )
    if( debug ){
    	print("Initial Values:")
    	print(starting.points)
    }
  }#end prior init points
  else if ( !is.null( init.point ) && init.point$type == "user's choice" )
  {
    read.init.point <- 1

    starting.points <- validate.initial.points.bhlm( n.groups, ncol( X ), ncol( M ), ncol( Z ), init.point$values[[ 1 ]], random.effects, fixed.effects, second.effects, error.var.common, error.var.type, random.coef.type, random.var.type )

    if ( !is.null( starting.points$random.var ) )    
    {
      s.random.var <- vector( "list", number.chains )
      s.random.var[[1]] <- starting.points$random.var
    }

    if ( number.chains > 1 )
    {
      for ( i in seq( 2, number.chains ) )
      {
        s.points <- validate.initial.points.bhlm( n.groups, ncol( X ), ncol( M ), ncol( Z ), init.point$values[[ i ]], random.effects, fixed.effects, second.effects, error.var.common, error.var.type, random.coef.type, random.var.type )
        if ( !is.null( starting.points$error.var ) )
        {
          starting.points$error.var <- cbind( starting.points$error.var, s.points$error.var )
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

    }#end if chains > 1
    else
    {
      if (  !is.null( starting.points$random.var ) )
      {
        starting.points$random.var <- s.random.var
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
      starting.points$error.var <-  matrix( 0, ncol = number.chains ) 
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


  bhlmodel <- vector( "list", number.chains )

  for ( i in seq( 1, number.chains ) )
  {
    simulation.seed <- .Random.seed
    fit <- fit.bayeshlm( n.groups, idx, n.responses, dim.random, dim.fixed, 
      dim.level2, response.names, random.names, fixed.names, level2.names, 
      group.names, X, M, Z, Y, total.missing.response, missing.response.counts, 
      missing.response.components, total.missing.random, missing.random.counts, 
      missing.random.components, total.missing.fixed, missing.fixed.counts, 
      missing.fixed.components, unique.lklhd.Cov, dim.error.Cov, 
      likelihood$errorCov, likelihood$df, likelihood$type, random.coef.mean, 
      random.coef.Cov, random.coef.df, random.coef.type, fixed.coef.mean, 
      fixed.coef.Cov, fixed.coef.df, fixed.coef.type, level2.coef.mean, 
      level2.coef.Cov, level2.coef.df, level2.coef.type, error.var.nu, 
      error.var.scale, error.var.power, error.var.common, error.var.type, 
      random.var.nu, random.var.scale, random.var.power, random.var.type, 
      read.init.point, starting.points$random.coef[,i], 
      starting.points$fixed.coef[,i], starting.points$level2.coef[,i], 
      starting.points$error.var[,i], starting.points$random.var[[i]], 
      sampler.type, burnInLength, simulationsToPerform, sampleFrequency, 
      print.stats, dimnames( Y )[[2]], debug = debug )
    bhlmodel[[i]] <- fit

  }#end for chains

  #return(posterior(mcmc.list(bhlmodel)))
  post <- mcmc.list(bhlmodel)
  oldClass(post) <- c("posterior", class(post))
  post

}#end






############################################################################################
##
## FIT BAYES HIERARCHICAL LINEAR MODEL
##
############################################################################################

fit.bayeshlm <- function( n.groups, sorted.by.group.idx, n.responses, dim.random, dim.fixed, dim.level2, response.names,
                          random.names, fixed.names, level2.names, group.names, X, M, Z, Y, 
                          total.missing.response, missing.response.counts, missing.response.components,
                          total.missing.random, missing.random.counts, missing.random.components,
                          total.missing.fixed, missing.fixed.counts, missing.fixed.components,
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
                          print.stats, observation.names = NULL, debug = F )
{

  number.data <- sum( n.responses )

  #now count how many variables there are.
  # output.dim is the total number of parameters
  # and n.vars appears to be the number of parameters 
  # for which sampling information (acceptance rates?) is 
  # generated.
  n.vars <- 0
  output.dim <- 0


  if ( error.var.type != 4 )
  {
    if ( error.var.common == 0 || error.var.common == 1 )
    {
      #different error variance for each group
      n.vars <- n.groups
    }
    else
    {
      #same error variance
      n.vars <- 1
    }
    output.dim <- n.vars
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
        # changed by DBW
        #dim.random.var <- 1
        dim.random.var <- dim.random
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

  #account for missing random predictors
  if ( total.missing.random > 0 )
  {
    n.vars <- n.vars + sum( missing.random.counts > 0 ) 
    output.dim <- output.dim + total.missing.random
  }

  #account for missing fixed predictors
  if ( total.missing.fixed > 0 )
  {
    n.vars <- n.vars + sum( missing.fixed.counts > 0 ) 
    output.dim <- output.dim + total.missing.fixed
  }
  
  if( debug ){
    print( paste( "number of parameters: ", output.dim ) )
  }

  #now create output array with the right dimensions
  output.samples <- matrix( 0, nrow = simulationsToPerform, ncol = output.dim )

  if ( print.stats == 1 )
  {
    if ( sampler.type == 0 )
    {
      #the Gibbs sampler case
      gibbs.stats <- rep( 0, n.vars )
    }
    else
    {
      gibbs.stats <- 0
    }
  }
  else
  {
    gibbs.stats <- 0
  }

  if(is.na(starting.random.var[1])) starting.random.var <- NULL

  #call the function
  if (debug){
    cat("n.groups:", n.groups, "\n")
    cat("n.responses", n.responses, "\n")
    cat("dim.random:", dim.random, "\n")
    cat("dim.level2:", dim.level2, "\n")
    cat("dim.fixed:", dim.fixed, "\n")
    cat("X", X, "\n")
    cat("M", M, "\n")
    cat("Z", Z, "\n")
    cat("Y", Y, "\n")
    cat("fixed.coef.mean", fixed.coef.mean, "\n")
    cat("fixed.coef.Cov", fixed.coef.Cov, "\n")
    cat("fixed.coef.type", fixed.coef.type, "\n")
    cat("fixed.coef.df", fixed.coef.df, "\n")
    cat("random.coef.mean", random.coef.mean, "\n")
    cat("random.coef.Cov", random.coef.Cov, "\n")
    cat("random.coef.type", random.coef.type, "\n")
    cat("random.var.nu", random.var.nu, "\n")
    cat("random.var.scale", random.var.scale, "\n")
    cat("random.var.power", random.var.power, "\n")
    cat("random.var.type", random.var.type, "\n")
    cat("error.var.nu", error.var.nu, "\n")
    cat("error.var.scale", error.var.scale, "\n")
    cat("error.var.power", error.var.power, "\n")
    cat("error.var.common", error.var.common, "\n")
    cat("output.dim:", output.dim, "\n")
    cat("read.init.point:", read.init.point, "\n")
    cat("starting.random.coef:", starting.random.coef, "\n")
    cat("starting.fixed.coef:", starting.fixed.coef, "\n")
    cat("starting.level2.coef:", starting.level2.coef, "\n")
    cat("starting.random.var:", starting.random.var, "\n")
    cat("starting.error.var:", starting.error.var, "\n")
    cat("burnInLength:", burnInLength, "\n")
    cat("simulationsToPerform:", simulationsToPerform, "\n")
    cat("sampleFrequency:", sampleFrequency, "\n")
  }

  fit <- .C( "fitBayesianHLM",
             as.integer( n.groups ),
             as.integer( n.responses ),
             as.integer( dim.random ), 
             as.integer( dim.fixed ), 
             as.integer( dim.level2 ), 
             #
             as.double( t(X) ), 
             as.double( t(M) ), 
             as.double( t(Z) ), 
             as.double( Y ), 
             #
             as.integer( total.missing.response ),
             as.integer( missing.response.counts ),
             as.integer( missing.response.components ),
             #
             as.integer( total.missing.random ),
             as.integer( missing.random.counts ),
             as.integer( missing.random.components ),
             #
             as.integer( total.missing.fixed ),
             as.integer( missing.fixed.counts ),
             as.integer( missing.fixed.components ),
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
             as.integer(burnInLength ), 
             as.integer(simulationsToPerform ), 
             as.integer(sampleFrequency ), 
             #
             as.integer( print.stats ),
             output.samples = as.double( output.samples ),
             gibbs.stats = as.double( gibbs.stats ),
             #
             PACKAGE = "FlexBayes")


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
  X.imputed <- t( X )
  M.imputed <- t( M )

  original.random.names <- random.names
  group.names <- as.character( group.names )

  end.index <- 0
  if ( dim.random > 0 )
  {
    random.coef <- matrix( out.samples[ , seq( 1, n.groups * dim.random ) ], nrow = simulationsToPerform )
    random.names <- as.vector( outer( outer( random.names, response.names, paste, sep = ":" ), group.names, paste, sep = ":" ) )   
    end.index <- n.groups * dim.random
  }

  if ( dim.fixed > 0 )
  {
    fixed.coef <- matrix( out.samples[ , end.index + seq( 1, dim.fixed ) ], nrow = simulationsToPerform )
    end.index <- end.index + dim.fixed
  }

  if ( dim.random > 0 && dim.level2 > 0 )
  {
    level2.coef <- matrix( out.samples[ , end.index + seq( 1, dim.level2 ) ], nrow = simulationsToPerform )
    level2.names <- as.vector( outer( level2.names, original.random.names, paste, sep = ":" ) )
    end.index <- end.index + dim.level2
  }

  if ( error.var.type != 4 )
  {
    if ( error.var.common == 0 || error.var.common == 1 )
    {
      #different error variance for each group
      error.var <- matrix( out.samples[ , end.index + seq( 1, n.groups ) ], nrow = simulationsToPerform )
      error.var.names <- as.vector( outer( "SIGMA:", group.names, paste, sep = "" ) )
      end.index <- end.index + n.groups
    }
    else
    {
      #same error variance prior
      error.var <- matrix( out.samples[ , end.index + 1 ], nrow = simulationsToPerform )
      error.var.names <- "SIGMA"
      end.index <- end.index + 1
    }
  }

  #if beta has an informative prior
  if ( dim.random > 0 && random.coef.type != 3 )
  {
    if ( random.var.type == 4 )
    {
    	# inverse Wishart prior case
      random.var <- matrix( out.samples[ , end.index + seq( 1, dim.random * dim.random ) ], nrow = simulationsToPerform )
      random.var.names <- as.vector( t( outer( paste( "RANDOM:TAU:", seq( 1, dim.random), sep = "" ), seq( 1, dim.random ), paste, sep = "." ) ) )
      end.index <- end.index + dim.random * dim.random
    }
    else
    {
      # case of random variances a priori iid
      # code changed by DBW.  Original preserved in comments
      #random.var <- matrix( out.samples[ , end.index + 1 ], nrow = simulationsToPerform )
      random.var <- matrix( out.samples[ , end.index + seq( 1, dim.random ) ], nrow = simulationsToPerform )
      #random.var.names <- "RANDOM:TAU"
      random.var.names <- paste("RANDOM:TAU:", original.random.names, sep = "")
      #end.index <- end.index + 1
      end.index <- end.index + dim.random
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
  imputed.values <- NULL
  imputed.names <- NULL
  if ( total.missing.response > 0 )
  {
    imputed.values <- matrix( out.samples[ , end.index + seq( 1, total.missing.response ) ], nrow = simulationsToPerform )
    end.index <- end.index + total.missing.response

    obs.missing <- seq( 1, length( missing.response.counts ) )[ missing.response.counts > 0 ]
    k <- 1
    for ( i in seq( 1, length( obs.missing ) ) )
    {
      for ( j in seq( 1, missing.response.counts[ obs.missing[i] ] ) )
      {
        if ( !is.null( observation.names ) )
        {
          imputed.names <- c( imputed.names, paste( "IMPUTED", observation.names[ obs.missing[i] ], response.names[ missing.response.components[ k ] ], sep = ":" ) )
        }
        else
        {
          imputed.names <- c( imputed.names, paste( "IMPUTED", obs.missing[i], response.names[ missing.response.components[ k ] ], sep = ":" ) )
        }

        #complete data
        Y.imputed[ response.names[ missing.response.components[ k ] ], obs.missing[i] ] <- mean( imputed.values[ , k ] )

        k <- k + 1
      }
    }
  }


  if ( total.missing.random > 0 )
  {
    if ( is.null( imputed.values ) )
    {
      imputed.values <- matrix( out.samples[ , end.index + seq( 1, total.missing.random ) ], nrow = simulationsToPerform )
      k <- 1
    }
    else
    {
      imputed.values <- cbind( imputed.values, matrix( out.samples[ , end.index + seq( 1, total.missing.random ) ], nrow = simulationsToPerform ) )
    }
    
    end.index <- end.index + total.missing.random
    obs.missing <- seq( 1, length( missing.random.counts ) )[ missing.random.counts > 0 ]
    k.rd <- 1
    for ( i in seq( 1, length( obs.missing ) ) )
    {
      for ( j in seq( 1, missing.random.counts[ obs.missing[i] ] ) )
      {
        if ( !is.null( observation.names ) )
        {
          imputed.names <- c( imputed.names, paste( "IMPUTED", observation.names[ obs.missing[i] ], original.random.names[ missing.random.components[ k.rd ] ], sep = ":" ) )
        }
        else
        {
          imputed.names <- c( imputed.names, paste( "IMPUTED", obs.missing[i], original.random.names[ missing.random.components[ k.rd ] ], sep = ":" ) )
        }

        #complete data
        X.imputed[ original.random.names[ missing.random.components[ k.rd ] ], obs.missing[i] ] <- mean( imputed.values[ , k ] )

        k.rd <- k.rd + 1
        k <- k + 1
      }
    }
  }


  if ( total.missing.fixed > 0 )
  {
    if ( is.null( imputed.values ) )
    {
      imputed.values <- matrix( out.samples[ , end.index + seq( 1, total.missing.fixed ) ], nrow = simulationsToPerform )
      k <- 1
    }
    else
    {
      imputed.values <- cbind( imputed.values, matrix( out.samples[ , end.index + seq( 1, total.missing.fixed ) ], nrow = simulationsToPerform ) )
    }
    
    end.index <- end.index + total.missing.fixed
    obs.missing <- seq( 1, length( missing.fixed.counts ) )[ missing.fixed.counts > 0 ]
    k.fx <- 1
    for ( i in seq( 1, length( obs.missing ) ) )
    {
      for ( j in seq( 1, missing.fixed.counts[ obs.missing[i] ] ) )
      {
        if ( !is.null( observation.names ) )
        {
          imputed.names <- c( imputed.names, paste( "IMPUTED", observation.names[ obs.missing[i] ], fixed.names[ missing.fixed.components[ k.fx ] ], sep = ":" ) )
        }
        else
        {
          imputed.names <- c( imputed.names, paste( "IMPUTED", obs.missing[i], fixed.names[ missing.fixed.components[ k.fx ] ], sep = ":" ) )
        }

        #complete data
        M.imputed[ fixed.names[ missing.fixed.components[ k.fx ] ], obs.missing[i] ] <- mean( imputed.values[ , k ] )

        k.fx <- k.fx + 1
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

    random.effects <- list( dim = dim.random, names = random.names, orig.names = original.random.names, coef = random.coef, scale = random.var, scale.names = random.var.names, tau = random.tau, level2 = level2.effects )
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

  if ( error.var.type != 4 )
  {
    scale <- list( scale = error.var, names = error.var.names, common = error.var.common, type = error.var.type )
  }
  else
  {
    scale <- NULL
  }

  if ( total.missing.response > 0 || total.missing.random > 0 || total.missing.fixed > 0  )
  {
    data.imputed <- list( Y = t( Y.imputed ), X = t( X.imputed ), M = t( M.imputed ) )
    imputed.components <- list( values = imputed.values, names = imputed.names, data = data.imputed )
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

  # The "scale" parameter
  if ( error.var.type != 4 )
  {
  	dimnames(error.var) <- list(NULL, error.var.names)
  	post.samples <- cbind(post.samples, error.var)
  }

	# the random effects
	if ( dim.random > 0 ){
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
		if (!is.null(random.var)){
		  dimnames(random.var) <- list(NULL, random.var.names)
		  post.samples <- cbind(post.samples, random.var)
		}
	  dimnames(random.coef) <- list(NULL, random.names)
		post.samples <- cbind(post.samples, random.coef)
		if ( random.coef.df > 0 )
    {
      dimnames(random.coef.tau) <- list(NULL, paste( "TAU:RANDOM:", seq( 1, n.groups ), sep = "" ))
      post.samples <- cbind(post.samples, random.coef.tau)
    }
	}

  # Imputed values
  if ( total.missing.response > 0 || total.missing.random > 0 || total.missing.fixed > 0  )
  {
    warning("Imputed values are not currently returned.")
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

  return( mcmc( data = post.samples, thin = sampleFrequency) )
  #, burnin = burnInLength ) )
}#end



