##################################################
# bhbm function in the FlexBayes package          #
##################################################
## Kjell: 02/13/2014 moved zero and identity functions to bayes.distribution.R


#####################################################################################
##
##  BAYESIAN HIERARCHICAL BINOMIAL MODEL
##
#####################################################################################

#trials.formula, random.formula, fixed.formula, level2.formula and group.formula are of type "formula"
bhbm <- function( fixed.formula = NULL, 
                  trials.formula = NULL,
                  random.formula = NULL, 
                  level2.formula = NULL, 
                  group.formula = NULL, 
                  data,
                  overdispersion = "none",
                  prior = NULL, 
      	          sampler = bhpm.sampler(), 
                  random.seed = .Random.seed, 
                  na.action = NULL, contrasts = NULL,
                  print.stats = F, debug = F )
{

  binomial.model <- overdispersion
  binomial.model <- match.arg( binomial.model, 
    choices = c("beta-conj", "logit-normal", "none") )

  if( is.null( prior ) ){
    if( binomial.model == "beta-conj" )
      prior <- bhpm.prior( xi = bayes.uniformShrinkage(0.5), common.glm = 2 )
    else if( binomial.model == "logit-normal" )
      prior <- bhpm.prior( sigma2 = bayes.uniformShrinkage(0.5), common.glm = 2 )
    else   # no overdispersion case
      prior <- bhpm.prior()
  }

  if( ( binomial.model == "beta-conj" ) && !prior$conjModel ){
    stop("For the beta-conjugate model, a prior must be specified for xi")
  } else if( ( binomial.model == "logit-normal" ) && prior$conjModel ){
    stop("A prior was specified for the xi parameters, which do not exist in the logit-normal model")
  } else if( ( binomial.model == "none" ) && ( prior$common.glm != 3 ) ){
    stop("For models with no overdispersion parameter, the prior cannot be specified with common.glm = 0, 1, or 2")
  }

  ####-------------------------------------------------------
  ## Get Data 
  ##

  if ( is.null( data ) )
  {
    stop( "bhbm: data must be provided." )
  }

  data <- as.data.frame( data )

  if ( !is.null( na.action ) )
  {
    #get rid of rows with missing data if na.omit, otherwise stop and fail (na.fail)

    data <- na.action( data )
    if ( nrow( data ) == 0 )
    {
      stop( "bhbm: every observation (row) of the data set contains at least a missing value.\n" )
    }
  }
  else
  {
    na.action <- na.fail
  }


  #name of variables must be preserved for output
  response.names <- NULL
  trials.names <- NULL
  random.names <- NULL
  fixed.names <- NULL
  level2.names <- NULL

  used.contrasts <- NULL

  X <- 0
  M <- 0
  Z <- 0
  NTrials <- 0

  random.effects <- F
  fixed.effects <- F
  random.vars <- -1
  fixed.vars <- -1

  if ( !is.null( random.formula ) )
  {
    if ( length( terms( random.formula )@term.labels ) == 0 && terms( random.formula )@response == 0 )
    {
      #may contain only the intercept
      if ( terms( random.formula )@intercept != 1 )
      {
        stop( "bhbm: formula for random effects must contain at least one variable or intercept." )
      }
      else
      {
        #formula contains the intercept
        random.vars <- 0
      }
    }
    else if ( length( terms( random.formula )@term.labels ) > 0 && terms( random.formula )@response == 1 )
    {
      #formula contains the response variable, and predictors other than intercept
      random.vars <- 1
    }
    else if ( length( terms( random.formula )@term.labels ) == 0 )
    {
      #formula contains response and perhaps intercept
      if ( terms( random.formula )@intercept != 1 )
      {
        stop( "bhbm: formula for random effects must contain at least one variable or intercept." )
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

  if ( !is.null( fixed.formula ) )
  {
    if ( length( terms( fixed.formula )@term.labels ) == 0 && terms( fixed.formula )@response == 0 )
    {
      #may contain only the intercept
      if ( terms( fixed.formula )@intercept != 1 )
      {
        stop( "bhbm: formula for fixed effects must contain at least one variable or intercept." )
      }
      else
      {
        #formula contains the intercept
        fixed.vars <- 0
      }
    }
    else if ( length( terms( fixed.formula )@term.labels ) > 0 && terms( fixed.formula )@response == 1 )
    {
      #formula contains the response variable (and predictors other than just intercept)
      fixed.vars <- 1
    }
    else if ( length( terms( fixed.formula )@term.labels ) == 0 )
    {
      #formula contains response and perhaps intercept
      if ( terms( fixed.formula )@intercept != 1 )
      {
        stop( "bhbm: formula for fixed effects must contain at least one variable or intercept." )
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
    
  if ( !is.element( random.vars, c( 1, 2 ) ) && !is.element( fixed.vars, c( 1, 2 ) ) )
  {
    stop( "bhbm: model specification is not valid. You need to provide a response variable." )  
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
        used.contrasts <- contrasts( X )
        random.names <- dimnames( X )[[2]]

      }
      else
      {
        #only intercept as predictor
        X <- matrix( rep( 1, nrow( data ) ), nrow = nrow( data ) )
        random.names <- "(Intercept)"
      }

      if ( is.element( random.vars, c( 1, 2 ) ) )
      {
        Y <- model.extract( model.random, response )
        response.names <- dimnames( attr( Terms, "factors" ) )[[1]][1]
      }

      random.effects <- T
    }
    else #random effects contain only the intercept
    {
      X <- matrix( rep( 1, nrow( data ) ), nrow = nrow( data ) )
      random.effects <- T
      random.names <- "(Intercept)"
    }

    dimnames( X ) <- list( dimnames( data )[[1]], random.names )

  }#end if random formula


  if ( fixed.vars >= 0 )
  {
    if ( fixed.vars > 0 )
    {
      model.fixed <- call( "model.frame", formula = fixed.formula, data = data, na.action = na.action )
      model.fixed <- eval( model.fixed, sys.parent() )
      Terms <- attr( model.fixed, "terms" )

      if ( fixed.vars != 2 )
      {
        M <- model.matrix( Terms, model.fixed, contrasts )
        used.contrasts <- c( used.contrasts, contrasts( M ) )
        fixed.names <- dimnames( M )[[2]]

      }
      else
      {
        M <- matrix( rep( 1, nrow( data ) ), nrow = nrow( data ) )
        fixed.names <- "(Intercept)"
      }

      if ( is.element( fixed.vars, c( 1, 2 ) ) && !is.element( random.vars, c( 1, 2 ) ) )
      {
        Y <- model.extract( model.fixed, response ) 
        response.names <- dimnames( attr( Terms, "factors" ) )[[1]][1]
      }

      fixed.effects <- T
    }    
    else #fixed effects contain only the intercept
    {
      M <- matrix( rep( 1, nrow( data ) ), nrow = nrow( data ) )
      fixed.effects <- T
      fixed.names <- "(Intercept)"
    }

    dimnames( M ) <- list( dimnames( data )[[1]], fixed.names )

  }#end if fixed formula


  second.effects <- F
  if ( !is.null( level2.formula ) )
  {
    level2.vars <- -1
    if ( length( terms( level2.formula )@term.labels ) == 0 )
    {
      #may contain only the intercept
      if ( terms( level2.formula )@intercept != 1 )
      {
        stop( "bhbm: formula for second level effects must contain at least one variable or intercept." )
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
      model.level2 <- call( "model.frame", formula = level2.formula, data = data, na.action = na.action )
      model.level2 <- eval( model.level2, sys.parent() )
      Terms <- attr( model.level2, "terms" )
  
      Z <- model.matrix( Terms, model.level2, contrasts )
      used.contrasts <- c( used.contrasts, contrasts( Z ) )
      level2.names <- dimnames( Z )[[2]]

      second.effects <- T

    }
    else if ( level2.vars == 0 )
    {
      #only the intercept
      Z <- matrix( rep( 1, nrow( data ) ), nrow = nrow( data ) )

      second.effects <- T
      level2.names <- "(Intercept)"
    }

  }#end if level2 effects


  #now get trials
  if ( !is.null( trials.formula ) )
  {
    model.NTrials <- call( "model.frame", formula = trials.formula, data = data, na.action = na.action )
    model.NTrials <- eval( model.NTrials, sys.parent() )
    Terms <- attr( model.NTrials, "terms" )
    Terms@intercept <- 0

    NTrials <- model.matrix( Terms, model.NTrials )
    trial.names <- dimnames( NTrials )[[2]]
  }
  else
  {
    stop( "bhbm: The number of trials associated to each count must be provided.\n" )
  }



  #sort data by group
  if ( is.null( Y ) )
  {
    stop( "bhbm: response variable must be provided.\n" )
  }

  if ( !is.null( group.formula ) )
  {
    group.var <- attr( terms( group.formula ), "term.labels" )
    idx <- order( data[[ group.var ]] )

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

    NTrials <- as.matrix( NTrials[ idx, ] )

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
      # Bug fix by Dawn: the old code had:
      # Z <- Z[ 1, ]
      Z <- t( matrix( Z[ 1, ] %o% diag( ncol( X ) ), nrow = ncol( Z ) * ncol( X ) ) )
    }
  }

  number.data <- length( Y )

  Y <- t( Y )
  if ( length( Y ) == sum( n.responses ) )
  {
    #a vector
    Y <- matrix( Y, nrow = 1 )
    dimnames( Y ) <- list( response.names, NULL )
  }

  dim.response <- nrow( Y )

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

####end Data-------------------------------------------------------




####----------------------------------------------------------------------
  ## Validate/interpret the arguments specifying the prior
  ##
  
  # check that the user has not specified priors for non-existent parameters
  if( !random.effects && prior$priorSpec$random.var )
    stop( "bhbm: a prior has been specified for the random effect variance, but there are no random effects in the model." )
  if( !fixed.effects && prior$priorSpec$fixed.coef )
    stop( "bhbm: a prior has been specified for fixed effects, but there are no fixed effects in the model." )
  if( !second.effects && prior$priorSpec$level2.coef )
    stop( "bhbm: a prior has been specified for level 2 effects, but there are no level 2 effects in the model." )

  if( binomial.model == "beta-conj" )
    valid.prior <- bhpm.prior( xi = prior$glm, fixed.coef = prior$fixed.coef, level2.coef = prior$level2.coef, 
      random.var = prior$random.var, common.glm = prior$common.glm )
  else if( binomial.model == "logit-normal" )
    valid.prior <- bhpm.prior( sigma2 = prior$glm, fixed.coef = prior$fixed.coef, level2.coef = prior$level2.coef, 
      random.var = prior$random.var, common.glm = prior$common.glm )
  else
    valid.prior <- bhpm.prior( fixed.coef = prior$fixed.coef, level2.coef = prior$level2.coef, 
      random.var = prior$random.var )
  
  glm.pars <- 1
  glm.type <- 0
  level2.coef.mean <- 0
  level2.coef.Cov <- 1
  level2.coef.df <- 0
  level2.coef.type <- 1 #non-informative

  fixed.coef.mean <- 0
  fixed.coef.Cov <- 1
  fixed.coef.df <- 0
  fixed.coef.type <- 1 #non-informative

  random.coef.mean <- 0
  random.coef.Cov <- 1
  random.coef.df <- 0
  random.coef.type <- 0

  random.var.scale <- 1
  random.var.nu <- 0
  random.var.power <- -0.5 #Haar measure
  random.var.type <- 0  #invChisq

  if ( binomial.model == "beta-conj" )
  {
    lambda.out <- F

    # the only prior implemented for xi in the beta-conj model is the 
    # uniform shrinkage distribution
    if ( valid.prior$glm@name == "uniform shrinkage" )
    {
      glm.pars <- valid.prior$glm@parameters[["median"]]
      glm.type <- 1
      glm.names <- "XI"
      common.glm <- valid.prior$common.glm
    }
    else
    {
      stop( "bhbm: Invalid prior specification for xi parameter.\n" )
    }
  }
  else if( binomial.model == "logit-normal" )
  {
    #assume logit-normal model

    lambda.out <- T

    ##sigma prior
    error.var.scale <- 1
    error.var.nu <- 0
    error.var.power <- -0.5 #haar measure
    error.var.common <- valid.prior$common.glm
    error.var.type <- 0

    if ( error.var.common == 2 )
    {
      #a common sigma
      glm.names <- "sigma"
    }
    else
    {
      glm.names <- paste( "sigma", seq( 1, n.groups ), sep = ":" )
    }

    if ( is.element( error.var.common, c( 1, 2 ) ) )
    {
      # common prior parameters for error variance

      error.var.type <- switch( valid.prior$glm@name,
      "invChisq" =
      {
        error.var.scale <- valid.prior$glm@parameters[["sigma0.sq"]]
        error.var.nu <- valid.prior$glm@parameters[["df"]]
        0
      },

      "non-informative power" =
      {
        error.var.power <- valid.prior$glm@parameters[["power"]]
        1
      },

      "uniform shrinkage" =
      {
        error.var.scale <- valid.prior$glm@parameters[["median"]]
        2
      },

      "du Mouchel" =
      {
        error.var.scale <- valid.prior$glm@parameters[["dispersion"]]
        3
      },

      "mass point" =
      {
        error.var.scale <- valid.prior$glm@parameters[["value"]]
        4
      },

      stop( "bhbm: distribution for error variance must be either \"invChisq\", \"non-informative power\", \"uniform shrinkage\", \"du Mouchel\", or \"mass point\" \n" )

      )#end switch
    }
    else
    {
      # different prior parameters for error variance

      if ( length( valid.prior$glm ) != n.groups )
      {
        stop( "bhbm: list of prior parameters for error variance is not of the appropriate size." )
      }

      if ( is.list( valid.prior$glm ) )
      { 
        error.var.scale <- rep( 0, length( valid.prior$glm ) )
        error.var.nu <- rep( 0, length( valid.prior$glm ) )
        error.var.power <- rep( 0, length( valid.prior$glm ) )

        for ( i in seq( 1, length( valid.prior$glm ) ) )
        {
          error.var.type <- switch( valid.prior$glm[[i]]@name,
          "invChisq" =
          {
            error.var.scale[i] <- valid.prior$glm[[i]]@parameters[["sigma0.sq"]]
            error.var.nu[i] <- valid.prior$glm[[i]]@parameters[["df"]]
            0
          },

          "non-informative power" =
          {
            error.var.power[i] <- valid.prior$glm[[i]]@parameters[["power"]]
            1
          },

          "uniform shrinkage" =
          {
            error.var.scale[i] <- valid.prior$glm[[i]]@parameters[["median"]]
            2
          },

          "du Mouchel" =
          {
            error.var.scale[i] <- valid.prior$glm[[i]]@parameters[["dispersion"]]
            3
          },

          "mass point" =
          {
            error.var.scale[i] <- valid.prior$glm[[i]]@parameters[["value"]]
            4
          },

          stop( "bhbm: distribution for error variance must be either \"invChisq\", \"non-informative power\", \"uniform shrinkage\", \"du Mouchel\", or \"mass point\" \n" )
          )#end switch

          if ( i == 1 )
          {
            prev.type <- error.var.type
          }
          else if ( prev.type != error.var.type )
          {
            stop( "bhbm: all priors for error variances should be of the same type." )
          }
        
        }#end for loop
      }#end list
      else
      {
        stop( "bhbm: prior for error variance is not valid." )
      }
    }#end different prior parameters  
  }#end sigma
  else {   # no overdispersion parameter
    glm.names <- NULL
    common.glm <- 3
    glm.type <- 1
  }
  
  if ( random.effects )
  {
    if ( second.effects )
    {
      if ( valid.prior$level2.coef@name == "non-informative" )
      {
        level2.coef.mean <- rep( 0.0, ncol( Z ) )
        level2.coef.Cov <- diag( ncol( Z ) )
        level2.coef.type <- 1
      }
      else
      {
        level2.coef.mean <- valid.prior$level2.coef@parameters[["mean vector"]]
        level2.coef.Cov <- valid.prior$level2.coef@parameters[["covariance matrix"]]

        level2.coef.mean <- valid.mean.specification.bhpm( level2.coef.mean, ncol( Z ) )
        level2.coef.Cov <- valid.Cov.specification.bhpm( level2.coef.Cov, ncol( Z ) )

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
        random.coef.mean <- rep( 0.0, ncol( X ) )
        random.coef.Cov <- diag( ncol( X ) )
      }
      else
      {
        #ignore prior mean in this case

        random.coef.mean <- rep( 0.0, ncol( X ) )
        random.coef.Cov <- valid.prior$random.coef@parameters[["covariance matrix"]]
        random.coef.Cov <- valid.Cov.specification.bhpm( random.coef.Cov, ncol( X ) )

        if ( valid.prior$random.coef@name == "t" )
        {
          random.coef.df <- valid.prior$random.coef@parameters[["t degrees of freedom"]]
        }
      }
    }#end second effects
    else
    {
      ## beta for the mixed effects model
      random.coef.type <- 1
      
      if ( class( valid.prior$random.coef ) == "bayes.distribution" )
      {
        if ( valid.prior$random.coef@name == "non-informative" )
        {
          random.coef.mean <- rep( 0, ncol( X ) )
          random.coef.Cov <- diag( ncol( X ) )
          random.coef.type <- 3
        }
        else
        {
          random.coef.mean <- valid.prior$random.coef@parameters[["mean"]]
          random.coef.Cov <- valid.prior$random.coef@parameters[["covariance matrix"]]

          random.coef.mean <- valid.mean.specification.bhpm( random.coef.mean, ncol( X ) )
          random.coef.Cov <- valid.Cov.specification.bhpm( random.coef.Cov, ncol( X ) )

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
          stop( "bhbm: list of prior parameters for random coefficients is not of the appropriate size." )
        }

        ##all covariances are assume equal

        if ( valid.prior$random.coef[[1]]@name == "non-informative" )
        {
          random.coef.Cov <- diag( ncol( X ) )
          random.coef.mean <- rep ( 0, ncol( X ) )
          random.coef.type <- 1
        }
        else
        {
          random.coef.Cov <- valid.prior$random.coef[[1]]@parameters[["covariance matrix"]]
          random.coef.Cov <- valid.Cov.specification.bhpm( random.coef.Cov, ncol( X ) )         

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
                                         r.mean <- valid.mean.specification.bhpm( r.mean, dim )
                                       }
                                       r.mean
                                     },
                                     y = valid.prior$random.coef,
                                     dim = ncol( X ) )

        }
 
        ## all t's are assume to have the same df's
        if ( valid.prior$random.coef[[1]]@name == "t" )
        {
          random.coef.df <- valid.prior$random.coef[[1]]@parameters[["t degrees of freedom"]]
        }

      }#end several prior parameters

    }#end no second effects


    ##check prior for random coefficients

    random.var.type <- switch( valid.prior$random.var@name,
      "invChisq" =
      {
        random.var.scale <- valid.prior$random.var@parameters[["sigma0.sq"]]
        random.var.nu <- valid.prior$random.var@parameters[["df"]]
        if( length( random.var.scale ) != dim.random )
        	random.var.scale <- rep( random.var.scale[1], dim.random )
        if( length( random.var.nu ) != dim.random )
        	random.var.nu <- rep( random.var.nu[1], dim.random )
        0
      },

      "non-informative power" =
      {
        random.var.power <- valid.prior$random.var@parameters[["power"]]
        if( length( random.var.power ) != dim.random )
        	random.var.power <- rep( random.var.power[1], dim.random )
        1
      },

      "uniform shrinkage" =
      {
        random.var.scale <- valid.prior$random.var@parameters[["median"]]
        if( length( random.var.scale ) != dim.random )
        	random.var.scale <- rep( random.var.scale[1], dim.random )
        2
      },

      "du Mouchel" =
      {
        random.var.scale <- valid.prior$random.var@parameters[["dispersion"]]
        if( length( random.var.scale ) != dim.random )
        	random.var.scale <- rep( random.var.scale[1], dim.random )
        3
      },

      "invWishart" =
      {
        random.var.scale <- valid.Cov.specification.bhpm( valid.prior$random.var@parameters[["scale"]], ncol( X ) )
        random.var.nu <- valid.prior$random.var@parameters[["degrees of freedom"]]
        4
      },

      stop( "bhbm: distribution for random coefficients variance must be either \"invChisq\", \"non-informative power\", \"uniform shrinkage\", \"du Mouchel\", or \"invWishart\" \n" )
    )#end switch

  }#end random effects prior



  ##fixed effects prior
  if ( fixed.effects )
  {
    if ( valid.prior$fixed.coef@name == "non-informative" )
    {
      fixed.coef.mean <- rep( 0.0, ncol( M ) )
      fixed.coef.Cov <- diag( ncol( M ) )
      fixed.coef.type <- 1
    }
    else
    {
      fixed.coef.mean <- valid.prior$fixed.coef@parameters[["mean vector"]]
      fixed.coef.Cov <- valid.prior$fixed.coef@parameters[["covariance matrix"]]

      fixed.coef.mean <- valid.mean.specification.bhpm( fixed.coef.mean, ncol( M ) )
      fixed.coef.Cov <- valid.Cov.specification.bhpm( fixed.coef.Cov, ncol( M ) )

      if ( valid.prior$fixed.coef@name == "t" )
      {
        fixed.coef.df <- valid.prior$fixed.coef@parameters[["t degrees of freedom"]]
      }
      fixed.coef.type <- 0
    }
  }#end fixed effects

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

  #Metropolis Hastings proposal flag:
  #if non-zero, then covariance for proposal depends on proposal
  #otherwise covariance is kep constant during simulation
  update.cov <- sampler$update.cov


  #one common value for all groups (glm paremeters)
  starting.glm.type <- 0

  ## flag for default starting points
  number.chains <- sampler$number.chains
  read.init.point <- 0
  init.point <- sampler$init.point
  if ( is.null( init.point ) || init.point$type == "prior" )
  {
    ## flag for specify/mixture starting points
    read.init.point = 1

    #add this line for compatibility with Poisson models
    trials.in.model <- 1
    # draw number.chains initial values at random from "prior"
    if ( ( binomial.model == "beta-conj" ) || ( binomial.model == "none" ) )
    {
      starting.points <- generateInitPoints.bhpm( number.draws = number.chains,
                                                  n.responses = n.responses, 
                                                  Y = Y, Expos = NTrials, X = X, 
                                                  M = M, Z = Z,
                                                  dim.random = dim.random, 
                                                  dim.fixed = dim.fixed, 
                                                  dim.level2 = dim.level2, 
                                                  exposure.in.model = trials.in.model,
                                                  glm.pars = glm.pars, glm.type = glm.type,
                                                  common.glm = common.glm,
                                                  random.coef.df = random.coef.df, 
                                                  random.coef.mean = random.coef.mean, 
                                                  random.coef.Cov = random.coef.Cov, 
                                                  random.coef.type = random.coef.type,
                                                  fixed.coef.df = fixed.coef.df, 
                                                  fixed.coef.mean = fixed.coef.mean, 
                                                  fixed.coef.Cov = fixed.coef.Cov,
                                                  random.var.nu = random.var.nu, 
                                                  random.var.scale = random.var.scale, 
                                                  random.var.type = random.var.type,
                                                  level2.coef.df = level2.coef.df, 
                                                  level2.coef.mean = level2.coef.mean, 
                                                  level2.coef.Cov = level2.coef.Cov,
                                                  mix.with.MLE = 0, debug = debug )
    }
    else if( binomial.model == "logit-normal" )
    {
      starting.points <- generateInitPoints.bhpmil( number.chains, n.responses, Y, NTrials, X, M, Z,
                                                dim.random, dim.fixed, dim.level2, trials.in.model,
                                                error.var.nu, error.var.scale, error.var.type, error.var.common,
                                                random.coef.df, random.coef.mean, random.coef.Cov, random.coef.type,
                                                fixed.coef.df, fixed.coef.mean, fixed.coef.Cov,
                                                random.var.nu, random.var.scale, random.var.type,
                                                level2.coef.df, level2.coef.mean, level2.coef.Cov,
                                                mix.with.MLE = 0, debug = debug )         
    }

    starting.glm.type <- 1

    if (debug){
      cat( "starting points are \n")
      print(starting.points)
    }

  }#end prior init points
  else if ( !is.null( init.point ) && init.point$type == "user's choice" )
  {
    read.init.point <- 1

    starting.points <- validate.initial.points.bhpm( n.groups, ncol( X ), ncol( M ), ncol( Z ), 
      init.point$values[[ 1 ]], random.effects, fixed.effects, second.effects, 
      valid.prior$common.glm, random.coef.type, random.var.type, 
      (binomial.model == "beta-conj") )

    if ( !is.null( starting.points$random.var ) )    
    {
      s.random.var <- vector( "list", number.chains )
      s.random.var[[1]] <- starting.points$random.var
    }

    if ( number.chains > 1 )
    {
      for ( i in seq( 2, number.chains ) )
      {
        s.points <- validate.initial.points.bhpm( n.groups, ncol( X ), ncol( M ), ncol( Z ), 
          init.point$values[[ i ]], random.effects, fixed.effects, second.effects, 
          valid.prior$common.glm, random.coef.type, random.var.type, 
          (binomial.model == "beta-conj") )

        if ( !is.null( starting.points$glm ) )
        {
          starting.points$glm <- cbind( starting.points$glm, s.points$glm )
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
        if (length(s.random.var[[1]]) != ncol(X)) {
          stop("length(random.var) is ", length(s.random.var[[1]]), " does not math number of random terms (", ncol(X), ")")
        }
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

    if ( is.null( starting.points$glm ) )
    {
      starting.points$glm <-  matrix( 1, ncol = number.chains ) 
    }


    if ( nrow( starting.points$glm ) > 1 )
    {
      #different starting points
      starting.glm.type <- 1
    }

  }#end user's choice
  else if ( !is.null( init.point ) )
  {
    stop( "bhbm: initialization procedure [", init.point$type, "] has not been implemented yet." )
  }


####end Initial Points------------------------------------------------------------------------------------------


####-------------------------------------------------------------------------------------------
  ##Fit the model
  ##

  if ( is.null( sampler$sampler ) )
  {
    sampler.type <- 0
  }
  else if ( sampler$sampler == "Metropolis" )
  {
    sampler.type <- 0
  }
  else
  {
    stop( "bhbm: sampler type [", sampler$sampler, "] has not been implemented." )
  }  


  bhpmodel <- vector( "list", number.chains )

  for ( i in seq( 1, number.chains ) )
  {
    simulation.seed <- .Random.seed

    if ( (binomial.model == "beta-conj") || (binomial.model == "none") )
    {
      fit <- fit.bayeshbm( n.groups, n.responses, dim.random, dim.fixed, dim.level2, 
        response.names, trials.names, glm.names, random.names, 
        fixed.names, level2.names, group.names, X, M, Z, Y, NTrials, random.coef.mean, 
        random.coef.Cov, random.coef.df, random.coef.type, fixed.coef.mean, 
        fixed.coef.Cov, fixed.coef.df, fixed.coef.type, level2.coef.mean, 
        level2.coef.Cov, level2.coef.df, level2.coef.type, random.var.nu, 
        random.var.scale, random.var.power, random.var.type, glm.pars, glm.type, common.glm,
        read.init.point, starting.points$random.coef[,i], 
        starting.points$fixed.coef[,i], starting.points$level2.coef[,i], 
        starting.points$random.var[[i]], starting.points$glm[,i], starting.glm.type, 
        sampler.type, burnInLength, simulationsToPerform, sampleFrequency, 
        update.cov, print.stats, dimnames( Y )[[2]], debug = debug )
    }
    else
    {
      fit <- fit.bayeshbmil( n.groups, n.responses, dim.random, dim.fixed, dim.level2, 
        response.names, trials.names, glm.names, random.names, 
        fixed.names, level2.names, group.names, X, M, Z, Y, NTrials, random.coef.mean, 
        random.coef.Cov, random.coef.df, random.coef.type, fixed.coef.mean, 
        fixed.coef.Cov, fixed.coef.df, fixed.coef.type, level2.coef.mean, 
        level2.coef.Cov, level2.coef.df, level2.coef.type, error.var.nu, 
        error.var.scale, error.var.power, error.var.common, error.var.type, 
        random.var.nu, random.var.scale, random.var.power, random.var.type, 
        read.init.point, starting.points$random.coef[,i], 
        starting.points$fixed.coef[,i], starting.points$level2.coef[,i], 
        starting.points$random.var[[i]], starting.points$glm[,i], sampler.type, 
        burnInLength, simulationsToPerform, sampleFrequency, print.stats, 
        dimnames( Y )[[2]], debug = debug )
    }
    bhpmodel[[i]] <- fit

  }#end for chains

  return( posterior( sims = mcmc.list(bhpmodel), call = match.call() ) )
  
}#end






############################################################################################
##
## FIT BAYES HIERARCHICAL BINOMIAL REGRESSION MODEL ( LOGIT-NORMAL )
##
############################################################################################

fit.bayeshbmil <- function( n.groups, n.responses, dim.random, dim.fixed, dim.level2,
                          response.names, trials.names, glm.names.orig,
                          random.names, fixed.names, level2.names, group.names, X, M, Z, Y, NTrials,
                          random.coef.mean, random.coef.Cov, random.coef.df, random.coef.type, 
                          fixed.coef.mean, fixed.coef.Cov, fixed.coef.df, fixed.coef.type, 
                          level2.coef.mean, level2.coef.Cov, level2.coef.df, level2.coef.type, 
                          error.var.nu, error.var.scale, error.var.power, error.var.common, error.var.type, 
                          random.var.nu, random.var.scale, random.var.power, random.var.type,
                          read.init.point, starting.random.coef, 
                          starting.fixed.coef, starting.level2.coef, 
                          starting.random.var, starting.glm,
                          sampler.type = 0, burnInLength, simulationsToPerform, sampleFrequency, 
                          print.stats, observation.names = NULL, debug = F )
{

  number.data <- sum( n.responses )

  #now count how many variables there are
  n.vars <- 0
  output.dim <- 0

  #check sigma is not known
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
  }
  output.dim <- n.vars


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

    #account for random coef variance
    if ( random.coef.type != 3 )
    {
      n.vars <- n.vars + 1
      if ( random.var.type == 4 )
      {
        dim.random.var <- dim.random * dim.random
      }
      else
      {
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
  # Fix up initial estimate arguments to have proper length
  if (length(starting.random.coef) != n.groups * dim.random) {
    if (length(starting.random.coef) == 1) {
       starting.random.coef <- rep(starting.random.coef, length=n.groups*dim.random)
    } else {
       stop("length(starting.random.coef) (=", length(starting.random.coef), ") is not compatible with n.groups*dim.random (", deparse(n.groups), "*", deparse(dim.random), ")")
    }
  }

  #now add space for lambdas
  output.dim <- output.dim + number.data
  n.vars <- n.vars + number.data

  #now create output array with the right dimensions
  output.samples <- matrix( 0, nrow = simulationsToPerform, ncol = output.dim )

  if ( (sampler.type == 0) && print.stats )
  {
    #the Metropolis Hastings sampler case
    mh.stats <- rep( 0, n.vars )
  }
  else
  {
    mh.stats <- 0
  }

  #call the function
  trials.flag <- 1

  #Binomial model
  poisson.fit <- 0

  if (debug){
    cat("n.groups:", n.groups, "\n")
    cat("n.responses", n.responses, "\n")
    cat("dim.random:", dim.random, "\n")
    cat("dim.fixed:", dim.fixed, "\n")
    cat("dim.level2:", dim.level2, "\n")
    cat("trials.flag", trials.flag, "\n")
    cat("X", X, "\n")
    cat("M", M, "\n")
    cat("Z", Z, "\n")
    cat("Y", Y, "\n")
    cat("NTrials", NTrials, "\n")
    cat("fixed.coef.mean", fixed.coef.mean, "\n")
    cat("fixed.coef.Cov", fixed.coef.Cov, "\n")
    cat("fixed.coef.df", fixed.coef.df, "\n")
    cat("fixed.coef.type", fixed.coef.type, "\n")
    cat("random.coef.mean", random.coef.mean, "\n")
    cat("random.coef.Cov", random.coef.Cov, "\n")
    cat("random.coef.df", random.coef.df, "\n")
    cat("random.coef.type", random.coef.type, "\n")
    cat("random.var.nu", random.var.nu, "\n")
    cat("random.var.scale", random.var.scale, "\n")
    cat("random.var.power", random.var.power, "\n")
    cat("random.var.type", random.var.type, "\n")
    cat("error.var.nu", error.var.nu, "\n")
    cat("error.var.scale", error.var.scale, "\n")
    cat("error.var.power", error.var.power, "\n")
    cat("error.var.common", error.var.common, "\n")
    cat("error.var.type", error.var.type, "\n")
    cat("output.dim:", output.dim, "\n")
  }
  
  fit <- .C( "fitBayesianHPMIL",
             as.integer( n.groups ),
             as.integer( n.responses ),
             as.integer( dim.random ), 
             as.integer( dim.fixed ), 
             as.integer( dim.level2 ), 
             as.integer( trials.flag ),
             #
             as.double( t(X) ), 
             as.double( t(M) ), 
             as.double( t(Z) ), 
             as.double( Y ), 
             as.double( t( NTrials ) ),
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
             as.double( random.var.nu ), 
             as.double( random.var.scale ),
             as.double( random.var.power ),
             as.integer( random.var.type ),
             #
             as.double( error.var.nu ), 
             as.double( error.var.scale ), 
             as.double( error.var.power ),
             as.integer( error.var.common ),
             as.integer( error.var.type ),
             #
             as.integer( read.init.point ), 
             as.double( starting.random.coef ),
             as.double( starting.fixed.coef ),
             as.double( starting.level2.coef ),
             as.double( starting.random.var ),
             #
             as.double( starting.glm ),
             #
             as.integer(burnInLength ), 
             as.integer(simulationsToPerform ), 
             as.integer(sampleFrequency ), 
             #
             as.integer( poisson.fit ),
             as.integer( print.stats ),
             output.samples = as.double( output.samples ),
             mh.stats = as.double( mh.stats ) )


  # OUTPUT
  # ------
  # The simulationsToPerform output simulations are returned in the array output_simulations.
  # This should be an array of length 
  #   J * [p] + [r] + q + < 1 | q*q > + [J] + [1] + [1] + < J | 1 > + [n.obs] ) * simulationsToPerform
  # where 
  # J: number of groups
  # p: dimension of beta vector
  # r: dimension of gamma vector
  # q: dimension of alpha vector
  # < 1 | q*q >: dimension of tau2 (either a random variable or a random matrix )
  # J: augmented variables for t-distributed betas
  # 1: augmented variable for t-distributed gamma
  # 1: augmented variable for t-distributed alpha
  # J: number of group variables for glm parameters
  #
  # < .|. > denotes a choice.
  # [.] denotes an optional argument.
  #
  # The above dimensions specification also denotes the order in which the simulations are returned:
  #       beta (by group), gamma , alpha, tau | tau2, 
  #             t_beta (by group), t_gamma, t_alpha, glm parameters (by group)
  #
  #
  # Note: Simulations for tau (not tau2) are returned if tau2 is a scalar,
  #       otherwise, the matrices tau2 are returned.
  #


  out.samples <- matrix( fit$output.samples, nrow = simulationsToPerform, ncol = output.dim )

  original.random.names <- random.names
  group.names <- as.character( group.names )

  end.index <- 0
  if ( dim.random > 0 )
  {
    random.coef <- matrix( out.samples[ , seq( 1, n.groups * dim.random ) ], nrow = simulationsToPerform )
    random.names <- as.vector( outer( random.names, group.names, paste, sep = ":" ) )
   
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

  if ( dim.random > 0 )
  {
    if ( random.coef.type != 3 )
    {
      if ( random.var.type == 4 )
      {
        random.var <- matrix( out.samples[ , end.index + seq( 1, dim.random * dim.random ) ], nrow = simulationsToPerform )
        random.var.names <- as.vector( t( outer( paste( "random:tau:", seq( 1, dim.random), sep = "" ), seq( 1, dim.random ), paste, sep = "." ) ) )
        end.index <- end.index + dim.random * dim.random
      }
      else
      {
        random.var <- matrix( out.samples[ , end.index + seq( 1, dim.random ) ], nrow = simulationsToPerform )
        random.var.names <- paste("RANDOM:TAU:", original.random.names, sep = "")
        end.index <- end.index + dim.random
      }
    }
    else
    {
      random.var <- NULL
      random.var.names <- NULL
    }
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

  #now the glm parameters
  glm.params <- NULL
  glm.names <- NULL
  if (  error.var.type != 4 )
  {
    if ( error.var.common == 0 || error.var.common == 1 )
    {
      #different error variance for each group
      glm.params <- matrix( out.samples[ , end.index + seq( 1, n.groups ) ], nrow = simulationsToPerform )
      glm.names <- as.vector( outer( glm.names.orig, group.names, paste, sep = ":" ) )
      end.index <- end.index + n.groups
    }
    else
    {
      #same error variance prior
      glm.params <- matrix( out.samples[ , end.index + 1 ], nrow = simulationsToPerform )
      glm.names <- glm.names.orig
      end.index <- end.index + 1
    }
  }

  #include lambda's
  lambdas <- matrix( out.samples[ , end.index + seq( 1, number.data ) ], nrow = simulationsToPerform )
  #Binomial fit
  lambda.names <- as.vector( outer( "THETA", seq( 1, number.data ), paste, sep = ":" ) )

  glm.parameters <- NULL
  level2.effects <- NULL
  random.effects <- NULL
  fixed.effects <- NULL

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

    random.effects <- list( dim = dim.random, names = random.names, orig.names = original.random.names, coef = random.coef, 
      scale = random.var, scale.names = random.var.names, tau = random.tau, level2 = level2.effects )
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

  glm.parameters <- list( glm = glm.params, names = glm.names, lambda = list( values = lambdas, names = lambda.names, par.names = "PROPORTION" ) )

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

  # The overdispersion parameter
  dimnames(glm.params) <- list(NULL, glm.names)
  post.samples <- cbind(post.samples, glm.params)

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

  # The lambdas
  dimnames(lambdas) <- list(NULL, lambda.names)
  post.samples <- cbind(post.samples, lambdas)

  return( mcmc( data = post.samples, 
    thin = sampleFrequency, burnin = burnInLength ) )
}#end







############################################################################################
##
## FIT BAYES HIERARCHICAL BINOMIAL REGRESSION MODEL (beta-conj)
##
############################################################################################

fit.bayeshbm <- function( n.groups, n.responses, dim.random, dim.fixed, dim.level2,
                          response.names, trials.names, glm.names.orig,
                          random.names, fixed.names, level2.names, group.names, X, M, Z, Y, NTrials,
                          random.coef.mean, random.coef.Cov, random.coef.df, random.coef.type, 
                          fixed.coef.mean, fixed.coef.Cov, fixed.coef.df, fixed.coef.type, 
                          level2.coef.mean, level2.coef.Cov, level2.coef.df, level2.coef.type, 
                          random.var.nu, random.var.scale, random.var.power, random.var.type,
                          glm.pars, glm.type, common.glm,
                          read.init.point, starting.random.coef, 
                          starting.fixed.coef, starting.level2.coef, 
                          starting.random.var, starting.glm, starting.glm.type,
                          sampler.type = 0, burnInLength, simulationsToPerform, sampleFrequency, 
                          update.cov, print.stats, observation.names = NULL, debug = F )
{

  number.data <- sum( n.responses )

  #now count how many variables there are
  n.vars <- 0
  output.dim <- 0

  if ( glm.type == 1 )
  {
    if ( common.glm == 0 || common.glm == 1 )
    {
      #different error variance for each group
      n.vars <- n.groups
    }
    else if ( common.glm == 2 )
    {
      #same error variance
      n.vars <- 1
    }
  }
  else
  {
    stop( "fit.bhbm: Wrong prior type [", glm.type, "] for glm parameters.\n" )
  }
  output.dim <- n.vars


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

    #account for random coef variance
    if ( random.coef.type != 3 )
    {
      n.vars <- n.vars + 1
      if ( random.var.type == 4 )
      {
        dim.random.var <- dim.random * dim.random
      }
      else
      {
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

  #now add space for lambdas, if there is overdispersion
  if( common.glm != 3 ) {     # check for overdispersion
    output.dim <- output.dim + number.data
  }
  n.vars <- n.vars + number.data

  #now create output array with the right dimensions
  output.samples <- matrix( 0, nrow = simulationsToPerform, ncol = output.dim )

  if ( (sampler.type == 0) && print.stats )
  {
    #the Metropolis Hastings sampler case
    mh.stats <- rep( 0, n.vars )
  }
  else
  {
    mh.stats <- 0
  }

  # Fix up initial estimate arguments to have proper length
  if (length(starting.random.coef) != n.groups * dim.random) {
    if (length(starting.random.coef) == 1) {
       starting.random.coef <- rep(starting.random.coef, length=n.groups*dim.random)
    } else {
       stop("length(starting.random.coef) (=", length(starting.random.coef), ") is not compatible with n.groups*dim.random (", deparse(n.groups), "*", deparse(dim.random), ")")
    }
  }

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
    cat("NTrials", NTrials, "\n")
    cat("random.coef.mean", random.coef.mean, "\n")
    cat("random.coef.Cov", random.coef.Cov, "\n")
    cat("random.coef.type", random.coef.type, "\n")
    cat("random.var.nu", random.var.nu, "\n")
    cat("random.var.scale", random.var.scale, "\n")
    cat("random.var.power", random.var.power, "\n")
    cat("random.var.type", random.var.type, "\n")
    cat("glm.pars", glm.pars, "\n")
    cat("common.glm", common.glm, "\n")
    cat("update.cov", update.cov, "\n")
    cat("sampler.type", sampler.type, "\n")
    cat("output.dim:", output.dim, "\n")
    cat("read.init.point:", read.init.point, "\n")
    cat("starting.random.coef:", starting.random.coef, "\n")
    cat("starting.fixed.coef:", starting.fixed.coef, "\n")
    cat("starting.level2.coef:", starting.level2.coef, "\n")
    cat("starting.random.var:", starting.random.var, "\n")
    cat("starting.glm.type:", starting.glm.type, "\n")
    cat("starting.glm:", starting.glm, "\n")
    cat("burnInLength:", burnInLength, "\n")
    cat("simulationsToPerform:", simulationsToPerform, "\n")
    cat("sampleFrequency:", sampleFrequency, "\n")
  }

  fit <- .C( "fitBayesianHBM",
             as.integer( n.groups ), # arg 1 (number_groups)
             as.integer( n.responses ), # arg 2 (number_data)
             as.integer( dim.random ),  # arg 3 (dim_beta)
             as.integer( dim.fixed ), 
             as.integer( dim.level2 ), 
             #
             as.double( t(X) ), # arg 6
             as.double( t(M) ), 
             as.double( t(Z) ), 
             as.double( Y ), 
             as.double( t( NTrials ) ),
             #
             as.double( fixed.coef.mean ), # arg 11A (gamma mean)
             as.double( fixed.coef.Cov ), 
             as.double( fixed.coef.df ), 
             as.integer( fixed.coef.type ),
             #
             as.double( random.coef.mean ), # arg 15 (betamean)
             as.double( random.coef.Cov ), 
             as.double( random.coef.df ), 
             as.integer( random.coef.type ), 
             #
             as.double( level2.coef.mean ), # arg 19 (alphamean)
             as.double( level2.coef.Cov ), 
             as.double( level2.coef.df ), 
             as.integer( level2.coef.type ),
             #
             as.double( random.var.nu ), # arg 23 (tauDF)
             as.double( random.var.scale ),
             as.double( random.var.power ),
             as.integer( random.var.type ),
             #
             as.double( glm.pars ), # arg 27 (xi_z0)
             as.integer( common.glm ),
             #
             as.integer( read.init.point ), # arg 29 (read_init_point)
             as.double( starting.random.coef ), # arg 30 (betaInit)
             as.double( starting.fixed.coef ),
             as.double( starting.level2.coef ),
             as.double( starting.random.var ),
             #
             as.integer( starting.glm.type ),
             as.double( starting.glm ),
             #
             as.integer(burnInLength ), 
             as.integer(simulationsToPerform ), 
             as.integer(sampleFrequency ), 
             #
             as.integer( update.cov ),
             as.integer( print.stats ),
             output.samples = as.double( output.samples ),
             mh.stats = as.double( mh.stats ) )


  # OUTPUT
  # ------
  # The simulationsToPerform output simulations are returned in the array output_simulations.
  # This should be an array of length 
  #   J * [p] + [r] + q + < 1 | q*q > + [J] + [1] + [1] + [J] ) * simulationsToPerform
  # where 
  # J: number of groups
  # p: dimension of beta vector
  # r: dimension of gamma vector
  # q: dimension of alpha vector
  # < 1 | q*q >: dimension of tau2 (either a random variable or a random matrix )
  # J: augmented variables for t-distributed betas
  # 1: augmented variable for t-distributed gamma
  # 1: augmented variable for t-distributed alpha
  # J: number of group variables for glm parameters
  #
  # < .|. > denotes a choice.
  # [.] denotes an optional argument.
  #
  # The above dimensions specification also denotes the order in which the simulations are returned:
  #       beta (by group), gamma , alpha, tau | tau2, 
  #             t_beta (by group), t_gamma, t_alpha, glm parameters (by group)
  #
  #
  # Note: Simulations for tau (not tau2) are returned if tau2 is a scalar,
  #       otherwise, the matrices tau2 are returned.
  #

  out.samples <- matrix( fit$output.samples, nrow = simulationsToPerform, ncol = output.dim )

  original.random.names <- random.names
  group.names <- as.character( group.names )

  end.index <- 0
  if ( dim.random > 0 )
  {
    random.coef <- matrix( out.samples[ , seq( 1, n.groups * dim.random ) ], nrow = simulationsToPerform )
    random.names <- as.vector( outer( random.names, group.names, paste, sep = ":" ) )
   
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

  if ( dim.random > 0 )
  {
    if ( random.coef.type != 3 )
    {
      if ( random.var.type == 4 )
      {
        random.var <- matrix( out.samples[ , end.index + seq( 1, dim.random * dim.random ) ], nrow = simulationsToPerform )
        random.var.names <- as.vector( t( outer( paste( "Random:tau:", seq( 1, dim.random), sep = "" ), seq( 1, dim.random ), paste, sep = "." ) ) )
        end.index <- end.index + dim.random * dim.random
      }
      else
      {
        random.var <- matrix( out.samples[ , end.index + seq( 1, dim.random ) ], nrow = simulationsToPerform )
        random.var.names <- paste("RANDOM:TAU:", original.random.names, sep = "")
        end.index <- end.index + dim.random
      }
    }
    else
    {
      random.var <- NULL
      random.var.names <- NULL
    }
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

  #now the glm parameters
  if ( common.glm == 0 || common.glm == 1 )
  {
    #different error variance for each group
    glm.params <- matrix( out.samples[ , end.index + seq( 1, n.groups ) ], nrow = simulationsToPerform )
    glm.names <- as.vector( outer( glm.names.orig, group.names, paste, sep = ":" ) )
    end.index <- end.index + n.groups 
  }
  else if( common.glm == 2 )
  {
    #same error variance
    glm.params <- matrix( out.samples[ , end.index + 1 ], nrow = simulationsToPerform )
    glm.names <- glm.names.orig
    end.index <- end.index + 1
  }

  #include lambda's
  if( common.glm != 3 ){
    lambdas <- matrix( out.samples[ , end.index + seq( 1, number.data ) ], nrow = simulationsToPerform )
    #Binomial fit
    lambda.names <- as.vector( outer( "THETA", seq( 1, number.data ), paste, sep = ":" ) )
  }

  glm.parameters <- NULL
  level2.effects <- NULL
  random.effects <- NULL
  fixed.effects <- NULL

  if ( dim.random > 0 )
  {
    if ( dim.level2 > 0 )
    {
      if ( level2.coef.df > 0 )
      {
        level2.tau <- list( tau = level2.coef.tau, names = "tau:level2" )
      }
      else
      {
        level2.tau <- NULL
      }
      level2.effects <- list( dim = dim.level2, names = level2.names, coef = level2.coef, tau = level2.tau )
    }

    if ( random.coef.df > 0 )
    {
      random.tau <- list( tau = random.coef.tau, names = paste( "tau:random:", seq( 1, n.groups ), sep = "" ) )
    }
    else
    {
      random.tau <- NULL
    }

    random.effects <- list( dim = dim.random, names = random.names, 
      orig.names = original.random.names, coef = random.coef, 
      scale = random.var, scale.names = random.var.names, 
      tau = random.tau, level2 = level2.effects )
  }

  if ( dim.fixed > 0 )
  {
    if ( fixed.coef.df > 0 )
    {
      fixed.tau <- list( tau = fixed.coef.tau, names = "tau:fixed" )
    }
    else
    {
      fixed.tau <- NULL
    }
    fixed.effects <- list( dim = dim.fixed, names = fixed.names, coef = fixed.coef, tau = fixed.tau )
  }

  if( common.glm != 3 ){
    glm.parameters <- list( glm = glm.params, names = glm.names, 
      lambda = list( values = lambdas, names = lambda.names, 
      par.names = "PROPORTION" ) )
  }
  
  ######## Create a named matrix of the posterior samples.
	post.samples <- matrix(nrow=simulationsToPerform, ncol=0)
  
  # First the fixed effects
  if ( dim.fixed > 0 ){
  	dimnames(fixed.coef) <- list(NULL, fixed.names)
		post.samples <- cbind(post.samples, fixed.coef)
		if ( fixed.coef.df > 0 ){
      fixed.tau.mat <- matrix( fixed.coef.tau, ncol=1, dimnames = list(NULL, c("tau:fixed")) )
      post.samples <- cbind(post.samples, fixed.tau.mat)
    }
	}

  # The overdispersion parameter
  if( common.glm != 3 ){
    dimnames(glm.params) <- list(NULL, glm.names)
    post.samples <- cbind(post.samples, glm.params)
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

  # The lambdas
  if( common.glm != 3 ){
    dimnames(lambdas) <- list(NULL, lambda.names)
    post.samples <- cbind(post.samples, lambdas)
  }

  return( mcmc( data = post.samples, 
    thin = sampleFrequency, burnin = burnInLength ) )

}#end




