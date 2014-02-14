##################################################
# bhpm class in the Bayes Module                 #
# Insightful: NIH Bayes II                       #
# Alejandro Murua 08/03                           #
##################################################




#####################################################################################
##
##  PRIOR
##
#####################################################################################

bhpm.prior <- function( sigma2 = NULL,
                        xi = NULL,
                        fixed.coef = "non-informative",
                        level2.coef = "non-informative",
                        random.var = bayes.nonInfoPower( -1.0 ),
                        common.glm = 3 )
{

  # create a list of indicators of whether the prior has been specified 
  # by the user for each parameter.  This should be used in bhpm to check that the user 
  # has not specified a prior for a parameter that is not in the model.
  priorSpec <- list( sigma2 = !missing( sigma2 ),
    xi = !missing( xi ),
    fixed.coef = !missing( fixed.coef ),
    level2.coef = !missing( level2.coef ),
    random.var = !missing( random.var ) )

  random.coef <- bayes.normal( mean = zero, cov = identity )
  
  # for common.glm:
  # 0 ==> independent overdispersion params; different hyper prior parameters
  # 1 ==> independent overdispersion params; same hyper prior parameters
  # 2 ==> common overdispersion parameter
  # 3 ==> no overdispersion parameter

  # xi is the overdispersion parameter in the Poisson-gamma and beta-binomial 
  # models, while sigma2 is the overdispersion parameter in the Poisson 
  # log-normal and binomial logit-normal models
  if( !is.null( sigma2 ) && is.null( xi ) ){
    glm <- sigma2
    conjModel <- F
    # conjModel indicates whether the model is conjugate,
    # i.e. whether this is the Poisson-gamma or beta-binomial model as 
    # opposed to the Poisson log-normal or binomial logit-normal model
  } else if( !is.null( xi ) && is.null( sigma2 ) ) {
    glm <- xi
    conjModel <- T
  } else {
    # if there is an overdispersion parameter, its prior must be specified
    if( common.glm == 3 ){
      glm <- NULL
      conjModel <- F
    } else {
      stop( "Must specify a prior for exactly one of sigma2 or xi, or specify no overdispersion parameter by setting 'common.glm = 3'." )
    } 
  }

  if ( !is.null( glm ) && class( glm ) == "bayes.distribution" )
  {
    if ( common.glm == 0 ){
      stop( "bhpm.prior: need a list of overdispersion prior parameters if common.glm = 0 in the Poisson Model." )
    } else if ( common.glm == 3 ){
      stop( "bhpm.prior: If common.glm = 3 then no prior should be specified for xi or sigma2." )
    }

    if ( !is.element( glm@name, c( "invChisq", "uniform shrinkage", "du Mouchel", "non-informative power", "mass point" ) ) )
    {
      stop( "bhpm.prior: Invalid prior distribution for Poisson Model parameter." )
    }
  }
  else if ( !is.null( glm ) && is.list( glm ) )
  {
    if ( common.glm != 0 )
    {
      stop( "bhpm.prior: need a single overdispersion prior when common.glm != 0 for the Poisson model." )
    }

    for ( i in seq( 1, length( glm ) ) )
    {
      if ( class( glm[[i]] ) == "bayes.distribution" )
      {
        if ( !is.element( glm[[i]]@name, c( "invChisq", "uniform shrinkage", "du Mouchel", "non-informative power", "mass point" ) ) )
        {
          stop( "bhpm.prior: Invalid prior distribution for Poisson model parameters." )
        }
      }
      else
      {
        stop( "bhpm.prior: Invalid prior distribution for Poisson model parameters." )
      } 
    }#end for loop
  }
  else if( common.glm != 3 )   
  {
    stop( "bhpm.prior: Invalid prior distribution for Poisson model parameters." )
  }  

  
  ##
  ## check the variance for the random coefficients
  ##

  if ( !is.null( random.var ) && class( random.var ) == "bayes.distribution" )
  {
    if ( !is.element( random.var@name, c( "invChisq", "uniform shrinkage", "du Mouchel", "invWishart", "non-informative power" ) ) )
    {
      stop( "bhpm.prior: Invalid prior distribution for random effects coefficients variance." )
    }
  }
  else if ( !is.null( random.var ) && is.list( random.var ) )
  {
    for ( i in seq( 1, length( random.var ) ) )
    {
      if ( class( random.var[[i]] ) == "bayes.distribution" )
      {
        if ( !is.element( random.var[[i]]@name, c( "invChisq", "uniform shrinkage", "du Mouchel", "non-informative power" ) ) )
        {
          stop( "bhpm.prior: Invalid prior distribution for random effects coefficients." )
        }        
      } else {
        stop( "bhpm.prior: Invalid prior distribution for random effects coefficients." )
      }
    }
  } else if ( !is.null( random.var ) )
  {
    stop( "bhpm.prior: Invalid prior distribution for random effects coefficients variance." )
  }



  ##
  ## check the prior for random coefficients
  ## if a list then check each element of the list
  ##

  if ( !is.null( random.coef ) && class( random.coef ) == "bayes.distribution" )
  {
    if ( !is.element( random.coef@name, c( "normal", "t", "non-informative" ) ) )
    {
      stop( "bhpm.prior: Invalid prior distribution for random effects coefficients." )
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
          stop( "bhpm.prior: Invalid prior distribution for random effects coefficients." )
        }        
      }
      else if ( random.coef[[i]] != "non-informative" )
      {
        stop( "bhpm.prior: Invalid prior distribution for random effects coefficients." )
      }
      else if ( random.coef[[i]] == "non-informative" )  
      {
        random.coef[[i]] <- bayes.nonInformative()
      }
    }#end for loop
  }
  else if ( !is.null( random.coef )  && random.coef != "non-informative" )
  {
    stop( "bhpm.prior: Invalid prior distribution for random effects coefficients." )
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
      stop( "bhpm.prior: Invalid prior distribution for fixed effects coefficients." )
    }
  }
  else if ( !is.null( fixed.coef ) && fixed.coef != "non-informative" )
  {  
    stop( "bhpm.prior: Invalid prior distribution for fixed effects coefficients." )
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
      stop( "bhpm.prior: Invalid prior distribution for second level coefficients." )
    }
  }
  else if ( !is.null( level2.coef ) && level2.coef != "non-informative" )
  {  
    stop( "bhpm.prior: Invalid prior distribution for second level coefficients." )
  }
  else if ( !is.null( random.coef )  && level2.coef == "non-informative" )  
  {
    level2.coef <- bayes.nonInformative()
  }



  ##
  ## create prior object for hierarchical model
  ##

  list( glm = glm, 
        random.coef = random.coef,
        fixed.coef = fixed.coef,
        level2.coef = level2.coef,
        random.var = random.var,
        common.glm = common.glm,
        conjModel = conjModel,
        priorSpec = priorSpec )

}#end



#####################################################################################
##
##  SAMPLER
##
#####################################################################################


bhpm.sampler <- function( nBurnin=1000, nSamples=1000, nThin=1,
                          nChains = 1,
                          update.cov = 1,
                          init.point = "prior",
                          params.init = NULL )
{
  random.effects = T
  fixed.effects = T
                          
  control <- sampler.control( nBurnin = nBurnin, nSamples = nSamples, nThin = nThin )
  number.chains <- nChains

  #assume params.init is a list of lists with initial points for each chain

  if ( number.chains <= 0 )
  {
    warning( "bhpm.sampler: zero or negative number of chains to simulate. Assuming 1 chain." )
    number.chains <- 1
  }

  #just in case this is not an integer as it should be  
  if ( ceiling( number.chains ) - number.chains > 0 )
  {
    new.number.chains <- ceiling( number.chains )
    warning( "bhpm.sampler: Non-integer number of chains specified [", number.chains, "]. Set number of chains to [", new.number.chains, "]." )
    number.chains <- new.number.chains
  }

  if ( update.cov != 0 )
  {
    update.cov <- 1
  }

  init.point <- match.arg( init.point, choices = c( "prior", "user's choice" ) )
  if ( init.point == "user's choice" )
  {
    if ( is.null( params.init ) )
    {
      stop( "bhpm.sampler: initial points must be specified." )
    }
    else
    { 
      if ( number.chains == 1 && !is.list( params.init[[1]] ) )
      {
        params.init <- list( params.init )
      }
      if ( length( params.init ) < number.chains )
      {
        stop( "bhpm.sampler: initial values for all parameters in all chains must be specified." )
      }
      else
      {
        if( !is.null( params.init[[1]]$xi ) && is.null( params.init[[1]]$sigma2 ) )
          conjModel <- T
        else if( !is.null( params.init[[1]]$sigma2 ) && is.null( params.init[[1]]$xi ) )
          conjModel <- F
        else
          # Case of no overdispersion parameter
          conjModel <- F
          
        for ( i in seq( 1, number.chains ) )
        {
          #cat("init values are\n")
          #print( params.init[[i]] )

          if( conjModel )
            glm <- params.init[[ i ]]$xi
          else if( !is.null( params.init[[ i ]]$sigma2 ) )
            # if there is an initial value specified for sigma2,
            glm <- params.init[[ i ]]$sigma2
          random.coef <- params.init[[ i ]]$random.coef
          fixed.coef <- params.init[[ i ]]$fixed.coef
          level2.coef <-  params.init[[ i ]]$level2.coef
          random.var <- params.init[[ i ]]$random.var

          #check glm parameters
          if ( !is.null( glm ) && !is.vector( glm ) && !is.matrix( glm ) && !is.list( glm ) )
          {
            stop( "bhpm.sampler: if initial values are specified for overdispersion parameters, these must be vector, matrix, or list form." )
          }
          else if ( is.list( glm ) )
          {
            if ( any ( diag( glm[[1]] < 0 ) ) )
            {
              stop( "bhlm.sampler: initial value for overdispersion parameters must be positive." )
            }

            this.glm <- glm[[1]]
            if ( length( glm ) > 1 )
            {
              for (  j in seq( 2, length( glm ) ) )
              {
                if ( any ( diag( glm[[j]] < 0 ) ) )
                {
                  stop( "bhpm.sampler: initial value for Poisson Model parameters must be positive." )
                }
                this.glm <- cbind( this.glm, glm[[j]] )
              }
            }
            params.init[[ i ]]$glm <- this.glm
          }

          
          if ( random.effects && is.null( level2.coef ) )
          {     
            if ( is.null( random.coef ) )
            {
              stop( "bhpm.sampler: initial value for random effects coefficients must be provided." )
            }
          }

          if ( random.effects && !is.null( random.coef ) && !is.null( level2.coef ) )
          {
            #warning( "bhpm.sampler: initial values for random effects coefficients are not allowed." )
          }


          if ( fixed.effects && is.null( fixed.coef ) )
          {
            stop( "bhpm.sampler: initial values for the fixed effects coefficients must be provided." )
          }

          #check parameters variance
          if ( is.null( random.var ) )
          {
            stop( "bhpm.sampler: initial value for random effects coefficients scale must be provided." )
          }    
          else if ( is.vector( random.var ) )
          {
            if ( any( random.var ) <= 0 )
            {
              stop( "bhpm.sampler: Initial value for random effects coefficients scale must be positive. " )
            }
          }
          else if ( is.matrix( random.var ) )
          {
            if ( any( diag( random.var ) ) <= 0 )
            {
              stop( "bhpm.sampler: Initial value for random effects coefficients variances must be positive. " )
            }

            if ( any( random.var - t( random.var ) != 0 ) )
            {
              stop( "bhpm.sampler: Initial value for random effects coefficients variances must be a symmetric matrix." )
            }
          }  
          else if ( is.list( random.var ) ){
            for( j in 1:length( random.var ) ){
              if ( !is.numeric( random.var ) || (random.var <= 0) ){
                stop( "bhpm.sampler: When initial values for random effects coefficients variances are specified as a list, the elements must be numeric and positive." )
              }
            }            
          } else
          {
            stop( "bhpm.sampler: Initial value for random effects coefficients scale are not valid." )
          }
        }#end for i
      }
    }
  }#end user's choice

  type <- "Metropolis"
  
  list( control = control, 
        number.chains = number.chains,
        update.cov = update.cov,
        init.point = list( values = params.init, type = init.point ), 
        sampler = type )
}#end




#####################################################################################
##
##  BINOMIAL MODEL PRIOR
##
#####################################################################################

bhbm.prior <- bhpm.prior


#####################################################################################
##
##  BINOMIAL MODEL SAMPLER
##
#####################################################################################


bhbm.sampler <- bhpm.sampler


#####################################################################################
##
##  VALIDATION FOR MEAN AND COVARIANCE
##
#####################################################################################

valid.mean.specification.bhpm <- function( x, dim )
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
    stop( "bhpm: the prior mean for one of the regression coefficients has the wrong dimension." )
  }

  x
}#end


valid.Cov.specification.bhpm <- function( x, dim )
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
      stop( "bhpm: covariance main diagonal must be positive." )
    }
    
    if ( length( x ) == 1 && dim > 1 )
    {
      x <- x * diag( dim )
    }
    else if ( length( x ) != dim )
    {
      stop( "bhpm: covariance specification has wrong dimension." )
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
      stop( "bhpm: covariance specification has wrong dimensions." )
    }
    else if ( any( x - t(x) != 0 ) )
    {
      stop( "bhpm: covariance specification is not symmetric." )
    }
    else if ( any( diag( x ) <= 0 ) )
    {
      stop( "bhpm: covariance specification contains negative or zero variances." )
    }    
  }
  else
  {
    stop( "bhpm: covariance specification is not a matrix!." )
  }

  x
}#end



##############################################################################################
##
##VALIDATION OF INITIAL POINTS
##
##############################################################################################

validate.initial.points.bhpm <- function( n.groups, dimX, dimM, dimZ, init.point, 
  random.effects, fixed.effects, second.effects, common.glm, random.coef.type, 
  random.var.type, conjModel )
{
  s.glm <- NULL
  s.random.coef <- NULL
  s.fixed.coef <- NULL
  s.random.var <- NULL
  s.level2.coef <- NULL
  
  if( conjModel ){    # if this is a Poisson-gamma or beta-binomial model
    if( is.null( init.point$xi ) || !is.null( init.point$sigma2 ) )
      stop( "Must specify an initial value for xi but not for sigma2 for this model" )
    glm <- init.point$xi
  } else if( common.glm != 3 ){    # if this is a Poisson log-normal or binomial logit-normal model
    if( is.null( init.point$sigma2 ) || !is.null( init.point$xi ) )
      stop( "Must specify an initial value for exactly one of xi, sigma2, depending on the model" )
    glm <- init.point$sigma2
  } else {    # this model has no overdispersion parameter
    if( !is.null( init.point$sigma2 ) || !is.null( init.point$xi ) )
      stop( "For the model without overdispersion, cannot specify a prior for xi or sigma2." )
  }
  
  if( is.element( common.glm, c(1,2) ) )
  {
    if ( length( glm ) > 1 )
    {
      stop( "bhpm: initial value for Poisson Model parameters should not be an array." )
    }
    else
    {
      s.glm <- matrix( glm, nrow = 1 )
    }
  }
  else if( common.glm == 0 )
  {
    if ( length( glm ) != n.groups )
    {
      stop( "bhpm: initial values for Poisson Model parameters should be an array of number of groups elements." )
    }
    else
    {
      s.glm <- matrix( glm, nrow = n.groups )  
    }
  } else {  # no overdispersion parameter
    s.glm <- NULL
  }


  if ( random.effects )
  {
    if ( second.effects )
    {
      if ( is.null( init.point$level2.coef ) )
      {
        stop( "bhpm: initial values for second stage regression coefficients must be provided." )
      }
      else if ( !is.vector( init.point$level2.coef ) )
      {
        stop( "bhpm: invalid initial values for second stage regression coefficients." )
      }
      else if ( length( init.point$level2.coef ) != dimZ )
      {
        stop( "bhpm: initial values for second stage regression coefficients have wrong dimension." )
      }       
      else
      {
        s.level2.coef <- matrix( init.point$level2.coef, nrow = dimZ )
      }
    }

    #if no non-informative
    if ( random.coef.type != 3 )
    {
      if ( random.var.type == 4 )
      {
        #inv wishart case
        if ( is.null( init.point$random.var ) )
        {
          stop( "bhpm: initial value for the covariance matrix associated to the random coefficients must be provided." )
        }

        s.random.var <- valid.Cov.specification.bhpm( init.point$random.var, dimX )
      }
      else
      {
        if ( is.null( init.point$random.var ) )
        {
          stop( "bhpm: initial value for the variance associated to the random coefficients must be provided." )
        }

        if ( is.vector( init.point$random.var ) )
        {
          if ( any( init.point$random.var <= 0 ) )
          {
            stop( "bhpm: initial value for the variance associated to the random coefficients must be positive." )
          }

          s.random.var <- matrix( init.point$random.var, nrow = 1, ncol=dimX )
        }
        else
        {
          stop( "bhpm: wrong specification for initial value for the variance associated to the random coefficients." ) 
        }

      }#end scalar / vector variance specification
    }
    else
    {
      s.random.var <- matrix( 0, nrow = 1, ncol=dimX )
    }


    if ( random.coef.type != 0 )
    {
      if ( is.null( init.point$random.coef ) )
      {
        stop( "bhpm: initial values for random coefficients must be specified." )
      }

      if ( random.coef.type == 1 )
      {  # single prior for all groups
        if ( is.vector( init.point$random.coef ) && length( init.point$random.coef ) == dimX )
        {
          s.random.coef <- matrix( init.point$random.coef, nrow = dimX )
        }
        else if ( is.matrix( init.point$random.coef ) && ( nrow( init.point$random.coef ) == dimX ) && 
          ( ncol( init.point$random.coef ) == n.groups ) ) {}
        else
        {
          stop( "bhpm: wrong specification for initial values for random coefficients." )
        }
      }
      else #random.coef.type == 2 or 3
      {
        if ( is.list( init.point$random.coef ) && length( init.point$random.coef ) == n.groups )
        {
          s.random.coef <- apply( matrix( seq( 1, n.groups ), nrow = 1 ), 
                                   MARGIN = 2, 
                                   FUN = function( x, y, dim )
                                   {
                                     r.mean <- y[[ x ]]
                                     r.mean <- valid.mean.specification.bhpm( r.mean, dim )
                                     r.mean
                                   },
                                   y = init.point$random.coef
                                   dim = dimX )            

        }
        else if ( is.matrix( init.point$random.coef ) && ncol( init.point$random.coef ) == n.groups && nrow( init.point$random.coef ) == dimX )
        {
          s.random.coef <- init.point$random.coef
        }
        else if ( is.vector( init.point$random.coef ) && length( init.point$random.coef ) == dimX && n.groups == 1 )
        {
          s.random.coef <- matrix( init.point$random.coef, ncol = 1 )
        }
        else
        {
          stop( "bhpm: wrong specification for initial values of random effects coefficients." )
        }
      }
    }#end no second level regression

  }#end if random effects

  if ( fixed.effects )
  {
    if ( is.null( init.point$fixed.coef ) )
    {
      stop( "bhpm: initial values of fixed effects coefficients must be provided." )
    }

    if ( is.vector( init.point$fixed.coef ) && length( init.point$fixed.coef ) == dimM )
    {
      s.fixed.coef <- matrix( init.point$fixed.coef, nrow = dimM )
    }
    else
    {
      stop( "bhpm: wrong specification for initial values of fixed effects coefficients." )
    }
  }#end if fixed effects

  list( glm = s.glm, random.coef = s.random.coef, fixed.coef = s.fixed.coef, random.var = s.random.var, level2.coef = s.level2.coef )

}#end validate initial point




#####################################################################################
##
##  BAYESIAN HIERARCHICAL POISSON MODEL
##
#####################################################################################

# exposure.formula, random.formula, fixed.formula, level2.formula and group.formula 
# are of type "formula"
bhpm <- function( fixed.formula = NULL, 
                  exposure.formula = NULL,
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

  poisson.model <- overdispersion
  poisson.model <- match.arg( poisson.model, 
    choices = c("gamma-conj", "log-normal", "none") )

  if( is.null( prior ) ){
    if( poisson.model == "gamma-conj" )
      prior <- bhpm.prior( xi = bayes.uniformShrinkage(0.5), common.glm = 2 )
    else if( poisson.model == "log-normal" )
      prior <- bhpm.prior( sigma2 = bayes.uniformShrinkage(0.5), common.glm = 2 )
    else   # no overdispersion case
      prior <- bhpm.prior()
  }

  if( ( poisson.model == "gamma-conj" ) && !prior$conjModel ){
    stop("For the gamma-conjugate model, a prior must be specified for xi")
  } else if( ( poisson.model == "log-normal" ) && prior$conjModel ){
    stop("A prior was specified for the xi parameters, which do not exist in the log-normal model")
  } else if( ( poisson.model == "none" ) && ( prior$common.glm != 3 ) ){
    stop("For models with no overdispersion parameter, the prior cannot be specified with common.glm = 0, 1, or 2")
  }

  ####-------------------------------------------------------
  ## Get Data 
  ##

  if ( is.null( data ) )
  {
    stop( "bhpm: data must be provided." )
  }

  data <- as.data.frame( data )

  if ( !is.null( na.action ) )
  {
    #get rid of rows with missing data if na.omit, otherwise stop and fail (na.fail)

    data <- na.action( data )
    if ( nrow( data ) == 0 )
    {
      stop( "bhpm: every observation (row) of the data set contains at least a missing value.\n" )
    }
  }
  else
  {
    na.action <- na.fail
  }


  #name of variables must be preserved for output
  response.names <- NULL
  exposure.names <- NULL
  random.names <- NULL
  fixed.names <- NULL
  level2.names <- NULL

  used.contrasts <- NULL

  X <- 0
  M <- 0
  Z <- 0
  Expos <- 0

  random.effects <- F
  fixed.effects <- F
  exposure.in.model <- F
  random.vars <- -1
  fixed.vars <- -1

  if ( !is.null( random.formula ) )
  {
    if ( length( terms( random.formula )@term.labels ) == 0 && terms( random.formula )@response == 0 )
    {
      #may contain only the intercept
      if ( terms( random.formula )@intercept != 1 )
      {
        stop( "bhpm: formula for random effects must contain at least one variable or intercept." )
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
        stop( "bhpm: formula for random effects must contain at least one variable or intercept." )
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
        stop( "bhpm: formula for fixed effects must contain at least one variable or intercept." )
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
        stop( "bhpm: formula for fixed effects must contain at least one variable or intercept." )
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
    stop( "bhpm: model specification is not valid. You need to provide a response variable." )  
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
        stop( "bhpm: formula for second level effects must contain at least one variable or intercept." )
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


  #now get exposure
  if ( !is.null( exposure.formula ) )
  {
    model.expos <- call( "model.frame", formula = exposure.formula, data = data, na.action = na.action )
    model.expos <- eval( model.expos, sys.parent() )
    Terms <- attr( model.expos, "terms" )
    Terms@intercept <- 0

    Expos <- model.matrix( Terms, model.expos )
    exposure.in.model <- T
    exposure.names <- dimnames( Expos )[[2]]
  }
  else
  {
    exposure.names <- "(Exposure)"
  }



  #sort data by group
  if ( is.null( Y ) )
  {
    stop( "bhpm: response variable must be provided.\n" )
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

    if ( exposure.in.model )
    {
      Expos <- as.matrix( Expos[ idx, ] )
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
    stop( "bhpm: a prior has been specified for the random effect variance, but there are no random effects in the model." )
  if( !fixed.effects && prior$priorSpec$fixed.coef )
    stop( "bhpm: a prior has been specified for fixed effects, but there are no fixed effects in the model." )
  if( !second.effects && prior$priorSpec$level2.coef )
    stop( "bhpm: a prior has been specified for level 2 effects, but there are no level 2 effects in the model." )

  if( poisson.model == "gamma-conj" )
    valid.prior <- bhpm.prior( xi = prior$glm, fixed.coef = prior$fixed.coef, level2.coef = prior$level2.coef, 
      random.var = prior$random.var, common.glm = prior$common.glm )
  else if( poisson.model == "log-normal" )
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

  if( poisson.model == "gamma-conj" )
  {
    lambda.out <- F

    # the only prior implemented for xi in the gamma-conj model is the 
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
      stop( "bhpm: Invalid prior specification for xi parameter.\n" )
    }
  }
  else if( poisson.model == "log-normal" )
  {
    #assume log-normal model

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

      stop( "bhpm: distribution for error variance must be either \"invChisq\", \"non-informative power\", \"uniform shrinkage\", \"du Mouchel\", or \"mass point\" \n" )

      )#end switch
    }
    else
    {
      # different prior parameters for error variance

      if ( length( valid.prior$glm ) != n.groups )
      {
        stop( "bhpm: list of prior parameters for error variance is not of the appropriate size." )
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

          stop( "bhpm: distribution for error variance must be either \"invChisq\", \"non-informative power\", \"uniform shrinkage\", \"du Mouchel\", or \"mass point\" \n" )
          )#end switch

          if ( i == 1 )
          {
            prev.type <- error.var.type
          }
          else if ( prev.type != error.var.type )
          {
            stop( "bhpm: all priors for error variances should be of the same type." )
          }
        
        }#end for loop
      }#end list
      else
      {
        stop( "bhpm: prior for error variance is not valid." )
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
          stop( "bhpm: list of prior parameters for random coefficients is not of the appropriate size." )
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


    ##check prior for random coef. variance

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

      stop( "bhpm: distribution for random coefficients variance must be either \"invChisq\", \"non-informative power\", \"uniform shrinkage\", \"du Mouchel\", or \"invWishart\" \n" )
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

  poisson.fit <- 1

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

    # draw number.chains initial values at random from "prior"
    if ( ( poisson.model == "gamma-conj" ) || ( poisson.model == "none" ) )
    {
      starting.points <- generateInitPoints.bhpm( number.draws = number.chains,
                                                  n.responses = n.responses, 
                                                  Y = Y, Expos = Expos, X = X, 
                                                  M = M, Z = Z,
                                                  dim.random = dim.random, 
                                                  dim.fixed = dim.fixed, 
                                                  dim.level2 = dim.level2, 
                                                  exposure.in.model = exposure.in.model,
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
                                                  mix.with.MLE = 0 )
    }
    else if( poisson.model == "log-normal" )
    {
      starting.points <- generateInitPoints.bhpmil( number.chains, n.responses, Y, Expos, X, M, Z,
                                                dim.random, dim.fixed, dim.level2, exposure.in.model,
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
      (poisson.model == "gamma-conj") )

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
          (poisson.model == "gamma-conj") )

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
    stop( "bhpm: initialization procedure [", init.point$type, "] has not been implemented yet." )
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
    stop( "bhpm: sampler type [", sampler$sampler, "] has not been implemented." )
  }  


  bhpmodel <- vector( "list", number.chains )

  for ( i in seq( 1, number.chains ) )
  {
    simulation.seed <- .Random.seed

    if ( (poisson.model == "gamma-conj") || (poisson.model == "none") )
    {
      fit <- fit.bayeshpm( n.groups, n.responses, dim.random, dim.fixed, dim.level2, 
        exposure.in.model, response.names, exposure.names, glm.names, random.names, 
        fixed.names, level2.names, group.names, X, M, Z, Y, Expos, random.coef.mean, 
        random.coef.Cov, random.coef.df, random.coef.type, fixed.coef.mean, 
        fixed.coef.Cov, fixed.coef.df, fixed.coef.type, level2.coef.mean, 
        level2.coef.Cov, level2.coef.df, level2.coef.type, random.var.nu, 
        random.var.scale, random.var.power, random.var.type, glm.pars, glm.type, common.glm,
        read.init.point, starting.points$random.coef[,i], 
        starting.points$fixed.coef[,i], starting.points$level2.coef[,i], 
        starting.points$random.var[[i]], starting.points$glm[,i], starting.glm.type, 
        sampler.type, burnInLength, simulationsToPerform, sampleFrequency, 
        update.cov, print.stats, dimnames( Y )[[2]], debug )
    }
    else 
    {
      fit <- fit.bayeshpmil( n.groups, n.responses, dim.random, dim.fixed, dim.level2, 
        exposure.in.model, response.names, exposure.names, glm.names, random.names, 
        fixed.names, level2.names, group.names, X, M, Z, Y, Expos, random.coef.mean, 
        random.coef.Cov, random.coef.df, random.coef.type, fixed.coef.mean, 
        fixed.coef.Cov, fixed.coef.df, fixed.coef.type, level2.coef.mean, 
        level2.coef.Cov, level2.coef.df, level2.coef.type, error.var.nu, 
        error.var.scale, error.var.power, error.var.common, error.var.type, 
        random.var.nu, random.var.scale, random.var.power, random.var.type, 
        read.init.point, starting.points$random.coef[,i], 
        starting.points$fixed.coef[,i], starting.points$level2.coef[,i], 
        starting.points$random.var[[i]], starting.points$glm[,i], sampler.type, 
        burnInLength, simulationsToPerform, sampleFrequency, poisson.fit, print.stats, 
        dimnames( Y )[[2]], debug = debug )
    }
    bhpmodel[[i]] <- fit

  }#end for chains

  return( posterior( sims = mcmc.list(bhpmodel), call = match.call() ) )
        
}#end






############################################################################################
##
## FIT BAYES HIERARCHICAL POISSON REGRESSION MODEL ( LOG-NORMAL )
##
############################################################################################

fit.bayeshpmil <- function( n.groups, n.responses, dim.random, dim.fixed, dim.level2, exposure.in.model,
                          response.names, exposure.names, glm.names.orig,
                          random.names, fixed.names, level2.names, group.names, X, M, Z, Y, Expos,
                          random.coef.mean, random.coef.Cov, random.coef.df, random.coef.type, 
                          fixed.coef.mean, fixed.coef.Cov, fixed.coef.df, fixed.coef.type, 
                          level2.coef.mean, level2.coef.Cov, level2.coef.df, level2.coef.type, 
                          error.var.nu, error.var.scale, error.var.power, error.var.common, error.var.type, 
                          random.var.nu, random.var.scale, random.var.power, random.var.type,
                          read.init.point, starting.random.coef, 
                          starting.fixed.coef, starting.level2.coef, 
                          starting.random.var, starting.glm,
                          sampler.type = 0, burnInLength, simulationsToPerform, sampleFrequency, 
                          poisson.fit, print.stats, observation.names = NULL, debug = F )
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
  if ( exposure.in.model )
  {
    exposure.flag <- 1
  }
  else
  {
    exposure.flag <- 0
  }

  if (debug) {
    cat("n.groups:", n.groups, "\n")
    cat("n.responses", n.responses, "\n")
    cat("dim.random:", dim.random, "\n")
    cat("dim.fixed:", dim.fixed, "\n")
    cat("dim.level2:", dim.level2, "\n")
    cat("exposure.flag", exposure.flag, "\n")
    cat("X", X, "\n")
    cat("M", M, "\n")
    cat("Z", Z, "\n")
    cat("Y", Y, "\n")
    cat("Expos", Expos, "\n")
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
    cat("poisson.fit:", poisson.fit, "\n")
  }
  fit <- .C( "fitBayesianHPMIL",
             as.integer( n.groups ),
             as.integer( n.responses ),
             as.integer( dim.random ), 
             as.integer( dim.fixed ), 
             as.integer( dim.level2 ), 
             as.integer( exposure.flag ),
             #
             as.double( t(X) ), 
             as.double( t(M) ), 
             as.double( t(Z) ), 
             as.double( Y ), 
             as.double( t( Expos ) ),
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
      #same error variance
      glm.params <- matrix( out.samples[ , end.index + 1 ], nrow = simulationsToPerform )
      glm.names <- glm.names.orig
      end.index <- end.index + 1
    }
  }

  #include lambda's
  lambdas <- matrix( out.samples[ , end.index + seq( 1, number.data ) ], nrow = simulationsToPerform )
  end.index <- end.index + number.data

  if ( poisson.fit == 1 )
  {
    lambda.names <- as.vector( outer( "LAMBDA", seq( 1, number.data ), paste, sep = ":" ) )
  }
  else
  {
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

  glm.parameters <- list( glm = glm.params, names = glm.names, lambda = list( values = lambdas, names = lambda.names, par.names = "INTENSITY" ) )

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

  # The "scale" parameter xi
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

  return(mcmc(post.samples, 
    thin = sampleFrequency, burnin = burnInLength))

}#end






############################################################################################
##
## FIT BAYES HIERARCHICAL POISSON REGRESSION MODEL
##
############################################################################################

fit.bayeshpm <- function( n.groups, n.responses, dim.random, dim.fixed, dim.level2, exposure.in.model,
                          response.names, exposure.names, glm.names.orig,
                          random.names, fixed.names, level2.names, group.names, X, M, Z, Y, Expos,
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
    stop( "fit.bhpm: Wrong prior type [", glm.type, "] for glm parameters.\n" )
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

  #call the function
  if ( exposure.in.model )
  {
    exposure.flag <- 1
  }
  else
  {
    exposure.flag <- 0
  }

  if (debug) {
    cat("n.groups:", n.groups, "\n")
    cat("n.responses", n.responses, "\n")
    cat("dim.random:", dim.random, "\n")
    cat("dim.level2:", dim.level2, "\n")
    cat("dim.fixed:", dim.fixed, "\n")
    cat("exposure.flag", exposure.flag, "\n")
    cat("X", X, "\n")
    cat("M", M, "\n")
    cat("Z", Z, "\n")
    cat("Y", Y, "\n")
    cat("Expos", Expos, "\n")
    cat("random.coef.mean", random.coef.mean, "\n")
    cat("random.coef.Cov", random.coef.Cov, "\n")
    cat("random.coef.type", random.coef.type, "\n")
    cat("random.var.nu", random.var.nu, "\n")
    cat("random.var.scale", random.var.scale, "\n")
    cat("random.var.power", random.var.power, "\n")
    cat("random.var.type", random.coef.type, "\n")
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
  
  fit <- .C( "fitBayesianHPM",
             as.integer( n.groups ),
             as.integer( n.responses ),
             as.integer( dim.random ), 
             as.integer( dim.fixed ), 
             as.integer( dim.level2 ), 
             as.integer( exposure.flag ),
             #
             as.double( t(X) ), 
             as.double( t(M) ), 
             as.double( t(Z) ), 
             as.double( Y ), 
             as.double( t( Expos ) ),
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
             as.double( glm.pars ),
             as.integer( common.glm ),
             #
             as.integer( read.init.point ), 
             as.double( starting.random.coef ),
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
    #Poisson fit
    lambda.names <- as.vector( outer( "LAMBDA", seq( 1, number.data ), paste, sep = ":" ) )
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
      par.names = "INTENSITY" ) )
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




