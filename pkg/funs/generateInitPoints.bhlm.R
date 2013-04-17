
#
#return: draws from random and fixed effects coefficients,
#        as well as for level2 coefficients, error variance
#        and random effects variance.
#        Draws are based on prior.
#        (useful for creating different chains starting at the draws)
#

generateInitPoints.bhlm <- function( number.draws, n.responses, Y, X, M, Z,  
  dim.random, dim.fixed, dim.level2, error.var.nu, error.var.scale, 
  error.var.type, common.error.var,  random.coef.df, random.coef.mean, 
  random.coef.Cov, random.coef.type, fixed.coef.df, fixed.coef.mean, 
  fixed.coef.Cov, random.var.nu, random.var.scale, random.var.type, 
  level2.coef.df, level2.coef.mean, level2.coef.Cov,  mix.with.MLE = 0 )
{

  n.groups <- length( n.responses )
  number.data <- sum( n.responses )

  #now count how many variables there are
  output.dim <- 0

  if ( dim.random )
  {
    output.dim <- output.dim + dim.random * n.groups

    if ( dim.level2 > 0 )
    {
      output.dim <- output.dim + dim.level2
    }

    #if no non-informative prior for random coefs
    if ( random.coef.type != 3 )
    {
      if ( is.matrix( random.var.scale ) )
      {
        output.dim <- output.dim + ncol( random.var.scale ) * nrow( random.var.scale ) 
      }
      else
      {
        # output.dim <- output.dim + 1
        output.dim <- output.dim + dim.random
        
        # if random.var.scale and random.var.nu are scalars, make them vectors of 
        # length dim.random
        if( dim.random > 1 ){
        	if( length( random.var.scale ) != dim.random )
        		random.var.scale <- rep( random.var.scale[1], dim.random )
        	if( length( random.var.nu ) != dim.random )
        		random.var.nu <- rep( random.var.nu[1], dim.random )        
        }
      }
    }
  }

  if ( dim.fixed > 0 )
  {
    output.dim <- output.dim + dim.fixed
  }

  if ( common.error.var == 0 )
  {
    output.dim <- output.dim + n.groups
  }
  else
  {
    output.dim <- output.dim + 1
  }
  n.scales.error.var <- length( error.var.scale )

  init.points <- matrix( 0.0, ncol = number.draws, nrow = output.dim )

  #call the function

  fit <- .C( "getInitialPointsBayesHLM",
             as.integer( number.draws ),
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
             as.double( random.coef.mean ), 
             as.double( random.coef.Cov ), 
             as.double( random.coef.df ), 
             as.integer( random.coef.type ),
             #
             as.double( fixed.coef.mean ), 
             as.double( fixed.coef.Cov ), 
             as.double( fixed.coef.df ), 
             #
             as.double( level2.coef.mean ), 
             as.double( level2.coef.Cov ), 
             as.double( level2.coef.df ), 
             #
             as.double( error.var.nu ), 
             as.double( error.var.scale ), 
             as.integer( n.scales.error.var ),
             as.integer( error.var.type ),
             as.integer( common.error.var ),
             #
             as.double( random.var.nu ), 
             as.double( random.var.scale ),
             as.integer( random.var.type ),
             #
             init.points = as.double( init.points ) )

  init.points <- matrix( fit$init.points, ncol = number.draws, nrow = output.dim )


  #####check number draws###########################
  init.random.coef <- matrix( 0, ncol = number.draws )
  init.random.var <-  matrix( 0, ncol = number.draws )
  init.level2.coef <-  matrix( 0, ncol = number.draws )
  init.fixed.coef <-  matrix( 0, ncol = number.draws )
  init.error.var <-  matrix( 0, ncol = number.draws )

  end.index <- 0
  if ( dim.random > 0 )
  {
    init.random.coef <- matrix( init.points[ seq( 1, n.groups * dim.random ), ], ncol = number.draws )
    end.index <- n.groups * dim.random

  }


  if ( dim.fixed > 0 )
  {
    init.fixed.coef <- matrix( init.points[ end.index + seq( 1, dim.fixed ), ], ncol = number.draws )
    end.index <- end.index + dim.fixed
  }

  if ( dim.random > 0 && dim.level2 > 0 )
  {
    init.level2.coef <- matrix( init.points[ end.index + seq( 1, dim.level2 ), ], ncol = number.draws )
    end.index <- end.index + dim.level2
  }


  if ( common.error.var == 0 )
  {
    #different error variance for each group
    init.error.var <- matrix( init.points[ end.index + seq( 1, n.groups ), ], ncol = number.draws )
    end.index <- end.index + n.groups
  }
  else
  {
    #same error variance prior
    init.error.var <- matrix( init.points[ end.index + 1, ], ncol = number.draws )

    end.index <- end.index + 1
  }


  #if non non-informative
  if ( dim.random > 0 && random.coef.type != 3 )
  {
    init.random.var <- vector( "list", number.draws )
    for ( i in seq( 1, number.draws ) )
    {
      if ( random.var.type == 4 )
      {
        init.random.var[[ i ]] <- matrix( init.points[ end.index + seq( 1, dim.random * dim.random ), i ], nrow = dim.random, ncol = dim.random )
      }
      else
      {
        #init.random.var[[ i ]]  <- matrix( init.points[ end.index + 1, i ], nrow = 1, ncol = 1 )
        init.random.var[[ i ]]  <- matrix( init.points[ end.index + 1:dim.random, i ], nrow = 1, ncol = dim.random )
      }
    }#end for

    if ( random.var.type == 4 )
    {
      end.index <- end.index + dim.random * dim.random
    }
    else
    {
      #end.index <- end.index + 1
      end.index <- end.index + dim.random
    }
  }
  else
  {
    init.random.var <- vector( "list", number.draws )
    for ( i in seq( 1, number.draws ) )
    {
      init.random.var[[i]] <- NULL
    }
  }

  list( error.var = init.error.var, random.coef = init.random.coef, fixed.coef = init.fixed.coef, level2.coef = init.level2.coef, random.var =  init.random.var )

}#end generate init points