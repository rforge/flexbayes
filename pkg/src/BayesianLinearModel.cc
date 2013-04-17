#include <math.h>
#include <string.h>
#include "R.h"

#include "DistributionParameter.h"
#include "BayesianLinearModel.h"



/* constructor.
   returns: an empty Distribution of type BayesianLinearModel
*/
BayesianLinearModel::BayesianLinearModel()
{

  emptyModel();

}//end



/* returns: an Distribution of type BayesianLinearModel.
   y is the response
   x is the matrix of predictors
   Assumes normal beta prior Normal( mean = p_beta, covariance = p_betaCov );
           InvChisquare prior for sigma2 degrees of freedom = p_nuSigma and scale = p_sigma2 
           (the scale is sigma0^2 (not sigma0) )

   Important: Note that this Distribution only accept errors with scale (covariance) sigma2 * identity.
              In order to handle sigma2 * Sigma, transform y and x before creating this Distribution,
              e.g. new_y = L^{-1}y,   new_x = L^{-1} x, 
              where Sigma = L L^t, is the Cholesky decomposition of Sigma 
              (one could also use L = Sigma^{1/2}).

*/
void BayesianLinearModel::initialize( CVector * y, CMatrix * x, 
                                      CVector * p_beta, CMatrix * p_betaCov,
                                      double p_sigma2, double p_nuSigma )
{  

  if ( y->Len() > 0 && y->Len() == x->Row() )
  {
    //load data
    response = new CVector( y->Len() );
    (*response) = (*y);
  
    predictors = new CMatrix( x->Row(), x->Col() );
    (*predictors) = (*x);

    //prior information on sigma2 (Inv-Xi2)
    sigma2 = new InvChisqDistribution( p_nuSigma );
    if ( p_sigma2 > 0 )
    {
      sigma2->setScale( p_sigma2 );
    }
    else
    {
      printf( "BayesianLinearModel::nonInformative: Initial scale for sigma is zero or negative. Replacing initial scale to one.\n\n" );
      sigma2->setScale( 1.0 );
    }

    //prior information on beta (assumed for the moment Normal)
    beta = new NormalDistribution( p_beta->Len() );
    beta->setMean( p_beta );
    beta->setCovariance( p_betaCov );

    number_of_variables = 2;

    //mixture parameters
    betas = NULL;
    betas_t = NULL;
    beta_mix = NULL;
    sigma2s = NULL;
    sigma2_mix = NULL;
    mixture_prior = false;

    //set default parameters
    tau2 = NULL;
    tau2_error = NULL;
    t_lkhd = false;
    t_beta = false;
    t_lkhd_beta = false;
    keep_tau2_error = false;
    conjugate_model = false;
    non_informative = false;
    simulations_performed = false;

    simulated_beta = NULL;
    simulated_sigma2 = NULL;

    xTx = NULL;
    xTy = NULL;
    tau2_error_weights = NULL;
    residuals = NULL;

    distr_map = NULL;
  }//end if dimensions are right
  else
  {
    printf("data has wrong dimensions. len(y) = %d, nrow(x) = %d.", y->Len(), x->Row() );
    emptyModel();
  }
    

}//end


/* for normal mixture prior: gamma is the prob( mixture components ) 
*/
void BayesianLinearModel::initialize( CVector * y, CMatrix * x, 
                                      int number_mixtures, CVector ** p_beta, CMatrix ** p_betaCov,
                                      double * gamma, double * d_freedom, char * mixture_type,
                                      double p_sigma2, double p_nuSigma, int save_stats )
{
  int i;

  if ( y->Len() > 0 && y->Len() == x->Row() )
  {
    //load data
    response = new CVector( y->Len() );
    (*response) = (*y);
  
    predictors = new CMatrix( x->Row(), x->Col() );
    (*predictors) = (*x);

    //prior information on sigma2 (Inv-Xi2)
    sigma2 = new InvChisqDistribution( p_nuSigma );
    if ( p_sigma2 > 0 )
    {
      sigma2->setScale( p_sigma2 );
    }
    else
    {
      printf( "BayesianLinearModel::nonInformative: Initial scale for sigma is zero or negative. Replacing initial scale to one.\n\n" );
      sigma2->setScale( 1.0 );
    }

    //prior information on beta (assumed for the moment Normal mixture)

    mixture_prior = true;
    beta_mix = new DistributionMixture( number_mixtures );

    if ( !strcmp( mixture_type, "normal" ) )
    {
      t_beta = false;
      betas = new NormalDistribution * [ number_mixtures ];
      for ( i = 0; i < number_mixtures; i++ )
      {
        betas[i] = new NormalDistribution( p_beta[i]->Len() );
        betas[i]->setMean( p_beta[i] );
        betas[i]->setCovariance( p_betaCov[i] );

        beta_mix->set( i, betas[i], gamma[i] );
      }
      //sigma and beta_mix
      number_of_variables = 2;
    }
    else if ( !strcmp( mixture_type, "t" ) )
    {
      t_beta = true;
      betas_t = new StudentTDistribution * [ number_mixtures ];
      for ( i = 0; i < number_mixtures; i++ )
      {
        betas_t[i] = new StudentTDistribution( p_beta[i]->Len(), d_freedom[i] );
        betas_t[i]->setLocation( p_beta[i] );
        betas_t[i]->setScaleMatrix( p_betaCov[i] );

        beta_mix->set( i, betas_t[i], gamma[i] );
      }      
      //sigma and beta_mix
      number_of_variables = 2;
    }

    if ( save_stats )
    {
      beta_mix->saveStatistics( true );
    }

    //normal distribution
    beta = NULL;
    //sigma2 mixture full conditional
    sigma2s = NULL;
    sigma2_mix = NULL;

    //set default parameters
    tau2 = NULL;
    tau2_error = NULL;
    t_lkhd = false;
    t_lkhd_beta = false;
    keep_tau2_error = false;
    conjugate_model = false;
    non_informative = false;
    simulations_performed = false;

    simulated_beta = NULL;
    simulated_sigma2 = NULL;

    xTx = NULL;
    xTy = NULL;
    tau2_error_weights = NULL;
    residuals = NULL;

    distr_map = NULL;
  }//end if dimensions are right
  else
  {
    printf("data has wrong dimensions. len(y) = %d, nrow(x) = %d.", y->Len(), x->Row() );
    emptyModel();
  }
    

}//end



/* returns: an Distribution of type BayesianLinearModel.
   y is the response
   x is the matrix of predictors
   Assumes a non-informative prior for beta. This can be seen as 
           a normal beta prior Normal( mean = zero , covariance = +infinity );
           InvChisquare prior for sigma2: degrees of freedom = p_nuSigma and scale = p_sigma2 
           (the scale is sigma0^2 (not sigma0) )
           If prior for sigma2 is also non-informative, then prior is InvChisq( df = 0, scale = 1 )
*/
void BayesianLinearModel::nonInformative( CVector * y, CMatrix * x,
                                          double p_sigma2, double p_nuSigma )
{

  if ( y->Len() > 0 && y->Len() == x->Row() )
  {
    //load data
    response = new CVector( y->Len() );
    (*response) = (*y);
  
    predictors = new CMatrix( x->Row(), x->Col() );
    (*predictors) = (*x);

    //non-informative prior on sigma2 (Inv-Xi2 with 0 df)
    sigma2 = new InvChisqDistribution( p_nuSigma );
    if ( p_sigma2 > 0 )
    {
      sigma2->setScale( p_sigma2 );
    }
    else
    {
      //scale should be any positive number
      //printf( "BayesianLinearModel::nonInformative: Initial scale for sigma is zero or negative. Replacing initial scale to one.\n\n" );
      sigma2->setScale( 1.0 );
    }

    //non-information prior on beta (assumed for the moment Normal)
    beta = new NormalDistribution( predictors->Col() );
    beta->setNonInformative();

    number_of_variables = 2;

    //mixture parameters
    betas = NULL;
    betas_t = NULL;
    beta_mix = NULL;
    sigma2s = NULL;
    sigma2_mix = NULL;
    mixture_prior = false;

    //set default parameters
    tau2 = NULL;
    tau2_error = NULL;
    t_lkhd = false;
    t_beta = false;
    t_lkhd_beta = false;
    keep_tau2_error = false;
    conjugate_model = false;

    non_informative = true;
    simulations_performed = false;

    simulated_beta = NULL;
    simulated_sigma2 = NULL;

    xTx = NULL;
    xTy = NULL;
    tau2_error_weights = NULL;
    residuals = NULL;

    distr_map = NULL;
  }//end if dimensions are right
  else
  {
    printf( "data has wrong dimensions. len(y) = %d, nrow(x) = %d.", y->Len(), x->Row() );
    emptyModel();
  }
    
}//end


/* destructor: of an BayesianLinearModel Distribution
 */
BayesianLinearModel::~BayesianLinearModel()
{
  int i;

  if ( response )
  {
    if ( tau2 )
    {
      delete tau2;
    }
    if ( tau2_error )
    {
      for ( i = 0; i < response->Len(); i++ )
      {
        delete tau2_error[i];
      }
      delete [] tau2_error;
    }

    delete response;
    delete  predictors;
    
    if ( sigma2 != NULL )
    {
      delete sigma2;
    }

    if ( mixture_prior )
    {

      //beta_mix->printDrawingStats();

      if ( t_beta )
      {
        for ( i = 0; i < beta_mix->numberOfMixtures(); i++ )
        {
          delete betas_t[i];
        }
        delete [] betas_t;
      }
      else
      {
        for ( i = 0; i < beta_mix->numberOfMixtures(); i++ )
        {
          delete betas[i];
        }
        delete [] betas;

        if ( conjugate_model )
	{
          for ( i = 0; i < sigma2_mix->numberOfMixtures(); i++ )
          {
            delete sigma2s[i];
          }
          delete [] sigma2s;
          delete sigma2_mix;
        }
      }
      delete beta_mix;
      mixture_prior = false;
    }
    else
    {
      delete beta;
    }

    if ( simulated_beta )
    {
      delete simulated_beta;
      delete simulated_sigma2;
      if ( keep_tau2_error )
      {
        delete simulated_tau2_errors;
      }
    }

    if ( xTx )
    {
      delete xTx;
      delete xTy;
    }

    if ( tau2_error_weights )
    {
      delete tau2_error_weights;
    }

    if ( residuals )
    {
      delete residuals;
    }

    if ( distr_map != NULL )
    {
      for ( i = 0; i < number_of_variables; i++ )
      {
        delete [] distr_map[i];
      }
      delete [] distr_map;
    }
  }//end if 

  number_of_variables = 0;

}//end


/* creates an empty BayesianLinearModel Distribution
 */
void BayesianLinearModel::emptyModel()
{
  response = NULL;
  predictors = NULL;
  sigma2 = NULL;
  beta = NULL;

  //mixture parameters
  betas = NULL;
  betas_t = NULL;
  beta_mix = NULL;
  sigma2s = NULL;  
  sigma2_mix = NULL;
  mixture_prior = false;

  tau2 = NULL;
  tau2_error = NULL;
  t_lkhd = false;
  t_beta = false;
  t_lkhd_beta = false;
  keep_tau2_error = false;
  conjugate_model = false;
  non_informative = false;
  simulations_performed = false;

  number_of_variables = 0;

  simulated_beta = NULL;
  simulated_sigma2 = NULL;

  xTx = NULL;
  xTy = NULL;
  tau2_error_weights = NULL;
  residuals = NULL;

  distr_map = NULL;
}//end
  

/* set likelihood distribution to t with p_nuError degrees of freedom.
   For Gibbs sampler purposes, this forces data augmentation by modeling
   each error as Normal( mean = 0, covariance = tau_i^2 ), with 
   tau_i^2 InvChisq( df = p_nuError, scale = sigma2 )
*/
void BayesianLinearModel::t_likelihood( double p_nuError )
{
  int i;

  keep_tau2_error = true;

  if ( !conjugate_model )
  {
    t_lkhd = true;
    tau2_error = new InvChisqDistribution * [ response->Len() ];
    for ( i = 0; i < response->Len(); i++ )
    {
      tau2_error[i] = new InvChisqDistribution( p_nuError );
      tau2_error[i]->setScale( 1 );
    }

    number_of_variables += response->Len();
  }
  else
  {
    printf( "Likelihood must be normally-distributed for conjugate prior specification." );
  }
}//end


/* set prior distribution for beta to a t with df = p_nuBeta, location = p_beta,
   and scale = p_betaCov.
   For Gibbs sampler purposes, this forces data augmentation by modeling
   beta as Normal( mean = p_beta, covariance = tau^2 p_betaCov ) with
   tau^2 InvChisq( df = p_nuBeta, scale = 1 )
*/
void BayesianLinearModel::t_beta_prior( double p_nuBeta )
{
  t_beta = true;
  tau2 = new InvChisqDistribution( p_nuBeta );
  tau2->setScale( 1 );

  number_of_variables++;

}//end


/* set likelihood distribution to t with p_nuError degrees of freedom,
   and prior distribution for beta to a t with df = p_nuBeta, location = p_beta,
   and scale = p_betaCov.
   See t_likelihood() and t_beta_prior() above.
*/
void BayesianLinearModel::t_likelihood_and_beta_prior( double p_nuError, double p_nuBeta )
{
  t_lkhd_beta = true;
  t_likelihood( p_nuError );
  t_beta_prior( p_nuBeta );

}//end


void BayesianLinearModel::conjugatePrior()
{
  int i;

  conjugate_model = true;

  //create mixture for sigma
  if ( mixture_prior )
  {
    //prior information on sigma2 (Inv-Xi2)
    sigma2s = new InvChisqDistribution * [ beta_mix->numberOfMixtures() ];
    sigma2_mix = new DistributionMixture( beta_mix->numberOfMixtures() );

    for ( i = 0; i < beta_mix->numberOfMixtures(); i++ )
    {
      sigma2s[i] = new InvChisqDistribution( sigma2->initialDegreesOfFreedom() );
      sigma2s[i]->setScale( sigma2->initialScale() );

      sigma2_mix->set( i, sigma2s[i], beta_mix->initialProportion( i ) );
    }

    delete sigma2;
    sigma2 = NULL;
  }

}//end



/* set default starting values for beta and sigma2 for Gibbs sampler
 */
void BayesianLinearModel::samplerDefaultInitialPoint()
{
  //set default starting points for simulation
  if ( mixture_prior )
  {
#ifdef FIX1
    DistributionParameter tmpdist = beta_mix->mean();
    beta_mix->setLastDraw( tmpdist );
#else
    beta_mix->setLastDraw( beta_mix->mean() );
#endif
//    beta_mix->setLastDraw( beta_mix->mean() );
    if ( conjugate_model )
    {
      DistributionParameter val( sigma2s[0]->scale() );
      sigma2_mix->setLastDraw( val );
    }
    else
    {
      sigma2->setLastDraw( sigma2->scale() );
    }
  }
  else
  {
#ifdef FIX1
    DistributionParameter tmpdist1 = beta->mean();
    DistributionParameter tmpdist2 = sigma2->scale();
    beta->setLastDraw( tmpdist1 );
    sigma2->setLastDraw( tmpdist2 );
#else
    beta->setLastDraw( beta->mean() );
    sigma2->setLastDraw( sigma2->scale() );
#endif
//    beta->setLastDraw( beta->mean() );
//    sigma2->setLastDraw( sigma2->scale() );
  }



}//end


/* set a user-supplied starting point for beta
 */
void BayesianLinearModel::samplerBetaInitialPoint( CVector & init_beta )
{
  if ( mixture_prior )
  {
    DistributionParameter beta0( init_beta );
    beta_mix->setLastDraw( beta0 );
  }
  else
  {
    beta->setLastDraw( init_beta );
  }
}//end
  

/* set a user-supplied starting point for sigma2
 */
void BayesianLinearModel::samplerSigma2InitialPoint( double init_sigma2 )
{
  if ( mixture_prior && conjugate_model )
  {
    DistributionParameter val( init_sigma2 );
    sigma2_mix->setLastDraw( val );
  }
  else
  {
    sigma2->setLastDraw( init_sigma2 );
  }
}//end


void BayesianLinearModel::createOutput( int simulations_to_keep )
{

  if ( simulations_to_keep > 0 )
  {
    simulated_beta = new CMatrix (simulations_to_keep, predictors->Col() );
    simulated_sigma2 = new CVector ( simulations_to_keep );

    if ( keep_tau2_error )
    {
      simulated_tau2_errors = new CMatrix( simulations_to_keep, response->Len() );
    }
  }

}//end


void BayesianLinearModel::dataAugmentationInitialDraws()
{
  int i;

  //initial draws from tau2 and tau2_error
  if ( t_lkhd_beta || t_beta )
  {
    if ( !mixture_prior )
    {
      tau2->draw();
    }
  }
  if ( t_lkhd_beta || t_lkhd )
  {
    for ( i = 0; i < response->Len(); i++ )
    {
      tau2_error[i]->draw();
      //keep it in an array of weights for regression
      updateRegressionWeight( i );
    }
  }

}//end


void BayesianLinearModel::updateRegressionWeight( int index )
{
  tau2_error_weights->Val( index ) = 1 / tau2_error[ index ]->lastItemDrawn();
}//end



void BayesianLinearModel::initializeTemporaryStructures()
{
  int i, index;

  //create working matrices
  xTx = new CMatrix( predictors->Col(), predictors->Col() );
  xTy = new CVector( predictors->Col() );

  (*xTx) = predictors->xTransposedX();
  (*xTy) = predictors->T() * (*response);

  //create indices for variables
  distr_map = new char * [ number_of_variables ];

  distr_map[0] = new char [ 5 ];
  sprintf( distr_map[0], "beta" );
  distr_map[1] = new char [ 6 ];
  sprintf( distr_map[1], "sigma" );

  index = 2;
  if ( t_beta && !mixture_prior )
  {
    distr_map[ index ] = new char [ 5 ];
    sprintf( distr_map[ index ], "tau2" );
    index++;
  }
  
  //printf( "Bayesian: t-likelihood is %d\n", t_lkhd );

  if ( tLikelihood() )
  {
    for ( i = 0; i < response->Len(); i++ )
    {
      //six characters for "error:", 10 characters for index
      distr_map[ index ] = new char [ 6 + 10 + 1 ];
      sprintf( distr_map[ index ], "tau2e:%d", i );
      index++;
    }
    tau2_error_weights = new CVector( response->Len() );
  }

  //printf(" distr map is: n vars = %d\n", number_of_variables);
  //for ( i = 0; i < number_of_variables; i++ )
  //{
  // printf("[%d] = %s\n", i, distr_map[ i ] );
  //}

  //initialize residuals
  residuals = new CVector( response->Len() );

}//end


void BayesianLinearModel::keepSimulation( int simulations_kept )
{
  if ( simulations_kept < simulated_sigma2->Len() )
  {
    if ( mixture_prior )
    {
#ifdef FIX1
      CVector tmpvec = beta_mix->lastDraw().getVector();
      simulated_beta->setRow( simulations_kept, tmpvec );
#else
      simulated_beta->setRow( simulations_kept, beta_mix->lastDraw().getVector() );
#endif
//      simulated_beta->setRow( simulations_kept, beta_mix->lastDraw().getVector() );
      if ( conjugate_model )
      {
        simulated_sigma2->Val( simulations_kept ) = sigma2_mix->lastDraw().getScalar();
      }
      else
      {
        simulated_sigma2->Val( simulations_kept ) = sigma2->lastItemDrawn();
      }
    }
    else
    {
#ifdef FIX1
      CVector tmpvec = beta->lastItemDrawn();
      simulated_beta->setRow( simulations_kept, tmpvec );
#else
      simulated_beta->setRow( simulations_kept, beta->lastItemDrawn() );
#endif
//      simulated_beta->setRow( simulations_kept, beta->lastItemDrawn() );
      simulated_sigma2->Val( simulations_kept ) = sigma2->lastItemDrawn();
    }

    if ( keep_tau2_error )
    {
      int i;
      for ( i = 0; i < response->Len(); i++ )
      {
        simulated_tau2_errors->Val( simulations_kept, i ) = tau2_error[i]->lastItemDrawn();
        //printf(" tau2[%d] = %f ", i, simulated_tau2_errors->Val( simulations_kept, i ) );
      }
      //printf("\n");
    }
  }
  else
  {
    printf( "BayesianLinearModel::keepSimulation: Index [%d] out of range. Maximum is [%d].\n", simulations_kept, simulated_sigma2->Len() );
  }

}//end


void BayesianLinearModel::drawVariable( int index )
{
  if ( index < number_of_variables )
  {
    if ( !strcmp( distr_map[ index ], "beta" ) )
    {
      if ( mixture_prior )
      {
        beta_mix->draw();

        if ( conjugate_model )
        {
          betas[ beta_mix->lastComponentDrawn() ]->scaleLastDraw( sqrt( sigma2_mix->lastDraw().getScalar() ) );
        }

        //printf( "Bayesian: mix draw = \n");
        //beta_mix->lastDraw().Print(); fflush(stdout);

      }
      else
      {
        beta->draw();
        if ( conjugate_model )
        {
          beta->scaleLastDraw( sqrt( sigma2->lastItemDrawn() ) );
        }
      }
    }
    else if ( !strcmp( distr_map[ index ], "sigma" ) )
    {
      if ( mixture_prior && conjugate_model )
      {
        sigma2_mix->draw();
      }
      else
      {
        sigma2->draw();
        //printf("sigma draw = %f\n", sigma2->lastItemDrawn() ); fflush(stdout);
      }
    }
    else if ( !strcmp( distr_map[ index ], "tau2" ) )
    {
      if ( !mixture_prior )
      {
        tau2->draw();
      }
    }
    else if ( !strncmp( distr_map[ index ], "tau2e", 5 ) )
    {
      int error_index;
      char * ptr, * temp_distr;

      temp_distr = new char [ strlen( distr_map[ index ] ) + 1 ];
      sprintf( temp_distr, "%s", distr_map[ index ] );
      ptr = strtok( temp_distr, ":" );
      if ( ptr != NULL )
      {
        ptr = strtok( NULL, ":" );
        if ( ptr != NULL )
	{
          error_index = atoi( ptr );
          if ( error_index >= 0 && error_index < response->Len() )
	  {
            tau2_error[ error_index ]->draw();
            updateRegressionWeight( error_index );
          }
          else
          {
            printf( " BayesianLinearModel::drawVariable: Wrong argument in [%s]. Index out of range.\n", distr_map[ index ] );
          }
        }
        else
        {
          printf( " BayesianLinearModel::drawVariable: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }
      else
      {
        printf( " BayesianLinearModel::drawVariable: Wrong argument in [%s]. Number expected.\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }
    else
    {
      printf( " BayesianLinearModel::drawVariable: Unknown variable [%s].\n", distr_map[ index ] );
    }
  }
  else
  {
    printf( "BayesianLinearModel::drawVariable: Variable index [%d] does not exist.\n", index );
  }

}//end


void BayesianLinearModel::fullConditionalUpdateVariable( int index )
{
  if ( index < number_of_variables )
  {
    if ( !strcmp( distr_map[ index ], "beta" ) )
    {
      gibbsUpdateBeta();
    }
    else if ( !strcmp( distr_map[ index ], "sigma" ) )
    {
      gibbsUpdateSigma2();
    }
    else if ( !strcmp( distr_map[ index ], "tau2" ) )
    {
      gibbsUpdateTau2();
    }
    else if ( !strncmp( distr_map[ index ], "tau2e", 5 ) )
    {
      int error_index;
      char * ptr, * temp_distr;

      temp_distr = new char [ strlen( distr_map[ index ] ) + 1 ];
      sprintf( temp_distr, "%s", distr_map[ index ] );
      ptr = strtok( temp_distr, ":" );
      if ( ptr != NULL )
      {
        ptr = strtok( NULL, ":" );
        if ( ptr != NULL )
	{
          error_index = atoi( ptr );
          if ( error_index >= 0 && error_index < response->Len() )
	  {
            gibbsUpdateTau2Error( error_index );
          }
          else
          {
            printf( " BayesianLinearModel::fullConditionalUpdateVariable: Wrong argument in [%s]. Index out of range.\n", distr_map[ index ] );
          }
        }
        else
        {
          printf( " BayesianLinearModel::fullConditionalUpdateVariable: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }
      else
      {
        printf( " BayesianLinearModel::fullConditionalUpdateVariable: Wrong argument in [%s]. Number expected.\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }
    else
    {
      printf( " BayesianLinearModel::fullConditionalUpdateVariable: Unknown variable [%s].\n", distr_map[ index ] );
    }
  }
  else
  {
    printf( "BayesianLinearModel::fullConditionalUpdateVariable: Variable index [%d] does not exist.\n", index );
  }

}//end




/* returns: beta's and sigmas into one array for output
 */
double * BayesianLinearModel::simulationsToArray( int simulations_to_keep )
{
  int i, j;
  double * simul_output;

  simul_output = new double  [ ( simulated_beta->Col() + 1 ) * simulations_to_keep ];

  for ( i = 0; i < simulations_to_keep; i++ ) 
  {
    for ( j = 0; j < simulated_beta->Col(); j++ )
    {
      simul_output[ j * simulations_to_keep + i ] = simulated_beta->Val( i, j );
    }

    simul_output[ simulated_beta->Col() * simulations_to_keep + i ] = sqrt( simulated_sigma2->Val( i ) );
  }//end for i

  return ( simul_output );
}//end


/* returns: beta's and sigmas into one array for output.
            array must be supplied by user.
 */
void BayesianLinearModel::simulationsToArray( double * simul_output, int simulations_to_keep )
{
  int i, j;

  for ( i = 0; i < simulations_to_keep; i++ ) 
  {
    for ( j = 0; j < simulated_beta->Col(); j++ )
    {
      simul_output[ j * simulations_to_keep + i ] = simulated_beta->Val( i, j );
    }

    simul_output[ simulated_beta->Col() * simulations_to_keep + i ] = sqrt( simulated_sigma2->Val( i ) );

    if ( keep_tau2_error )
    {
      for ( j = 0; j < simulated_tau2_errors->Col(); j++ )
      {
        simul_output[ ( simulated_beta->Col() + j + 1 ) *  simulations_to_keep + i ] = sqrt( simulated_tau2_errors->Val( i, j ) );
      }
    }
  }//end for i

}//end



void BayesianLinearModel::saveMixtureDrawingStatistics( double * stats )
{
  if ( mixture_prior && beta_mix->statisticsAreSaved() )
  {
    beta_mix->printDrawingStats( stats );
    
  }
}//end


/* update parameters of conditional distribution for tau^2.
   This is done when beta prior is a t distribution.
   updated scale is based on M_D, the Mahalonobis distance from the current beta to the initial one.
   actual update of parameters is done by the call to  InvChisqDistribution::update()
*/
void BayesianLinearModel::gibbsUpdateTau2()
{
  //update
  double s_beta;

  CVector beta_diff( beta->dimension() );
  beta_diff = beta->lastItemDrawn() - beta->initialMean();

  s_beta = M_D( beta->initialInverseCovariance(), beta_diff );

  tau2->update( beta->dimension(), s_beta ); 

}//end


/* updates the data matrices transposed(x) * x and transposed(x) * y
   x is the matrix of predictors, and y is the response vector.
   When the error priors are modeled as t-distibutions, then the data must be
   transformed by the diagonal matrix diag( 1/tau_i ).
*/
void BayesianLinearModel::gibbsUpdateWorkingMatrices()
{
  if ( t_lkhd_beta || t_lkhd )
  {
    //weighted (or transformed) predictors. 
    //update must be made before updating beta
    (*xTx) = predictors->xTransposedX( (*tau2_error_weights) );
    (*xTy) = ( predictors->T() ) * ( response->weighted( (*tau2_error_weights) ) );
  }

  //printf( "BLM: working matrices: xTx and xTy:\n" );
  //xTx->Print();
  //xTy->Print();

}//end


/* updates the parameters of the conditional distribution of beta.
   When beta prior is modeled as a t-distribution, then the precision
   of 1/tau^2 enters into the updating equations.
   actual update is made with a call to NormalDistribution::update().
*/
void BayesianLinearModel::gibbsUpdateBeta()
{
  //update beta

  if ( !mixture_prior )
  {
    if ( !conjugate_model )
    {  
      CMatrix sigma_beta( beta->dimension(), beta->dimension() );
      CVector mu_beta( beta->dimension() );

      //this update must be done first
      gibbsUpdateWorkingMatrices();

      mu_beta = (*xTy) * (1/sigma2->lastItemDrawn());
      sigma_beta = (1/sigma2->lastItemDrawn()) * (*xTx);

      if ( t_beta || t_lkhd_beta )
      {
        double precision_tau2 = 1/tau2->lastItemDrawn();
        beta->update( precision_tau2, mu_beta, sigma_beta );
      }
      else
      {
        beta->update( mu_beta, sigma_beta );
      }
    }
    else // conjugate_model
    {
      if ( !simulations_performed )
      {
        beta->update( (*xTy), (*xTx) );
        simulations_performed = true;
      }
    }

    //printf( "BLM: beta mean is \n" );
    //beta->mean().Print();
    //printf( "BLM: beta: initial inv Cov beta is \n" );
    //beta->initialInvCovBeta().Print();
    //printf( "BLM: chol Inv Cov is \n");
    //beta->choleskyInvCov().Print();

  }//no mixture
  else
  {
    int i;

    if ( !conjugate_model )
    {  
      CMatrix sigma_beta( predictors->Col(), predictors->Col() );
      CVector mu_beta( predictors->Col() );

      //this update must be done first
      gibbsUpdateWorkingMatrices();

      mu_beta = (*xTy) * (1/sigma2->lastItemDrawn());
      sigma_beta = (1/sigma2->lastItemDrawn()) * (*xTx);

      for ( i = 0; i < beta_mix->numberOfMixtures(); i++ )
      {
        if ( t_beta || t_lkhd_beta )
        {
          betas_t[i]->update( mu_beta, sigma_beta );
          updateMixtureProportion( i );
        }
        else
        {
          betas[i]->update( mu_beta, sigma_beta );
          updateMixtureProportion( i );
        }
      }//end for i
      //printf("bayesian: finish update beta mixture\n"); fflush(stdout);

      beta_mix->normalizeLogProportions();
      //beta_mix->viewCurrentProportions(); fflush(stdout);
      //beta_mix->viewInitialProportions(); fflush( stdout);
    }
    else
    {
      if ( !simulations_performed )
      {
        for ( i = 0; i < beta_mix->numberOfMixtures(); i++ )
        {
          betas[i]->update( (*xTy), (*xTx) );
        }
        simulations_performed = true;      
      }

      for ( i = 0; i < beta_mix->numberOfMixtures(); i++ )
      {
        updateMixtureProportion( i );
      } 
      beta_mix->normalizeLogProportions();
    }
  }//end mixture case
}//end


void BayesianLinearModel::updateMixtureProportion( int i )
{
  double tau_factor, beta0_factor, beta_factor, beta_diff_factor, det_factor, init_det_factor, prop_factor, updated_prop;

  if ( t_beta || t_lkhd_beta )
  {
    //now update proportions
    CVector vbeta( betas_t[i]->dimension() );

#ifdef FIX1
    DistributionParameter tmpdist = betas_t[i]->tau2()->lastDraw();
    tau_factor = betas_t[i]->tau2()->logDensity( tmpdist );
#else
    tau_factor = betas_t[i]->tau2()->logDensity( betas_t[i]->tau2()->lastDraw() );
#endif
//    tau_factor = betas_t[i]->tau2()->logDensity( betas_t[i]->tau2()->lastDraw() );

    beta0_factor = betas_t[i]->mu()->initialMean() * betas_t[i]->mu()->initialInvCovBeta();
    beta0_factor /= betas_t[i]->tau2()->lastItemDrawn();

#ifdef FIX1
    CVector tmpvec = betas_t[i]->mean().getVector();
    vbeta = betas_t[i]->mu()->choleskyInvCov().asUpperTriangularMultiplyRight( tmpvec );
#else
    vbeta = betas_t[i]->mu()->choleskyInvCov().asUpperTriangularMultiplyRight( betas_t[i]->mean().getVector() );
#endif
//    vbeta = betas_t[i]->mu()->choleskyInvCov().asUpperTriangularMultiplyRight( betas_t[i]->mean().getVector() );
    beta_factor = vbeta * vbeta;

    det_factor = betas_t[i]->mu()->invCovDeterminant();
    init_det_factor = betas_t[i]->mu()->initialInvDeterminant();
   
    prop_factor = beta_mix->initialProportion( i );

    //updated_prop = prop_factor * exp( -0.5 * ( beta0_factor - beta_factor ) + tau_factor );
    //updated_prop /= det_factor;
    updated_prop = log( prop_factor) - 0.5 * ( beta0_factor - beta_factor ) + tau_factor - 0.5 * log( det_factor ) + 0.5 * log( init_det_factor ) - ( ( double) betas_t[i]->dimension() / 2.0 ) * log( betas_t[i]->tau2()->lastItemDrawn() );

    beta_mix->updateLogProportion( i, updated_prop );
  }
  else
  {
    CVector vbeta( betas[i]->dimension() );

    beta0_factor = betas[i]->initialMean() * betas[i]->initialInvCovBeta();
    
#ifdef FIX1
    CVector tmpvec = betas[i]->mean().getVector();
    vbeta = betas[i]->choleskyInvCov().asUpperTriangularMultiplyRight( tmpvec );
#else
    vbeta = betas[i]->choleskyInvCov().asUpperTriangularMultiplyRight( betas[i]->mean().getVector() );
#endif
//    vbeta = betas[i]->choleskyInvCov().asUpperTriangularMultiplyRight( betas[i]->mean().getVector() );
    beta_factor = vbeta * vbeta;

    beta_diff_factor = 0.5 * ( beta0_factor - beta_factor );
    if ( conjugate_model )
    {
      beta_diff_factor /= sigma2_mix->lastDraw().getScalar();
    }

    //printf("vbeta is \n"); vbeta->Print();
    //printf("beta mean is \n"); betas[i]->mean().Print();
    

    //det_factor = sqrt( betas[i]->invCovDeterminant() );
    det_factor = betas[i]->invCovDeterminant();
    init_det_factor = betas[i]->initialInvDeterminant();

    prop_factor = beta_mix->initialProportion( i );

    //printf("Bayesian:updateMixtureProps:[%d] beta0_factor - beta_factor = %f\n", i, beta0_factor - beta_factor );
    //fflush(stdout);

    //updated_prop = prop_factor * exp( -0.5 * ( beta0_factor - beta_factor ) );
    //updated_prop /= det_factor;
    updated_prop = log( prop_factor) - beta_diff_factor - 0.5 * log( det_factor ) + 0.5 * log( init_det_factor );

    //printf("Bayesian:updateMixtureProps:[%d]  beta0_factor = %f, beta_factor = %f, \ndet_factor = %f,  prop_factor = %f, \nupdated_prop = %f\n", i, beta0_factor, beta_factor, det_factor, prop_factor, updated_prop); fflush(stdout);

    beta_mix->updateLogProportion( i, updated_prop );
  }

}//end



void BayesianLinearModel::updateSigma2MixtureProportion( int i, double s_beta, double add_df, bool include_Cov )
{
  double beta_factor, det_factor, init_det_factor, prop_factor, power_factor, updated_prop;

  beta_factor = s_beta + sigma2s[i]->initialDegreesOfFreedom() * sigma2s[i]->initialScale();
  init_det_factor = betas[i]->initialInvDeterminant();
  prop_factor = beta_mix->initialProportion( i );
  power_factor = ( add_df + sigma2s[i]->initialDegreesOfFreedom() ) / 2.0;

  updated_prop = log( prop_factor) + 0.5 * log( init_det_factor ) - power_factor * log( beta_factor );

  if ( include_Cov )
  {
    det_factor = betas[i]->invCovDeterminant();
    updated_prop -= ( 0.5 * log( det_factor ) );
  }


  sigma2_mix->updateLogProportion( i, updated_prop );

}//end




/* update parameters of conditional distribution for sigma2.
   sigma2 update is based on a scale derived from (perhaps weighted) the sum of squares of residuals.
   The residuals need to be weighted according to the precision of tau_i^2 when errors are modeled as
   t-distributed.
   actual update of parameters is done by the call to  InvChisqDistribution::update()
*/
void BayesianLinearModel::gibbsUpdateSigma2()
{
  double s_beta, add_df;

  //update
  if ( !mixture_prior )
  {
    (*residuals) = (*response) - ( (*predictors) * beta->lastItemDrawn() );
  }
  else
  {
    (*residuals) = (*response) - ( (*predictors) * beta_mix->lastDraw().getVector() );
  }

  if ( t_lkhd || t_lkhd_beta )
  {
    s_beta = (*residuals) * ( residuals->weighted( (*tau2_error_weights) ) );
    add_df = response->Len();
    sigma2->update( add_df, s_beta ); 
  }
  else if ( conjugate_model && !mixture_prior )
  {
    double s_beta_dist;
    CVector beta_diff( beta->dimension() );
    beta_diff = beta->lastItemDrawn() - beta->initialMean();

    s_beta_dist = M_D( beta->initialInverseCovariance(), beta_diff );
    s_beta  = (*residuals) * (*residuals);
    s_beta += s_beta_dist;
    add_df = response->Len() + beta->dimension();
    sigma2->update( add_df, s_beta ); 
  }
  else if ( conjugate_model && mixture_prior )
  {
    double s_beta_dist;
    int i, index;
    bool include_Cov;

    include_Cov = false;
    index = beta_mix->lastComponentDrawn();
    CVector beta_diff( betas[ 0 ]->dimension() );

    add_df = response->Len() + betas[0]->dimension();
    for ( i = 0; i < sigma2_mix->numberOfMixtures(); i++ )
    {
      beta_diff = betas[ index ]->lastItemDrawn() - betas[i]->initialMean();
      s_beta_dist = M_D( betas[i]->initialInverseCovariance(), beta_diff );
      s_beta  = (*residuals) * (*residuals);
      s_beta += s_beta_dist;
      sigma2s[i]->update( add_df, s_beta );
      updateSigma2MixtureProportion( i, s_beta, add_df, include_Cov );
    }
    sigma2_mix->normalizeLogProportions();
  }
  else 
  {
    s_beta = (*residuals) * (*residuals);
    add_df = response->Len();
    sigma2->update( add_df, s_beta ); 
  }

}//end


/* update parameters of conditional distribution for tau_i^2, i=1,...,n.
   This is needed when the error priors are modeled as t-distributions.
   tau_i^2 update is based on a scale derived from a weighted i-th residual.
   The residuals need to be weighted according to the precision of sigma2.
   actual update of parameters is done by the call to  InvChisqDistribution::update()
*/
void BayesianLinearModel::gibbsUpdateTau2Error( int i )
{
  double s_residual;

  s_residual = ( residuals->Val(i) * residuals->Val(i) ) / sigma2->lastItemDrawn();
  tau2_error[i]->update( 1, s_residual );
}//end





/* Get simul_to_keep samples from the joint posterior of (beta, sigma2) | response
   This is done only for the conjugate and non-informative normal-likelihood case
   (non-informative means non-informative prior for beta, and either InvChisq prior on sigma2
    or InvChisq( df = 0 ) on sigma2, i.e. non-informative in log(sigma2) ).

    actual sampling is done in exactSamplerStep() (see below).
*/
void BayesianLinearModel::exactSampler( int simul_to_keep )
{
  int simulations_to_keep, number_simulations_performed;

  //printf( "BayesianLinearModel:: started exact sampling for this model.\n" );

  simulations_to_keep = simul_to_keep;

  if ( simulations_to_keep > 0 )
  {
    number_simulations_performed = 0;

    //get memory for output
    createOutput( simulations_to_keep );
    //temporary storage
    initializeTemporaryStructures();
    //update beta (get correct covariance and mean)
    gibbsUpdateBeta();       
    //get the posterior parameters for sigma2
    exactSamplerPosteriorSigma2();

    while ( number_simulations_performed < simulations_to_keep )
    {
      exactSamplerStep();
      keepSimulation( number_simulations_performed );
      number_simulations_performed++;
    }//end loop
  }//end if
}//end



/* performs a draw from the joint posterior of (beta, sigma2) | response.
   See exactSampler() above for details on priors restrictions.
*/
void BayesianLinearModel::exactSamplerStep()
{
  //draw sigma2
  if ( mixture_prior && conjugate_model )
  {
    sigma2_mix->draw();
  }
  else
  {
    sigma2->draw();
  }

  //update and draw beta
  gibbsUpdateBeta();
  if ( mixture_prior )
  {
    beta_mix->draw();
    if ( conjugate_model )
    {
      betas[ beta_mix->lastComponentDrawn() ]->scaleLastDraw( sqrt( sigma2s[ sigma2_mix->lastComponentDrawn() ]->lastItemDrawn() ) );
    }
  }
  else
  {
    beta->draw();
    if ( conjugate_model )
    {
      beta->scaleLastDraw( sqrt( sigma2->lastItemDrawn() ) );
    }
  }
}//end


/* update parameters for the posterior of sigma2 | response
   The posterior of sigma2 is assumed to be an InvChisq distribution.
   See exactSampler() above for details on priors restrictions.
*/ 
void BayesianLinearModel::exactSamplerPosteriorSigma2()
{
  double scale, ss_response, ss_beta, ss_beta_mean, add_df;

  if ( !mixture_prior )
  {
    //compute additional scale: y^T y + beta0^T betaCov0^{-1} beta0 - mu_beta^T Sigma_beta^{-1} mu_beta
    CVector LTMean( beta->dimension() );
#ifdef FIX1
    CVector tmpvec = beta->mean().getVector();
    LTMean = beta->choleskyInvCov().asUpperTriangularMultiplyRight( tmpvec );
#else
    LTMean = beta->choleskyInvCov().asUpperTriangularMultiplyRight( beta->mean().getVector() );
#endif
//    LTMean = beta->choleskyInvCov().asUpperTriangularMultiplyRight( beta->mean().getVector() );

    ss_response = (*response) * (*response);
    ss_beta_mean = LTMean * LTMean;

    if ( conjugate_model )
    {
      ss_beta = beta->initialMean() * beta->initialInvCovBeta();
      ss_response += ss_beta;

      add_df = response->Len();
    }
    else //is non-informative
    {
      ss_beta_mean *= sigma2->lastItemDrawn();
      add_df = response->Len() - beta->dimension();
    }

    scale = ss_response - ss_beta_mean;

    sigma2->update( add_df, scale ); 
  }
  else if ( conjugate_model && mixture_prior )
  {
    int i;
    bool include_Cov;
    //compute additional scale: y^T y + beta0^T betaCov0^{-1} beta0 - mu_beta^T Sigma_beta^{-1} mu_beta
    CVector LTMean( betas[0]->dimension() );

    include_Cov = true;
    ss_response = (*response) * (*response);
    for ( i = 0; i < sigma2_mix->numberOfMixtures(); i++ )
    {
#ifdef FIX1
      CVector tmpvec = betas[i]->mean().getVector();
      LTMean = betas[i]->choleskyInvCov().asUpperTriangularMultiplyRight( tmpvec );
#else
      LTMean = betas[i]->choleskyInvCov().asUpperTriangularMultiplyRight( betas[i]->mean().getVector() );
#endif
//      LTMean = betas[i]->choleskyInvCov().asUpperTriangularMultiplyRight( betas[i]->mean().getVector() );
      ss_beta_mean = LTMean * LTMean;
      ss_beta = betas[i]->initialMean() * betas[i]->initialInvCovBeta();
      scale = ss_response + ss_beta - ss_beta_mean;

      add_df = response->Len();
      
      sigma2s[i]->update( add_df, scale );

      //now update proportions
      updateSigma2MixtureProportion( i, scale, add_df, include_Cov );      

    }//end for i
    sigma2_mix->normalizeLogProportions();

    //printf("Bayesian:exact sigma2: proportions are: \n" );
    //for ( i = 0; i < sigma2_mix->numberOfMixtures(); i++ )
    //{
    //  printf("prop[%d] = %f  scale = %f  nu = %f \n", i, sigma2_mix->proportion(i), sigma2s[i]->scale(), sigma2s[i]->degreesOfFreedom() );
    //}
    //printf("\n\n");

  }
}//end
