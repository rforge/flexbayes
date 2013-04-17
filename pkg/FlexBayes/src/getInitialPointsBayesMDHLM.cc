#include "R.h"

#include "rtErr.h"
#include "Vector.h"
#include "Matrix.h"
#include "Const.h"
#include "DistributionExtended.h"
#include "DistributionMixture.h"

extern "C" {

/* creates a set of initial points for starting
   a simulation (Gibbs sampling)
*/


void getInitialPointsBayesMDHLM( 
                     long * number_draws,
                     long * number_groups,
                     long * number_data,
                     long * dim_response,
                     long * dim_beta,
                     long * dim_gamma,
                     long * dim_alpha,

                     double * data_response,

                     double * betamean,
                     double * betaCov,
                     double * betaDF,
                     long * prior_beta_type,

                     double * gammamean,
                     double * gammaCov,
                     double * gammaDF,

                     double * alphamean,
                     double * alphaCov,
                     double * alphaDF,

                     double * sigmaDF,
                     double * sigmaScale,
                     long * dim_sigmaScale,
                     long * prior_sigma_type,
                     long * common_sigma,

                     double * tauDF,
                     double * tauScale,
                     long * prior_tau_type,

                     double * output_points ) throw( rtErr )
{
  bool random_effects, fixed_effects, second_effects;
  int i, j, k, l, q, q2, dim, start_index;

  CMatrix ** response;

  /*set random seed*/
  GetRNGstate();

  random_effects = false;
  fixed_effects = false;
  second_effects = false;

  //the data  
  //the responses
  response = new CMatrix * [ (int) (*number_groups) ];
  start_index = 0;
  for ( i = 0; i < (int) (*number_groups); i++ )
  {
    response[i] = new CMatrix( (int) (*dim_response), (int) number_data[i] );
    for ( j = 0; j < (int) number_data[i]; j++ )
    {
      for ( k = 0; k < (int) (*dim_response); k++ )
      {
        response[i]->Val( k, j ) = data_response[ start_index + j * ((int) (*dim_response)) + k ]; 
      }
    }
    start_index += ((int) number_data[i]) * ((int) (*dim_response));

  }//end for i


  dim = (int) (*dim_beta);
  if ( dim > 0 )
  {
    random_effects = true;
  }

  dim = (int) (*dim_gamma);
  if ( dim > 0 )
  {
    fixed_effects = true;
  }

  if ( !random_effects && !fixed_effects )
  {
    printf( "getInitialPointsBhlm: No fixed effects nor random effects provided. This model is not valid.\n" );
    char the_error[] = "getInitialPointsBhlm: No fixed effects nor random effects provided. This model is not valid.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

  dim = (int) (*dim_alpha);
  if ( dim > 0 )
  {
    second_effects = true;
  }



  bool multiple_hyper;
  int n_generate;
  double beta_df, gamma_df,  alpha_df, sigma_df, sigma_scale, tau_df, tau_scale;
  
  CMatrix * gammaInit;
  CMatrix ** betaInit;
  CMatrix * alphaInit;
  CMatrix *** sigmaInit;
  CMatrix * sigma_cov;
  CMatrix ** tauInit;

  gammaInit = NULL;
  betaInit = NULL;
  alphaInit = NULL;
  sigmaInit = NULL;
  sigma_cov = NULL;
  tauInit = NULL;


  if ( (*gammaDF) > 0 && (*gammaDF) < 3 )
  {
    gamma_df = (*gammaDF);
  }
  else 
  {
    gamma_df = 3;
  }


  //generate fixed effects coefficients
  if ( fixed_effects )
  {
    StudentTDistribution gamma( ((int) (*dim_gamma) ), gamma_df );

    CVector p_gamma0( gammamean, ((int) (*dim_gamma) ) );
    CMatrix p_gammaCov( gammaCov, ((int) (*dim_gamma) ), ((int) (*dim_gamma) ) );

    gamma.setScaleMatrix( &p_gammaCov );
    gamma.setLocation( &p_gamma0 );

    gammaInit = new CMatrix( ((int) (*dim_gamma) ), ((int) (*number_draws) ) );
    for ( j = 0; j < ((int) (*number_draws) ); j++ )
    {
#ifdef FIX1
      CVector tmpvec = gamma.draw().getVector();
      gammaInit->setColumn( j , tmpvec );
#else
      gammaInit->setColumn( j, gamma.draw().getVector() );
#endif
//      gammaInit->setColumn( j, gamma.draw().getVector() );  
    }

  }//end fixed effects



  //generate error variance
  if ( ((int) (*prior_sigma_type)) != 4 && ((int) (*prior_sigma_type)) != 5 )
  {
    if ( (*sigmaDF) > 0 && (*sigmaDF) < 3 )
    {
      sigma_df = (*sigmaDF);
    }
    else
    {
      sigma_df = 3;
    }
  }
  else
  {
    sigma_df = ( (double) (*dim_response) );
  }


  // for common sigma
  // 0 ==> independent sigmas; different hyper prior parameters
  // 1 ==> independent sigmas; same hyper prior parameters
  // 2 ==> common sigma; same hyper prior parameters
  if ( ( (int) (*common_sigma) ) == 0 )
  {
    n_generate = ((int) (*number_groups) );

    if ( ((int) (*dim_sigmaScale) ) > 1 && ((int) (*dim_sigmaScale) ) == ((int) (*number_groups) )  )
    {
      multiple_hyper = true;
    }
    else if ( ((int) (*dim_sigmaScale) ) > 1 )
    {
      printf( "getInitialPointsBhlm: Wrong dimensions in error variance hyperparameter array.\n" );
      char the_error[] = "getInitialPointsBhlm: Wrong dimensions in error variance hyperparameter array.";
      rtErr runtime_error( the_error );
      throw runtime_error;
    }
    else
    {
      multiple_hyper = false;
    }
  }
  else
  {
    n_generate = 1;
    multiple_hyper = false;
  }   


  sigmaInit = new CMatrix ** [ ((int) (*number_draws) ) ];
  for ( j = 0; j < ((int) (*number_draws) ); j++ )
  {
    sigmaInit[j] = new CMatrix * [ n_generate ];
  }


  if ( ( (int) (*prior_sigma_type) ) == 4 || ( (int) (*prior_sigma_type) ) == 5 )
  {
    sigma_cov = new CMatrix( ((int) (*dim_response)), ((int) (*dim_response)) );
  }


  start_index = 0;
  for ( i = 0; i < n_generate; i++ )
  {
    if ( ( (int) (*prior_sigma_type) ) != 4 && ( (int) (*prior_sigma_type) ) != 5 )
    {
      if ( multiple_hyper )
      {
        if ( sigmaScale[i] > 0 )
        {
          sigma_scale = sigmaScale[i];
        }
        else
        {
          sigma_scale = 1.0;
        }
      }
      else if ( i == 0 )
      {
        if ( (*sigmaScale) > 0 )
        {
          sigma_scale = (*sigmaScale);
        }
        else
        {
          sigma_scale = 1.0;
        }
      }
    }
    else
    {
      if ( multiple_hyper )
      {
        for ( j = 0; j < ((int) (*dim_response)); j++ )
        {
          for ( k = 0; k < ((int) (*dim_response)); k++ )
	  {
            sigma_cov->Val( j, k ) = sigmaScale[ start_index + j * ((int) (*dim_response)) + k ];
          }
        }
        start_index += ((int) (*dim_response)) * ((int) (*dim_response));
      }
      else if ( i == 0 )
      {
        for ( j = 0; j < ((int) (*dim_response)); j++ )
        {
          for ( k = 0; k < ((int) (*dim_response)); k++ )
	  {
            sigma_cov->Val( j, k ) = sigmaScale[ j * ((int) (*dim_response)) + k ];
          }
        }
      }
    }


    //for prior sigma type
    // 0 ==> InvChisq
    // 1 ==> non informative
    // 2 ==> uniform shrinkage
    // 3 ==> du Mouche
    // 4 ==> invWishart
    // 5 ==> mass Point
    if ( ( (int) (*prior_sigma_type) ) == 0 || ( (int) (*prior_sigma_type) ) == 1 )
    {
      InvChisqDistribution sigma2( sigma_df );
      sigma2.setScale( sigma_scale );
      for ( j = 0; j < ((int) (*number_draws) ); j++ )
      {
        sigmaInit[j][i] = new CMatrix( sigma2.draw().getMatrix() );
      }
    }
    else if ( ( (int) (*prior_sigma_type) ) == 2 )
    {
      UniformShrinkageDistribution sigma2( sigma_scale );
      for ( j = 0; j < ((int) (*number_draws) ); j++ )
      {
        sigmaInit[j][i] = new CMatrix( sigma2.draw().getMatrix() );
      }
    }        
    else if ( ( (int) (*prior_sigma_type) ) == 3 )
    {
      DuMouchelDistribution sigma2( sigma_scale );
      for ( j = 0; j < ((int) (*number_draws) ); j++ )
      {
        sigmaInit[j][i] = new CMatrix( sigma2.draw().getMatrix() );
      }
    }
    else if ( ( (int) (*prior_sigma_type) ) == 4 )
    {        

      InvWishartDistribution sigma2( sigma_df );
      sigma2.setScaleMatrix( sigma_cov );

      for ( j = 0; j < ((int) (*number_draws) ); j++ )
      {
        sigmaInit[j][i] = new CMatrix( sigma2.draw().getMatrix() );
      }      
    }
    else if ( ( (int) (*prior_sigma_type) ) == 5 )
    {        
      for ( j = 0; j < ((int) (*number_draws) ); j++ )
      {
        sigmaInit[j][i] = new CMatrix( (*sigma_cov) );
      }      
    }
  }//end for i


  //generate alpha
  if ( random_effects && second_effects )
  {
    if ( (*alphaDF) > 0 && (*alphaDF) < 3 )
    {
      alpha_df = (*alphaDF);
    }
    else
    {
      alpha_df = 3;
    }

    CVector p_alpha0( alphamean, ((int) (*dim_alpha) ) );
    CMatrix p_alphaCov( alphaCov, ((int) (*dim_alpha) ), ((int) (*dim_alpha) ) );


    StudentTDistribution alpha( ((int) (*dim_alpha) ), alpha_df );
    alpha.setScaleMatrix( &p_alphaCov );
    alpha.setLocation( &p_alpha0 );
    alphaInit = new CMatrix( ((int) (*dim_alpha) ), ((int) (*number_draws) ) );
    for ( j = 0; j < ((int) (*number_draws) ); j++ )
    {
#ifdef FIX1
      CVector tmpvec = alpha.draw().getVector();
      alphaInit->setColumn( j , tmpvec );
#else
      alphaInit->setColumn( j, alpha.draw().getVector() );
#endif
//      alphaInit->setColumn( j, alpha.draw().getVector() );
    }
  }//end level2


  bool by_column = false;//indicate not to read by row

  //generate random coefficients variance
  if ( random_effects && ((int) (*prior_beta_type)) != 3 )
  {
    tauInit = new CMatrix * [ ((int) (*number_draws) ) ];
    if ( ( (int) (*prior_tau_type) ) != 4 )
    {
      if ( (*tauScale) > 0 )
      {
        tau_scale = (*tauScale);
      }
      else
      {
        tau_scale = 1.0;
      }
    }


    if ( ((int) (*prior_tau_type)) != 4 )
    {
      if ( (*tauDF) > 0 && (*tauDF) < 3 )
      {
        tau_df = (*tauDF);
      }
      else
      {
        tau_df = 3;
      }
    }
    else
    {
      tau_df = ((double) (*dim_beta));
    }



    // 0 ==> InvChisq
    // 1 ==> non informative
    // 2 ==> uniform shrinkage
    // 3 ==> du Mouchel
    // 4 ==> InvWishart
    if ( ( (int) (*prior_tau_type) ) == 0 || ( (int) (*prior_tau_type) ) == 1 )
    {
      InvChisqDistribution tau2( tau_df );
      tau2.setScale( tau_scale );
      for ( j = 0; j < ((int) (*number_draws) ); j++ )
      {
        tauInit[j] = new CMatrix( tau2.draw().getMatrix() );
      }
    }
    else if ( ( (int) (*prior_tau_type) ) == 2 )
    {
      UniformShrinkageDistribution tau2( tau_scale );
      for ( j = 0; j < ((int) (*number_draws) ); j++ )
      {
        tauInit[j] = new CMatrix( tau2.draw().getMatrix() );
      }
    }        
    else if ( ( (int) (*prior_tau_type) ) == 3 )
    {
      DuMouchelDistribution tau2( tau_scale );
      for ( j = 0; j < ((int) (*number_draws) ); j++ )
      {
        tauInit[j] = new CMatrix( tau2.draw().getMatrix() );
      }
    }   
    else if ( ( (int) (*prior_tau_type) ) == 4 )
    {
      CMatrix tScale( tauScale, ((int) (*dim_beta)), ((int) (*dim_beta) ), by_column );
      InvWishartDistribution tau2( tau_df, &tScale );
      for ( j = 0; j < ((int) (*number_draws) ); j++ )
      {
        tauInit[j] = new CMatrix( tau2.draw().getMatrix() );
      }
    }

    //generate beta
    betaInit = new CMatrix * [ ((int) (*number_draws) ) ];
    CMatrix p_betaCov( betaCov, ((int) (*dim_beta) ), ((int) (*dim_beta) ), by_column );

    CVector * p_beta0;
    CMatrix * p_beta0_m;

    p_beta0 = NULL;
    p_beta0_m = NULL;

    if ( (*betaDF) > 0 && (*betaDF) < 3 )
    {
      beta_df = (*betaDF);
    }
    else 
    {
      beta_df = 3;
    }


    if ( ( (int) (*prior_beta_type) ) == 1 )
    {
      p_beta0 = new CVector( betamean, ((int) (*dim_beta) ) );
    }
    else if ( ( (int) (*prior_beta_type) ) == 2 )
    {
      p_beta0_m = new CMatrix( betamean, ((int) (*dim_beta) ), ((int) (*number_draws) ), by_column );
    }
       
    for ( j = 0; j < ((int) (*number_draws) ); j++ )
    {
      StudentTDistribution beta( ((int) (*dim_beta) ), beta_df );

      if ( ( (int) (*prior_tau_type) ) == 4 )
      {
        CMatrix W( (*tauInit[j]) );
        beta.setScaleMatrix( &W );
      }
      else
      {
        CMatrix W( ((int) (*dim_beta) ), ((int) (*dim_beta) ) );
        W = tauInit[j]->Val( 0, 0 ) * p_betaCov;
        beta.setScaleMatrix( &W );
      }

      // 0 ==> usual full model
      // 1 ==> no alpha parameter in model. A single beta0
      // 2 ==> no alpha parameter in model. beta0 for each group
      if ( ( (int) (*prior_beta_type) ) == 0 )
      {
        CVector bMean( alphaInit->getColumn( j ) );
        beta.setLocation( &bMean );

      }
      else if ( ( (int) (*prior_beta_type) ) == 1 )
      {
        beta.setLocation( p_beta0 );
      }
      else if ( ( (int) (*prior_beta_type) ) == 2 )
      {
        CVector bMean = p_beta0_m->getColumn( j );
        beta.setLocation( &bMean );
      }

      betaInit[j] = new CMatrix( ((int) (*dim_beta) ), ((int) (*number_groups) ) );
      for ( i = 0; i < ((int) (*number_groups) ); i++ )
      {
#ifdef FIX1
        CVector tmpvec = beta.draw().getVector();
        betaInit[j]->setColumn( i , tmpvec );
#else
        betaInit[j]->setColumn( i, beta.draw().getVector() );
#endif
//        betaInit[j]->setColumn( i, beta.draw().getVector() );
      }
    }//end for j


    if ( p_beta0 != NULL )
    {
      delete p_beta0;
    }

    if ( p_beta0_m != NULL )
    {
      delete p_beta0_m;
    }
  }//end if random effects if random effects (informative prior)
  else if ( random_effects )
  {
    //generate beta
    betaInit = new CMatrix * [ ((int) (*number_draws) ) ];

    CVector p_betaMean( ((int) (*dim_beta)) );
    CMatrix p_betaCov( ((int) (*dim_beta)), ((int) (*dim_beta)) );

    p_betaMean.setToZero();
    p_betaCov.setDiagonal( 1000.0 );

    if ( (*betaDF) > 0 && (*betaDF) < 3 )
    {
      beta_df = (*betaDF);
    }
    else 
    {
      beta_df = 3;
    }

    StudentTDistribution beta( ((int) (*dim_beta)), beta_df );
    beta.setScaleMatrix( &p_betaCov );
    beta.setLocation( &p_betaMean );
      
    for ( j = 0; j < ((int) (*number_draws) ); j++ )
    {
      betaInit[j] = new CMatrix( ((int) (*dim_beta)), ((int) (*number_groups) ) );
      for ( i = 0; i < ((int) (*number_groups) ); i++ )
      {
#ifdef FIX1
        CVector tmpvec = beta.draw().getVector();
        betaInit[j]->setColumn( i , tmpvec );
#else
        betaInit[j]->setColumn( i, beta.draw().getVector() );
#endif
//        betaInit[j]->setColumn( i, beta.draw().getVector() );
      }
    }//end for j
  }

  //now save results
  int index;

  index = 0;
  for ( j = 0; j < ((int) (*number_draws) ); j++ )
  {
    if ( random_effects )
    {
      for ( i = 0; i < ((int) (*number_groups) ); i++ )
      {
        //save vector for group i draw j
        for ( k = 0; k < ((int) (*dim_beta) ); k++ )
	{
          output_points[ index + k ] = betaInit[j]->Val( k, i );
        } 
        index += ((int) (*dim_beta) );
      }//end for i
      delete betaInit[j];
    }

    if ( fixed_effects )
    {
      for ( k = 0; k < ((int) (*dim_gamma) ); k++ )
      {
        output_points[ index + k ] = gammaInit->Val( k, j );
      } 
      index += ((int) (*dim_gamma) );
    }

    if ( random_effects && second_effects )
    {
      for ( k = 0; k < ((int) (*dim_alpha) ); k++ )
      {
        output_points[ index + k ] = alphaInit->Val( k, j );
      }
      index += ((int) (*dim_alpha) );
    }


    q = sigmaInit[j][0]->Row();
    q2 = sigmaInit[j][0]->Row() * sigmaInit[j][0]->Col();
    for ( i = 0; i < n_generate; i++ )
    {
      for ( k = 0; k < sigmaInit[j][i]->Col(); k++ )
      {
        for ( l = 0; l < sigmaInit[j][i]->Row(); l++ )
	{            
          output_points[ index + i * q2 + k * q + l ] = sigmaInit[j][i]->Val( l, k );
        }
      }
      delete sigmaInit[j][i];
    }//end for i
    index += ( n_generate * q2 );

    delete [] sigmaInit[j];


    if ( random_effects && ((int) (*prior_beta_type)) != 3 )
    {
      q = tauInit[j]->Row();
      q2 = tauInit[j]->Row() * tauInit[j]->Col();
      for ( k = 0; k < tauInit[j]->Col(); k++ )
      {
        for ( i = 0; i < tauInit[j]->Row(); i++ )
	{
          output_points[ index + k * q + i ] = tauInit[j]->Val( i, k );
        }
      }//end for k
      index += q2;

      delete tauInit[j];
    }
  }//end for j    


  delete [] sigmaInit;
  if ( fixed_effects )
  {
    delete gammaInit;
  }
  if ( random_effects )
  {
    delete [] betaInit;

    if ( ((int) (*prior_beta_type)) != 3 )
    {
      delete [] tauInit;
    }

    if ( second_effects )
    {
      delete alphaInit;
    }
  }

  for ( i = 0; i < ((int) (*number_groups) ); i++ )
  {
    delete response[i];
  }
  delete [] response;

  /* update random seed */
  PutRNGstate();

  
}//end

}//end extern C
