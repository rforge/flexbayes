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


void getInitialPointsBayesHPMIL( 
                     int * number_draws,
                     int * number_groups,
                     int * number_data,
                     int * dim_beta,
                     int * dim_gamma,
                     int * dim_alpha,
                     int * exposure_flag,

                     double * random_data,
                     double * fixed_data,
                     double * second_data,
                     double * data_response,
                     double * data_exposure,

                     double * betamean,
                     double * betaCov,
                     double * betaDF,
                     int * prior_beta_type,

                     double * gammamean,
                     double * gammaCov,
                     double * gammaDF,

                     double * alphamean,
                     double * alphaCov,
                     double * alphaDF,

                     double * tauDF,
                     double * tauScale,
                     int * prior_tau_type,

                     double * sigmaDF,
                     double * sigmaScale,
                     int * dim_sigmaScale,
                     int * prior_sigma_type,
                     int * common_sigma,

                     double * output_points ) throw( rtErr )
{
  bool random_effects, fixed_effects, second_effects;
  int i, j, k, dim, start_index;

  CVector ** response;
  CVector ** exposures;
  CMatrix ** random_predictors;
  CMatrix ** fixed_predictors;
  CMatrix ** second_predictors;

  /*set random seed*/
  GetRNGstate();

  random_effects = false;
  fixed_effects = false;
  second_effects = false;

  //the data  
  //the responses
  response = new CVector * [ (int) (*number_groups) ];
  start_index = 0;
  for ( i = 0; i < (int) (*number_groups); i++ )
  {
    response[i] = new CVector( (int) number_data[i] );
    for ( j = 0; j < (int) number_data[i]; j++ )
    {
      response[i]->Val(j) = data_response[ start_index + j ]; 
    }
    start_index += (int) number_data[i];

    //printf( "getInitPoints: response[%d] = \n",i );
    //response[i]->Print();
  }//end for i



  exposures = new CVector * [ (int) (*number_groups) ];
  start_index = 0;
  for ( i = 0; i < (int) (*number_groups); i++ )
  {
    exposures[i] = new CVector( (int) number_data[i] );
    if ( ((int) (*exposure_flag)) == 1 )
    {
      for ( j = 0; j < (int) number_data[i]; j++ )
      {
        exposures[i]->Val(j) = data_exposure[ start_index + j ]; 
      }
     
      start_index += (int) number_data[i];
    }
    else
    {
      exposures[i]->setTo( 1.0 );
    }

    //printf( "getInitPoints: exposures[%d] = \n",i );
    //exposures[i]->Print();

  }//end for i

  dim = (int) (*dim_beta);
  if ( dim > 0 )
  {
    random_effects = true;
    
    random_predictors = new CMatrix * [ (int) (*number_groups) ];
    start_index = 0;
    
    for ( i = 0; i < (int) (*number_groups); i++ )
    {
      random_predictors[i] = new CMatrix ( response[i]->Len(), dim );
      for ( k = 0; k < response[i]->Len(); k++ )
      {
        for ( j = 0; j < dim; j++ )
        {
          random_predictors[i]->Val( k, j ) = random_data[ start_index + k * dim + j  ]; 
        } 
      }//end for k
      start_index += dim * response[i]->Len();

      //printf( "getInitPoints: random predictors[%d] are\n", i );
      //random_predictors[i]->Print();


    }//end for i
  }

  dim = (int) (*dim_gamma);
  if ( dim > 0 )
  {
    fixed_effects = true;
    
    fixed_predictors = new CMatrix * [ (int) (*number_groups) ];
    start_index = 0;
    
    for ( i = 0; i < (int) (*number_groups); i++ )
    {
      fixed_predictors[i] = new CMatrix ( response[i]->Len(), dim );
      for ( k = 0; k < response[i]->Len(); k++ )
      {
        for ( j = 0; j < dim; j++ )
        {
          fixed_predictors[i]->Val( k, j ) = fixed_data[ start_index + k * dim + j  ]; 
        } 
      }//end for k
      start_index += dim * response[i]->Len();

      //printf( "getInitPoints: fixed predictors[%d] are\n", i );
      //fixed_predictors[i]->Print();
    }//end for i
  }

  if ( !random_effects && !fixed_effects )
  {
    printf( "getInitialPointsBhpm: No fixed effects nor random effects provided. This model is not valid.\n" );
    char the_error[] = "getInitialPointsBhpm: No fixed effects nor random effects provided. This model is not valid.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

  dim = (int) (*dim_alpha);
  if ( dim > 0 )
  {
    second_effects = true;

    second_predictors = new CMatrix * [ (int) (*number_groups) ];
    start_index = 0;

    for ( i = 0; i < (int) (*number_groups); i++ )
    {
      second_predictors[i] = new CMatrix ( ( (int) (*dim_beta) ), dim );
      for ( k = 0; k < ( (int) (*dim_beta) ); k++ )
      {
        for ( j = 0; j < dim; j++ )
        {
          second_predictors[i]->Val( k, j ) = second_data[ start_index + k * dim + j ]; 
        } 
      }//end for j
      start_index += dim * ( (int) (*dim_beta) );

      //printf( "getInitPoints: level2 predictors[%d] are\n", i );
      //second_predictors[i]->Print();
    }//end for i
  }



  bool multiple_hyper;
  int n_generate;
  double beta_df, gamma_df,  alpha_df, sigma_df, sigma_scale;
  CVector tau_df (*dim_beta), tau_scale (*dim_beta);
  
  CMatrix * gammaInit;
  CMatrix ** betaInit;
  CMatrix * alphaInit;
  CMatrix ** tauInit;
  CMatrix * sigmaInit;

  gammaInit = NULL;
  betaInit = NULL;
  alphaInit = NULL;
  tauInit = NULL;
  sigmaInit = NULL;


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

    //printf("gamma dim is %d.  df = %f \n",gamma.dimension(), gamma.degreesOfFreedom() );
    //printf("gamma mean is \n");
    //gamma.mu()->initialMean().Print();
    //printf("gamma scale is \n" );
    //gamma.scaleMatrix().Print();

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

    //printf("gamma Init = \n" );
    //gammaInit->Print();
  }//end fixed effects


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

    //printf("getInitPoints: alpha mean is \n" );
    //p_alpha0.Print();
    //printf("getInitPoints: alpha scale is \n" );
    //p_alphaCov.Print();


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

  //generate random coefficients variance (if beta prior is no non-informative)
  if ( random_effects && ((int) (*prior_beta_type) ) != 3 )
  {
    tauInit = new CMatrix * [ ((int) (*number_draws) ) ];
    if ( ( (int) (*prior_tau_type) ) != 4 )
    {
    	for( k=0; k< ((int)(*dim_beta)); k++ ){
        if ( tauScale[k] > 0 )
          tau_scale.Val(k) = tauScale[k];
        else
        	tau_scale.Val(k) = 1.0;
        if( (tauDF[k] > 0) && (tauDF[k] < 3) ) 
        	tau_df.Val(k) = tauDF[k];
        else
        	tau_df.Val(k) = 3;
      }
    }

    if ( ( (int) (*prior_tau_type) ) == 4 )
    {
    	// Inverse Wishart case
      CMatrix tScale( tauScale, ((int) (*dim_beta)), ((int) (*dim_beta) ), by_column );
      InvWishartDistribution tau2( tauDF[0], &tScale );
      for ( j = 0; j < ((int) (*number_draws) ); j++ )
      {
        tauInit[j] = new CMatrix( tau2.draw().getMatrix() );
      }
    } else {
      for ( j = 0; j < ((int) (*number_draws) ); j++ ){
      	tauInit[j] = new CMatrix( 1, ((int)(*dim_beta)) );

        // 0 ==> InvChisq
        // 1 ==> non informative
        // 2 ==> uniform shrinkage
        // 3 ==> du Mouchel
        if ( ( (int) (*prior_tau_type) ) == 0 || ( (int) (*prior_tau_type) ) == 1 )
        {
          InvChisqDistribution tau2;
          for( k=0; k< ((int)(*dim_beta)); k++ ){
            tau2.setDegreesOfFreedom( tau_df.Val(k) );
            tau2.setScale( tau_scale.Val(k) );
            tauInit[j]->Val( 0, k ) = tau2.draw().getScalar();
          }
        } else if ( ( (int) (*prior_tau_type) ) == 2 ) {
          UniformShrinkageDistribution tau2;
          for( k=0; k< ((int)(*dim_beta)); k++ ){
            tau2.setParameter( tau_scale.Val(k) );
            tauInit[j]->Val( 0, k ) = tau2.draw().getScalar();
          }
        } else if ( ( (int) (*prior_tau_type) ) == 3 ) {
          DuMouchelDistribution tau2;
          for( k=0; k< ((int)(*dim_beta)); k++ ){
            tau2.setParameter( tau_scale.Val(k) );
            tauInit[j]->Val( 0, k ) = tau2.draw().getScalar();
          }
        }
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
      if ( ( (int) (*prior_beta_type) ) == 1 )
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
        if ( ( (int) (*prior_beta_type) ) == 0 )
        {
          CVector bMean( (*second_predictors[i]) * alphaInit->getColumn( j ) );
          beta.setLocation( &bMean );
        }

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

  }//end if random effects
  else if ( random_effects )
  {
    //non informative case
    //generate beta
    beta_df = 3;
    betaInit = new CMatrix * [ ((int) (*number_draws) ) ];
    CMatrix p_betaCov( ((int) (*dim_beta) ), ((int) (*dim_beta) ) );
    CVector p_beta0( ((int) (*dim_beta) ) );

    p_beta0.setToZero();
    p_betaCov.setDiagonal( 1.0 );

    StudentTDistribution beta( ((int) (*dim_beta) ), beta_df );
    beta.setLocation( &p_beta0 );
    beta.setScaleMatrix( &p_betaCov );

    for ( j = 0; j < ((int) (*number_draws) ); j++ )
    {
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
    }
  }

  //generate glm parameters
  //generate error variance
  if ( (*sigmaDF) > 0 && (*sigmaDF) < 3 )
  {
    sigma_df = (*sigmaDF);
  }
  else
  {
    sigma_df = 3;
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
      printf( "getInitialPointsBayesHPMIL: Wrong dimensions in error variance hyperparameter array.\n" );
      char the_error[] = "getInitialPointsBayesHPMIL: Wrong dimensions in error variance hyperparameter array.";
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
    multiple_hyper = false ;
  }   

  sigmaInit = new CMatrix( n_generate, ((int) (*number_draws) ) );
  for ( i = 0; i < n_generate; i++ )
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

    //for prior sigma type
    // 0 ==> InvChisq
    // 1 ==> non informative
    // 2 ==> uniform shrinkage
    // 3 ==> du Mouchel
    if ( ( (int) (*prior_sigma_type) ) == 0 || ( (int) (*prior_sigma_type) ) == 1 )
    {
      InvChisqDistribution sigma2( sigma_df );
      sigma2.setScale( sigma_scale );
      for ( j = 0; j < ((int) (*number_draws) ); j++ )
      {
        sigmaInit->Val( i, j ) = sigma2.drawOneItem();
      }
    }
    else if ( ( (int) (*prior_sigma_type) ) == 2 )
    {
      UniformShrinkageDistribution sigma2( sigma_scale );
      for ( j = 0; j < ((int) (*number_draws) ); j++ )
      {
        sigmaInit->Val( i, j ) = sigma2.draw().getScalar();
      }
    }        
    else if ( ( (int) (*prior_sigma_type) ) == 3 )
    {
      DuMouchelDistribution sigma2( sigma_scale );
      for ( j = 0; j < ((int) (*number_draws) ); j++ )
      {
        sigmaInit->Val( i, j ) = sigma2.draw().getScalar();
      }
    }
    #ifdef DEBUG_SIGMA
      printf( "getInitialPointsBHPMIL: sigmaInit = \n" );  
      sigmaInit->Print();  fflush(stdout);
    #endif        
  }//end for i


  //now save results
  int index;

  //printf("getInitPoints: will save results\n"); fflush(stdout);

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


    if ( random_effects && ((int) (*prior_beta_type)) != 3 )
    {
      for ( k = 0; k < tauInit[j]->Col(); k++ )
      {
        for ( i = 0; i < tauInit[j]->Row(); i++ )
	{
          output_points[ index ] = tauInit[j]->Val( i, k );
          index++;
        }
      }//end for k

      delete tauInit[j];
    }

    //glm parameters
    if ( ( (int) (*common_sigma) ) == 0 )      
    {
      for ( i = 0; i < ((int) (*number_groups) ); i++ )
      {
        output_points[ index + i ] = sigmaInit->Val( i, j );
      }
      index += ((int) (*number_groups) );
    }
    else
    {
      output_points[ index ] = sigmaInit->Val( 0, j );
      index++;
    }

  }//end for j    


  delete sigmaInit;
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
    delete exposures[i];
    if ( random_effects )
    {
      delete random_predictors[i];
    }

    if ( fixed_effects )
    {
      delete fixed_predictors[i];
    }

    if ( second_effects )
    {
      delete second_predictors[i];
    }
  }

  delete [] response;
  delete [] exposures;

  if ( random_effects )
  {
    delete [] random_predictors;
  }

  if ( fixed_effects )
  {
    delete [] fixed_predictors;
  }

  if ( second_effects )
  {
    delete [] second_predictors;
  }


  /* update random seed */
  PutRNGstate();

  //printf("getInitPoints: totally done\n" ); fflush(stdout);

}//end

}//end extern C
