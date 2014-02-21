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


void getInitialPointsBayesMVHLM( 
                     int * number_draws,
                     int * number_groups,
                     int * number_data,
                     int * dim_response,
                     int * dim_beta,
                     int * dim_gamma,
                     int * dim_alpha,

                     double * data_random,
                     double * data_fixed,
                     double * data_level2,
                     double * data_response,

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

                     double * sigmaDF,
                     double * sigmaScale,
                     int * dim_sigmaScale,
                     int * prior_sigma_type,
                     int * common_sigma,

                     double * tauDF,
                     double * tauScale,
                     int * prior_tau_type,

                     double * output_points ) throw( rtErr )
{
  bool random_effects, fixed_effects, second_effects;
  int i, j, k, l, j1, j2, q, q2, dim, beta_rxp, gamma_rxp, alpha_rxp, number_obs, start_index;

  CMatrix ** response;
  CMatrix *** random_predictors;
  CMatrix *** fixed_predictors;
  CMatrix ** second_predictors;

  /*set random seed*/
  GetRNGstate();

  random_effects = false;
  fixed_effects = false;
  second_effects = false;

  //the data  
  //the responses
  response = new CMatrix * [ (int) (*number_groups) ];
  start_index = 0;
  number_obs = 0;
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
    number_obs += ((int) number_data[i]);

    //printf( "fit: response[%d] = \n", i );
    //response[i]->Print();

  }//end for i


  dim = (int) (*dim_beta);
  if ( dim > 0 )
  {
    random_effects = true;
    start_index = 0;
    beta_rxp = ((int) (*dim_response)) * dim;    
    random_predictors = new CMatrix ** [ ((int) (*number_groups)) ];
    for ( i = 0; i < (int) (*number_groups); i++ )
    {
      random_predictors[i] = new CMatrix * [ (int) number_data[i] ];
      for ( j = 0; j < ((int) number_data[i]); j++ )
      {
        random_predictors[i][j] = new CMatrix ( ((int) (*dim_response)), beta_rxp );
        random_predictors[i][j]->setToZero();
        for ( k = 0; k < (int) (*dim_response); k++ )
        {  
          for ( l = 0; l < dim; l++ )
	  {
            random_predictors[i][j]->Val( k, k * dim + l ) = data_random[ start_index + l ];
          }
        } 
        start_index += dim;

        //printf( "fit: random predictors[%d][%d] = \n", i, j );
        //random_predictors[i][j]->Print();
      }//end for j
    }//end for i
  }//end if random effects


  dim = (int) (*dim_gamma);
  if ( dim > 0 )
  {
    fixed_effects = true;

    start_index = 0;
    gamma_rxp = ((int) (*dim_response)) * dim;    
    fixed_predictors = new CMatrix ** [ ((int) (*number_groups)) ];
    for ( i = 0; i < (int) (*number_groups); i++ )
    {
      fixed_predictors[i] = new CMatrix * [ (int) number_data[i] ];
      for ( j = 0; j < ((int) number_data[i]); j++ )
      {
        fixed_predictors[i][j] = new CMatrix ( ((int) (*dim_response)), gamma_rxp );
        fixed_predictors[i][j]->setToZero();
        for ( k = 0; k < (int) (*dim_response); k++ )
        {  
          for ( l = 0; l < dim; l++ )
	  {
            fixed_predictors[i][j]->Val( k, k * dim + l ) = data_fixed[ start_index + l ];
          }
        } 
        start_index += dim;

        //printf( "fit: fixed predictors[%d][%d] = \n", i, j );
        //fixed_predictors[i][j]->Print();
      }//end for j
    }//end for i
  }//end if fixed effects


  if ( !random_effects && !fixed_effects )
  {
    Rprintf( "getInitialPointsBhlm: No fixed effects nor random effects provided. This model is not valid.\n" );
    char the_error[] = "getInitialPointsBhlm: No fixed effects nor random effects provided. This model is not valid.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

  dim = (int) (*dim_alpha);
  if ( dim > 0 )
  {
    second_effects = true;

    start_index = 0;
    alpha_rxp = ((int) (*dim_response)) * dim;

    second_predictors = new CMatrix * [ ((int) (*number_groups)) ];
    for ( i = 0; i < (int) (*number_groups); i++ )
    {
      second_predictors[i] = new CMatrix ( beta_rxp, alpha_rxp );
      second_predictors[i]->setToZero();
      for ( j = 0; j < (int) (*dim_response); j++ )
      {
        j1 = j * ((int) (*dim_beta));
        j2 = j * dim;

        for ( k = 0; k < (int) (*dim_beta); k++ )
        {  
          for ( l = 0; l < dim; l++ )
	  {
            second_predictors[i]->Val( j1 + k, j2 + l ) = data_level2[ start_index + k * dim + l ];
          }
        } 
      }
      start_index += dim * ((int) (*dim_beta));

      //printf( "fit: second predictors[%d] = \n", i );
      //second_predictors[i]->Print();
    }//end for i
  }//end if second level



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
    StudentTDistribution gamma( gamma_rxp, gamma_df );

    CVector p_gamma0( gammamean, gamma_rxp );
    CMatrix p_gammaCov( gammaCov, gamma_rxp, gamma_rxp );

    gamma.setScaleMatrix( &p_gammaCov );
    gamma.setLocation( &p_gamma0 );

    //printf("gamma dim is %d.  df = %f \n",gamma.dimension(), gamma.degreesOfFreedom() );
    //printf("gamma mean is \n");
    //gamma.mu()->initialMean().Print();
    //printf("gamma scale is \n" );
    //gamma.scaleMatrix().Print();

    gammaInit = new CMatrix( gamma_rxp, ((int) (*number_draws) ) );
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


  //printf("**** getpoints: common sigma is %d, and dim sigma scale is %d\n", ((int)(*common_sigma)), ((int)(*dim_sigmaScale)) );
 

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
      Rprintf( "getInitialPointsBhlm: Wrong dimensions in error variance hyperparameter array.\n" );
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


  //printf(" gen: sigma type is %d.  n generate is %d\n", ( (int) (*prior_sigma_type) ), n_generate );

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
      //printf( "gen: sigma type 4or5: multiple_hyper is %d\n", multiple_hyper );

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
        //printf( "gen: will init sigma_cov\n" );

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

      //printf( "gen: sigma df is %d,  sigma cov is \n", sigma_df );
      //sigma_cov->Print();
      //fflush(stdout);

      InvWishartDistribution sigma2( sigma_df );
      sigma2.setScaleMatrix( sigma_cov );

      //printf( "gen: ready to draw sigma\n"); fflush(stdout);
      for ( j = 0; j < ((int) (*number_draws) ); j++ )
      {
        sigmaInit[j][i] = new CMatrix( sigma2.draw().getMatrix() );

        //sigmaInit[j][i]->Print(); fflush(stdout);
      }      
    }
    else if ( ( (int) (*prior_sigma_type) ) == 5 )
    {        
      //printf("sigma cov[%d] = \n", i ); fflush(stdout);
      //sigma_cov->Print();
   
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

    CVector p_alpha0( alphamean, alpha_rxp );
    CMatrix p_alphaCov( alphaCov, alpha_rxp, alpha_rxp );

    //printf("getInitPoints: alpha mean is \n" );
    //p_alpha0.Print();
    //printf("getInitPoints: alpha scale is \n" );
    //p_alphaCov.Print();


    StudentTDistribution alpha( alpha_rxp, alpha_df );
    alpha.setScaleMatrix( &p_alphaCov );
    alpha.setLocation( &p_alpha0 );
    alphaInit = new CMatrix( alpha_rxp, ((int) (*number_draws) ) );
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
      tau_df = ((double) beta_rxp);
    }


    //printf( "gen: prior tau type is %d. dim beta is %d\n", ((int) (*prior_tau_type)), beta_rxp );

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
      CMatrix tScale( tauScale, beta_rxp, beta_rxp, by_column );
      InvWishartDistribution tau2( tau_df, &tScale );
      for ( j = 0; j < ((int) (*number_draws) ); j++ )
      {
        tauInit[j] = new CMatrix( tau2.draw().getMatrix() );
      }
    }

    //generate beta
    betaInit = new CMatrix * [ ((int) (*number_draws) ) ];
    CMatrix p_betaCov( betaCov, beta_rxp, beta_rxp, by_column );

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
      p_beta0 = new CVector( betamean, beta_rxp );
    }
    else if ( ( (int) (*prior_beta_type) ) == 2 )
    {
      p_beta0_m = new CMatrix( betamean, beta_rxp, ((int) (*number_draws) ), by_column );
    }
       
    for ( j = 0; j < ((int) (*number_draws) ); j++ )
    {
      StudentTDistribution beta( beta_rxp, beta_df );

      if ( ( (int) (*prior_tau_type) ) == 4 )
      {
        CMatrix W( (*tauInit[j]) );
        beta.setScaleMatrix( &W );
      }
      else
      {
        CMatrix W( beta_rxp, beta_rxp );
        W = tauInit[j]->Val( 0, 0 ) * p_betaCov;
        beta.setScaleMatrix( &W );
      }

      // 0 ==> usual full model (see below)
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

      betaInit[j] = new CMatrix( beta_rxp, ((int) (*number_groups) ) );


      for ( i = 0; i < ((int) (*number_groups) ); i++ )
      {
        if ( ( (int) (*prior_beta_type) ) == 0 )
        {
          //printf( "getInitPoints: full model \n" );

          //printf("second predictors are \n" );
          //second_predictors[i]->Print();
          //printf("alpha Init[%d] is \n", j );
          //alphaInit->getColumn( j ).Print();
          //fflush(stdout);

          CVector bMean( (*second_predictors[i]) * alphaInit->getColumn( j ) );
          beta.setLocation( &bMean );

          //printf( "getInitPoints: beta mean is \n" );
          //bMean.Print();
        }

#ifdef FIX1
        CVector tmpvec = beta.draw().getVector();
        betaInit[j]->setColumn( i , tmpvec );
#else
        betaInit[j]->setColumn( i, beta.draw().getVector() );
#endif
//        betaInit[j]->setColumn( i, beta.draw().getVector() );
      }//end for i
    }//end for j


    //for ( j = 0; j < ((int) (*number_draws) ); j++ )
    //{
    //  printf("beta[%d] = \n", j );
    //  betaInit[j]->Print();
    //}


    if ( p_beta0 != NULL )
    {
      delete p_beta0;
    }

    if ( p_beta0_m != NULL )
    {
      delete p_beta0_m;
    }
  }//end if random effects (informative prior)
  else if ( random_effects )
  {
    //generate beta
    betaInit = new CMatrix * [ ((int) (*number_draws) ) ];

    CVector p_betaMean( beta_rxp );
    CMatrix p_betaCov( beta_rxp, beta_rxp );

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

    StudentTDistribution beta( beta_rxp, beta_df );
    beta.setScaleMatrix( &p_betaCov );
    beta.setLocation( &p_betaMean );
      
    for ( j = 0; j < ((int) (*number_draws) ); j++ )
    {
      betaInit[j] = new CMatrix( beta_rxp, ((int) (*number_groups) ) );
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

  //printf("getInitPoints: will save results\n"); fflush(stdout);

  index = 0;
  for ( j = 0; j < ((int) (*number_draws) ); j++ )
  {
    if ( random_effects )
    {
      for ( i = 0; i < ((int) (*number_groups) ); i++ )
      {
        //save vector for group i draw j
        for ( k = 0; k < beta_rxp; k++ )
	{
          output_points[ index + k ] = betaInit[j]->Val( k, i );
        } 
        index += beta_rxp;
      }//end for i
      delete betaInit[j];
    }

    if ( fixed_effects )
    {
      for ( k = 0; k < gamma_rxp; k++ )
      {
        output_points[ index + k ] = gammaInit->Val( k, j );
      } 
      index += gamma_rxp;
    }

    if ( random_effects && second_effects )
    {
      for ( k = 0; k < alpha_rxp; k++ )
      {
        output_points[ index + k ] = alphaInit->Val( k, j );
      }
      index += alpha_rxp;
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


  if ( random_effects )
  {
    for ( i = 0; i < ((int) (*number_groups) ); i++ )
    {
      for ( j = 0; j < ((int) number_data[i]); j++ )
      {
        delete random_predictors[i][j];
      }
      delete [] random_predictors[i];
    }
    delete [] random_predictors;
  }

  if ( fixed_effects )
  {
    for ( i = 0; i < ((int) (*number_groups) ); i++ )
    {
      for ( j = 0; j < ((int) number_data[i]); j++ )
      {
        delete fixed_predictors[i][j];
      }
      delete [] fixed_predictors[i];
    }
    delete [] fixed_predictors;
  }

  if ( second_effects )
  {
    for ( i = 0; i < ((int) (*number_groups) ); i++ )
    {
      delete second_predictors[i];
    }
    delete [] second_predictors;
  }


  /* update random seed */
  PutRNGstate();

  //printf("getInitPoints: totally done\n" ); fflush(stdout);

}//end

}//end extern C
