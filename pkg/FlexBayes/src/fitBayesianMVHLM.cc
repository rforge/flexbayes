#include "R.h"
#include "rtErr.h"
#include "Vector.h"
#include "Matrix.h"

#include "BayesianMVHLM.h"
#include "Gibbs.h"

//#include <stdio.h>

extern "C" {

void  transformResponseAndPredictorsMV( int * number_groups, 
                                        int * number_data,
                                        bool random_effects,
                                        bool fixed_effects,
                                        CMatrix ** response, 
                                        CMatrix *** random_predictors, 
                                        CMatrix *** fixed_predictors,   
                                        int * unique_error_Cov, 
                                        int * dim_error_Cov, 
		 			double * error_Cov ) throw( rtErr );

/*
  Computes a Bayes Hierarchical linear model:

  |  y_j = X_j beta_j + M_j gamma + error_j
  |
  |  beta_j = Z_j alpha + beta_error :the random effects
  |
  |  gamma ~ flat prior == N( 0, +infty ) :the fixed effects
  |        or ~ Normal( gamma_0, Gamma_0 )
  |
  |  alpha ~ N( alpha_0, V_0 )
  |  beta_error ~ N( 0, tau2 V )
  |  error_j ~ Normal( 0, sigma2 Sigma_j )  ||   error_j ~ Normal( 0, sigma2_j Sigma_j )  
  |
  |  sigma2 ~ p(sigma2) || sigma2_j ~ p(sigma2_j)
  |  tau2 ~ p(tau2) or tau2V ~ InvWishart( nu, S )

    y_j: data response for the j-th group
    X_j: data predictors for the random effects of the j-th group
    M_j: data predictors for the fixed effects
    Z_j: data predictors for the "population" or secodn level predictors

    j = 1, ..., J.

  p(sigma2): could be InvChisq( nuSigma0, sigma02 )
                     non Informative: proper:   Uniform Shrinkage ( sigma02 )
                                                Du Mouchel( sigma02 )
                                      improper: proportional to sigma2^(1 + a ), with a <= 0.

  p(tau2): could be InvChisq( nuTau0, tau02 )
                     non Informative: proper:   Uniform Shrinkage ( tau02 )
                                                Du Mouchel( tau02 )
                                      improper: proportional to tau2^(1 + b ), with b <= 0.


   Robustness: distributions for gamma, beta and alpha can be made to follow Multivariate Student's t's.
               errors could be distributed as Student's t's at the observation level or at the
               group level.

   Computations are based on Gibbs sampler drawings from the full conditionals.
   burnInLength: number of initial Gibbs sampler iterations to discard.
   simulationsToPerform: number of Gibbs sampler drawings to keep as output.
   sampleFrequency: only iterations that are multiples of sampleFrequency are to be kept.


   OUTPUT
   ------
   The simulationsToPerform output simulations are returned in the array output_simulations.
   This should be an array of length 
     ( J * [p] + [r]  + q + < 1 | J > + < 1 | q*q > + [J] + +[1] + [1] + [J] + [n] ) * simulationsToPerform
   where 
   J: number of groups
   p: dimension of beta vector
   r: dimension of gamma vector
   q: dimension of alpha vector
   < 1 | J >: sigma2 (either a common sigma2 or J independent sigma2's)
   < 1 | q*q >: dimension of tau2 (either a random variable or a random matrix )
   J: augmented variables for t-distributed betas
   1: augmented variable for t-distributed gamma
   1: augmented variable for t-distributed alpha
   J: augmented variables for t-distributed group errors
   n: number of observations. augmented variables for t-distributed errors.

   < .|. > denotes a choice.
   [.] denotes an optional argument.

   The above dimensions specification also denotes the order in which the simulations are returned:
         beta (by group), gamma, alpha, sigma, tau | tau2, 
               t_beta (by group), t_gamma, t_alpha, t_groups, t_errors.
                                                       k = 0,..., dimCov - 1.

   Note: Simulations for sigma (not sigma2) are returned.
         Simulations for tau (not tau2) are returned if tau2 is a scalar,
         otherwise, the matrices tau2 are returned.

*/
void fitBayesianMVHLM( int * number_groups,
                     int * number_data,
                     int * dim_response,
                     int * dim_beta,
                     int * dim_gamma,
                     int * dim_alpha,

                     double * data_response,
                     double * data_random,
                     double * data_fixed,
                     double * data_level2,

                     int * total_missingR_data,
                     int * number_response_missing,
                     int * response_missing,

                     int * total_missingRP_data,
                     int * number_random_predictors_missing,
                     int * random_predictors_missing,

                     int * total_missingFP_data,
                     int * number_fixed_predictors_missing,
                     int * fixed_predictors_missing,

                     int * unique_error_Cov,
                     int * dim_error_Cov,
                     double * error_Cov,

                     double * degreesOfFreedom_likelihood,
                     int * likelihood_type,

                     double * gammamean,
                     double * gammaCov,
                     double * gammaDF,
                     int * prior_gamma_type,

                     double * betamean,
                     double * betaCov,
                     double * betaDF,
                     int * prior_beta_type,

                     double * alphamean,
                     double * alphaCov,
                     double * alphaDF,
                     int * prior_alpha_type,

                     double * sigmaDF,
                     double * sigmaScale,
                     double * sigmaPower,
                     int * common_sigma,
                     int * prior_sigma_type,

                     double * tauDF,
                     double * tauScale,
                     double * tauPower,
                     int * prior_tau_type,

                     int * read_init_point,
                     double * betaInit,
                     double * gammaInit,
                     double * alphaInit,
                     double * sigma2Init,
                     double * tau2Init,

                     int * burnInLength,
                     int * simulationsToPerform,
                     int * sampleFrequency,

                     int * print_statistics,
                     double * output_simulations,
		       double * gibbs_drawing_stats ) throw( rtErr )
{
  bool random_effects, fixed_effects, second_effects;
  int i, j, k, l, j1, j2, dim, beta_rxp, gamma_rxp, alpha_rxp, total_dim, number_obs, start_index, simulations_kept;

  CMatrix ** response;
  CMatrix *** random_predictors;
  CMatrix *** fixed_predictors;
  CMatrix ** second_predictors;

  

  /*set random seed*/
  GetRNGstate();

  random_effects = false; 
  fixed_effects = false;
  second_effects = false;

  beta_rxp = 0;
  gamma_rxp = 0;
  alpha_rxp = 0;

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

        //printf("data[%d]  ",  start_index + j * ((int) (*dim_response)) + k );
      }
    }
    start_index += ((int) number_data[i]) * ((int) (*dim_response));
    number_obs += ((int) number_data[i]);

    //printf( "\nfit: number data[%d] = %d,  response[%d] = \n", i, ((int) number_data[i]), i );
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
    Rprintf( "fitBayesianHLM: No fixed effects nor random effects provided. This model is not valid.\n" );
    char the_error[] = "fitBayesianHLM: No fixed effects nor random effects provided. This model is not valid.";
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


  //construct appropriate error covariance matrix (if needed)
  //check if error_Cov is not the identity
  transformResponseAndPredictorsMV( number_groups, number_data, random_effects, fixed_effects,
                                    response, random_predictors, fixed_predictors,
                                    unique_error_Cov, dim_error_Cov, error_Cov );



  BayesianMVHLM bayes_MVHLinear;
  
  bayes_MVHLinear.initialize( response, ((int) (*number_groups) ) );

  if ( random_effects )
  {
    bayes_MVHLinear.randomEffects( random_predictors );
  }

  if ( fixed_effects )
  {
    bayes_MVHLinear.fixedEffects( fixed_predictors );
  }

  if ( second_effects )
  {
    bayes_MVHLinear.secondStageRandomEffects( second_predictors );
  }


  if ( (int) (*total_missingR_data) > 0 )
  {
    k = 0;
    bayes_MVHLinear.initializeResponseMissingData( true );
    for ( i = 0; i < number_obs; i++ )
    { 
      if ( (int) number_response_missing[i] > 0 )
      {
        CVector m_comp( (int) number_response_missing[i] );
        for ( j = 0; j < (int) number_response_missing[i]; j++ )
	{
          //counts starts from 0 in C++, but 1 in SPlus
          m_comp.Val(j) = ((int) response_missing[ k ]) - 1 ;
          k++;
        }

        bayes_MVHLinear.responseMissingData( i, m_comp );
      } 
    }
  }//end if missing responses


  if ( (int) (*total_missingRP_data) > 0 && random_effects )
  {
    k = 0;
    bayes_MVHLinear.initializeRandomPredictorsMissingData( true );
    for ( i = 0; i < number_obs; i++ )
    { 
      if ( (int) number_random_predictors_missing[i] > 0 )
      {
        CVector m_comp( (int) number_random_predictors_missing[i] );
        for ( j = 0; j < (int) number_random_predictors_missing[i]; j++ )
	{
          //counts starts from 0 in C++, but 1 in SPlus
          m_comp.Val(j) = ((int) random_predictors_missing[ k ]) - 1 ;
          k++;
        }

        bayes_MVHLinear.randomPredictorsMissingData( i, m_comp );
      } 
    }
  }//end if missing random_predictors


  if ( (int) (*total_missingFP_data) > 0 && fixed_effects )
  {
    k = 0;
    bayes_MVHLinear.initializeFixedPredictorsMissingData( true );
    for ( i = 0; i < number_obs; i++ )
    { 
      if ( (int) number_fixed_predictors_missing[i] > 0 )
      {
        CVector m_comp( (int) number_fixed_predictors_missing[i] );
        for ( j = 0; j < (int) number_fixed_predictors_missing[i]; j++ )
	{
          //counts starts from 0 in C++, but 1 in SPlus
          m_comp.Val(j) = ((int) fixed_predictors_missing[ k ]) - 1 ;
          k++;
        }

        bayes_MVHLinear.fixedPredictorsMissingData( i, m_comp );
      } 
    }
  }//end if missing fixed_predictors


  //printf("fit: number_obs = %d,  beta_rxp = %d,  gamma_rxp = %d,  alpha_rxp = %d. \n", number_obs, beta_rxp, gamma_rxp, alpha_rxp );

  bool by_column = false; //indicate not to read by row

  //the beta prior
  // 0 ==> usual full model
  // 1 ==> no alpha parameter in model. A single beta0
  // 2 ==> no alpha parameter in model. beta0 for each group
  if ( random_effects )
  {
    CMatrix p_betaCov( betaCov, beta_rxp, beta_rxp, by_column );
  
    if ( ( (int) (*prior_beta_type) ) == 0 )
    {
      bayes_MVHLinear.betaPrior( &p_betaCov );
    }
    else if ( ( (int) (*prior_beta_type) ) == 1 )
    {
      CVector p_beta0( betamean, beta_rxp );
      bayes_MVHLinear.betaPrior( &p_beta0, &p_betaCov );
    }    
    else if ( ( (int) (*prior_beta_type) ) == 2 )
    {
      CMatrix p_beta0( betamean, beta_rxp, ((int) (*number_groups) ), by_column );
      bayes_MVHLinear.betaPrior( &p_beta0, &p_betaCov );
    }    
  }

  if ( fixed_effects )
  {
    // 0 ==> usual Normal distribution
    // 1 ==> non informative
    if ( ( (int) (*prior_gamma_type) ) == 0 )
    {
      CVector p_gamma0( gammamean, gamma_rxp );
      CMatrix p_gammaCov( gammaCov, gamma_rxp, gamma_rxp, by_column );
      bayes_MVHLinear.gammaPrior( &p_gamma0, &p_gammaCov );
    }
    else
    {
      bayes_MVHLinear.gammaPriorNonInformative( gamma_rxp );
    }
  }

  if ( random_effects )
  {
    if ( second_effects )
    {
      // 0 ==> usual Normal distribution
      // 1 ==> non informative
      if ( ( (int) (*prior_alpha_type) ) == 0 )
      {
        CVector p_alpha0( alphamean, alpha_rxp );
        CMatrix p_alphaCov( alphaCov, alpha_rxp, alpha_rxp, by_column );
        bayes_MVHLinear.alphaPrior( &p_alpha0, &p_alphaCov );
      }
      else if ( ( (int) (*prior_alpha_type) ) == 1 )
      {
        bayes_MVHLinear.alphaPriorNonInformative( alpha_rxp );
      }
    }
  }

  // for common sigma
  // 0 ==> independent sigmas; different hyper prior parameters
  // 1 ==> independent sigmas; same hyper prior parameters
  // 2 ==> common sigma; same hyper prior parameters
  if ( ( (int) (*common_sigma) ) == 0 || ( (int) (*common_sigma) ) == 1 )
  {
    bayes_MVHLinear.sigma2CommonPrior( false );

    //for prior sigma type
    // 0 ==> InvChisq
    // 1 ==> non informative
    // 2 ==> uniform shrinkage
    // 3 ==> du Mouchel
    // 4 ==> InvWishart
    // 5 ==> mass Point (known)
    if ( ( (int) (*common_sigma) ) == 0 )
    {
      if ( ( (int) (*prior_sigma_type) ) == 0 )
      {
        bayes_MVHLinear.sigma2PriorInvChisq( sigmaDF, sigmaScale );
      }
      else if ( ( (int) (*prior_sigma_type) ) == 1 )
      {
        bayes_MVHLinear.sigma2PriorNonInformative( sigmaPower );
      }
      else if ( ( (int) (*prior_sigma_type) ) == 2 )
      {
        bayes_MVHLinear.sigma2PriorUniformShrinkage( sigmaScale );
      }
      else if ( ( (int) (*prior_sigma_type) ) == 3 )
      {
        bayes_MVHLinear.sigma2PriorDuMouchel( sigmaScale );
      }
      else if (  ( (int) (*prior_sigma_type) ) == 4  )
      {
        total_dim = ((int)  ( (*dim_response) * (*dim_response) ) );
        CMatrix ** p_sigmaCov;
        p_sigmaCov = new CMatrix * [ (int) (*number_groups) ];

        for ( i = 0; i < (int) (*number_groups); i++ )
	{
          p_sigmaCov[i] = new CMatrix( ((int) (*dim_response) ), ((int) (*dim_response) ) );
          for ( j = 0; j < (int) (*dim_response); j++ )
	  {
            for ( k = 0; k < (int) (*dim_response); k++ )
	    {
              p_sigmaCov[i]->Val( j, k ) = sigmaScale[ i * total_dim + j * ((int) (*dim_response)) + k ];
            }
          }

          //printf( "fit: p sigmaCov[ %d ] is \n", i);
          //p_sigmaCov[i]->Print();

        }//end for i


        bayes_MVHLinear.sigma2PriorInvWishart( sigmaDF, p_sigmaCov );

        for ( i = 0; i < (int) (*number_groups); i++ )
	{
          delete p_sigmaCov[i];
        }
        delete [] p_sigmaCov;
        
      }
      else if (  ( (int) (*prior_sigma_type) ) == 5  )
      {
        total_dim = ((int)  ( (*dim_response) * (*dim_response) ) );
        CMatrix ** p_sigmaCov;
        p_sigmaCov = new CMatrix * [ (int) (*number_groups) ];

        for ( i = 0; i < (int) (*number_groups); i++ )
	{
          p_sigmaCov[i] = new CMatrix( ((int) (*dim_response) ), ((int) (*dim_response) ) );
          for ( j = 0; j < (int) (*dim_response); j++ )
	  {
            for ( k = 0; k < (int) (*dim_response); k++ )
	    {
              p_sigmaCov[i]->Val( j, k ) = sigmaScale[ i * total_dim + j * ((int) (*dim_response)) + k ];
            }
          }

        }//end for i


        bayes_MVHLinear.sigma2Known( p_sigmaCov );

        for ( i = 0; i < (int) (*number_groups); i++ )
	{
          delete p_sigmaCov[i];
        }
        delete [] p_sigmaCov;
        
      }
    }
  }
  else
  {
    bayes_MVHLinear.sigma2CommonPrior( true );
  }

  if ( ( (int) (*common_sigma) ) != 0 )
  {
    if ( ( (int) (*prior_sigma_type) ) == 0 )
    {
      bayes_MVHLinear.sigma2PriorInvChisq( (*sigmaDF), (*sigmaScale) );
    }
    else if ( ( (int) (*prior_sigma_type) ) == 1 )
    {
      bayes_MVHLinear.sigma2PriorNonInformative( (*sigmaPower) );
    }
    else if ( ( (int) (*prior_sigma_type) ) == 2 )
    {
      bayes_MVHLinear.sigma2PriorUniformShrinkage( (*sigmaScale) );
    }
    else if ( ( (int) (*prior_sigma_type) ) == 3 )
    {
      bayes_MVHLinear.sigma2PriorDuMouchel( (*sigmaScale) );
    }
    else if (  ( (int) (*prior_sigma_type) ) == 4 )
    {
      CMatrix p_sigmaCov( sigmaScale, ((int) (*dim_response) ), ((int) (*dim_response) ), by_column );


      //printf( "fit: p sigmaCov is \n", i);
      //p_sigmaCov.Print();

      bayes_MVHLinear.sigma2PriorInvWishart( (*sigmaDF), &p_sigmaCov );
    }
    else if (  ( (int) (*prior_sigma_type) ) == 5 )
    {
      CMatrix p_sigmaCov( sigmaScale, ((int) (*dim_response) ), ((int) (*dim_response) ), by_column );

      bayes_MVHLinear.sigma2Known( &p_sigmaCov );
    }
  }

  if ( random_effects )
  {
    // 0 ==> InvChisq
    // 1 ==> non informative
    // 2 ==> uniform shrinkage
    // 3 ==> du Mouchel
    // 4 ==> InvWishart
    if ( ( (int) (*prior_tau_type) ) == 0 )
    {
      bayes_MVHLinear.tau2PriorInvChisq( (*tauDF), (*tauScale) );
    }
    else if ( ( (int) (*prior_tau_type) ) == 1 )
    {
      bayes_MVHLinear.tau2PriorNonInformative( (*tauPower) );
    }
    else if ( ( (int) (*prior_tau_type) ) == 2 )
    {
      bayes_MVHLinear.tau2PriorUniformShrinkage( (*tauScale) );
    }
    else if ( ( (int) (*prior_tau_type) ) == 3 )
    {
      bayes_MVHLinear.tau2PriorDuMouchel( (*tauScale) );
    }
    else if ( ( (int) (*prior_tau_type) ) == 4 )
    {
      CMatrix p_betaCov( tauScale, beta_rxp, beta_rxp, by_column );

      //printf( "fit: tau2 betaCovis \n");
      //p_betaCov.Print();

      bayes_MVHLinear.betaCovPriorInvWishart( (*tauDF), &p_betaCov );
    }
  }
  
  if ( random_effects )
  {
    if ( (*betaDF) > 0 )
    {
      bayes_MVHLinear.betaTPrior( (*betaDF) );
    }

    if ( second_effects && (*alphaDF) > 0 )
    {
      bayes_MVHLinear.alphaTPrior( (*alphaDF) );
    }
  }

  if ( fixed_effects )
  {
    if ( (*gammaDF) > 0 )
    {
      bayes_MVHLinear.gammaTPrior( (*gammaDF) );
    }
  }

  // 0 ==> normal
  // 1 ==> t errors
  // 2 ==> group errors
  if ( (*likelihood_type) == 1 )
  {
    if ( (*degreesOfFreedom_likelihood) > 0 )
    {
      bayes_MVHLinear.tLikelihood( (*degreesOfFreedom_likelihood) );
    }
    else
    {
      bayes_MVHLinear.tLikelihood( 3.0 );
    }
  }
  else if ( (*likelihood_type) == 2 )
  {
    if ( (*degreesOfFreedom_likelihood) > 0 )
    {
      bayes_MVHLinear.groupTLikelihood( (*degreesOfFreedom_likelihood) );
    }
    else
    {
      bayes_MVHLinear.groupTLikelihood( 3.0 );
    }
  }


  //get ready for sampling
  if ( (int) (*read_init_point) )
  {
    if ( random_effects )
    {
      if ( second_effects )
      {
        CVector init_alpha( alphaInit, alpha_rxp );

        //printf( "BHL: init alpha \n" );
        //init_alpha.Print();

        bayes_MVHLinear.samplerAlphaInitialPoint( init_alpha );

        //full model initial point for beta
        //this is only valid for prior parameters, not for initial points
        //CMatrix init_beta( beta_rxp, ((int) (*number_groups)) );
        //for ( i = 0; i < ((int) (*number_groups)); i++ )
	//{
	//  init_beta.setColumn( i, (*second_predictors[i]) * init_alpha );
        //}
        //bayes_MVHLinear.samplerBetaInitialPoint( init_beta );
       
      }
      //this is only valid for prior parameters, not for initial points
      //else if ( ((int) (*prior_beta_type) ) == 1 )
      //{
      //  CVector init_beta( betaInit, ((int) (*dim_beta) ) );
      //  bayes_MVHLinear.samplerBetaInitialPoint( init_beta );
      //}
      else //if ( ((int) (*prior_beta_type) ) == 2 )
      {
        CMatrix init_beta( betaInit, beta_rxp, ((int) (*number_groups) ), by_column );

        //printf( "BHL: init beta \n" );
        //init_beta.Print();

        bayes_MVHLinear.samplerBetaInitialPoint( init_beta );
      }
    
      if ( ((int) (*prior_tau_type) ) != 4 )
      {
        bayes_MVHLinear.samplerTau2InitialPoint( (*tau2Init) );

        //printf( "BHL: init random var \n %f\n", (*tau2Init) );
        

      }
      else
      {
        CMatrix init_tau2( tau2Init, beta_rxp, beta_rxp, by_column );

        //printf( "BHL: init random var \n" );
        //init_tau2.Print();

        bayes_MVHLinear.samplerTau2InitialPoint( init_tau2 );
      }
    }

    if ( fixed_effects )
    {
      CVector init_gamma( gammaInit, gamma_rxp );

      //printf( "BHL: init gamma \n" );
      //init_gamma.Print();

      bayes_MVHLinear.samplerGammaInitialPoint( init_gamma );
    }


    if ( ((int) (*common_sigma) ) != 0 )
    {
      if ( ((int) (*prior_sigma_type)) != 4 && ((int) (*prior_sigma_type)) != 5 )
      {
        bayes_MVHLinear.samplerSigma2InitialPoint( (*sigma2Init) );
        //printf( "BHL: init error var \n %f\n", (*sigma2Init) );
      }
      else
      {
        CMatrix init_sigma2( sigma2Init, ((int) (*dim_response)), ((int) (*dim_response)) );

        //printf( "BHL: init error var \n" );
        //init_sigma2.Print();

        bayes_MVHLinear.samplerSigma2InitialPoint( init_sigma2 );
      }
    }
    else
    {
      if ( ((int) (*prior_sigma_type)) != 4 && ((int) (*prior_sigma_type)) != 5 )
      {
        CVector init_sigma2( sigma2Init, ((int) (*number_groups) ) );

        //printf( "BHL: init error var \n" );
        //init_sigma2.Print();

        bayes_MVHLinear.samplerSigma2InitialPoint( init_sigma2 );
      }
      else
      {
        CMatrix sigma2_init_pts( sigma2Init, ((int) (*dim_response)), ((int) (*dim_response)) * ((int) (*number_groups)), by_column );         
        CMatrix ** sigma2_init;
        sigma2_init = new CMatrix * [ ((int) (*number_groups)) ];
        start_index = 0;
        for ( i = 0; i < ((int) (*number_groups)); i++ )
	{
          sigma2_init[i] = new CMatrix( ((int) (*dim_response)), ((int) (*dim_response)) );
          for ( j = 0; j < ((int) (*dim_response)); j++ )
	  {
#ifdef FIX1
            CVector tmpvec = sigma2_init_pts.getColumn( start_index + j );
            sigma2_init[i]->setColumn( j , tmpvec );
#else
            sigma2_init[i]->setColumn( j, sigma2_init_pts.getColumn( start_index + j ) );
#endif
//            sigma2_init[i]->setColumn( j, sigma2_init_pts.getColumn( start_index + j ) );
          }
          start_index += ((int) (*dim_response));
        }//end for i
        bayes_MVHLinear.samplerSigma2InitialPoint( sigma2_init );

        for ( i = 0; i < ((int) (*number_groups)); i++ )
	{
          delete sigma2_init[i];
        }
        delete [] sigma2_init;
      }
    }

  }//end if initial points provided
  else
  {
    bayes_MVHLinear.samplerDefaultInitialPoint();
  }

  //initialize missing data points
  bayes_MVHLinear.samplerMissingVariablesInitialPoint();


  //printf("fit: will do Gibbs sampler\n"); fflush(stdout);

  //get ready to start Gibbs sampler
  Gibbs sampler( (int) (*burnInLength), 
                 (int) (*simulationsToPerform),
                 (int) (*sampleFrequency) );

  if ( (*print_statistics) )
  {
    sampler.saveStatistics( true );
  }


  //start Gibbs sampler
  sampler.doGibbs( &bayes_MVHLinear );

  simulations_kept = sampler.simulationsKept();

  //printf("fit: simulations_kept = %d\n", simulations_kept );
  //fflush(stdout);

  if ( (*print_statistics) )
  {
    sampler.printDrawingStats( gibbs_drawing_stats );   
  }


  //get the output simulation samples 
  bayes_MVHLinear.simulationsToArray( output_simulations, simulations_kept );

  /* update random seed */
  PutRNGstate();


  for ( i = 0; i < ((int) (*number_groups) ); i++ )
  {
    delete response[i];  
  }//end for i
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

  //printf("fitBayesianMD: Done\n"); fflush(stdout);

}//end

    

void  transformResponseAndPredictorsMV( int * number_groups, 
                                        int * number_data,
                                        bool random_effects,
                                        bool fixed_effects,
                                        CMatrix ** response, 
                                        CMatrix *** random_predictors, 
                                        CMatrix *** fixed_predictors,   
                                        int * unique_error_Cov, 
                                        int * dim_error_Cov, 
                                        double * error_Cov ) throw( rtErr )
{
  //construct appropriate error covariance matrix (if needed)
  //check if error_Cov is not the identity

  bool identityCov;
  int i, j, k, start_index;
  double weight;

  if ( ( (int) (*unique_error_Cov) ) == 1 )
  {
    if ( ((int) (*dim_error_Cov) ) == 0 )
    {
      //a scalar times the identity covariance
      if ( (*error_Cov) != 1.0 )
      {
        //transform predictors and response
        weight = 1.0 / sqrt( (*error_Cov) );
        for ( i = 0; i < (int) (*number_groups); i++ )
        {
          response[i]->multiplyByScalar( weight );

          if ( random_effects )
	  {
            for ( j = 0; j < ((int) number_data[i]); j++ )
	    {
              random_predictors[i][j]->multiplyByScalar( weight );
            }
          }

          if ( fixed_effects )
	  {
            for ( j = 0; j < ((int) number_data[i]); j++ )
	    {
              fixed_predictors[i][j]->multiplyByScalar( weight );
            }
          }
        }//end for i
      }
    }
    else
    {
      Rprintf( "fitBayesianHierarchicalLinearModel: wrong error covariance argument.\n" );
      char the_error[] = "fitBayesianHierarchicalLinearModel: wrong error covariance argument.";
      rtErr runtime_error( the_error );
      throw runtime_error;
    }
  }//end if unique
  else
  {
    if ( ((int) (*dim_error_Cov) ) == 0 )
    {
      //a scalar times the identity covariance
      for ( i = 0; i < (int) (*number_groups); i++ )
      {
        if ( error_Cov[i] != 1.0 )
        {
          //transform predictors and response
          weight = 1.0 / sqrt( error_Cov[i] );
          response[i]->multiplyByScalar( weight );

          if ( random_effects )
	  {
            for ( j = 0; j < ((int) number_data[i]); j++ )
	    {
              random_predictors[i][j]->multiplyByScalar( weight );
            }
          }

          if ( fixed_effects )
	  {
            for ( j = 0; j < ((int) number_data[i]); j++ )
	    {
              fixed_predictors[i][j]->multiplyByScalar( weight );
            }
          }
        }
      }//end for i
    }
    else if ( ((int) (*dim_error_Cov) ) == 1 )
    {
      //a diagonal covariance
      start_index = 0;
      for ( i = 0; i < (int) (*number_groups); i++ )
      {
        identityCov = true;
        j = 0;
        while ( identityCov && j < response[i]->Row() )
        {
          if ( error_Cov[ start_index + j ] != 1.0 )
          {
            identityCov = false;
          }
          j++;
        }//end loop

        if ( !identityCov )
        {
          //transform predictors
          CVector errorCov( response[i]->Row() );
          for ( j = 0; j < response[i]->Row(); j++ )
          {
            errorCov.Val(j) = 1.0 / sqrt( error_Cov[ start_index + j ] );
          }

          for ( j = 0; j < response[i]->Col(); j++ )
	  {
            response[i]->setToWeightedColumn( j, errorCov );
          }

          if ( random_effects )
	  {
            for ( j = 0; j < ((int) number_data[i]); j++ )
	    {
              for ( k = 0; k < random_predictors[i][j]->Col(); k++ )
              {
                random_predictors[i][j]->setToWeightedColumn( k, errorCov );
              }
            }
          }

          if ( fixed_effects )
	  {
            for ( j = 0; j < ((int) number_data[i]); j++ )
	    {
              for ( k = 0; k < fixed_predictors[i][j]->Col(); k++ )
              {
                fixed_predictors[i][j]->setToWeightedColumn( k, errorCov );
              }
            }
          }

        }
        start_index += response[i]->Row();
      }//end for i
    }
    else
    {
      //a matrix for covariance: transform predictors
      CMatrix * working_Cov;
      CVector * cholLT;

      start_index = 0;
      for ( i = 0; i < (int) (*number_groups); i++ )
      {
        working_Cov = new CMatrix( response[i]->Row(), response[i]->Row() );
        cholLT = new CVector( ( response[i]->Row() * ( response[i]->Row() + 1 ) ) / 2 );
        for ( j = 0; j < response[i]->Row(); j++ )
	{
          for ( k = 0; k < response[i]->Row(); k++ )
	  {
            working_Cov->Val( k, j ) = error_Cov[ start_index + j * response[i]->Row() + k ];
          }
	}

        (*cholLT) = working_Cov->choleskyDecomposition();
        working_Cov->assignInverseOfLowerTriangular( (*cholLT) );
        delete cholLT;

        CVector weighted_vector( response[i]->Row() );

        //do the transformation
        for ( j = 0; j < response[i]->Col(); j++ )
	{
          weighted_vector = (*working_Cov) * response[i]->getColumn( j );
          //assign the transformation
          response[i]->setColumn( j, weighted_vector );
        }

        if ( random_effects )
	{
          for ( j = 0; j < ((int) number_data[i]); j++ )
	  {
            for ( k = 0; k < random_predictors[i][j]->Col(); k++ )
            {
              //do the transformation
              weighted_vector = (*working_Cov) * ( random_predictors[i][j]->getColumn(k) );
              //assign the transformation
              random_predictors[i][j]->setColumn( k, weighted_vector );
            }
          }
        }

        if ( fixed_effects )
	{
          for ( j = 0; j < ((int) number_data[i]); j++ )
	  {
            for ( k = 0; k < fixed_predictors[i][j]->Col(); k++ )
            {
              //do the transformation
              weighted_vector = (*working_Cov) * ( fixed_predictors[i][j]->getColumn(k) );
              //assign the transformation
              fixed_predictors[i][j]->setColumn( k, weighted_vector );
            }
          }
        }

        start_index += response[i]->Row() * response[i]->Row();
        delete working_Cov;
      }//end for i
    }//end matrix
  }//end cov is not a unique scalar

}//end

}
