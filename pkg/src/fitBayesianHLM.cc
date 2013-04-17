#include "R.h"
#include "rtErr.h"
#include "Vector.h"
#include "Matrix.h"

#include "BayesianHierarchicalLinearModel.h"
#include "Gibbs.h"


extern "C" {

void  transformResponseAndPredictors( long * number_groups, 
                                      bool random_effects, 
                                      bool fixed_effects,
                                      CVector ** response, 
                                      CMatrix ** random_predictors, 
                                      CMatrix ** fixed_predictors, 
                                      long * dim_beta, 
                                      long * dim_gamma, 
                                      long * unique_error_cov, 
                                      long * dim_error_Cov, 
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
void fitBayesianHLM( long * number_groups,
                     long * number_data,
                     long * dim_beta,
                     long * dim_gamma,
                     long * dim_alpha,

                     double * random_data,
                     double * fixed_data,
                     double * second_data,
                     double * data_response,

                     long * total_missingR_data,
                     long * number_response_missing,
                     long * response_missing,

                     long * total_missingRP_data,
                     long * number_random_predictors_missing,
                     long * random_predictors_missing,

                     long * total_missingFP_data,
                     long * number_fixed_predictors_missing,
                     long * fixed_predictors_missing,

                     long * unique_error_Cov,
                     long * dim_error_Cov,
                     double * error_Cov,

                     double * degreesOfFreedom_likelihood,
                     long * likelihood_type,


                     double * gammamean,
                     double * gammaCov,
                     double * gammaDF,
                     long * prior_gamma_type,

                     double * betamean,
                     double * betaCov,
                     double * betaDF,
                     long * prior_beta_type,

                     double * alphamean,
                     double * alphaCov,
                     double * alphaDF,
                     long * prior_alpha_type,

                     double * sigmaDF,
                     double * sigmaScale,
                     double * sigmaPower,
                     long * common_sigma,
                     long * prior_sigma_type,

                     double * tauDF,
                     double * tauScale,
                     double * tauPower,
                     long * prior_tau_type,

                     long * read_init_point,
                     double * betaInit,
                     double * gammaInit,
                     double * alphaInit,
                     double * sigma2Init,
                     double * tau2Init,

                     long * burnInLength,
                     long * simulationsToPerform,
                     long * sampleFrequency,

                     long * print_statistics,
                     double * output_simulations,
                     double * gibbs_drawing_stats ) throw( rtErr ) 
{
  bool random_effects, fixed_effects, second_effects;
  int i, j, k, dim, number_obs, start_index, simulations_kept;

  CVector ** response;
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
  number_obs = 0;
  for ( i = 0; i < (int) (*number_groups); i++ )
  {
    response[i] = new CVector( (int) number_data[i] );
    for ( j = 0; j < (int) number_data[i]; j++ )
    {
      response[i]->Val(j) = data_response[ start_index + j ]; 
    }
    start_index += ((int) number_data[i]);
    number_obs += ((int) number_data[i]);

    //printf("fit: response[%d] = \n", i );
    //response[i]->Print();

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
          random_predictors[i]->Val( k, j ) = random_data[ start_index + k * dim + j ];
        } 
      }//end for j
      start_index += dim * response[i]->Len();


      //printf("fit: random predictors[%d] = \n", i );
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
          fixed_predictors[i]->Val( k, j ) = fixed_data[ start_index + k * dim + j ];
        } 
      }//end for j
      start_index += dim * response[i]->Len();

      //printf("fit: fixed predictors[%d] = \n", i );
      //fixed_predictors[i]->Print();
    }//end for i
  }

  if ( !random_effects && !fixed_effects )
  {
    printf( "fitBayesianHLM: No fixed effects nor random effects provided. This model is not valid.\n" );
    char the_error[] = "fitBayesianHLM: No fixed effects nor random effects provided. This model is not valid.";
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

      //printf("fit: level2 predictors[%d] = \n", i );
      //second_predictors[i]->Print();
    }//end for i
  }


  //construct appropriate error covariance matrix (if needed)
  //check if error_Cov is not the identity
  transformResponseAndPredictors( number_groups, random_effects, fixed_effects,
                                  response, random_predictors, fixed_predictors, 
                                  dim_beta, dim_gamma, unique_error_Cov, dim_error_Cov, error_Cov );




  BayesianHierarchicalLinearModel bayes_HLinear;
  
  bayes_HLinear.initialize( response, ((int) (*number_groups) ) );

  if ( random_effects )
  {
    bayes_HLinear.randomEffects( random_predictors );
  }

  if ( fixed_effects )
  {
    bayes_HLinear.fixedEffects( fixed_predictors );
  }

  if ( second_effects )
  {
    bayes_HLinear.secondStageRandomEffects( second_predictors );
  }



  //missing data handling
  if ( (int) (*total_missingR_data) > 0 )
  {
    k = 0;
    bayes_HLinear.initializeResponseMissingData( true );
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

        bayes_HLinear.responseMissingData( i, m_comp );
      } 
    }
  }//end if missing responses


  if ( (int) (*total_missingRP_data) > 0 && random_effects )
  {
    k = 0;
    bayes_HLinear.initializeRandomPredictorsMissingData( true );
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

        bayes_HLinear.randomPredictorsMissingData( i, m_comp );
      } 
    }
  }//end if missing random_predictors


  if ( (int) (*total_missingFP_data) > 0 && fixed_effects )
  {
    k = 0;
    bayes_HLinear.initializeFixedPredictorsMissingData( true );
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

        bayes_HLinear.fixedPredictorsMissingData( i, m_comp );
      } 
    }
  }//end if missing fixed_predictors



  for ( i = 0; i < ((int) (*number_groups) ); i++ )
  {
    delete response[i];
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


  bool by_column = false; //indicate not to read by row

  //the beta prior
  // 0 ==> usual full model
  // 1 ==> no alpha parameter in model. A single beta0
  // 2 ==> no alpha parameter in model. beta0 for each group
  // 3 ==> non-informative prior
  if ( random_effects )
  {
    CMatrix p_betaCov( betaCov, ((int) (*dim_beta) ), ((int) (*dim_beta) ), by_column );
  
    if ( ( (int) (*prior_beta_type) ) == 0 )
    {
      bayes_HLinear.betaPrior( &p_betaCov );
    }
    else if ( ( (int) (*prior_beta_type) ) == 1 )
    {
      CVector p_beta0( betamean, ((int) (*dim_beta) ) );
      bayes_HLinear.betaPrior( &p_beta0, &p_betaCov );
    }    
    else if ( ( (int) (*prior_beta_type) ) == 2 )
    {
      CMatrix p_beta0( betamean, ((int) (*dim_beta) ), ((int) (*number_groups) ), by_column );
      bayes_HLinear.betaPrior( &p_beta0, &p_betaCov );
    }    
    else if ( ( (int) (*prior_beta_type) ) == 3 )
    {
      bayes_HLinear.betaPriorNonInformative( ((int) (*dim_beta)) );
    }
  }

  if ( fixed_effects )
  {
    // 0 ==> usual Normal distribution
    // 1 ==> non informative
    if ( ( (int) (*prior_gamma_type) ) == 0 )
    {
      CVector p_gamma0( gammamean, ((int) (*dim_gamma) ) );
      CMatrix p_gammaCov( gammaCov, ((int) (*dim_gamma) ), ((int) (*dim_gamma) ), by_column );
      bayes_HLinear.gammaPrior( &p_gamma0, &p_gammaCov );
    }
    else
    {
      bayes_HLinear.gammaPriorNonInformative( ((int) (*dim_gamma) ) );
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
        CVector p_alpha0( alphamean, ((int) (*dim_alpha) ) );
        CMatrix p_alphaCov( alphaCov, ((int) (*dim_alpha) ), ((int) (*dim_alpha) ), by_column );
        bayes_HLinear.alphaPrior( &p_alpha0, &p_alphaCov );
      }
      else if ( ( (int) (*prior_alpha_type) ) == 1 )
      {
        bayes_HLinear.alphaPriorNonInformative( ((int) (*dim_alpha) ) );
      }
    }
  }

  // for common sigma
  // 0 ==> independent sigmas; different hyper prior parameters
  // 1 ==> independent sigmas; same hyper prior parameters
  // 2 ==> common sigma; same hyper prior parameters
  if ( ( (int) (*common_sigma) ) == 0 || ( (int) (*common_sigma) ) == 1 )
  {
    bayes_HLinear.sigma2CommonPrior( false );

    //for prior sigma type
    // 0 ==> InvChisq
    // 1 ==> non informative
    // 2 ==> uniform shrinkage
    // 3 ==> du Mouchel
    // 4 ==> known
    if ( ( (int) (*common_sigma) ) == 0 )
    {
      if ( ( (int) (*prior_sigma_type) ) == 0 )
      {
        bayes_HLinear.sigma2PriorInvChisq( sigmaDF, sigmaScale );
      }
      else if ( ( (int) (*prior_sigma_type) ) == 1 )
      {
        bayes_HLinear.sigma2PriorNonInformative( sigmaPower );
      }
      else if ( ( (int) (*prior_sigma_type) ) == 2 )
      {
        bayes_HLinear.sigma2PriorUniformShrinkage( sigmaScale );
      }
      else if ( ( (int) (*prior_sigma_type) ) == 3 )
      {
        bayes_HLinear.sigma2PriorDuMouchel( sigmaScale );
      }
      else if ( ( (int) (*prior_sigma_type) ) == 4 )
      {
        bayes_HLinear.sigma2Known( sigmaScale );
      }
    }
  }
  else
  {
    bayes_HLinear.sigma2CommonPrior( true );
  }

  if ( ( (int) (*common_sigma) ) != 0 )
  {
    if ( ( (int) (*prior_sigma_type) ) == 0 )
    {
      bayes_HLinear.sigma2PriorInvChisq( (*sigmaDF), (*sigmaScale) );
    }
    else if ( ( (int) (*prior_sigma_type) ) == 1 )
    {
      bayes_HLinear.sigma2PriorNonInformative( (*sigmaPower) );
    }
    else if ( ( (int) (*prior_sigma_type) ) == 2 )
    {
      bayes_HLinear.sigma2PriorUniformShrinkage( (*sigmaScale) );
    }
    else if ( ( (int) (*prior_sigma_type) ) == 3 )
    {
      bayes_HLinear.sigma2PriorDuMouchel( (*sigmaScale) );
    }
    else if ( ( (int) (*prior_sigma_type) ) == 4 )
    {
      bayes_HLinear.sigma2Known( (*sigmaScale) );
    }
  }


  //make sure beta has an informative prior
  if ( random_effects && ((int) (*prior_beta_type)) != 3 )
  {
    // 0 ==> InvChisq
    // 1 ==> non informative
    // 2 ==> uniform shrinkage
    // 3 ==> du Mouchel
    // 4 ==> InvWishart
    if ( ( (int) (*prior_tau_type) ) == 0 )
    {
      bayes_HLinear.tau2PriorInvChisq( tauDF, tauScale );
    }
    else if ( ( (int) (*prior_tau_type) ) == 1 )
    {
      bayes_HLinear.tau2PriorNonInformative( tauPower );
    }
    else if ( ( (int) (*prior_tau_type) ) == 2 )
    {
      bayes_HLinear.tau2PriorUniformShrinkage( tauScale );
    }
    else if ( ( (int) (*prior_tau_type) ) == 3 )
    {
      bayes_HLinear.tau2PriorDuMouchel( tauScale );
    }
    else if ( ( (int) (*prior_tau_type) ) == 4 )
    {
      CMatrix p_betaCov( tauScale, ((int) (*dim_beta) ), ((int) (*dim_beta) ), by_column );

      //printf("fit: setting invWishart. original df = %f. original cov dim is %d, cov is\n",  (*tauDF), ((int) (*dim_beta) ) );
      //p_betaCov.Print(); fflush(stdout);
      
      bayes_HLinear.betaCovPriorInvWishart( (*tauDF), &p_betaCov );
      #ifdef DEBUG1
        printf("Inverse Wishart prior set\n"); fflush(stdout);
      #endif
    }
  }
  
  if ( random_effects )
  {
    if ( (*betaDF) > 0 )
    {
      bayes_HLinear.betaTPrior( (*betaDF) );
    }

    if ( second_effects && (*alphaDF) > 0 )
    {
      bayes_HLinear.alphaTPrior( (*alphaDF) );
    }
  }

  if ( fixed_effects )
  {
    if ( (*gammaDF) > 0 )
    {
      bayes_HLinear.gammaTPrior( (*gammaDF) );
    }
  }

  // 0 ==> normal
  // 1 ==> t errors
  // 2 ==> group errors
  if ( (*likelihood_type) == 1 )
  {
    if ( (*degreesOfFreedom_likelihood) > 0 )
    {
      bayes_HLinear.tLikelihood( (*degreesOfFreedom_likelihood) );
    }
    else
    {
      bayes_HLinear.tLikelihood( 3.0 );
    }
  }
  else if ( (*likelihood_type) == 2 )
  {
    if ( (*degreesOfFreedom_likelihood) > 0 )
    {
      bayes_HLinear.groupTLikelihood( (*degreesOfFreedom_likelihood) );
    }
    else
    {
      bayes_HLinear.groupTLikelihood( 3.0 );
    }
  }

  //get ready for sampling
  if ( (int) (*read_init_point) )
  {
    if ( random_effects )
    {
      if ( second_effects )
      {
        CVector init_alpha( alphaInit, ((int) (*dim_alpha) ) );

        //printf( "BHL: init alpha \n" );
        //init_alpha.Print();

        bayes_HLinear.samplerAlphaInitialPoint( init_alpha );

      }
      //this is only valid for prior parameters, not for initial points
      //else if ( ((int) (*prior_beta_type) ) == 1 )
      //{
      //  CVector init_beta( betaInit, ((int) (*dim_beta) ) );
      //  bayes_HLinear.samplerBetaInitialPoint( init_beta );
      //}
      else //if ( ((int) (*prior_beta_type) ) == 2 )
      {
        CMatrix init_beta( betaInit, ((int) (*dim_beta) ), ((int) (*number_groups) ), by_column );


        //printf( "BHL: init beta \n" );
        //init_beta.Print();

        bayes_HLinear.samplerBetaInitialPoint( init_beta );
      }

      //there is an informative prior for beta
      if ( ((int) (*prior_beta_type)) != 3 )    
      {
        if ( ((int) (*prior_tau_type) ) != 4 )
        {
          bayes_HLinear.samplerTau2InitialPoint( tau2Init );

          //printf( "BHL: init random var \n %f\n", (*tau2Init) );
        
        }
        else
        {
          CMatrix init_tau2( tau2Init, ((int) (*dim_beta) ), ((int) (*dim_beta) ), by_column );

          //printf( "BHL: init random var \n" );
          //init_tau2.Print();

          bayes_HLinear.samplerTau2InitialPoint( init_tau2 );
        }
      }
    }

    if ( fixed_effects )
    {
      CVector init_gamma( gammaInit, ((int) (*dim_gamma) ) );

      //printf( "BHL: init gamma \n" );
      //init_gamma.Print();

      bayes_HLinear.samplerGammaInitialPoint( init_gamma );
    }


    if ( ((int) (*prior_sigma_type)) != 4 )
    {
      if ( ((int) (*common_sigma) ) != 0 )
      {
        bayes_HLinear.samplerSigma2InitialPoint( (*sigma2Init) );
 
        //printf( "BHL: init error var \n %f\n", (*sigma2Init) );
      }
      else
      {
        CVector init_sigma2( sigma2Init, ((int) (*number_groups) ) );

        //printf( "BHL: init error var \n" );
        //init_sigma2.Print();

        bayes_HLinear.samplerSigma2InitialPoint( init_sigma2 );
      }
    }

  }//end if initial points provided
  else
  {
    bayes_HLinear.samplerDefaultInitialPoint();
  }

  //initialize missing data points
  bayes_HLinear.samplerMissingVariablesInitialPoint();
  #ifdef DEBUG1
     printf("Initial values set\n"); fflush(stdout);
  #endif
   
  //get ready to start Gibbs sampler
  Gibbs sampler( (int) (*burnInLength), 
                 (int) (*simulationsToPerform),
                 (int) (*sampleFrequency) );

  if ( (*print_statistics) )
  {
    sampler.saveStatistics( true );
  }

  //start Gibbs sampler
  sampler.doGibbs( &bayes_HLinear );

  simulations_kept = sampler.simulationsKept();

  if ( (*print_statistics) )
  {
    sampler.printDrawingStats( gibbs_drawing_stats );   
  }



  //get the output simulation samples 
  bayes_HLinear.simulationsToArray( output_simulations, simulations_kept );

  /* update random seed */
  PutRNGstate();

  //printf("fit: out of here\n\n" ); fflush(stdout);
}//end


}
    

void  transformResponseAndPredictors( long * number_groups, 
                                      bool random_effects, 
                                      bool fixed_effects,
                                      CVector ** response, 
                                      CMatrix ** random_predictors, 
                                      CMatrix ** fixed_predictors, 
                                      long * dim_beta, 
                                      long * dim_gamma, 
                                      long * unique_error_Cov, 
                                      long * dim_error_Cov, 
                                      double * error_Cov ) throw( rtErr )
{
  //construct appropriate error covariance matrix (if needed)
  //check if error_Cov is not the identity

  bool identityCov;
  int i, j, k, start_index;

  if ( ( (int) (*unique_error_Cov) ) == 1 )
  {
    if ( ((int) (*dim_error_Cov) ) == 0 )
    {
      //a scalar times the identity covariance
      if ( (*error_Cov) != 1.0 && (*error_Cov) > 0 )
      {
        //transform predictors and response
        for ( i = 0; i < (int) (*number_groups); i++ )
        {
          response[i]->multiplyByScalar( 1.0 / sqrt( (*error_Cov) ) );
          if ( random_effects )
	  {
            random_predictors[i]->multiplyByScalar( 1.0 / sqrt( (*error_Cov) ) );
          }

          if ( fixed_effects )
	  {
            fixed_predictors[i]->multiplyByScalar( 1.0 / sqrt( (*error_Cov) ) );
          }
        }//end for i
      }
      else if ( (*error_Cov) <= 0 )
      {
        printf( "fitBayesianHierarchicalLinearModel: Negative or zero error covariance argument.\n" );
        char the_error[] = "fitBayesianHierarchicalLinearModel: Negative or zero error covariance argument.";
        rtErr runtime_error( the_error );
        throw runtime_error;        
      }
    }
    else
    {
      printf( "fitBayesianHierarchicalLinearModel: wrong error covariance argument.\n" );
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
          response[i]->multiplyByScalar( 1.0 / sqrt( error_Cov[i] ) );

          if ( random_effects )
	  {
            random_predictors[i]->multiplyByScalar( 1.0 / sqrt( error_Cov[i] ) );
          }

          if ( fixed_effects )
	  {
            fixed_predictors[i]->multiplyByScalar( 1.0 / sqrt( error_Cov[i] ) );
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
        while ( identityCov && j < response[i]->Len() )
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
          CVector errorCov( response[i]->Len() );
          for ( j = 0; j < response[i]->Len(); j++ )
          {
            errorCov.Val(j) = 1.0 / sqrt( error_Cov[ start_index + j ] );
          }

          response[i]->setToWeighted( errorCov );
          if ( random_effects )
	  {
            for ( j = 0; j < random_predictors[i]->Col(); j++ )
            {
              random_predictors[i]->setToWeightedColumn( j, errorCov );
            }
          }

          if ( fixed_effects )
	  {
            for ( j = 0; j < fixed_predictors[i]->Col(); j++ )
            {
              fixed_predictors[i]->setToWeightedColumn( j, errorCov );
            }
          }
        }

        start_index += response[i]->Len();
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
        working_Cov = new CMatrix( response[i]->Len(), response[i]->Len() );
        cholLT = new CVector( ( response[i]->Len() * ( response[i]->Len() + 1 ) ) / 2 );
        for ( j = 0; j < response[i]->Len(); j++ )
	{
          for ( k = 0; k < response[i]->Len(); k++ )
	  {
            working_Cov->Val( k, j ) = error_Cov[ start_index + j * response[i]->Len() + k ];
          }
	}

        (*cholLT) = working_Cov->choleskyDecomposition();
        working_Cov->assignInverseOfLowerTriangular( (*cholLT) );
        delete cholLT;

        CVector weighted_vector( response[i]->Len() );

        //do the transformation
        weighted_vector = (*working_Cov) * ( *(response[i]) );
        //assign the transformation
        ( *(response[i]) ) = weighted_vector;

        if ( random_effects )
	{
          for ( j = 0; j < ( (int) (*dim_beta) ); j++ )
          {
            //do the transformation
            weighted_vector = (*working_Cov) * ( random_predictors[i]->getColumn(j) );
            //assign the transformation
            random_predictors[i]->setColumn( j, weighted_vector );
          }
        }

        if ( fixed_effects )
	{
          for ( j = 0; j < ( (int) (*dim_gamma) ); j++ )
          {
            //do the transformation
            weighted_vector = (*working_Cov) * ( fixed_predictors[i]->getColumn(j) );
            //assign the transformation
            fixed_predictors[i]->setColumn( j, weighted_vector );
          }
        }

        start_index += response[i]->Len() * response[i]->Len();
        delete working_Cov;
      }//end for i
    }//end matrix
  }//end cov is not a unique scalar

}//end

