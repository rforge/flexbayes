#include "R.h"
#include "rtErr.h"
#include "Vector.h"
#include "Matrix.h"

#include "BayesianBinomialLogitLinkModel.h"
#include "BayesianHierarchicalGlmModel.h"
#include "MetropolisHastings.h"


extern "C" {
/*
  Computes a Bayes Hierarchical Binomial Regression model:

  |  counts_ij ~ Binomial( n_ij, theta_ij )
  |
  |  theta_ij ~ Beta( xi_j mu_ij, xi_j ( 1 -  mu_ij ) )
  |
  |  logit( mu_ij ) = X_ij beta_j + M_ij gamma
  |
  |  beta_j = Z_j alpha + beta_error :the random effects
  |
  |  gamma ~ flat prior == N( 0, +infty ) :the fixed effects
  |        or ~ Normal( gamma_0, Gamma_0 )
  |
  |  alpha ~ N( alpha_0, V_0 )
  |  beta_error ~ N( 0, tau2 V )
  |
  |  tau2 ~ p(tau2) or tau2V ~ InvWishart( nu, S )
  |
  |  xi_j ~ Uniform Shrinkage( xi_0 )

    counts_ij: data response for i-th case on j-th group (counts)
    e_ij: data trials 
    X_j: data predictors for the random effects of the j-th group
    M_j: data predictors for the fixed effects
    Z_j: data predictors for the "population" or second level predictors

    j = 1, ..., J.

  theta_ij are nuisance parameters, and are integrated out in the algorithm.

  p(tau2): could be InvChisq( nuTau0, tau02 )
                     non Informative: proper:   Uniform Shrinkage ( tau02 )
                                                Du Mouchel( tau02 )
                                      improper: proportional to tau2^(1 + b ), with b <= 0.


   Robustness: distributions for gamma, beta and alpha can be made to follow Multivariate Student's t's.

   Computations are based on Gibbs and Metropolis-Hastings sampler drawings.
   burnInLength: number of initial Gibbs sampler iterations to discard.
   simulationsToPerform: number of Gibbs sampler drawings to keep as output.
   sampleFrequency: only iterations that are multiples of sampleFrequency are to be kept.


   OUTPUT
   ------
   The simulationsToPerform output simulations are returned in the array output_simulations.
   This should be an array of length 
     ( J * [p] + [r]  + q  + < 1 | q*q > + [J] + +[1] + [1] + [J] ) * simulationsToPerform
   where 
   J: number of groups
   p: dimension of beta vector
   r: dimension of gamma vector
   q: dimension of alpha vector
   < 1 | q*q >: dimension of tau2 (either a random variable or a random matrix )
   J: augmented variables for t-distributed betas
   1: augmented variable for t-distributed gamma
   1: augmented variable for t-distributed alpha
   J: the xi variables


   < .|. > denotes a choice.
   [.] denotes an optional argument.

   The above dimensions specification also denotes the order in which the simulations are returned:
         beta (by group), gamma, alpha, tau | tau2, 
               t_beta (by group), t_gamma, t_alpha, xi (by group).

   Note: Simulations for tau (not tau2) are returned if tau2 is a scalar,
         otherwise, the matrices tau2 are returned.

*/
void fitBayesianHBM( int * number_groups,
                     int * number_data,
                     int * dim_beta,
                     int * dim_gamma,
                     int * dim_alpha,

                     double * random_data,
                     double * fixed_data,
                     double * second_data,
                     double * data_counts,
                     double * data_trials,

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

                     double * tauDF,
                     double * tauScale,
                     double * tauPower,
                     int * prior_tau_type,

                     double * xi_z0,
                     int * common_xi,

                     int * read_init_point,
                     double * betaInit,
                     double * gammaInit,
                     double * alphaInit,
                     double * tau2Init,

                     int * xiInit_type,
                     double * xiInit,

                     int * burnInLength,
                     int * simulationsToPerform,
                     int * sampleFrequency,

                     int * update_cov_flag,
                     int * print_statistics,
                     double * output_simulations,
                     double * mh_drawing_stats ) throw( rtErr ) 
{
  bool random_effects, fixed_effects, second_effects;
  int i, j, k, dim, number_obs, start_index, simulations_kept;

  CVector ** counts;
  CVector ** trials;
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
  counts = new CVector * [ (int) (*number_groups) ];
  start_index = 0;
  number_obs = 0;
  for ( i = 0; i < (int) (*number_groups); i++ )
  {
    counts[i] = new CVector( (int) number_data[i] );
    for ( j = 0; j < (int) number_data[i]; j++ )
    {
      counts[i]->Val(j) = data_counts[ start_index + j ]; 
    }
    start_index += ((int) number_data[i]);
    number_obs += ((int) number_data[i]);

#ifdef DEBUG1
    Rprintf( "fit: counts[%d] = \n",i );
    counts[i]->Print();
#endif

  }//end for i

  trials = new CVector * [ (int) (*number_groups) ];
  start_index = 0;
  for ( i = 0; i < (int) (*number_groups); i++ )
  {
    trials[i] = new CVector( (int) number_data[i] );
    for ( j = 0; j < (int) number_data[i]; j++ )
    {
      trials[i]->Val(j) = data_trials[ start_index + j ]; 
    }
    start_index += ((int) number_data[i]);

    //printf( "fit: trials[%d] = \n",i );
    //trials[i]->Print();
  }//end for i


  dim = (int) (*dim_beta);
  if ( dim > 0 )
  {
    random_effects = true;
    
    random_predictors = new CMatrix * [ (int) (*number_groups) ];
    start_index = 0;
    
    for ( i = 0; i < (int) (*number_groups); i++ )
    {
      random_predictors[i] = new CMatrix ( counts[i]->Len(), dim );
      for ( k = 0; k < counts[i]->Len(); k++ )
      {
        for ( j = 0; j < dim; j++ )
        {
          random_predictors[i]->Val( k, j ) = random_data[ start_index + k * dim + j ];
        } 
      }//end for j
      start_index += dim * counts[i]->Len();

#ifdef DEBUG1
      Rprintf( "fit: random[%d] = \n",i );
      random_predictors[i]->Print();
#endif

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
      fixed_predictors[i] = new CMatrix ( counts[i]->Len(), dim );
      for ( k = 0; k < counts[i]->Len(); k++ )
      {
        for ( j = 0; j < dim; j++ )
        {
          fixed_predictors[i]->Val( k, j ) = fixed_data[ start_index + k * dim + j ];
        } 
      }//end for j
      start_index += dim * counts[i]->Len();
#ifdef DEBUG1
      Rprintf( "fit: fixed[%d] = \n",i );
      fixed_predictors[i]->Print();
#endif
    }//end for i
  }

  if ( !random_effects && !fixed_effects )
  {
    Rprintf( "fitBayesianHBM: No fixed effects nor random effects provided. This model is not valid.\n" );
    char the_error[] = "fitBayesianHBM: No fixed effects nor random effects provided. This model is not valid.";
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

      //printf( "fit: second[%d] = \n",i );
      //second_predictors[i]->Print();
    }//end for i
  }

  BayesianHierarchicalGlmModel bayes_HGlm;
  
#ifdef DEBUG1
  Rprintf( "counts[0] = \n");
  counts[0]->Print();
#endif
  // Just sets the number of groups and number of obs in each group
  bayes_HGlm.initialize( counts, ((int) (*number_groups) ) );

  if ( random_effects )
  {
  	// Sets the random effects predictors matrix, increases the 
  	// number of parameters, and sets the random effects flag 
  	// to true
    bayes_HGlm.randomEffects( random_predictors );
  }

  if ( fixed_effects )
  {
  	// Sets the fixed effects predictors matrix, increases the 
  	// number of parameters, and sets the fixed effects flag 
  	// to true
    bayes_HGlm.fixedEffects( fixed_predictors );
  }

  if ( second_effects )
  {
    bayes_HGlm.secondStageRandomEffects( second_predictors );
  }


  for ( i = 0; i < ((int) (*number_groups) ); i++ )
  {
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
  // 3 ==> no alpha parameter in model. beta is non-informative
  if ( random_effects )
  {
    if ( ( (int) (*prior_beta_type) ) == 3 )
    {
      bayes_HGlm.betaPriorNonInformative( ((int) (*dim_beta) ) );
    }
    else
    {
      CMatrix p_betaCov( betaCov, ((int) (*dim_beta) ), ((int) (*dim_beta) ), by_column );
  
      if ( ( (int) (*prior_beta_type) ) == 0 )
      {
        bayes_HGlm.betaPrior( &p_betaCov );
      }
      else if ( ( (int) (*prior_beta_type) ) == 1 )
      {
        CVector p_beta0( betamean, ((int) (*dim_beta) ) );
        bayes_HGlm.betaPrior( &p_beta0, &p_betaCov );
      }    
      else if ( ( (int) (*prior_beta_type) ) == 2 )
      {
        CMatrix p_beta0( betamean, ((int) (*dim_beta) ), ((int) (*number_groups) ), by_column );
        bayes_HGlm.betaPrior( &p_beta0, &p_betaCov );
      }    
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
      bayes_HGlm.gammaPrior( &p_gamma0, &p_gammaCov );
    }
    else
    {
      bayes_HGlm.gammaPriorNonInformative( ((int) (*dim_gamma) ) );
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
        bayes_HGlm.alphaPrior( &p_alpha0, &p_alphaCov );
      }
      else if ( ( (int) (*prior_alpha_type) ) == 1 )
      {
        bayes_HGlm.alphaPriorNonInformative( ((int) (*dim_alpha) ) );
      }
    }
  }


  //no non-informative case
  if ( random_effects && ((int) (*prior_beta_type)) != 3 )
  {
    // 0 ==> InvChisq
    // 1 ==> non informative
    // 2 ==> uniform shrinkage
    // 3 ==> du Mouchel
    // 4 ==> InvWishart
    if ( ( (int) (*prior_tau_type) ) == 0 )
    {
      bayes_HGlm.tau2PriorInvChisq( tauDF, tauScale );
    }
    else if ( ( (int) (*prior_tau_type) ) == 1 )
    {
      bayes_HGlm.tau2PriorNonInformative( tauPower );
    }
    else if ( ( (int) (*prior_tau_type) ) == 2 )
    {
      bayes_HGlm.tau2PriorUniformShrinkage( tauScale );
    }
    else if ( ( (int) (*prior_tau_type) ) == 3 )
    {
      bayes_HGlm.tau2PriorDuMouchel( tauScale );
    }
    else if ( ( (int) (*prior_tau_type) ) == 4 )
    {
      CMatrix p_betaCov( tauScale, ((int) (*dim_beta) ), ((int) (*dim_beta) ), by_column );

#ifdef DEBUG1
      Rprintf("fit: setting invWishart. original df = %f. original cov dim is %d, cov is\n",  (*tauDF), ((int) (*dim_beta) ) );
      p_betaCov.Print();
#endif      
      bayes_HGlm.betaCovPriorInvWishart( (*tauDF), &p_betaCov );
    }
  }


  if ( random_effects )
  {
    if ( (*betaDF) > 0 && ((int) (*prior_beta_type)) != 3 )
    {
      bayes_HGlm.betaTPrior( (*betaDF) );
    }

    if ( second_effects && (*alphaDF) > 0 && ((int) (*prior_alpha_type)) != 1 )
    {
      bayes_HGlm.alphaTPrior( (*alphaDF) );
    }
  }

  if ( fixed_effects )
  {
    if ( (*gammaDF) > 0 && ((int) (*prior_gamma_type)) != 1 )
    {
      bayes_HGlm.gammaTPrior( (*gammaDF) );
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

#ifdef DEBUG1
        Rprintf( "BHL: init alpha \n" );
        init_alpha.Print();
#endif

        bayes_HGlm.samplerAlphaInitialPoint( init_alpha );

      }
      //this is only valid for prior parameters, not for initial points
      //else if ( ((int) (*prior_beta_type) ) == 1 )
      //{
      //  CVector init_beta( betaInit, ((int) (*dim_beta) ) );
      //  bayes_HGlm.samplerBetaInitialPoint( init_beta );
      //}
      else //if ( ((int) (*prior_beta_type) ) == 2 )
      {
        CMatrix init_beta( betaInit, ((int) (*dim_beta) ), ((int) (*number_groups) ), by_column );

#ifdef DEBUG1
        Rprintf( "BHL: init beta \n" );
        init_beta.Print();
#endif

        bayes_HGlm.samplerBetaInitialPoint( init_beta );
      }
    
      if ( ((int) (*prior_beta_type) ) != 3 )
      {
        if ( ((int) (*prior_tau_type) ) != 4 )
        {
          bayes_HGlm.samplerTau2InitialPoint( tau2Init );
        }
        else
        {
          CMatrix init_tau2( tau2Init, ((int) (*dim_beta) ), ((int) (*dim_beta) ), by_column );
#ifdef DEBUG1
          Rprintf( "BHL: init random var \n" );
          init_tau2.Print();
#endif

          bayes_HGlm.samplerTau2InitialPoint( init_tau2 );
        }
      }
    }

    if ( fixed_effects )
    {
      CVector init_gamma( gammaInit, ((int) (*dim_gamma) ) );

#ifdef DEBUG1
      Rprintf( "BHL: init gamma \n" );
      init_gamma.Print();
#endif

      bayes_HGlm.samplerGammaInitialPoint( init_gamma );
    }
  }//end if initial points provided
  else
  {
    bayes_HGlm.samplerDefaultInitialPoint();
  }


  //connect with the Binomial part
  BayesianBinomialLogitLinkModel bayes_BinomialLL;

  bayes_BinomialLL.initialize( counts, trials, ((int) (*number_groups) ) );

  for ( i = 0; i < ((int) (*number_groups)); i++ )
  {
    delete counts[i];
    delete trials[i];
  }  
  delete [] counts;
  delete [] trials;

  //prior for Xi in Poisson part
  // common_xi
  // 0 ==> independent xis; different hyper prior parameters
  // 1 ==> independent xis; same hyper prior parameters
  // 2 ==> common xi; same hyper prior parameters
  if ( ( (int) (*common_xi) ) == 0 || ( (int) (*common_xi) ) == 1 )
  {
    bayes_BinomialLL.xiType( 0 );   // group-specific overdispersion parameter case

    if ( ( (int) (*common_xi) ) == 0 )
    {
      bayes_BinomialLL.logXiPrior( xi_z0 );
    }
  }
  else if( ( (int) (*common_xi) ) == 2 )
  {
    bayes_BinomialLL.xiType( 1 );   // common overdispersion parameter case
  }
  else {
  	bayes_BinomialLL.xiType( 2 );  // no overdispersion case
  }

  if ( ( ( (int) (*common_xi) ) == 1 ) || 
  	   ( ( (int) (*common_xi) ) == 2 ) )
  {
    #ifdef DEBUG1
      Rprintf( "*xi_z0 = %f\n", *xi_z0 );
    #endif
    bayes_BinomialLL.logXiPrior( (*xi_z0) );
  }

  if ( ( (int) (*read_init_point) ) &&
  	   ( ( (int) (*common_xi) ) != 3 ) )  // overdispersion parameter case
  {
    // 0 ==> same starting point for all groups
    // 1 ==> different starting point for each group
    if ( ((int) (*xiInit_type)) == 0 )
    {
      bayes_BinomialLL.samplerXiInitialPoint( (*xiInit) );
    }
    else
    {
#ifdef DEBUG1
      Rprintf( "xiInit = %f\n", *xiInit );
#endif
      CVector init_xi( xiInit, ((int) (*number_groups)) );
      bayes_BinomialLL.samplerXiInitialPoint( init_xi );
    }  
  }  
  else
  {
    bayes_BinomialLL.samplerDefaultInitialPoint();
  }


  //connect both parts
  bayes_HGlm.glmModel( &bayes_BinomialLL );

  if ( ((int) (*update_cov_flag)) == 1 )
  {
    bayes_HGlm.setUpdateProposalCovariance( true );
  }
  else
  {
    bayes_HGlm.setUpdateProposalCovariance( false );
  }



  //get ready to start Gibbs sampler
  MetropolisHastings sampler( (int) (*burnInLength), 
                              (int) (*simulationsToPerform),
                              (int) (*sampleFrequency) );

  if ( (*print_statistics) )
  {
    sampler.saveStatistics( true );
  }
  

  //start MetropolisHastings sampler
  sampler.doMetropolisHastings( &bayes_HGlm );

  simulations_kept = sampler.simulationsKept();

  if ( (*print_statistics) )
  {
    sampler.printDrawingStats( mh_drawing_stats );   
    //sampler.printDrawingStats();   
  }

  #ifdef DEBUG1
    Rprintf( "fitBayesHPM: total simulations = %d\n", simulations_kept );
  #endif

  //get the output simulation samples 
  bayes_HGlm.simulationsToArrayWithMeans( output_simulations, simulations_kept );

  /* update random seed */
  PutRNGstate();

}//end


}
    

