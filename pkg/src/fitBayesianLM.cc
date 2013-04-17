#include "R.h"
#include "Vector.h"
#include "Matrix.h"

#include "BayesianLinearModel.h"
#include "Gibbs.h"

//#include <stdio.h>

extern "C" {
/*
  Computes a Bayes linear model:
    y = X beta + error

    (y = data_response is a vector of length number_data;
     X = data_predictors is a matrix of dimensions number_data * dimCov )

  with priors:
   1) the semi-conjugate case
      beta ~ Normal( mean = betamean, covariance = betaCov )
      sigma2 ~ InvChisq( df = prior_sigmaDF, scale = prior_sigmaScale * prior_sigmaScale )

   or
   2) the non-conjugate case
      beta ~ t( df = betaDF, mean = betamean, covariance = betaCov )
      sigma2 ~ InvChisq( df = prior_sigmaDF, scale = prior_sigmaScale * prior_sigmaScale )

   or
   3) the non-informative case
      (beta, sigma2) ~ InvChisq( df = prior_sigmaDF, scale = prior_sigmaScale * prior_sigmaScale )
                       if also prior on sigma2 is non-informative then
                       (beta, sigma2) ~ InvChisq( df = 0, scale = 1 )
                       (prior_sigmaDF should be set to 0 to signal a non-informative prior for sigma2)

   where the errors are either
   1) normally distributes, error_i ~ Normal( mean = 0, variance = sigma2 )
   or
   2) t-distributed: error_i ~ t( df = degreesOfFreedom_likelihood, location = 0, scale = sigma2 )

   Computations are based on Gibbs sampler drawings from the full conditionals.
   burnInLength: number of initial Gibbs sampler iterations to discard.
   simulationsToPerform: number of Gibbs sampler drawings to keep as output.
   sampleFrequency: only iterations that are multiples of sampleFrequency are to be kept.

   The simulationsToPerform output simulations are returned in the array output_simulations.
   This should be an array of length (dimCov + 1) * simulationsToPerform.
   This array is such that the simulations for the k-th coefficient are stored in
   output_simulations[ k * simulationsToPerform + i ], i = 0,..., simulationsToPerform - 1,
                                                       k = 0,..., dimCov - 1.

   Simulations for sigma (not sigma2) are returned in
   output_simulations[ dimCov * simulationsToPerform + i ], i = 0,..., simulationsToPerform - 1.

*/
void fitBayesianLM( long * number_data,
                    long * dim_Cov,
                    double * data_response,
                    double * data_predictors,
                    long * dim_error_Cov,
                    double * error_Cov,
                    double * degreesOfFreedom_likelihood,
                    double * prior_sigmaDF,
                    double * prior_sigmaScale,
                    long * prior_Type,
                    double * betamean,
                    double * betaCov,
                    double * betaDF,
		    double * proportions,
                    long * number_mixtures,
                    long * burnInLength,
                    long * simulationsToPerform,
                    long * sampleFrequency,
                    long * read_init_point,
                    double * betaInit,
                    double * sigmaInit,
                    long * sampler_type,
                    long * print_statistics,
                    double * output_simulations,
                    double * gibbs_drawing_stats,
                    double * mixture_drawing_stats )
{
  bool identityCov, mixture_prior;
  int i, simulations_kept;

  CVector * response;
  CMatrix * predictors;

  /*set random seed*/
  GetRNGstate();

  //the data  
  response = new CVector( data_response, (int) (*number_data) );
  predictors = new CMatrix( data_predictors, (int) (*number_data), (int) (*dim_Cov) );


  //construct appropriate error covariance matrix (if needed)
  //check if error_Cov is not the identity
  if ( ( (int) (*dim_error_Cov) ) == 1 )
  {
    //a scalar times the identity covariance
    if ( (*error_Cov) != 1.0 )
    {
      //transform predictors and response
      response->multiplyByScalar( 1.0 / sqrt( (*error_Cov) ) );
      predictors->multiplyByScalar( 1.0 / sqrt( (*error_Cov) ) );
    }
  }
  else if (  (*dim_error_Cov) == (*number_data) )
  {
    //a diagonal covariance
    identityCov = true;
    i = 0;
    while ( identityCov && i < ( (int) (*number_data) ) )
    {
      if ( error_Cov[i] != 1.0 )
      {
        identityCov = false;
      }
      i++;
    }

    if ( !identityCov )
    {
      //transform predictors
      for ( i = 0; i < ( (int) (*number_data) ); i++ )
      {
        error_Cov[i] = 1.0 / sqrt( error_Cov[i] );
      }

      CVector errorCov( error_Cov, (int) (*number_data) );

      response->setToWeighted( errorCov );
      for ( i = 0; i < ( (int) (*dim_Cov) ); i++ )
      {
        predictors->setToWeightedColumn( i, errorCov );
      }
    }
  }
  else
  {
    //a matrix for covariance: transform predictors
    CMatrix * working_Cov;
    CVector * cholLT;

    working_Cov = new CMatrix( error_Cov, (int) (*number_data), (int) (*number_data) );
    cholLT = new CVector( (int) ( (*number_data) * ( (*number_data) + 1) ) / 2 );

    (*cholLT) = working_Cov->choleskyDecomposition();
    working_Cov->assignInverseOfLowerTriangular( (*cholLT) );
    delete cholLT;

    CVector weighted_vector( (int) (*number_data) );

    //do the transformation
    weighted_vector = (*working_Cov) * (*response);
    //assign the transformation
    (*response) = weighted_vector;

    for ( i = 0; i < ( (int) (*dim_Cov) ); i++ )
    {
      //do the transformation
      weighted_vector = (*working_Cov) * ( predictors->getColumn(i) );
      //assign the transformation
      predictors->setColumn( i, weighted_vector );
    }
    delete working_Cov;
  }



  //find out what kind of model this is 
  mixture_prior = false;

  BayesianLinearModel bayes_linear;

  //the beta prior
  //mixture case: normal mixture (normal lklhd), or t-mixture (normal lklhd), 
  //or normal mixture (t-lklhd), or t mixture (t-lklhd)
  if ( ( (int) (*prior_Type) ) == 7 || ( (int) (*prior_Type) ) == 8 || 
       ( (int) (*prior_Type) ) == 9 || ( (int) (*prior_Type) ) == 10 ||
       ( (int) (*prior_Type) ) == 11 || ( (int) (*prior_Type) ) == 12 || 
       ( (int) (*prior_Type) ) == 13 || ( (int) (*prior_Type) ) == 14 ||
       ( (int) (*prior_Type) ) == 15 )
  {
    //the mixture case
    if ( ( (int) (*number_mixtures) ) > 1 )
    {
      CVector ** beta0s;
      CMatrix ** prior_betaCovs;
      double * beta_vector, * betaCov_matrix;
      int i, j, k, start_vector, start_matrix;
      char * mixture_type;

      beta0s = new CVector * [ (int) (*number_mixtures) ];
      prior_betaCovs = new CMatrix * [ (int) (*number_mixtures) ];
      beta_vector = new double [ (int) (*dim_Cov) ];
      betaCov_matrix = new double [ (int) ( (*dim_Cov) * (*dim_Cov) ) ];

      start_vector = 0;
      start_matrix = 0;
      for ( i = 0; i < (int) (*number_mixtures); i++ )
      {
        for ( j = 0; j < (int) (*dim_Cov); j++ )
	{
          beta_vector[j] = betamean[ start_vector + j ];
          for ( k = 0; k < (int) (*dim_Cov); k++ )
	  {
            betaCov_matrix[ j * ((int) (*dim_Cov) ) + k ] = betaCov[ start_matrix + j * ((int) (*dim_Cov) ) + k ];
          }//end for k
        }//end for j
        start_vector += ((int) (*dim_Cov) );
        start_matrix += ((int) ( (*dim_Cov) * (*dim_Cov) ) );

        beta0s[i] = new CVector( beta_vector, (int) (*dim_Cov) );
        prior_betaCovs[i] = new CMatrix( betaCov_matrix, (int) (*dim_Cov), (int) (*dim_Cov) );
      }//end for i

      if ( ((int) (*prior_Type) ) == 7 || ((int) (*prior_Type) ) == 8 ||
           ((int) (*prior_Type) ) == 11 || ((int) (*prior_Type) ) == 12 ||
           ((int) (*prior_Type) ) == 15 )
      {
        mixture_type = new char [7];
        sprintf( mixture_type, "normal" );
      }
      else
      {
        mixture_type = new char [2];
        sprintf( mixture_type, "t" );
      }

      bayes_linear.initialize( response, predictors, (int) (*number_mixtures), 
                               beta0s, prior_betaCovs, proportions, betaDF, mixture_type,
                               (*prior_sigmaScale), (*prior_sigmaDF),
                               (*print_statistics) );

      for ( i = 0; i < (int) (*number_mixtures); i++ )
      {
        delete beta0s[i];
        delete prior_betaCovs[i];
      }
      delete [] beta0s;
      delete [] prior_betaCovs;  

      delete [] beta_vector;
      delete [] betaCov_matrix;

      mixture_prior = true;
    }//end if number mixtures > 1

  }//end if mixture
  else
  {
    CVector * beta0;
    CMatrix * prior_betaCov;

    beta0 = new CVector( betamean, (int) (*dim_Cov) );
    prior_betaCov = new CMatrix( betaCov, (int) (*dim_Cov), (int) (*dim_Cov) );

    //printf("prior beta is \n"); prior_betaCov->Print();

    if ( ( (int) (*prior_Type) ) == 0 || ( (int) (*prior_Type) ) == 3 )
    {
      //the non-informative beta prior case
      bayes_linear.nonInformative( response, predictors,
                                   (*prior_sigmaScale),
                                   (*prior_sigmaDF) );
    }
    else
    {
      //the normal-likelihood normal-beta-prior case ( (*prior_Type) == 1 )
      bayes_linear.initialize( response, predictors,
                               beta0, prior_betaCov,
                               (*prior_sigmaScale),
                               (*prior_sigmaDF) );
    }

    delete beta0;
    delete prior_betaCov;
  }//end no mixture


  //the conjugate case
  if ( ( (int) (*prior_Type) ) == 6 || ((int) (*prior_Type) ) == 15 )
  {
    bayes_linear.conjugatePrior();
  }
  //the normal-likelihood t-beta-prior case
  else if ( ( (int) (*prior_Type) ) == 2 )
  {
    bayes_linear.t_beta_prior( (*betaDF) );
  }
  //the t-likelihood normal-beta-prior case or the t-likelihood non-informative prior case
  //or t-likelihood mixture prior case
  else if ( ( (int) (*prior_Type) ) == 4 || ( (int) (*prior_Type) ) == 3 ||
            ( (int) (*prior_Type) ) == 11 || ( (int) (*prior_Type) ) == 12 ||
            ( (int) (*prior_Type) ) == 13 || ( (int) (*prior_Type) ) == 14 ) 
  {
    bayes_linear.t_likelihood( (*degreesOfFreedom_likelihood) );
  }
  //the t-likelihood t-beta-prior case
  else if ( ( (int) (*prior_Type) ) == 5 )
  {
    bayes_linear.t_likelihood_and_beta_prior( (*degreesOfFreedom_likelihood), (*betaDF) );
  }


  //get ready for sampling
  if ( (int) (*read_init_point) )
  {
    CVector initBeta( betaInit, (int) (*dim_Cov) );
    bayes_linear.samplerBetaInitialPoint( initBeta );
    bayes_linear.samplerSigma2InitialPoint( (*sigmaInit) );
  }
  else
  {
    bayes_linear.samplerDefaultInitialPoint();
  }
  
  if ( ( bayes_linear.IsPriorConjugate() 
         || ( bayes_linear.IsPriorNonInformative() && !bayes_linear.tLikelihood() ) )
       && (*sampler_type) == 1 )
  {
    //get ready for exact sampling
    //run sampler
    bayes_linear.exactSampler( (int) (*simulationsToPerform) );
    simulations_kept = (int) (*simulationsToPerform);

  }
  else
  {
    //get ready to start Gibbs sampler
    Gibbs sampler( (int) (*burnInLength), 
                   (int) (*simulationsToPerform),
                   (int) (*sampleFrequency) );

    if ( (*print_statistics) )
    {
      sampler.saveStatistics( true );
    }

    //so far no particular initialization points for Gibbs sampler are given
    //start Gibbs sampler
    sampler.doGibbs( &bayes_linear );

    simulations_kept = sampler.simulationsKept();

    if ( (*print_statistics) )
    {
      sampler.printDrawingStats( gibbs_drawing_stats );   
    }
  }

  //get the output simulation samples 
  bayes_linear.simulationsToArray( output_simulations, simulations_kept );
  if ( mixture_prior && (*print_statistics) )
  {
    bayes_linear.saveMixtureDrawingStatistics( mixture_drawing_stats );
  }

  /* update random seed */
  PutRNGstate();

  delete response;
  delete predictors;

}//end


}
    
