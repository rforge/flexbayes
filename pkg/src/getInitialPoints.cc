#include "R.h"
#include "Vector.h"
#include "Matrix.h"

#include "Const.h"
#include "DistributionExtended.h"
#include "DistributionMixture.h"

extern "C" {

/* creates a set of initial points for starting
   a simulation (Gibbs or exact 
*/
void getInitialPointsBayesLM( long * number_data,
                              long * dim_Cov,
                              double * data_response,
                              double * data_predictors,
                              double * prior_sigmaDF,
                              double * prior_sigmaScale,
                              double * betamean,
                              double * betaCov,
                              double * betaDF,
                              double * beta_props,
                              long * number_mixtures,
                              long * number_draws,
                              long * mixture_with_MLE,
                              double * output_initial_points )
{
  int i;
  double local_sigmaDF, beta_dfreedom;
  CVector * response, * beta0;
  CMatrix * predictors, * prior_betaCov;
  InvChisqDistribution * sigma2_prior, * sigma2_mle;
  StudentTDistribution * beta_prior, * beta_mle;
  DistributionMixture * mixture;
  DistributionMixture * sigma_mixture;
  StudentTDistribution ** betas_t;
  DistributionMixture * beta_prior_mix;  

  /*set random seed*/
  GetRNGstate();

  //the data  
  response = new CVector( data_response, (int) (*number_data) );
  predictors = new CMatrix( data_predictors, (int) (*number_data), (int) (*dim_Cov) );


  //the beta prior
  if ( (*number_mixtures) <= 1 )
  {
    beta0 = new CVector( betamean, (int) (*dim_Cov) );
    prior_betaCov = new CMatrix( betaCov, (int) (*dim_Cov), (int) (*dim_Cov) );

    if ( (*betaDF) > 3 || (*betaDF) <= 0 )
    {
      beta_dfreedom = 3.0;
    }
    else
    {
      beta_dfreedom = (*betaDF);
    }
    //set over-disperse prior distributions
    //beta

    //printf( "getInitialPoints: beta mean and cov from prior:\n" );
    //beta0->Print();
    //prior_betaCov->Print();

    beta_prior = new StudentTDistribution( (int) (*dim_Cov), beta_dfreedom ); 
    beta_prior->setLocation( beta0 );
    beta_prior->setScaleMatrix( prior_betaCov );
  }
  else
  {
    CVector ** beta0s;
    CMatrix ** prior_betaCovs;
    double * beta_vector, * betaCov_matrix;
    int i, j, k, start_vector, start_matrix;

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

    beta_prior_mix = new DistributionMixture( (int) (*number_mixtures) );
    betas_t = new StudentTDistribution * [ (int) (*number_mixtures) ];
    for ( i = 0; i < (int) (*number_mixtures); i++ )
    {
      //printf(" betaDF[%d] = %f, props = %f \n", i, betaDF[i], beta_props[i] );
      //beta0s[i]->Print();
      //prior_betaCovs[i]->Print();

      if ( betaDF[i] > 3 || betaDF[i] <= 0 )
      {
        beta_dfreedom = 3.0;
      }
      else
      {
        beta_dfreedom = betaDF[i];
      }
      betas_t[i] = new StudentTDistribution( beta0s[i]->Len(), beta_dfreedom );
      betas_t[i]->setLocation( beta0s[i] );
      betas_t[i]->setScaleMatrix( prior_betaCovs[i] );

      beta_prior_mix->set( i, betas_t[i], beta_props[i] );
    }//end for i


    for ( i = 0; i < (int) (*number_mixtures); i++ )
    {
      delete beta0s[i];
      delete prior_betaCovs[i];
    }
    delete [] beta0s;
    delete [] prior_betaCovs;      

    delete [] beta_vector;
    delete [] betaCov_matrix;
  }//end if mixture


  //now check for degrees of freedom for InvChisq
  if ( (*prior_sigmaDF) == 0 )
  {
    //the non-informative prior case.
    //use a small number of degrees of freedom
    local_sigmaDF = NON_INFORMATIVE_SIGMA_DF;
  }
  else
  {
    local_sigmaDF = (*prior_sigmaDF);
  }

  //sigma2
  sigma2_prior = new InvChisqDistribution( local_sigmaDF );
  sigma2_prior->setScale( (*prior_sigmaScale) );

  //printf( "mixture = %d \n", (int) (*mixture_with_MLE) );
  //fflush(stdout);

  if ( ( (int) (*mixture_with_MLE) ) == 1 )
  {
    //now get MLE's for the parameters
    CMatrix xTx_inverse( predictors->xTransposedX().inverse() );
    CVector betaMLE( xTx_inverse * ( predictors->T() * (*response) ) );
    double sigma2MLE = ( (*response) - ( (*predictors) * betaMLE ) ) * ( (*response) - ( (*predictors) * betaMLE ) );
    sigma2MLE /= ( (double) (*number_data) );

    //printf( "sigma2MLE = %f\n", sigma2MLE );
    //printf( "xTx inverse is :\n" ); xTx_inverse.Print(); fflush( stdout );

    xTx_inverse.multiplyByScalar( sigma2MLE );//this is Sigma_MLE
    //beta
    beta_mle = new StudentTDistribution( (int) (*dim_Cov), beta_dfreedom ); 
    beta_mle->setLocation( &betaMLE );
    beta_mle->setScaleMatrix( &xTx_inverse );


    //printf( "betaMLE is :\n"); betaMLE.Print();
    //printf( "COvMLE is :\n" ); xTx_inverse.Print(); fflush( stdout );

    //sigma2
    sigma2_mle = new InvChisqDistribution( local_sigmaDF );
    sigma2_mle->setScale( sigma2MLE );

    //ready to sample (compute proportions first)

    //sample beta
    mixture = new DistributionMixture( 2 );
    if ( (*number_mixtures) <= 1 )
    {
      mixture->set( 0, beta_prior, 0.5 );
    }
    else
    {
      mixture->set( 0, beta_prior_mix, 0.5 );
    }
    mixture->set( 1, beta_mle, 0.5 );

    //sample sigma2
    sigma_mixture = new DistributionMixture( 2 );
    sigma_mixture->set( 0, sigma2_prior, 0.5 );
    sigma_mixture->set( 1, sigma2_mle, 0.5 );
  
  }
  else if ( ( (int) (*mixture_with_MLE) ) == 0 )
  {
    //sample beta
    mixture = new DistributionMixture( 1 );
    if ( (*number_mixtures) <= 1 )
    {
      mixture->set( 0, beta_prior, 1.0 ); 
    }
    else
    {
      mixture->set( 0, beta_prior_mix, 1.0 );
    }

    //sample sigma2
    sigma_mixture = new DistributionMixture( 1 );
    sigma_mixture->set( 0, sigma2_prior, 1.0 );

  }
  else if ( ( (int) (*mixture_with_MLE) ) == 2 )
  {
    //now get MLE's for the parameters
    CMatrix xTx_inverse( predictors->xTransposedX().inverse() );
    CVector betaMLE( xTx_inverse * ( predictors->T() * (*response) ) );
    double sigma2MLE = ( (*response) - ( (*predictors) * betaMLE ) ) * ( (*response) - ( (*predictors) * betaMLE ) );
    sigma2MLE /= ( (double) (*number_data) );

    xTx_inverse.multiplyByScalar( sigma2MLE );//this is Sigma_MLE
    //beta
    beta_mle = new StudentTDistribution( (int) (*dim_Cov), beta_dfreedom ); 
    beta_mle->setLocation( &betaMLE );
    beta_mle->setScaleMatrix( &xTx_inverse );

    //sigma2
    sigma2_mle = new InvChisqDistribution( local_sigmaDF );
    sigma2_mle->setScale( sigma2MLE );

    //ready to sample (compute proportions first)

    //sample beta
    mixture = new DistributionMixture( 1 );
    mixture->set( 0, beta_mle, 1.0 );

    //sample sigma2
    sigma_mixture = new DistributionMixture( 1 );
    sigma_mixture->set( 0, sigma2_mle, 1.0 );
  
  }  

  //printf("will go to stratified sampling\n"); fflush(stdout);
  DistributionParameter * samples = mixture->stratifiedDrawFromMixture( ( (int) (*number_draws) ) );
  int start_index = 0;
  //printf( "the beta samples:\n" ); fflush(stdout);
  for ( i = 0; i < (int) (*number_draws); i++ )
  {
    //samples[i].Print(); fflush(stdout);
    samples[i].toArray( output_initial_points, start_index );
    start_index += samples[i].length();
  }
  //fflush(stdout);
  delete [] samples;

  samples = sigma_mixture->stratifiedDrawFromMixture( ( (int) (*number_draws) ) );
  //printf( "the sigma2 samples:\n" ); fflush(stdout);
  for ( i = 0; i < (int) (*number_draws); i++ )
  {
    //samples[i].Print(); fflush(stdout);
    samples[i].toArray( output_initial_points, start_index );
    start_index += samples[i].length();
  }
  //fflush(stdout);
  delete [] samples;
  delete sigma_mixture;
  delete mixture;
  
  if ( ( ( (int) (*mixture_with_MLE) ) == 1 ) || ( ( (int) (*mixture_with_MLE) ) == 2 ) )
  {
    delete beta_mle;
    delete sigma2_mle;
  }

  if ( (*number_mixtures) <= 1 )
  {
    delete beta_prior;
    delete sigma2_prior;
  }
  else
  {
    for ( i = 0; i < beta_prior_mix->numberOfMixtures(); i++ )
    {
      delete betas_t[i];
    }
    delete [] betas_t;
    delete beta_prior_mix;
  }

  /* update random seed */
  PutRNGstate();

}//end
}

