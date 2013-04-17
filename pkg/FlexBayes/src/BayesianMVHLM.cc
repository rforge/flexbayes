#include <math.h>
#include <string.h>
#include "R.h"

#include "Const.h"
#include "DistributionParameter.h"
#include "BayesianMVHLM.h"


extern int randomStart( int n );

/* constructor.
   returns: an empty Distribution of type BayesianLinearModel
*/
BayesianMVHLM::BayesianMVHLM()
{

  emptyModel();

}//end


void BayesianMVHLM::emptyModel()
{
  response = NULL;
  random_predictors = NULL;
  fixed_predictors = NULL;
  beta_predictors = NULL;
  
  //distributions
  //first stage
  gamma = NULL;
  tau2_gamma = NULL;

  beta = NULL;
  //t prior
  tau2_betas = NULL;    

  sigma2 = NULL;
  sigma2ICS = NULL;
  sigma2PNIP = NULL;
  sigma2NIP = NULL;
  sigma2InvWishart = NULL;

  //second stage
  alpha = NULL;
  tau2 = NULL;
  tau2ICS = NULL;
  tau2PNIP = NULL;
  tau2NIP = 0;
  tau2InvWishart = NULL;
  tau2_alpha = NULL;

  sigma2_first_draw = NULL;
  tau2_first_draw = NULL;

  //t-likelihood
  tau2_groups = NULL;
  tau2_errors = NULL;

  //missing data 
  keep_missing_data = false;
  number_missing_response = NULL;
  missing_response_components = NULL;
  number_missing_random_predictors = NULL;
  missing_random_predictors_components = NULL;
  number_missing_fixed_predictors = NULL;
  missing_fixed_predictors_components = NULL;

  number_of_variables = 0;
  number_of_groups = 0;
  number_of_observations = 0;

  t_lkhd = false;
  t_beta = false;
  t_gamma = false;
  t_alpha = false;
  t_group_lkhd = false;
  non_informative_alpha = false;
  non_informative_gamma = false;
  invWishart_betaCov = false;
  keep_tau2_error = true;

  random_effects = false;
  fixed_effects = false;
  second_stage_effects = false;

  common_sigma = true;
  invChisq_sigma2 = false;
  duMouchel_sigma2 = false;
  uniformShrinkage_sigma2 = false;
  properNonInfoPrior_sigma2 = false;
  nonInformativePower_sigma2 = false;
  invWishart_sigma2 = false;
  known_sigma2 = false;

  invChisq_tau2 = false;
  duMouchel_tau2 = false;
  uniformShrinkage_tau2 = false;
  properNonInfoPrior_tau2 = false;
  nonInformativePower_tau2 = false;

  //simulation samples
  number_of_simulations = 0;
  simulated_beta = NULL;
  simulated_gamma = NULL;
  simulated_alpha = NULL;
  simulated_sigma2 = NULL;
  simulated_tau2 = NULL;
  simulated_tau2_beta = NULL;
  simulated_tau2_gamma = NULL;
  simulated_tau2_alpha = NULL;
  simulated_tau2_errors = NULL;
  simulated_tau2_groups = NULL;
  simulated_missingR = NULL;
  simulated_missingRP = NULL;
  simulated_missingFP = NULL;

  //working matrices for Gibbs sampler
  //working matrices for Gibbs sampler
  //random effects
  xTx = NULL;
  xTy = NULL;
  //fixed effects
  sMTm = NULL;
  sMTy = NULL;
  mTm = NULL;
  mTy = NULL;
  //second stage
  zTz = NULL;
  zTb = NULL;
  zTvIz = NULL;
  zTvIb = NULL;

  final_group_index = NULL;

  tau2_error_weights = NULL;
  tau2_group_weights = NULL;
  residuals = NULL;
  second_residuals = NULL;

  //array of distributions map (maps distributions names to actual objects)
  distr_map = NULL;


  //missing data priors
  randomEffects_means = NULL;
  randomEffects_vars = NULL;
  fixedEffects_means = NULL;
  fixedEffects_vars = NULL;

}//end



BayesianMVHLM::~BayesianMVHLM()
{
  int i, j;

  if ( number_missing_response != NULL )
  {
    delete number_missing_response;
  }

  if ( missing_response_components != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      for ( j = 0; j < response[i]->Col(); j++ )
      {
        delete missing_response_components[i][j];
      }
      delete [] missing_response_components[i];
    }
    delete [] missing_response_components;
  }


  if ( number_missing_random_predictors != NULL )
  {
    delete number_missing_random_predictors;
  }

  if ( missing_random_predictors_components != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      for ( j = 0; j < response[i]->Col(); j++ )
      {
        delete missing_random_predictors_components[i][j];
      }
      delete [] missing_random_predictors_components[i];
    }
    delete [] missing_random_predictors_components;
  }


  if ( number_missing_fixed_predictors != NULL )
  {
    delete number_missing_fixed_predictors;
  }

  if ( missing_fixed_predictors_components != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      for ( j = 0; j < response[i]->Col(); j++ )
      {
        delete missing_fixed_predictors_components[i][j];
      }
      delete [] missing_fixed_predictors_components[i];
    }
    delete [] missing_fixed_predictors_components;
  }



  if ( random_effects )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      for ( j = 0; j < response[i]->Col(); j++ )
      {
        delete random_predictors[i][j];
      }
      delete [] random_predictors[i];
    }
    delete [] random_predictors;
  }

  if ( fixed_effects )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      for ( j = 0; j < response[i]->Col(); j++ )
      {
        delete fixed_predictors[i][j];
      }
      delete [] fixed_predictors[i];
    }
    delete [] fixed_predictors;
  }

  if ( beta_predictors != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete beta_predictors[i];
    }
    delete [] beta_predictors;    
  }


  if ( response != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete response[i];
    }
    delete [] response;
  }

  
  //distributions
  //first stage
  if ( gamma != NULL )
  {
    delete gamma;
  }

  if ( tau2_gamma != NULL )
  {
    delete tau2_gamma;
  }


  if ( beta != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete beta[i];
    }
    delete [] beta;
  }


  //t prior
  if ( tau2_betas != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete tau2_betas[i];
    }
    delete [] tau2_betas;
  }

  if ( sigma2 != NULL )
  {
    delete sigma2;
  }

    if ( sigma2NIP != NULL )
  {
    delete [] sigma2NIP;
  }

  if ( sigma2ICS != NULL )
  {
    if ( common_sigma )
    {
      delete sigma2ICS[0];
    }
    else 
    {
      for ( i = 0; i < number_of_groups; i++ )
      {
        delete sigma2ICS[i];
      }
    }
    delete [] sigma2ICS;
  }

  if ( sigma2PNIP != NULL )
  {
    if ( common_sigma )
    {
      delete sigma2PNIP[0];
    }
    else 
    {
      for ( i = 0; i < number_of_groups; i++ )
      {
        delete sigma2PNIP[i];
      }
    }
    delete [] sigma2PNIP;
  }


  if ( sigma2InvWishart != NULL )
  {
    if ( common_sigma )
    {
      delete sigma2InvWishart[0];
    }
    else 
    {
      for ( i = 0; i < number_of_groups; i++ )
      {
        delete sigma2InvWishart[i];
      }
    }
    delete [] sigma2InvWishart;
  }


  //second stage
  if ( alpha != NULL )
  {
    delete alpha;
  }

  if ( tau2 != NULL )
  {
    delete tau2;
  }

  if ( tau2ICS != NULL )
  {
    delete tau2ICS;
  }

  if ( tau2PNIP != NULL )
  {
    delete tau2PNIP;
  }

  if ( tau2InvWishart != NULL )
  {
    delete tau2InvWishart;
  }

  if ( tau2_alpha != NULL )
  {
    delete tau2_alpha;
  }

  if ( sigma2_first_draw != NULL )
  {
    if ( common_sigma )
    {
      delete sigma2_first_draw[0];
    }
    else 
    {
      for ( i = 0; i < number_of_groups; i++ )
      {
        delete sigma2_first_draw[i];
      }
    }
    delete [] sigma2_first_draw;
  }

  if ( tau2_first_draw != NULL )
  {
    delete tau2_first_draw;
  }

  //t-likelihood
  if ( tau2_groups != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete tau2_groups[i];
    }
    delete [] tau2_groups;
  }

  if ( tau2_errors != NULL )
  {
    for ( i = 0; i < number_of_observations; i++ )
    {
      delete tau2_errors[i];
    }
    delete [] tau2_errors;
  }



  //simulation samples
  if ( simulated_beta != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete simulated_beta[i];
    }
    delete [] simulated_beta;
  }

  if ( simulated_gamma != NULL )
  {
    delete simulated_gamma;
  }

  if ( simulated_alpha != NULL )
  {
    delete simulated_alpha;
  }

  if ( simulated_sigma2 != NULL )
  {
    if ( common_sigma )
    {
      for ( i = 0; i < number_of_simulations; i++ )
      {
        delete simulated_sigma2[i][0];
        delete [] simulated_sigma2[i];
      }
    }
    else
    {
      for ( i = 0; i < number_of_simulations; i++ )
      {
        for ( j = 0; j < number_of_groups; j++ )
	{
          delete simulated_sigma2[i][j];
        }
        delete [] simulated_sigma2[i];
      }
    }
    delete [] simulated_sigma2;
  }

  if ( simulated_tau2 != NULL )
  {
    for ( i = 0; i < number_of_simulations; i++ )
    {
      delete simulated_tau2[i];
    }
    delete [] simulated_tau2;
  }

  if ( simulated_tau2_beta != NULL )
  {
    delete simulated_tau2_beta;
  }

  if ( simulated_tau2_gamma != NULL )
  {
    delete simulated_tau2_gamma;
  }

  if ( simulated_tau2_alpha != NULL )
  {
    delete simulated_tau2_alpha;
  }

  if ( simulated_tau2_errors != NULL )
  {
    delete simulated_tau2_errors;
  }

  if ( simulated_tau2_groups != NULL )
  {
    delete simulated_tau2_groups;
  }

  if ( simulated_missingR != NULL )
  {
    for ( i = 0; i < number_of_observations; i++ )
    {
      delete simulated_missingR[i];
    }
    delete [] simulated_missingR;
  }

  if ( simulated_missingRP != NULL )
  {
    for ( i = 0; i < number_of_observations; i++ )
    {
      delete simulated_missingRP[i];
    }
    delete [] simulated_missingRP;
  }

  if ( simulated_missingFP != NULL )
  {
    for ( i = 0; i < number_of_observations; i++ )
    {
      delete simulated_missingFP[i];
    }
    delete [] simulated_missingFP;
  }

  //working matrices for Gibbs sampler
  //working matrices for Gibbs sampler
  //random effects
  if ( xTx != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete xTx[i];
    }
    delete [] xTx;
  }

  if ( xTy != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete xTy[i];
    }
    delete [] xTy;
  }

  //fixed effects
  if ( mTm != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete mTm[i];
    }
    delete [] mTm;
  }

  if ( mTy != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete mTy[i];
    }
    delete [] mTy;
  }

  if ( sMTm != NULL )
  {
    delete sMTm;
  }

  if ( sMTy != NULL )
  {
    delete sMTy;
  }

  //second stage
  if ( zTz != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete zTz[i];
    }
    delete [] zTz;
  }

  if ( zTb != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete zTb[i];
    }
    delete [] zTb;
  }

  if ( zTvIz != NULL )
  {
    delete zTvIz;
  }

  if ( zTvIb != NULL )
  {
    delete zTvIb;
  }

  if ( final_group_index != NULL )
  {
    delete final_group_index;
  }

  if ( tau2_error_weights != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete tau2_error_weights[i];
    }
    delete [] tau2_error_weights;
  }

  if ( tau2_group_weights != NULL )
  {
    delete tau2_group_weights;
  }

  if ( residuals != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete residuals[i];
    }
    delete [] residuals;
  }

  if ( second_residuals != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete second_residuals[i];
    }
    delete [] second_residuals;
  }

  //array of distributions map (maps distributions names to actual objects)
  if ( distr_map != NULL )
  {
    for ( i = 0; i < number_of_variables; i++ )
    {
      delete [] distr_map[i];
    }
    delete [] distr_map;
  }


  if ( randomEffects_means != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete randomEffects_means[ i ];
      delete randomEffects_vars[ i ];
    }
    delete [] randomEffects_means;
    delete [] randomEffects_vars;
  }

  if ( fixedEffects_means != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete fixedEffects_means[ i ];
      delete fixedEffects_vars[ i ];
    }
    delete [] fixedEffects_means;
    delete [] fixedEffects_vars;
  }


}//end



void BayesianMVHLM::initialize( CMatrix ** y, int n_groups )
{
  //each column of y is an element in the group
  //so number of elements in the group is number of columns
  //and dimension of vector case is number of rows

  int i;

  number_of_groups = n_groups;
  response = new CMatrix * [ number_of_groups ];
  number_of_observations = 0;

  final_group_index = new CVector( number_of_groups );
  for ( i = 0; i < number_of_groups; i++ )
  {
    response[i] = new CMatrix( (*y[i]) );

    if ( i > 0 )
    {
      final_group_index->Val( i ) = final_group_index->Val( i - 1 ) + response[i]->Col();
    }
    else
    {
      final_group_index->Val( i ) = response[i]->Col() - 1;
    }

    number_of_observations += response[i]->Col();

    //printf("BMDHL: response[%d] = \n", i );
    //response[i]->Print();

  }


}//end


void BayesianMVHLM::initializeResponseMissingData( bool keep_them )
{
  int i, j;

  keep_missing_data = keep_them;

  number_missing_response = new CVector( number_of_observations ); 
  number_missing_response->setToZero();
  missing_response_components = new CVector ** [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    missing_response_components[i] = new CVector * [ response[i]->Col() ];
    for ( j = 0; j < response[i]->Col(); j++ )
    {
      missing_response_components[i][j] = new CVector( 1 );
      missing_response_components[i][j]->setToZero();
    }
  }
}//end


void BayesianMVHLM::responseMissingData( int index, CVector & m_comp )
{
  int group, obs;

  DistributionParameter grp( groupIndex( index ) );
  group = (int) grp.getVector().Val(0);
  obs = (int) grp.getVector().Val(1);     

  delete missing_response_components[ group ][ obs ];
  missing_response_components[ group ][ obs ] = new CVector( m_comp );

  number_missing_response->Val( index ) = m_comp.Len();
  number_of_variables++;
}//end


void BayesianMVHLM::initializeRandomPredictorsMissingData( bool keep_them )
{
  int i, j;

  keep_missing_data = keep_them;

  number_missing_random_predictors = new CVector( number_of_observations ); 
  number_missing_random_predictors->setToZero();
  missing_random_predictors_components = new CVector ** [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    missing_random_predictors_components[i] = new CVector * [ response[i]->Col() ];
    for ( j = 0; j < response[i]->Col(); j++ )
    {
      missing_random_predictors_components[i][j] = new CVector( 1 );
      missing_random_predictors_components[i][j]->setToZero();
    }
  }
}//end


void BayesianMVHLM::randomPredictorsMissingData( int index, CVector & m_comp )
{
  int group, obs;

  DistributionParameter grp( groupIndex( index ) );
  group = (int) grp.getVector().Val(0);
  obs = (int) grp.getVector().Val(1);     

  delete missing_random_predictors_components[ group ][ obs ];
  missing_random_predictors_components[ group ][ obs ] = new CVector( m_comp );

  number_missing_random_predictors->Val( index ) = m_comp.Len();
  number_of_variables++;
}//end


void BayesianMVHLM::initializeFixedPredictorsMissingData( bool keep_them )
{
  int i, j;

  keep_missing_data = keep_them;

  number_missing_fixed_predictors = new CVector( number_of_observations ); 
  number_missing_fixed_predictors->setToZero();
  missing_fixed_predictors_components = new CVector ** [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    missing_fixed_predictors_components[i] = new CVector * [ response[i]->Col() ];
    for ( j = 0; j < response[i]->Col(); j++ )
    {
      missing_fixed_predictors_components[i][j] = new CVector( 1 );
      missing_fixed_predictors_components[i][j]->setToZero();
    }
  }
}//end


void BayesianMVHLM::fixedPredictorsMissingData( int index, CVector & m_comp )
{
  int group, obs;

  DistributionParameter grp( groupIndex( index ) );
  group = (int) grp.getVector().Val(0);
  obs = (int) grp.getVector().Val(1);     

  delete missing_fixed_predictors_components[ group ][ obs ];
  missing_fixed_predictors_components[ group ][ obs ] = new CVector( m_comp );

  number_missing_fixed_predictors->Val( index ) = m_comp.Len();
  number_of_variables++;
}//end



void BayesianMVHLM::randomEffects( CMatrix *** x )
{
  int i, j;

  //count beta's
  number_of_variables += number_of_groups;

  random_predictors = new CMatrix ** [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    random_predictors[i] = new CMatrix * [ response[i]->Col() ];
    for ( j = 0; j < response[i]->Col(); j++ )
    {
      random_predictors[i][j] = new CMatrix( (*(x[i][j])) );
    
      //printf("BHL: random predictors[%d][%d] = \n", i, j );
      //random_predictors[i][j]->Print();
    }
  }

  random_effects = true;  
}//end


void BayesianMVHLM::fixedEffects( CMatrix *** m )
{
  int i, j;

  //count gamma
  number_of_variables++;

  fixed_predictors = new CMatrix ** [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    fixed_predictors[i] = new CMatrix * [ response[i]->Col() ];
    for ( j = 0; j < response[i]->Col(); j++ )
    {
      fixed_predictors[i][j] = new CMatrix( (*(m[i][j])) );
    
      //printf("BHL: fixed predictors[%d][%d] = \n", i, j );
      //fixed_predictors[i][j]->Print();
    }
  }

  fixed_effects = true;  
}//end


void BayesianMVHLM::secondStageRandomEffects( CMatrix ** z )
{
  int i;

  //count alpha
  number_of_variables++;

  beta_predictors = new CMatrix * [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    beta_predictors[i] = new CMatrix( (*(z[i])) );

    //printf("BHL: level2 predictors[%d] = \n", i );
    //beta_predictors[i]->Print();
  }
  //count alpha

  second_stage_effects = true;
}//end


void BayesianMVHLM::betaPrior( CMatrix * p_betaCov )
{
  int i;

  //printf(" BHL: full model for beta. p_betaCov is\n" );
  //p_betaCov->Print();

  beta = new NormalDistribution * [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    beta[i] = new NormalDistribution( p_betaCov->Col() );
    beta[i]->setCovariance( p_betaCov );
  }

}//end


void BayesianMVHLM::betaPrior( CVector * p_beta, CMatrix * p_betaCov )
{
  int i;
  beta = new NormalDistribution * [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    beta[i] = new NormalDistribution( p_betaCov->Col() );
    beta[i]->setCovariance( p_betaCov );
    beta[i]->setMean( p_beta );

    //printf("beta prior is \n");
    //beta[i]->initialMean().Print();
    //beta[i]->initialCovariance().Print();

  }

}//end


void BayesianMVHLM::betaPrior( CMatrix * p_beta, CMatrix * p_betaCov )
{
  int i;

  beta = new NormalDistribution * [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    beta[i] = new NormalDistribution( p_betaCov->Col() );
    beta[i]->setCovariance( p_betaCov );
    CVector mean_val( p_beta->getColumn( i ) );
    beta[i]->setMean( &mean_val );

    //printf( "BHL: beta[%d] prior\n", i );
    //beta[i]->initialMean().Print();
    //beta[i]->initialCovariance().Print();
  }

}//end



void BayesianMVHLM::gammaPriorNonInformative( int dim )
{
  gamma = new NormalDistribution( dim );
  gamma->setNonInformative();

  non_informative_gamma = true;

  //printf( "BHL: noninfo gamma prior\n" );
  //gamma->initialMean().Print();
  //gamma->initialCovariance().Print();

}//end


void BayesianMVHLM::gammaPrior( CVector * p_gamma, CMatrix * p_gammaCov )    
{
  gamma = new NormalDistribution( p_gamma->Len() );
  gamma->setMean( p_gamma );
  gamma->setCovariance( p_gammaCov );


  //printf( "BHL: gamma prior\n" );
  //gamma->initialMean().Print();
  //gamma->initialCovariance().Print();

}//end



void BayesianMVHLM::alphaPrior( CVector * p_alpha, CMatrix * p_alphaCov ) throw( rtErr ) 
{
  if ( second_stage_effects )
  {
    alpha = new NormalDistribution( p_alpha->Len() );
    alpha->setMean( p_alpha );
    alpha->setCovariance( p_alphaCov );
  }
  else
  {
    printf( "BayesianMVHLM::alphaPrior: alpha is not a parameter of this model.\n" );
    char the_error[] = "BayesianMVHLM::alphaPrior: alpha is not a parameter of this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

  //printf( "BHL: alpha prior\n" );
  //alpha->initialMean().Print();
  //alpha->initialCovariance().Print();

}//end


void BayesianMVHLM::alphaPriorNonInformative( int dim ) throw( rtErr )
{
  if ( second_stage_effects )
  {
    alpha = new NormalDistribution( dim );
    alpha->setNonInformative();
    non_informative_alpha = true;
  }
  else
  {
    printf( "BayesianMVHLM::alphaPriorNonInformative: alpha is not a parameter of this model.\n" );
    char the_error[] = "BayesianMVHLM::alphaPriorNonInformative: alpha is not a parameter of this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

  //printf( "BHL: alpha prior: non informative\n" );
  //alpha->initialMean().Print();
  //alpha->initialCovariance().Print();

}//end



void BayesianMVHLM::sigma2CommonPrior( bool val )
{
  common_sigma = val;
}//end



void BayesianMVHLM::sigma2PriorInvChisq( double p_nuSigma2, double p_sigma2 )
{
  invChisq_sigma2 = true;

  if ( common_sigma )
  {
    sigma2 = new DistributionMixture( 1 );
    sigma2ICS = new InvChisqDistribution * [ 1 ];
    sigma2ICS[0] = new InvChisqDistribution( p_nuSigma2 );
    sigma2ICS[0]->setScale( p_sigma2 );
    sigma2_first_draw = new DistributionParameter * [ 1 ];
    sigma2_first_draw[0] = new DistributionParameter( p_sigma2 );
    sigma2->set( 0, sigma2ICS[0], 1.0 );

    number_of_variables++;
  }
  else
  {
    int i;

    sigma2 = new DistributionMixture( number_of_groups );
    sigma2ICS = new InvChisqDistribution * [ number_of_groups ];
    sigma2_first_draw = new DistributionParameter * [ number_of_groups ];
    for ( i = 0; i < number_of_groups; i++ )
    {
      sigma2ICS[i] = new InvChisqDistribution( p_nuSigma2 );
      sigma2ICS[i]->setScale( p_sigma2 );
      sigma2_first_draw[i] = new DistributionParameter( p_sigma2 );
      sigma2->set( i, sigma2ICS[i], 1.0 );
      number_of_variables++;
    }
  }
}//end



void BayesianMVHLM::sigma2PriorInvChisq( double * p_nuSigma2, double * p_sigma2 )
{
  invChisq_sigma2 = true;

  if ( common_sigma )
  {
    sigma2 = new DistributionMixture( 1 );
    sigma2ICS = new InvChisqDistribution * [ 1 ];
    sigma2ICS[0] = new InvChisqDistribution( p_nuSigma2[0] );
    sigma2ICS[0]->setScale( p_sigma2[0] );
    sigma2_first_draw = new DistributionParameter * [ 1 ];
    sigma2_first_draw[0] = new DistributionParameter( p_sigma2[0] );
    sigma2->set( 0, sigma2ICS[0], 1.0 );

    number_of_variables++;
  }
  else
  {
    int i;

    sigma2 = new DistributionMixture( number_of_groups );
    sigma2ICS = new InvChisqDistribution * [ number_of_groups ];
    sigma2_first_draw = new DistributionParameter * [ number_of_groups ];
    for ( i = 0; i < number_of_groups; i++ )
    {
      sigma2ICS[i] = new InvChisqDistribution( p_nuSigma2[i] );
      sigma2ICS[i]->setScale( p_sigma2[i] );
      sigma2_first_draw[i] = new DistributionParameter( p_sigma2[i] );
      sigma2->set( i, sigma2ICS[i], 1.0 );
      number_of_variables++;
    }
  }
}


void BayesianMVHLM::sigma2PriorDuMouchel( double p_sigma2 )
{
  double df;

  df = number_of_observations * response[0]->Row() + 1;
  duMouchel_sigma2 = true;
  properNonInfoPrior_sigma2 = true;
  
  if ( common_sigma )
  {
    sigma2 = new DistributionMixture( 1 );
    sigma2PNIP = new ProperNonInfoPosteriorHLM * [ 1 ];
    sigma2PNIP[0] = new ProperNonInfoPosteriorHLM( p_sigma2, df );
    sigma2PNIP[0]->setPriorDistribution( "du_mouchel" );
    sigma2_first_draw = new DistributionParameter * [ 1 ];
    sigma2_first_draw[0] = new DistributionParameter( p_sigma2 );
    sigma2->set( 0, sigma2PNIP[0], 1.0 );

    number_of_variables++;
  }
  else
  {
    int i;

    sigma2 = new DistributionMixture( number_of_groups );
    sigma2PNIP = new ProperNonInfoPosteriorHLM * [ number_of_groups ];
    sigma2_first_draw = new DistributionParameter * [ number_of_groups ];
    for ( i = 0; i < number_of_groups; i++ )
    {
      sigma2PNIP[i] = new ProperNonInfoPosteriorHLM( p_sigma2, df );
      sigma2PNIP[i]->setPriorDistribution( "du_mouchel" );
      sigma2_first_draw[i] = new DistributionParameter( p_sigma2 );
      sigma2->set( i, sigma2PNIP[i], 1.0 );
      number_of_variables++;
    }
  }
}//end


void BayesianMVHLM::sigma2PriorDuMouchel( double * p_sigma2 )
{
  double df;

  df = number_of_observations * response[0]->Row() + 1;
  duMouchel_sigma2 = true;
  properNonInfoPrior_sigma2 = true;
  
  if ( common_sigma )
  {
    sigma2 = new DistributionMixture( 1 );
    sigma2PNIP = new ProperNonInfoPosteriorHLM * [ 1 ];
    sigma2PNIP[0] = new ProperNonInfoPosteriorHLM( p_sigma2[0], df );
    sigma2PNIP[0]->setPriorDistribution( "du_mouchel" );
    sigma2_first_draw = new DistributionParameter * [ 1 ];
    sigma2_first_draw[0] = new DistributionParameter( p_sigma2[0] );
    sigma2->set( 0, sigma2PNIP[0], 1.0 );

    number_of_variables++;
  }
  else
  {
    int i;

    sigma2 = new DistributionMixture( number_of_groups );
    sigma2PNIP = new ProperNonInfoPosteriorHLM * [ number_of_groups ];
    sigma2_first_draw = new DistributionParameter * [ number_of_groups ];
    for ( i = 0; i < number_of_groups; i++ )
    {
      sigma2PNIP[i] = new ProperNonInfoPosteriorHLM( p_sigma2[i], df );
      sigma2PNIP[i]->setPriorDistribution( "du_mouchel" );
      sigma2_first_draw[i] = new DistributionParameter( p_sigma2[i] );
      sigma2->set( i, sigma2PNIP[i], 1.0 );
      number_of_variables++;
    }
  }
}//end



void BayesianMVHLM::sigma2PriorUniformShrinkage( double p_sigma2 )
{
  double df;

  df = number_of_observations * response[0]->Row();
  uniformShrinkage_sigma2 = true;
  properNonInfoPrior_sigma2 = true;

  if ( common_sigma )
  {
    sigma2 = new DistributionMixture( 1 );
    sigma2PNIP = new ProperNonInfoPosteriorHLM * [ 1 ];
    sigma2_first_draw = new DistributionParameter * [ 1 ];

    sigma2PNIP[0] = new ProperNonInfoPosteriorHLM( p_sigma2, df );
    sigma2PNIP[0]->setPriorDistribution( "uniform_shrinkage" );

    sigma2_first_draw[0] = new DistributionParameter( p_sigma2 );
    sigma2->set( 0, sigma2PNIP[0], 1.0 );

    number_of_variables++;
  }
  else
  {
    int i;

    sigma2 = new DistributionMixture( number_of_groups );
    sigma2PNIP = new ProperNonInfoPosteriorHLM * [ number_of_groups ];
    sigma2_first_draw = new DistributionParameter * [ number_of_groups ];
    for ( i = 0; i < number_of_groups; i++ )
    {
      sigma2PNIP[i] = new ProperNonInfoPosteriorHLM( p_sigma2, df );
      sigma2PNIP[i]->setPriorDistribution( "uniform_shrinkage" );

      sigma2_first_draw[i] = new DistributionParameter( p_sigma2 );
      sigma2->set( i, sigma2PNIP[i], 1.0 );
      number_of_variables++;
    }
  }
}//end


void BayesianMVHLM::sigma2PriorUniformShrinkage( double * p_sigma2 )
{
  double df;

  df = number_of_observations * response[0]->Row();
  uniformShrinkage_sigma2 = true;
  properNonInfoPrior_sigma2 = true;

  if ( common_sigma )
  {
    sigma2 = new DistributionMixture( 1 );
    sigma2PNIP = new ProperNonInfoPosteriorHLM * [ 1 ];
    sigma2_first_draw = new DistributionParameter * [ 1 ];

    sigma2PNIP[0] = new ProperNonInfoPosteriorHLM( p_sigma2[0], df );
    sigma2PNIP[0]->setPriorDistribution( "uniform_shrinkage" );

    sigma2_first_draw[0] = new DistributionParameter( p_sigma2[0] );
    sigma2->set( 0, sigma2PNIP[0], 1.0 );

    number_of_variables++;
  }
  else
  {
    int i;

    sigma2 = new DistributionMixture( number_of_groups );
    sigma2PNIP = new ProperNonInfoPosteriorHLM * [ number_of_groups ];
    sigma2_first_draw = new DistributionParameter * [ number_of_groups ];
    for ( i = 0; i < number_of_groups; i++ )
    {
      sigma2PNIP[i] = new ProperNonInfoPosteriorHLM( p_sigma2[i], df );
      sigma2PNIP[i]->setPriorDistribution( "uniform_shrinkage" );

      sigma2_first_draw[i] = new DistributionParameter( p_sigma2[i] );
      sigma2->set( i, sigma2PNIP[i], 1.0 );
      number_of_variables++;
    }
  }
}//end



void BayesianMVHLM::sigma2PriorNonInformative( double p_power ) throw( rtErr )
{

  nonInformativePower_sigma2 = true;

  if ( common_sigma )
  {
    if ( p_power <= -0.5 * number_of_observations * response[0]->Row() )
    {
      printf( "BayesianMVHLM::sigma2PriorNonInformative: Power [%f] in noninformative prior is not valid.\n", p_power );  
      char the_error[] = "BayesianMVHLM::sigma2PriorNonInformative: Power in noninformative prior is not valid.";
      rtErr runtime_error( the_error );
      throw runtime_error;
    }

    sigma2NIP = new double [ 1 ];
    sigma2 = new DistributionMixture( 1 );
    sigma2ICS = new InvChisqDistribution * [ 1 ];
    sigma2_first_draw = new DistributionParameter * [ 1 ];

    sigma2ICS[0] = new InvChisqDistribution( 0 );
    sigma2ICS[0]->setScale( 1.0 );
    sigma2_first_draw[0] = new DistributionParameter( 1.0 );
    sigma2->set( 0, sigma2ICS[0], 1.0 );
    sigma2NIP[0] = p_power;

    number_of_variables++;
  }
  else
  {
    int i;

    sigma2NIP = new double [ number_of_groups ];
    sigma2 = new DistributionMixture( number_of_groups );
    sigma2ICS = new InvChisqDistribution * [ number_of_groups ];
    sigma2_first_draw = new DistributionParameter * [ number_of_groups ];

    for ( i = 0; i < number_of_groups; i++ )
    {
      if ( p_power <= -0.5 * response[i]->Row() * response[i]->Col() )
      {
        printf( "BayesianMVHLM::sigma2PriorNonInformative: Power [%f] in noninformative prior is not valid.\n", p_power );  
        char the_error[] = "BayesianMVHLM::sigma2PriorNonInformative: Power in noninformative prior is not valid.";
        rtErr runtime_error( the_error );
        throw runtime_error;
      }

      sigma2ICS[i] = new InvChisqDistribution( 0 );
      sigma2ICS[i]->setScale( 1.0 );
      sigma2_first_draw[i] = new DistributionParameter( 1.0 );
      sigma2->set( i, sigma2ICS[i], 1.0 );
      sigma2NIP[i] = p_power;
      number_of_variables++;
    }
  }
}//end


void BayesianMVHLM::sigma2PriorNonInformative( double * p_power ) throw( rtErr )
{

  nonInformativePower_sigma2 = true;

  if ( common_sigma )
  {
    if ( p_power[0] <= -0.5 * number_of_observations * response[0]->Row() )
    {
      printf( "BayesianMVHLM::sigma2PriorNonInformative: Power [%f] in noninformative prior is not valid.\n", p_power );  
      char the_error[] = "BayesianMVHLM::sigma2PriorNonInformative: Power in noninformative prior is not valid.";
      rtErr runtime_error( the_error );
      throw runtime_error;
    }

    sigma2NIP = new double [ 1 ];
    sigma2 = new DistributionMixture( 1 );
    sigma2ICS = new InvChisqDistribution * [ 1 ];
    sigma2_first_draw = new DistributionParameter * [ 1 ];

    sigma2ICS[0] = new InvChisqDistribution( 0 );
    sigma2ICS[0]->setScale( 1.0 );
    sigma2_first_draw[0] = new DistributionParameter( 1.0 );
    sigma2->set( 0, sigma2ICS[0], 1.0 );
    sigma2NIP[0] = p_power[0];

    number_of_variables++;
  }
  else
  {
    int i;

    sigma2NIP = new double [ number_of_groups ];
    sigma2 = new DistributionMixture( number_of_groups );
    sigma2ICS = new InvChisqDistribution * [ number_of_groups ];
    sigma2_first_draw = new DistributionParameter * [ number_of_groups ];

    for ( i = 0; i < number_of_groups; i++ )
    {
      if ( p_power[i] <= -0.5 * response[i]->Row() * response[i]->Col() )
      {
        printf( "BayesianMVHLM::sigma2PriorNonInformative: Power [%f] in noninformative prior is not valid.\n", p_power );  
        char the_error[] = "BayesianMVHLM::sigma2PriorNonInformative: Power in noninformative prior is not valid.";
        rtErr runtime_error( the_error );
        throw runtime_error;
      }

      sigma2ICS[i] = new InvChisqDistribution( 0 );
      sigma2ICS[i]->setScale( 1.0 );
      sigma2_first_draw[i] = new DistributionParameter( 1.0 );
      sigma2->set( i, sigma2ICS[i], 1.0 );
      sigma2NIP[i] = p_power[i];
      number_of_variables++;
    }
  }
}//end


void BayesianMVHLM::sigma2PriorInvWishart( double p_nu, CMatrix * p_Cov )
{

  invWishart_sigma2 = true;

  if ( common_sigma )
  {
    sigma2 = new DistributionMixture( 1 );    
    sigma2InvWishart = new InvWishartDistribution * [ 1 ];
    sigma2_first_draw = new DistributionParameter * [1 ];

    sigma2InvWishart[0] = new InvWishartDistribution( p_nu );
    sigma2InvWishart[0]->setScaleMatrix( p_Cov );
    sigma2_first_draw[0] = new DistributionParameter( (*p_Cov) );
    sigma2->set( 0, sigma2InvWishart[0], 1.0 );
  
    number_of_variables++;
  }
  else
  {
    int i;

    sigma2 = new DistributionMixture( number_of_groups );    
    sigma2InvWishart = new InvWishartDistribution * [ number_of_groups ];
    sigma2_first_draw = new DistributionParameter * [ number_of_groups ];

    for ( i = 0; i < number_of_groups; i++ )
    {
      sigma2InvWishart[i] = new InvWishartDistribution( p_nu );
      sigma2InvWishart[i]->setScaleMatrix( p_Cov );
      sigma2_first_draw[i] = new DistributionParameter( (*p_Cov) );
      sigma2->set( i, sigma2InvWishart[i], 1.0 );
      number_of_variables++;
    }
  }
}//end


void BayesianMVHLM::sigma2PriorInvWishart( double * p_nu, CMatrix ** p_Cov )
{

  invWishart_sigma2 = true;

  if ( common_sigma )
  {
    sigma2 = new DistributionMixture( 1 );    
    sigma2InvWishart = new InvWishartDistribution * [ 1 ];
    sigma2_first_draw = new DistributionParameter * [ 1 ];

    sigma2InvWishart[0] = new InvWishartDistribution( p_nu[0] );
    sigma2InvWishart[0]->setScaleMatrix( p_Cov[0] );
    sigma2_first_draw[0] = new DistributionParameter( (*(p_Cov[0])) );
    sigma2->set( 0, sigma2InvWishart[0], 1.0 );
  
    number_of_variables++;
  }
  else
  {
    int i;

    sigma2 = new DistributionMixture( number_of_groups );    
    sigma2InvWishart = new InvWishartDistribution * [ number_of_groups ];
    sigma2_first_draw = new DistributionParameter * [ number_of_groups ];

    for ( i = 0; i < number_of_groups; i++ )
    {
      sigma2InvWishart[i] = new InvWishartDistribution( p_nu[i] );
      sigma2InvWishart[i]->setScaleMatrix( p_Cov[i] );
      sigma2_first_draw[i] = new DistributionParameter( (*(p_Cov[i])) );
      sigma2->set( i, sigma2InvWishart[i], 1.0 );
      number_of_variables++;
    }
  }
}//end



void BayesianMVHLM::sigma2Known( CMatrix * p_Cov )
{
  known_sigma2 = true;
  invWishart_sigma2 = true;

  if ( common_sigma )
  {
    sigma2 = new DistributionMixture( 1 );    
    sigma2InvWishart = new InvWishartDistribution * [ 1 ];
    sigma2_first_draw = new DistributionParameter * [1 ];

    sigma2InvWishart[0] = new InvWishartDistribution( 1 );
    sigma2InvWishart[0]->setScaleMatrix( p_Cov );
    sigma2_first_draw[0] = new DistributionParameter( (*p_Cov) );
    sigma2->set( 0, sigma2InvWishart[0], 1.0 );
  
  }
  else
  {
    int i;

    sigma2 = new DistributionMixture( number_of_groups );    
    sigma2InvWishart = new InvWishartDistribution * [ number_of_groups ];
    sigma2_first_draw = new DistributionParameter * [ number_of_groups ];

    for ( i = 0; i < number_of_groups; i++ )
    {
      sigma2InvWishart[i] = new InvWishartDistribution( 1 );
      sigma2InvWishart[i]->setScaleMatrix( p_Cov );
      sigma2_first_draw[i] = new DistributionParameter( (*p_Cov) );
      sigma2->set( i, sigma2InvWishart[i], 1.0 );
    }
  }
}//end


void BayesianMVHLM::sigma2Known( CMatrix ** p_Cov )
{

  known_sigma2 = true;
  invWishart_sigma2 = true;

  if ( common_sigma )
  {
    sigma2 = new DistributionMixture( 1 );    
    sigma2InvWishart = new InvWishartDistribution * [ 1 ];
    sigma2_first_draw = new DistributionParameter * [ 1 ];

    sigma2InvWishart[0] = new InvWishartDistribution( 1 );
    sigma2InvWishart[0]->setScaleMatrix( p_Cov[0] );
    sigma2_first_draw[0] = new DistributionParameter( (*(p_Cov[0])) );
    sigma2->set( 0, sigma2InvWishart[0], 1.0 );
  
  }
  else
  {
    int i;

    sigma2 = new DistributionMixture( number_of_groups );    
    sigma2InvWishart = new InvWishartDistribution * [ number_of_groups ];
    sigma2_first_draw = new DistributionParameter * [ number_of_groups ];

    for ( i = 0; i < number_of_groups; i++ )
    {
      sigma2InvWishart[i] = new InvWishartDistribution( 1 );
      sigma2InvWishart[i]->setScaleMatrix( p_Cov[i] );
      sigma2_first_draw[i] = new DistributionParameter( (*(p_Cov[i])) );
      sigma2->set( i, sigma2InvWishart[i], 1.0 );
    }
  }
}//end



void BayesianMVHLM::tau2PriorInvChisq( double p_nuTau2, double p_tau2 ) throw( rtErr )
{
  if ( random_effects )
  {
    tau2 = new DistributionMixture( 1 );
    tau2ICS = new InvChisqDistribution( p_nuTau2 );
    tau2ICS->setScale( p_tau2 );
    tau2_first_draw = new DistributionParameter( p_tau2 );
    tau2->set( 0, tau2ICS, 1.0 );
    invChisq_tau2 = true;

    number_of_variables++;
  }
  else
  {
    printf( "BayesianMVHLM::tau2PriorInvChisq: There are no random effects in the model. Tau2 does not make sense in this model.\n" );
    char the_error[] = "BayesianMVHLM::tau2PriorInvChisq: There are no random effects in the model. Tau2 does not make sense in this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

}//end


void BayesianMVHLM::tau2PriorDuMouchel( double p_tau2 ) throw( rtErr )
{
  double df;

  if ( random_effects )
  {
    tau2 = new DistributionMixture( 1 );
    df = number_of_groups * beta[0]->dimension() + 1;

    tau2PNIP = new ProperNonInfoPosteriorHLM( p_tau2, df );
    tau2PNIP->setPriorDistribution( "du_mouchel" );

    tau2_first_draw = new DistributionParameter( p_tau2 );
    tau2->set( 0, tau2PNIP, 1.0 );
    duMouchel_tau2 = true;
    properNonInfoPrior_tau2 = true;

    number_of_variables++;
  }
  else
  {
    printf( "BayesianMVHLM::tau2PriorDuMouchel: There are no random effects in the model. Tau2 does not make sense in this model.\n" );
    char the_error[] = "BayesianMVHLM::tau2PriorDuMouchel: There are no random effects in the model. Tau2 does not make sense in this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}//end


void BayesianMVHLM::tau2PriorUniformShrinkage( double p_tau2 ) throw( rtErr )
{
  double df;

  if ( random_effects )
  {
    tau2 = new DistributionMixture( 1 );
    df = number_of_groups * beta[0]->dimension();

    tau2PNIP = new ProperNonInfoPosteriorHLM( p_tau2, df );
    tau2PNIP->setPriorDistribution( "uniform_shrinkage" );

    tau2_first_draw = new DistributionParameter( p_tau2 );
    tau2->set( 0, tau2PNIP, 1.0 );
    uniformShrinkage_tau2 = true;
    properNonInfoPrior_tau2 = true;

    number_of_variables++;
  }
  else
  {
    printf( "BayesianMVHLM::tau2PriorUniformShrinkage: There are no random effects in the model. Tau2 does not make sense in this model.\n" );
    char the_error[] = "BayesianMVHLM::tau2PriorUniformShrinkage: There are no random effects in the model. Tau2 does not make sense in this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}//end


void BayesianMVHLM::tau2PriorNonInformative( double p_power ) throw( rtErr )
{

  if ( random_effects )
  {
    //printf("BHLM: p_power = %f, n groups = %d,  p = %d\n", p_power, number_of_groups, random_predictors[0]->Col()); fflush(stdout);

    if ( p_power <= -0.5 * number_of_groups * beta[0]->dimension() )
    {
      printf( "BayesianMVHLM::tau2PriorNonInformative: Power [%f] in noninformative prior is not valid.\n", p_power );  
      char the_error[] = "BayesianMVHLM::tau2PriorNonInformative: Power in noninformative prior is not valid.";
      rtErr runtime_error( the_error );
      throw runtime_error;
    }

    tau2NIP = p_power;
    tau2 = new DistributionMixture( 1 );
    tau2ICS = new InvChisqDistribution( 0 );
    tau2ICS->setScale( 1.0 );
    tau2_first_draw = new DistributionParameter( 1.0 );
    tau2->set( 0, tau2ICS, 1.0 );
    nonInformativePower_tau2 = true;

    number_of_variables++;
  }
  else
  {
    printf( "BayesianMVHLM::tau2PriorNonInformative: There are no random effects in the model. Tau2 does not make sense in this model.\n" );
    char the_error[] = "BayesianMVHLM::tau2PriorNonInformative: There are no random effects in the model. Tau2 does not make sense in this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

}//end



void BayesianMVHLM::betaCovPriorInvWishart( double p_nuV, CMatrix * p_VCov ) throw( rtErr )
{
  if ( random_effects )
  {
    tau2 = new DistributionMixture( 1 );
    tau2InvWishart = new InvWishartDistribution( p_nuV );
    tau2InvWishart->setScaleMatrix( p_VCov );
    tau2_first_draw = new DistributionParameter( (*p_VCov) );
    tau2->set( 0, tau2InvWishart, 1.0 );
    invWishart_betaCov = true;
  
    number_of_variables++;
  }
  else
  {
    printf( "BayesianMVHLM::betaCovPriorInvWishart: There are no random effects in the model. Tau2 does not make sense in this model.\n" );
    char the_error[] = "BayesianMVHLM::betaCovPriorInvWishart: There are no random effects in the model. Tau2 does not make sense in this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

}//end



void BayesianMVHLM::tLikelihood( double p_nuError )
{
  int i;

  tau2_errors = new InvChisqDistribution * [ number_of_observations ];
  for ( i = 0; i < number_of_observations; i++ )
  {
    tau2_errors[i] = new InvChisqDistribution( p_nuError );
    tau2_errors[i]->setScale( 1.0 );
  }

  number_of_variables += number_of_observations;
  t_lkhd = true;
}//end



void BayesianMVHLM::groupTLikelihood( double p_nuGroup )
{
  int i;

  tau2_groups = new InvChisqDistribution * [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    tau2_groups[i] = new InvChisqDistribution( p_nuGroup );
    tau2_groups[i]->setScale( 1.0 );
  }

  number_of_variables += number_of_groups;
  t_group_lkhd = true;
}//end



void BayesianMVHLM::gammaTPrior( double p_nuGamma ) throw( rtErr )
{
  if ( fixed_effects )
  {
    tau2_gamma = new InvChisqDistribution( p_nuGamma );
    tau2_gamma->setScale( 1.0 );
    t_gamma = true;

    number_of_variables++;
  }
  else
  {
    printf( "BayesianMVHLM::gammaTPrior: gamma is not a parameter of this model.\n" );
    char the_error[] = "BayesianMVHLM::gammaTPrior: gamma is not a parameter of this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

}//end



void BayesianMVHLM::betaTPrior( double p_nuBeta ) throw( rtErr )
{
  int i;

  if ( random_effects )
  {
    tau2_betas = new InvChisqDistribution * [ number_of_groups ];
    for ( i = 0; i < number_of_groups; i++ )
    {
      tau2_betas[i] = new InvChisqDistribution( p_nuBeta );
      tau2_betas[i]->setScale( 1.0 );
    }
    t_beta = true;

    number_of_variables += number_of_groups;
  }
  else
  {
    printf( "BayesianMVHLM::betaTPrior: beta is not a parameter of this model.\n" );
    char the_error[] = "BayesianMVHLM::betaTPrior: beta is not a parameter of this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}//end


void BayesianMVHLM::alphaTPrior( double p_nuAlpha ) throw( rtErr )
{
  if ( second_stage_effects )
  {
    tau2_alpha = new InvChisqDistribution( p_nuAlpha );
    tau2_alpha->setScale( 1.0 );
    t_alpha = true;

    number_of_variables++;
  }
  else
  {
    printf( "BayesianMVHLM::alphaTPrior: alpha is not a parameter of this model.\n" );
    char the_error[] = "BayesianMVHLM::alphaTPrior: alpha is not a parameter of this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}//end



void BayesianMVHLM::samplerDefaultInitialPoint()
{
  int i;

  if ( random_effects )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
#ifdef FIX1
      DistributionParameter tmpdist = beta[i]->mean();
      beta[i]->setLastDraw( tmpdist );
#else
      beta[i]->setLastDraw( beta[i]->mean() );
#endif
//      beta[i]->setLastDraw( beta[i]->mean() );
    }

    tau2->setLastDraw( (*tau2_first_draw) );

    if ( second_stage_effects )
    {
#ifdef FIX1
      DistributionParameter tmpdist = alpha->mean();
      alpha->setLastDraw( tmpdist );
#else
      alpha->setLastDraw( alpha->mean() );
#endif
//      alpha->setLastDraw( alpha->mean() );
    }
  }

  if ( fixed_effects )
  {
#ifdef FIX1
      DistributionParameter tmpdist = gamma->mean();
      gamma->setLastDraw( tmpdist );
#else
      gamma->setLastDraw( gamma->mean() );
#endif
//    gamma->setLastDraw( gamma->mean() );
  }

  if ( common_sigma )
  {
    sigma2->setLastDraw( ( *(sigma2_first_draw[0]) ) );
  }
  else
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      sigma2->setLastDrawForComponent( i, ( *(sigma2_first_draw[i]) ) );
    }
  }
}//end


void BayesianMVHLM::samplerSigma2InitialPoint( double init_sigma2 )
{
  if ( common_sigma )
  {
    delete sigma2_first_draw[0];
    delete [] sigma2_first_draw;

    sigma2_first_draw = new DistributionParameter * [ 1 ];
    sigma2_first_draw[0] = new DistributionParameter( init_sigma2 );
    sigma2->setLastDraw( ( *(sigma2_first_draw[0]) ) );
  }
  else
  {
    int i;

    for ( i = 0; i < number_of_groups; i++ )
    {
      delete sigma2_first_draw[i];
    }
    delete [] sigma2_first_draw;

    sigma2_first_draw = new DistributionParameter * [ number_of_groups ];
    for ( i = 0; i < number_of_groups; i++ )
    {
      sigma2_first_draw[i] = new DistributionParameter( init_sigma2 );
      sigma2->setLastDrawForComponent( i, ( *(sigma2_first_draw[i]) ) );
    }
  }

}//end


void BayesianMVHLM::samplerSigma2InitialPoint( CVector & init_sigma2 )
{
  if ( common_sigma )
  {
    delete sigma2_first_draw[0];
    delete [] sigma2_first_draw;

    sigma2_first_draw = new DistributionParameter * [ 1 ];
    sigma2_first_draw[0] = new DistributionParameter( init_sigma2.Val(0) );
    sigma2->setLastDraw( ( *(sigma2_first_draw[0]) ) );
  }
  else
  {
    int i;

    for ( i = 0; i < number_of_groups; i++ )
    {
      delete sigma2_first_draw[i];
    }
    delete [] sigma2_first_draw;

    sigma2_first_draw = new DistributionParameter * [ number_of_groups ];
    for ( i = 0; i < number_of_groups; i++ )
    {
      sigma2_first_draw[i] = new DistributionParameter( init_sigma2.Val(i) );
      sigma2->setLastDrawForComponent( i, ( *(sigma2_first_draw[i]) ) );
    }
  }

}//end



void BayesianMVHLM::samplerSigma2InitialPoint( CMatrix & init_sigma2 )
{
  if ( common_sigma )
  {
    delete sigma2_first_draw[0];
    delete [] sigma2_first_draw;

    sigma2_first_draw = new DistributionParameter * [ 1 ];
    sigma2_first_draw[0] = new DistributionParameter( init_sigma2 );
    sigma2->setLastDraw( ( *(sigma2_first_draw[0]) ) );
  }
  else
  {
    int i;

    for ( i = 0; i < number_of_groups; i++ )
    {
      delete sigma2_first_draw[i];
    }
    delete [] sigma2_first_draw;

    sigma2_first_draw = new DistributionParameter * [ number_of_groups ];
    for ( i = 0; i < number_of_groups; i++ )
    {
      sigma2_first_draw[i] = new DistributionParameter( init_sigma2 );
      sigma2->setLastDrawForComponent( i, ( *(sigma2_first_draw[i]) ) );
    }
  }

}//end


void BayesianMVHLM::samplerSigma2InitialPoint( CMatrix ** init_sigma2 )
{
  if ( common_sigma )
  {
    delete sigma2_first_draw[0];
    delete [] sigma2_first_draw;

    sigma2_first_draw = new DistributionParameter * [ 1 ];
    sigma2_first_draw[0] = new DistributionParameter( (*(init_sigma2[0])) );
    sigma2->setLastDraw( ( *(sigma2_first_draw[0]) ) );
  }
  else
  {
    int i;

    for ( i = 0; i < number_of_groups; i++ )
    {
      delete sigma2_first_draw[i];
    }
    delete [] sigma2_first_draw;

    sigma2_first_draw = new DistributionParameter * [ number_of_groups ];
    for ( i = 0; i < number_of_groups; i++ )
    {
      sigma2_first_draw[i] = new DistributionParameter( (*(init_sigma2[i])) );
      sigma2->setLastDrawForComponent( i, ( *(sigma2_first_draw[i]) ) );
    }
  }

}//end




void BayesianMVHLM::samplerTau2InitialPoint( double init_tau2 )
{
  if ( !invWishart_betaCov )
  {
    delete tau2_first_draw;
    tau2_first_draw = new DistributionParameter( init_tau2 );
    tau2->setLastDraw( (*tau2_first_draw) );
  }
  else
  {
    printf( "BayesianMVHLM::samplerTau2InitialPoint: Initial point should be a matrix, not a scalar.\n" );
  }
}//end


void BayesianMVHLM::samplerTau2InitialPoint( CMatrix & init_tau2 )
{
  if ( invWishart_betaCov )
  {
    delete tau2_first_draw;
    tau2_first_draw = new DistributionParameter( init_tau2 );
    tau2->setLastDraw( (*tau2_first_draw) );
  }
  else
  {
    printf( "BayesianMVHLM::samplerTau2InitialPoint: Initial point should not be a matrix, but a scalar.\n" );
  }
}//end


void BayesianMVHLM::samplerAlphaInitialPoint( CVector & init_alpha )
{
  alpha->setLastDraw( init_alpha );
}//end


void BayesianMVHLM::samplerBetaInitialPoint( CVector & init_beta )
{
  int i;
  
  for ( i = 0; i < number_of_groups; i++ )
  {
    beta[ i ]->setLastDraw( init_beta );
  }
}//end


void BayesianMVHLM::samplerBetaInitialPoint( CMatrix & init_beta )
{
  int i;
  
  for ( i = 0; i < number_of_groups; i++ )
  {
    beta[ i ]->setLastDraw( init_beta.getColumn( i ) );
  }
}//end



void BayesianMVHLM::samplerGammaInitialPoint( CVector & init_gamma )
{
  gamma->setLastDraw( init_gamma );

}//end



void BayesianMVHLM::samplerMissingVariablesInitialPoint()
{
  if ( keep_missing_data )
  {
    if ( number_missing_response != NULL )    
    {
      samplerMissingResponseInitialPoint();
    }

    if ( number_missing_random_predictors != NULL )
    {
      samplerMissingRandomPredictorsInitialPoint();
    }

    if ( number_missing_fixed_predictors != NULL )
    {
      samplerMissingFixedPredictorsInitialPoint();
    }
  }
}//end



void BayesianMVHLM::samplerMissingResponseInitialPoint()
{
  int i, k, index, group, obs, nmss;

  //compute means for each group ignoring missing data

  DistributionParameter * grp;
  CVector ** response_means;
  CVector ** n_items;
  CVector missing_in_group( number_of_groups );
  
  missing_in_group.setToZero();
  response_means = new CVector * [ number_of_groups ];
  n_items = new CVector * [ number_of_groups ];
  for ( group = 0; group < number_of_groups; group++ )
  {
    response_means[ group ] = new CVector( response[ 0 ]->Row() );
    n_items[ group ] = new CVector( response[ 0 ]->Row() );
  
    response_means[ group]->setToZero();
    n_items[ group ]->setToZero();
  }

  index = 0;
  group = 0;
  obs = 0;
  nmss = (int) number_missing_response->Val( 0 );

  while ( index < number_of_observations )
  {
    for ( i = 0; i < response[ group ]->Row(); i++ )
    {
      if ( nmss > 0 )
      {
        if ( !missing_response_components[ group ][ obs ]->includes( ((double) i ) ) )
        {
          response_means[ group ]->Val( i ) += response[ group ]->Val( i, obs );
          n_items[ group ]->Val(i) = n_items[ group ]->Val(i) + 1;
        }
        else
	{
          missing_in_group.Val( group ) = missing_in_group.Val( group ) + 1;
        }
      }
      else
      {
        response_means[ group ]->Val( i ) += response[ group ]->Val( i, obs );
        n_items[ group ]->Val(i) = n_items[ group ]->Val(i) + 1;
      }
    }//end for i

    index++;
    grp = new DistributionParameter( groupIndex( index ) );
    group = (int) grp->getVector().Val(0);
    obs = (int) grp->getVector().Val(1);
    nmss = (int) number_missing_response->Val( index );

    delete grp;

  }//end while index

  //compute means
  for ( group = 0; group < number_of_groups; group++ )
  {
    if ( missing_in_group.Val( group ) > 0 )
    {
      for ( i = 0; i < response[ group ]->Row(); i++ )
      {
        if ( n_items[ group ]->Val( i ) > 0 )
	{
          response_means[ group ]->Val( i ) = response_means[ group ]->Val( i ) / n_items[ group ]->Val( i );
        }        
      }
    }
  }//end for



  //ready to impute by means
  index = 0;
  group = 0;
  obs = 0;
  nmss = (int) number_missing_response->Val( 0 );

  while ( index < number_of_observations )
  {
    if ( nmss > 0 )
    {   
      for ( i = 0; i < nmss; i++ )
      {
        k = (int) missing_response_components[ group ][ obs ]->Val( i );
        response[ group ]->Val( k, obs ) = response_means[ group ]->Val( k );
      }

      //printf("initial response[%d][%d] is \n", group, obs );
      //response[ group ]->Print();  

    }//end if

    index++;
    grp = new DistributionParameter( groupIndex( index ) );
    group = (int) grp->getVector().Val(0);
    obs = (int) grp->getVector().Val(1);
    nmss = (int) number_missing_response->Val( index );

    delete grp;
  }//end while index


  for ( group = 0; group < number_of_groups; group++ )
  {
    delete response_means[ group ];
    delete n_items[ group ];
  }
  delete [] response_means;
  delete [] n_items;


}//end



void BayesianMVHLM::samplerMissingRandomPredictorsInitialPoint()
{
  int i, j, k, index, group, obs, nmss, p_dim;
  double diff;

  //compute means for each group ignoring missing data
  p_dim = (int) ( ((double) random_predictors[ 0 ][ 0 ]->Col()) / ((double) response[ 0 ]->Row() ) );

  CVector x_pred( p_dim );
  DistributionParameter * grp;
  CVector ** n_items;
  CVector missing_in_group( number_of_groups );
  
  missing_in_group.setToZero();

  randomEffects_means = new CVector * [ number_of_groups ];
  randomEffects_vars = new CVector * [ number_of_groups ];

  n_items = new CVector * [ number_of_groups ];
  for ( group = 0; group < number_of_groups; group++ )
  {
    randomEffects_means[ group ] = new CVector( p_dim );
    randomEffects_vars[ group ] = new CVector( p_dim );
    n_items[ group ] = new CVector( p_dim );
  
    randomEffects_means[ group ]->setToZero();
    randomEffects_vars[ group ]->setToZero();
    n_items[ group ]->setToZero();
  }


  index = 0;
  group = 0;
  obs = 0;
  nmss = (int) number_missing_random_predictors->Val( 0 );

  while ( index < number_of_observations )
  {
    for ( i = 0; i < p_dim; i++ )
    {
      x_pred = random_predictors[ group ][ obs ]->getRow( 0 ).subVector( 0, p_dim - 1 );

      if ( nmss > 0 )
      {
        if ( !missing_random_predictors_components[ group ][ obs ]->includes( ((double) i ) ) )
        {
          randomEffects_means[ group ]->Val( i ) += x_pred.Val( i );
          n_items[ group ]->Val(i) += 1;
        }
        else
	{
          missing_in_group.Val( group ) += 1;
        }
      }
      else
      {
        randomEffects_means[ group ]->Val( i ) += x_pred.Val( i );
        n_items[ group ]->Val(i) += 1;
      }
    }//end for i

    index++;
    grp = new DistributionParameter( groupIndex( index ) );
    group = (int) grp->getVector().Val(0);
    obs = (int) grp->getVector().Val(1);
    nmss = (int) number_missing_random_predictors->Val( index );

    delete grp;

  }//end while index

  //compute means
  for ( group = 0; group < number_of_groups; group++ )
  {
    if ( missing_in_group.Val( group ) > 0 )
    {
      for ( i = 0; i < p_dim; i++ )
      {
        if ( n_items[ group ]->Val( i ) > 0 )
	{
          randomEffects_means[ group ]->Val( i ) = randomEffects_means[ group ]->Val( i ) / n_items[ group ]->Val( i );
        }        
      }
    }
  }//end for


  //compute variances
  index = 0;
  group = 0;
  obs = 0;
  while ( index < number_of_observations )
  {
    for ( i = 0; i < p_dim; i++ )
    {
      x_pred = random_predictors[ group ][ obs ]->getRow( 0 ).subVector( 0, p_dim - 1 );

      if ( nmss > 0 )
      {
        if ( !missing_random_predictors_components[ group ][ obs ]->includes( ((double) i ) ) )
        {
          diff = ( x_pred.Val( i ) - randomEffects_means[ group ]->Val( i ) );
          randomEffects_vars[ group ]->Val( i ) += ( diff * diff );
        }
      }
      else
      {
        diff = ( x_pred.Val( i ) - randomEffects_means[ group ]->Val( i ) );
        randomEffects_vars[ group ]->Val( i ) += ( diff * diff );
      }
    }//end for i

    index++;
    grp = new DistributionParameter( groupIndex( index ) );
    group = (int) grp->getVector().Val(0);
    obs = (int) grp->getVector().Val(1);
    nmss = (int) number_missing_random_predictors->Val( index );

    delete grp;

  }//end while index

  //compute vars (divide)
  for ( group = 0; group < number_of_groups; group++ )
  {
    if ( missing_in_group.Val( group ) > 0 )
    {
      for ( i = 0; i < p_dim; i++ )
      {
        if ( n_items[ group ]->Val( i ) > 0 )
	{
          randomEffects_vars[ group ]->Val( i ) = randomEffects_vars[ group ]->Val( i ) / n_items[ group ]->Val( i );
        }        
      }
    }
  }//end for




  //ready to impute by means
  index = 0;
  group = 0;
  obs = 0;
  nmss = (int) number_missing_random_predictors->Val( 0 );

  while ( index < number_of_observations )
  {
    if ( nmss > 0 )
    {   
      for ( i = 0; i < nmss; i++ )
      {
        k = (int) missing_random_predictors_components[ group ][ obs ]->Val( i );
        for ( j = 0; j < response[ group ]->Row(); j++ )
        {
          random_predictors[ group ][ obs ]->Val( j, j * p_dim + k ) = randomEffects_means[ group ]->Val( k );
        } 
 
        //printf("initial random predcitors[%d][%d] is \n", group, obs );
        //random_predictors[ group ][ obs ]->Print();  

      }
    }//end if

    index++;
    grp = new DistributionParameter( groupIndex( index ) );
    group = (int) grp->getVector().Val(0);
    obs = (int) grp->getVector().Val(1);
    nmss = (int) number_missing_random_predictors->Val( index );

    delete grp;
  }//end while index


  for ( group = 0; group < number_of_groups; group++ )
  {
    delete n_items[ group ];
  }
  delete [] n_items;


}//end



void BayesianMVHLM::samplerMissingFixedPredictorsInitialPoint()
{
  int i, j, k, index, group, obs, nmss, p_dim;
  double diff;

  //compute means for each group ignoring missing data
  p_dim = (int) ( ((double) fixed_predictors[ 0 ][ 0 ]->Col()) / ((double) response[ 0 ]->Row() ) );

  CVector x_pred( p_dim );
  DistributionParameter * grp;
  CVector ** n_items;
  CVector missing_in_group( number_of_groups );
  
  missing_in_group.setToZero();

  fixedEffects_means = new CVector * [ number_of_groups ];
  fixedEffects_vars = new CVector * [ number_of_groups ];

  n_items = new CVector * [ number_of_groups ];
  for ( group = 0; group < number_of_groups; group++ )
  {
    fixedEffects_means[ group ] = new CVector( p_dim );
    fixedEffects_vars[ group ] = new CVector( p_dim );
    n_items[ group ] = new CVector( p_dim );
  
    fixedEffects_means[ group ]->setToZero();
    fixedEffects_vars[ group ]->setToZero();
    n_items[ group ]->setToZero();
  }


  index = 0;
  group = 0;
  obs = 0;
  nmss = (int) number_missing_fixed_predictors->Val( 0 );

  while ( index < number_of_observations )
  {
    for ( i = 0; i < p_dim; i++ )
    {
      x_pred = fixed_predictors[ group ][ obs ]->getRow( 0 ).subVector( 0, p_dim - 1 );

      if ( nmss > 0 )
      {
        if ( !missing_fixed_predictors_components[ group ][ obs ]->includes( ((double) i ) ) )
        {
          fixedEffects_means[ group ]->Val( i ) += x_pred.Val( i );
          n_items[ group ]->Val(i) += 1;
        }
        else
	{
          missing_in_group.Val( group ) += 1;
        }
      }
      else
      {
        fixedEffects_means[ group ]->Val( i ) += x_pred.Val( i );
        n_items[ group ]->Val(i) += 1;
      }
    }//end for i

    index++;
    grp = new DistributionParameter( groupIndex( index ) );
    group = (int) grp->getVector().Val(0);
    obs = (int) grp->getVector().Val(1);
    nmss = (int) number_missing_fixed_predictors->Val( index );

    delete grp;

  }//end while index

  //compute means
  for ( group = 0; group < number_of_groups; group++ )
  {
    if ( missing_in_group.Val( group ) > 0 )
    {
      for ( i = 0; i < p_dim; i++ )
      {
        if ( n_items[ group ]->Val( i ) > 0 )
	{
          fixedEffects_means[ group ]->Val( i ) = fixedEffects_means[ group ]->Val( i ) / n_items[ group ]->Val( i );
        }        
      }
    }
  }//end for


  //compute variances
  index = 0;
  group = 0;
  obs = 0;
  while ( index < number_of_observations )
  {
    for ( i = 0; i < p_dim; i++ )
    {
      x_pred = fixed_predictors[ group ][ obs ]->getRow( 0 ).subVector( 0, p_dim - 1 );

      if ( nmss > 0 )
      {
        if ( !missing_fixed_predictors_components[ group ][ obs ]->includes( ((double) i ) ) )
        {
          diff = ( x_pred.Val( i ) - fixedEffects_means[ group ]->Val( i ) );
          fixedEffects_vars[ group ]->Val( i ) += ( diff * diff );
        }
      }
      else
      {
        diff = ( x_pred.Val( i ) - fixedEffects_means[ group ]->Val( i ) );
        fixedEffects_vars[ group ]->Val( i ) += ( diff * diff );
      }
    }//end for i

    index++;
    grp = new DistributionParameter( groupIndex( index ) );
    group = (int) grp->getVector().Val(0);
    obs = (int) grp->getVector().Val(1);
    nmss = (int) number_missing_fixed_predictors->Val( index );

    delete grp;

  }//end while index

  //compute vars (divide)
  for ( group = 0; group < number_of_groups; group++ )
  {
    if ( missing_in_group.Val( group ) > 0 )
    {
      for ( i = 0; i < p_dim; i++ )
      {
        if ( n_items[ group ]->Val( i ) > 0 )
	{
          fixedEffects_vars[ group ]->Val( i ) = fixedEffects_vars[ group ]->Val( i ) / n_items[ group ]->Val( i );
        }        
      }
    }
  }//end for



  //ready to impute by means
  index = 0;
  group = 0;
  obs = 0;
  nmss = (int) number_missing_fixed_predictors->Val( 0 );

  while ( index < number_of_observations )
  {
    if ( nmss > 0 )
    {   
      for ( i = 0; i < nmss; i++ )
      {
        k = (int) missing_fixed_predictors_components[ group ][ obs ]->Val( i );
        for ( j = 0; j < response[ group ]->Row(); j++ )
        {
          //printf( "p_dim = %d, k = %d,  j = %d, col = %d\n", p_dim, k, j, j * p_dim + k ); fflush(stdout);

          fixed_predictors[ group ][ obs ]->Val( j, j * p_dim + k ) = fixedEffects_means[ group ]->Val( k );
        }

        //printf("initial fixed predcitors[%d][%d] is \n", group, obs );
        //fixed_predictors[ group ][ obs ]->Print();  

      }
    }//end if

    index++;
    grp = new DistributionParameter( groupIndex( index ) );
    group = (int) grp->getVector().Val(0);
    obs = (int) grp->getVector().Val(1);
    nmss = (int) number_missing_fixed_predictors->Val( index );

    delete grp;
  }//end while index


  for ( group = 0; group < number_of_groups; group++ )
  {
    delete n_items[ group ];
  }
  delete [] n_items;


}//end




DistributionParameter BayesianMVHLM::groupIndex( int index )
{
  bool found;
  int i;

  CVector grp( 2 );

  //compute the corresponding group
  found = false;
  i = 0;
  grp.Val(0) = -1;
  grp.Val(1) = -1;



  while ( !found )
  {
    if ( index <= final_group_index->Val(i) )
    {
      found = true;
      grp.Val(0) = i;
    }
    else
    {
      i++;
    }
    if ( i == number_of_groups )
    {
      found = true;
    }
  }//end loop

  //find relative index
  if ( grp.Val(0) >= 0 )
  {
    if ( grp.Val(0) == 0 )
    {
      grp.Val(1) = index;
    }
    else
    {
      grp.Val(1) = index - 1 - final_group_index->Val( ( (int) grp.Val(0) - 1 ) );
    }

  }

  DistributionParameter grp_par( grp );
  return grp_par;
}//end



void BayesianMVHLM::updateRegressionWeight( int index ) throw( rtErr )
{
  int group, obs;

  DistributionParameter grp( groupIndex( index ) );

  group = (int) grp.getVector().Val(0);
  obs = (int) grp.getVector().Val(1);

  if ( group >= 0 )
  {
    tau2_error_weights[ group ]->Val( obs ) = 1 / tau2_errors[ index ]->lastItemDrawn();
  }
  else
  {
    printf( "BayesianMVHLM::updateRegressionWeight: Index [%d] out of range.\n", index );
    char the_error[] = "BayesianMVHLM::updateRegressionWeight: Index [%d] out of range.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}//end


void BayesianMVHLM::updateGroupRegressionWeight( int grp )
{
  int i;

  for ( i = 0; i < response[ grp ]->Col(); i++ )
  {
    tau2_error_weights[ grp ]->Val( i ) = 1 / tau2_groups[ grp ]->lastItemDrawn();
  }
}//end


void BayesianMVHLM::initializeTemporaryStructures()
{
  int i, j, k, len_digits_vars, len_digits_obs;

  //create working matrices
  if ( random_effects )
  {
    xTx = new CMatrix * [ number_of_groups ];
    xTy = new CVector * [ number_of_groups ];
    for ( i = 0; i < number_of_groups; i++ )
    {
      xTx[i] = new CMatrix( beta[0]->dimension(), beta[0]->dimension() );
      xTy[i] = new CVector( beta[0]->dimension() ); 

      xTx[i]->setToZero();
      for ( j = 0; j < response[i]->Col(); j++ )
      {
        if ( !invWishart_sigma2 )
        {
#ifdef FIX1
          CMatrix tmpmat = random_predictors[i][j]->xTransposedX();
          xTx[i]->add( tmpmat );
#else
          xTx[i]->add( random_predictors[i][j]->xTransposedX() );
#endif
//          xTx[i]->add( random_predictors[i][j]->xTransposedX() );
        }
        else if ( common_sigma )
        {
#ifdef FIX1
          CMatrix tmpmat1 = sigma2InvWishart[0]->inverseLastItemDrawn();
          CMatrix tmpmat2 = random_predictors[i][j]->xTransposedX( tmpmat1 );
          xTx[i]->add( tmpmat2 );
#else
          xTx[i]->add( random_predictors[i][j]->xTransposedX( sigma2InvWishart[0]->inverseLastItemDrawn() ) );
#endif
//          xTx[i]->add( random_predictors[i][j]->xTransposedX( sigma2InvWishart[0]->inverseLastItemDrawn() ) );
        }
        else
        { 
#ifdef FIX1
          CMatrix tmpmat1 = sigma2InvWishart[i]->inverseLastItemDrawn();
          CMatrix tmpmat2 = random_predictors[i][j]->xTransposedX( tmpmat1 );
          xTx[i]->add( tmpmat2 );
#else
          xTx[i]->add( random_predictors[i][j]->xTransposedX( sigma2InvWishart[i]->inverseLastItemDrawn() ) );
#endif
//          xTx[i]->add( random_predictors[i][j]->xTransposedX( sigma2InvWishart[i]->inverseLastItemDrawn() ) );
        }
      }//end for j


      xTy[i]->setToZero();
      for ( j = 0; j < response[i]->Col(); j++ )
      {
        if ( !invWishart_sigma2 )
        {
#ifdef FIX1
          CVector tmpvec = random_predictors[i][j]->T() * response[i]->getColumn( j );
          xTy[i]->add( tmpvec );
#else
          xTy[i]->add( random_predictors[i][j]->T() * response[i]->getColumn( j ) );
#endif
//          xTy[i]->add( random_predictors[i][j]->T() * response[i]->getColumn( j ) );
        }
        else if ( common_sigma )
    	{
#ifdef FIX1
          CVector tmpvec = random_predictors[i][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * response[i]->getColumn( j ) );
          xTy[i]->add( tmpvec );
#else
          xTy[i]->add( random_predictors[i][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * response[i]->getColumn( j ) ) );
#endif
//          xTy[i]->add( random_predictors[i][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * response[i]->getColumn( j ) ) );
        }
        else
	    {
#ifdef FIX1
          CVector tmpvec = random_predictors[i][j]->T() * ( sigma2InvWishart[i]->inverseLastItemDrawn() * response[i]->getColumn( j ) );
          xTy[i]->add( tmpvec );
#else
          xTy[i]->add( random_predictors[i][j]->T() * ( sigma2InvWishart[i]->inverseLastItemDrawn() * response[i]->getColumn( j ) ) );
#endif
//          xTy[i]->add( random_predictors[i][j]->T() * ( sigma2InvWishart[i]->inverseLastItemDrawn() * response[i]->getColumn( j ) ) );
        }

      }//end for j

    }//end for i
  }//end random effects

  if ( fixed_effects )
  {
    sMTm = new CMatrix( gamma->dimension(), gamma->dimension() );
    sMTy = new CVector( gamma->dimension() );
    sMTm->setToZero();
    sMTy->setToZero();

    mTm = new CMatrix * [ number_of_groups ];
    mTy = new CVector * [ number_of_groups ];
    for ( i = 0; i < number_of_groups; i++ )
    {
      mTm[i] = new CMatrix( gamma->dimension(), gamma->dimension() );
      mTy[i] = new CVector( gamma->dimension() );


      mTm[i]->setToZero();
      for ( j = 0; j < response[i]->Col(); j++ )
      {
        if ( !invWishart_sigma2 )
        {
#ifdef FIX1
          CMatrix tmpmat = fixed_predictors[i][j]->xTransposedX();
          mTm[i]->add( tmpmat );
#else
          mTm[i]->add( fixed_predictors[i][j]->xTransposedX() );
#endif
//          mTm[i]->add( fixed_predictors[i][j]->xTransposedX() );
        }
        else if ( common_sigma )
        {
#ifdef FIX1
          CMatrix tmpmat1 = sigma2InvWishart[0]->inverseLastItemDrawn();
          CMatrix tmpmat2 = fixed_predictors[i][j]->xTransposedX( tmpmat1 );
          mTm[i]->add( tmpmat2 );
#else
          mTm[i]->add( fixed_predictors[i][j]->xTransposedX( sigma2InvWishart[0]->inverseLastItemDrawn() ) );
#endif
//          mTm[i]->add( fixed_predictors[i][j]->xTransposedX( sigma2InvWishart[0]->inverseLastItemDrawn() ) );
        }
        else
        { 
#ifdef FIX1
          CMatrix tmpmat1 = sigma2InvWishart[i]->inverseLastItemDrawn();
          CMatrix tmpmat2 = fixed_predictors[i][j]->xTransposedX( tmpmat1 );
          mTm[i]->add( tmpmat2 );
#else
          mTm[i]->add( fixed_predictors[i][j]->xTransposedX( sigma2InvWishart[i]->inverseLastItemDrawn() ) );
#endif
//          mTm[i]->add( fixed_predictors[i][j]->xTransposedX( sigma2InvWishart[i]->inverseLastItemDrawn() ) );
        }
      }//end for j


      mTy[i]->setToZero();
      for ( j = 0; j < response[i]->Col(); j++ )
      {
        if ( !invWishart_sigma2 )
        {
#ifdef FIX1
          CVector tmpvec = fixed_predictors[i][j]->T() * response[i]->getColumn( j );
          mTy[i]->add( tmpvec );
#else
          mTy[i]->add( fixed_predictors[i][j]->T() * response[i]->getColumn( j ) );
#endif
//          mTy[i]->add( fixed_predictors[i][j]->T() * response[i]->getColumn( j ) );
        }
        else if ( common_sigma )
	    {
#ifdef FIX1
          CVector tmpvec = fixed_predictors[i][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * response[i]->getColumn( j ) );
          mTy[i]->add( tmpvec );
#else
          mTy[i]->add( fixed_predictors[i][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * response[i]->getColumn( j ) ) );
#endif
//          mTy[i]->add( fixed_predictors[i][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * response[i]->getColumn( j ) ) );
        }
        else
    	{
#ifdef FIX1
          CVector tmpvec = fixed_predictors[i][j]->T() * ( sigma2InvWishart[i]->inverseLastItemDrawn() * response[i]->getColumn( j ) );
          mTy[i]->add( tmpvec );
#else
          mTy[i]->add( fixed_predictors[i][j]->T() * ( sigma2InvWishart[i]->inverseLastItemDrawn() * response[i]->getColumn( j ) ) );
#endif
//          mTy[i]->add( fixed_predictors[i][j]->T() * ( sigma2InvWishart[i]->inverseLastItemDrawn() * response[i]->getColumn( j ) ) );
        }

      }//end for j

     sMTm->add( (*(mTm[i])) );
     sMTy->add( (*(mTy[i])) );

    }//end loop
  }//end fixed effects

  //second stage
  if ( second_stage_effects )
  {
    zTvIz = new CMatrix( alpha->dimension(), alpha->dimension() );
    zTvIb = new CVector( alpha->dimension() );
    zTvIz->setToZero();
    zTvIb->setToZero();

    zTz = new CMatrix * [ number_of_groups ];
    zTb = new CVector * [ number_of_groups ];
    for ( i = 0; i < number_of_groups; i++ )
    {
      zTz[i] = new CMatrix( alpha->dimension(), alpha->dimension() );
      zTb[i] = new CVector( alpha->dimension() );

      //note: beta cov has been updated according to whether invWishart and/or t_beta are in place
#ifdef FIX1
      CMatrix tmpmat = beta[i]->initialInverseCovariance();
      (*(zTz[i])) = beta_predictors[i]->xTransposedX( tmpmat );
#else
      (*(zTz[i])) = beta_predictors[i]->xTransposedX( beta[i]->initialInverseCovariance() );
#endif
//      (*(zTz[i])) = beta_predictors[i]->xTransposedX( beta[i]->initialInverseCovariance() );
      (*(zTb[i])) = beta_predictors[i]->T() * beta[i]->initialInverseCovariance() * beta[i]->lastDraw().getVector();

      zTvIz->add( (*(zTz[i])) );
      zTvIb->add( (*(zTb[i])) );
    }//end loop
  }

  //initialize residuals
  residuals = new CMatrix * [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    residuals[i] = new CMatrix( response[i]->Row(), response[i]->Col() );
    gibbsUpdateResiduals( i );
  }

  if ( random_effects )
  {
    second_residuals = new CVector * [ number_of_groups ];  
    for ( i = 0; i < number_of_groups; i++ )
    {
      second_residuals[i] = new CVector( beta[i]->dimension() );
      gibbsUpdateSecondStageResiduals( i );
    }
  }


  //create indices for variables
  len_digits_vars = ( (int) ( log( (double) number_of_variables ) / log( 10.0 ) ) ) + 1;

  distr_map = new char * [ number_of_variables ];
  k = 0;

  //printf("Bayes: n vars is %d\n", number_of_variables); fflush(stdout);

  if ( random_effects )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      distr_map[ k ] = new char [ 6 + len_digits_vars ];
      sprintf( distr_map[ k ], "beta:%d", i );
      k++;
    }
  }

  if ( fixed_effects )
  {
    distr_map[ k ] = new char [ 6 ];
    sprintf( distr_map[ k ], "gamma" );
    k++;
  }

  if ( second_stage_effects )
  {
    distr_map[k] = new char [ 6 ];
    sprintf( distr_map[ k ], "alpha" );
    k++;
  }

  if ( !known_sigma2 )
  {
    if ( common_sigma )
    {
      distr_map[k] = new char [ 6 ];
      sprintf( distr_map[k], "sigma" );
      k++;
    }
    else
    {
      for ( i = 0; i < number_of_groups; i++ )
      {
        distr_map[k] = new char [ 7 + len_digits_vars ];
        sprintf( distr_map[k], "sigma:%d", i );
        k++;
      }
    }
  }     

  if ( random_effects )
  {
    distr_map[k] = new char [ 5 ];
    sprintf( distr_map[k], "tau2" );
    k++;
  }


  len_digits_obs = ( (int) ( log( (double) number_of_observations ) / log( 10.0 ) ) ) + 1;

  if ( t_lkhd )
  {
    for ( i = 0; i < number_of_observations; i++ )
    {
      //six characters for "error:", remainder for characters for index
      distr_map[ k ] = new char [ 7 + len_digits_obs ];
      sprintf( distr_map[ k ], "tau2e:%d", i );
      k++;
    }
    tau2_error_weights = new CVector * [ number_of_groups ];
    for ( i = 0; i < number_of_groups; i++ )
    {
      tau2_error_weights[i] = new CVector( response[i]->Col() );
    }
  }

  if ( t_group_lkhd )
  {
    tau2_error_weights = new CVector * [ number_of_groups ];
    tau2_group_weights = new CVector( number_of_groups );
    for ( i = 0; i < number_of_groups; i++ )
    {
      //six characters for "error:", remainder for characters for index
      distr_map[ k ] = new char [ 7 + len_digits_vars ];
      sprintf( distr_map[ k ], "tau2g:%d", i );
      k++;
      tau2_error_weights[i] = new CVector( response[i]->Col() );
    }
  }


  if ( fixed_effects && t_gamma )
  {
    distr_map[ k ] = new char [ 7 ];
    sprintf( distr_map[ k ], "tau2gm" );
    k++;   
  }


  if ( random_effects && t_beta )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      distr_map[ k ] = new char [ 7 + len_digits_vars ];
      sprintf( distr_map[ k ], "tau2b:%d", i );
      k++;
    }
  }

  if ( second_stage_effects && t_alpha )
  {
    distr_map[ k ] = new char [ 6 ];
    sprintf( distr_map[ k ], "tau2a" );
    k++;   
  }


  if ( keep_missing_data )
  {
    if ( number_missing_response != NULL )
    {
      for ( i = 0; i < number_of_observations; i++ )
      {
        if ( number_missing_response->Val( i ) > 0 )
        {
          distr_map[ k ] = new char[ 10 + len_digits_obs ];
          sprintf( distr_map[ k ], "missingR:%d", i );
          k++;
        }
      }//end for i
    }//end if missing responses

    if ( number_missing_random_predictors != NULL )
    {
      for ( i = 0; i < number_of_observations; i++ )
      {
        if ( number_missing_random_predictors->Val( i ) > 0 )
        {
          distr_map[ k ] = new char[ 11 + len_digits_obs ];
          sprintf( distr_map[ k ], "missingRP:%d", i );
          k++;
        }
      }//end for i
    }//end if missing random predictors


    if ( number_missing_fixed_predictors != NULL )
    {
      for ( i = 0; i < number_of_observations; i++ )
      {
        if ( number_missing_fixed_predictors->Val( i ) > 0 )
        {
          distr_map[ k ] = new char[ 11 + len_digits_obs ];
          sprintf( distr_map[ k ], "missingFP:%d", i );
          k++;
        }
      }//end for i
    }//end if missing fixed predictors

  }//end if missing data


  //printf( "BHL: total vars = %d.  Number vars = %d. The map is:\n", k, number_of_variables  );
  //for ( i = 0; i < k; i++ )
  //{
  //  printf("map[%d] = %s\n", i, distr_map[ i ] );
  //}
  //fflush( stdout );


}//end



void BayesianMVHLM::gibbsUpdateWorkingMatrices()
{
  int i;

  //weighted (or transformed) predictors. 
  //update must be made before updating beta
  if ( random_effects )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      gibbsUpdateWorkingMatrices( i, "random_effects" );
    }
  }

  if ( fixed_effects )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      gibbsUpdateWorkingMatrices( i, "fixed_effects" );
    }
  }

}//end


void BayesianMVHLM::gibbsUpdateWorkingMatrices( int index, char * type )
{
  int j;
  double precision;

  if ( t_lkhd || t_group_lkhd )
  {
    //weighted (or transformed) predictors. 
    //update must be made before updating beta
    if ( !strcmp( type, "random_effects" ) )
    {
      xTx[ index ]->setToZero();
      xTy[ index ]->setToZero();


      xTx[ index ]->setToZero();
      for ( j = 0; j < response[ index ]->Col(); j++ )
      {
        if ( !invWishart_sigma2 )
        {
#ifdef FIX1
          CMatrix tmpmat = tau2_error_weights[ index ]->Val(j) * random_predictors[ index ][j]->xTransposedX();
          xTx[ index ]->add( tmpmat );
#else
          xTx[ index ]->add( tau2_error_weights[ index ]->Val(j) * random_predictors[ index ][j]->xTransposedX() );
#endif
//          xTx[ index ]->add( tau2_error_weights[ index ]->Val(j) * random_predictors[ index ][j]->xTransposedX() );
        }
        else if ( common_sigma )
        {
#ifdef FIX1
          CMatrix tmpmat1 = sigma2InvWishart[0]->inverseLastItemDrawn();
          CMatrix tmpmat2 = tau2_error_weights[ index ]->Val(j) * random_predictors[ index ][j]->xTransposedX( tmpmat1 );
          xTx[ index ]->add( tmpmat2 );
#else
          xTx[ index ]->add( tau2_error_weights[ index ]->Val(j) * random_predictors[ index ][j]->xTransposedX( sigma2InvWishart[0]->inverseLastItemDrawn() ) );
#endif
//          xTx[ index ]->add( tau2_error_weights[ index ]->Val(j) * random_predictors[ index ][j]->xTransposedX( sigma2InvWishart[0]->inverseLastItemDrawn() ) );
        }
        else
        { 
#ifdef FIX1
          CMatrix tmpmat1 = sigma2InvWishart[ index ]->inverseLastItemDrawn();
          CMatrix tmpmat2 = tau2_error_weights[ index ]->Val(j) * random_predictors[ index ][j]->xTransposedX( tmpmat1 );
          xTx[ index ]->add( tmpmat2 );
#else
          xTx[ index ]->add( tau2_error_weights[ index ]->Val(j) * random_predictors[ index ][j]->xTransposedX( sigma2InvWishart[ index ]->inverseLastItemDrawn() ) );
#endif
//          xTx[ index ]->add( tau2_error_weights[ index ]->Val(j) * random_predictors[ index ][j]->xTransposedX( sigma2InvWishart[ index ]->inverseLastItemDrawn() ) );
        }
      }//end for j


      xTy[ index ]->setToZero();
      for ( j = 0; j < response[ index ]->Col(); j++ )
      {
        if ( !invWishart_sigma2 )
        {
          if ( fixed_effects )
          {
#ifdef FIX1
            CVector tmpvec = ( random_predictors[ index ][j]->T() * ( response[ index ]->getColumn( j ) - ( (*(fixed_predictors[ index ][j])) * gamma->lastDraw().getVector() ) ) ) * tau2_error_weights[ index ]->Val(j);
            xTy[ index ]->add( tmpvec );
#else
            xTy[ index ]->add( ( random_predictors[ index ][j]->T() * ( response[ index ]->getColumn( j ) - ( (*(fixed_predictors[ index ][j])) * gamma->lastDraw().getVector() ) ) ) * tau2_error_weights[ index ]->Val(j) );
#endif
//            xTy[ index ]->add( ( random_predictors[ index ][j]->T() * ( response[ index ]->getColumn( j ) - ( (*(fixed_predictors[ index ][j])) * gamma->lastDraw().getVector() ) ) ) * tau2_error_weights[ index ]->Val(j) );
	      }
          else
	      {
#ifdef FIX1
            CVector tmpvec = ( random_predictors[ index ][j]->T() * response[ index ]->getColumn( j ) ) * tau2_error_weights[ index ]->Val(j);
            xTy[ index ]->add( tmpvec );
#else
            xTy[ index ]->add( ( random_predictors[ index ][j]->T() * response[ index ]->getColumn( j ) ) * tau2_error_weights[ index ]->Val(j) );
#endif
//            xTy[ index ]->add( ( random_predictors[ index ][j]->T() * response[ index ]->getColumn( j ) ) * tau2_error_weights[ index ]->Val(j) );
          }
        }
        else if ( common_sigma )
	    {
          if ( fixed_effects )
          {
#ifdef FIX1
            CVector tmpvec = ( random_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(fixed_predictors[ index ][j])) * gamma->lastDraw().getVector() ) ) ) ) * tau2_error_weights[ index ]->Val(j);
            xTy[ index ]->add( tmpvec );
#else
            xTy[ index ]->add( ( random_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(fixed_predictors[ index ][j])) * gamma->lastDraw().getVector() ) ) ) ) * tau2_error_weights[ index ]->Val(j) );
#endif
//            xTy[ index ]->add( ( random_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(fixed_predictors[ index ][j])) * gamma->lastDraw().getVector() ) ) ) ) * tau2_error_weights[ index ]->Val(j) );
          }
          else
	      {
#ifdef FIX1
            CVector tmpvec = ( random_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) ) ) * tau2_error_weights[ index ]->Val(j);
            xTy[ index ]->add( tmpvec );
#else
            xTy[ index ]->add( ( random_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) ) ) * tau2_error_weights[ index ]->Val(j) );
#endif
//            xTy[ index ]->add( ( random_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) ) ) * tau2_error_weights[ index ]->Val(j) );
          }
        }
        else
	    {
          if ( fixed_effects )
          {
#ifdef FIX1
            CVector tmpvec = ( random_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(fixed_predictors[ index ][j])) * gamma->lastDraw().getVector() ) ) ) ) * tau2_error_weights[ index ]->Val(j);
            xTy[ index ]->add( tmpvec );
#else
            xTy[ index ]->add( ( random_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(fixed_predictors[ index ][j])) * gamma->lastDraw().getVector() ) ) ) ) * tau2_error_weights[ index ]->Val(j) );
#endif
//            xTy[ index ]->add( ( random_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(fixed_predictors[ index ][j])) * gamma->lastDraw().getVector() ) ) ) ) * tau2_error_weights[ index ]->Val(j) );
          }
          else
	      {
#ifdef FIX1
            CVector tmpvec = ( random_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) ) ) * tau2_error_weights[ index ]->Val(j);
            xTy[ index ]->add( tmpvec );
#else
            xTy[ index ]->add( ( random_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) ) ) * tau2_error_weights[ index ]->Val(j) );
#endif
//            xTy[ index ]->add( ( random_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) ) ) * tau2_error_weights[ index ]->Val(j) );
          }
	    }
      }//end for j

    }//end if random effects

    if ( !strcmp( type, "fixed_effects" ) )
    {
      //update sum matrix: substract current value 
#ifdef FIX1
      CVector tmpvec = (*(mTy[ index ])) * ( -1.0 );
      sMTy->add( tmpvec );
#else
      sMTy->add( (*(mTy[ index ])) * ( -1.0 ) );
#endif
//      sMTy->add( (*(mTy[ index ])) * ( -1.0 ) );
#ifdef FIX1
      CMatrix tmpmat = ( -1.0 * ( *(mTm[ index ]) ) );
      sMTm->add( tmpmat );
#else
      sMTm->add( ( -1.0 * ( *(mTm[ index ]) ) ) );
#endif
//      sMTm->add( ( -1.0 * ( *(mTm[ index ]) ) ) );

      mTm[ index ]->setToZero();
      mTy[ index ]->setToZero();

      for ( j = 0; j < response[ index ]->Col(); j++ )
      {
        if ( common_sigma || invWishart_sigma2 )
	{
          precision = tau2_error_weights[ index ]->Val(j);
        }
        else
	{
          precision = tau2_error_weights[ index ]->Val(j)/ sigma2->lastDrawFromComponent( index ).getScalar();
        }

        if ( !invWishart_sigma2 )
        {
#ifdef FIX1
          CMatrix tmpmat = precision * fixed_predictors[ index ][j]->xTransposedX();
          mTm[ index ]->add( tmpmat );
#else
          mTm[ index ]->add( precision * fixed_predictors[ index ][j]->xTransposedX() );
#endif
//          mTm[ index ]->add( precision * fixed_predictors[ index ][j]->xTransposedX() );
        }
        else if ( common_sigma )
        {
#ifdef FIX1
          CMatrix tmpmat1 = sigma2InvWishart[0]->inverseLastItemDrawn();
          CMatrix tmpmat2 = precision * fixed_predictors[ index ][j]->xTransposedX( tmpmat1 );
          mTm[ index ]->add( tmpmat2 );
#else
          mTm[ index ]->add( precision * fixed_predictors[ index ][j]->xTransposedX( sigma2InvWishart[0]->inverseLastItemDrawn() ) );
#endif
//          mTm[ index ]->add( precision * fixed_predictors[ index ][j]->xTransposedX( sigma2InvWishart[0]->inverseLastItemDrawn() ) );
        }
        else
        { 
#ifdef FIX1
          CMatrix tmpmat1 = sigma2InvWishart[ index ]->inverseLastItemDrawn();
          CMatrix tmpmat2 = precision * fixed_predictors[ index ][j]->xTransposedX( tmpmat1 );
          mTm[ index ]->add( tmpmat2 );
#else
          mTm[ index ]->add( precision * fixed_predictors[ index ][j]->xTransposedX( sigma2InvWishart[ index ]->inverseLastItemDrawn() ) );
#endif
//          mTm[ index ]->add( precision * fixed_predictors[ index ][j]->xTransposedX( sigma2InvWishart[ index ]->inverseLastItemDrawn() ) );
        }
      }//end for j


      for ( j = 0; j < response[ index ]->Col(); j++ )
      {
        if ( common_sigma || invWishart_sigma2 )
	{
          precision = tau2_error_weights[ index ]->Val(j);
        }
        else
	{
          precision = tau2_error_weights[ index ]->Val(j) / sigma2->lastDrawFromComponent( index ).getScalar();
        }

        if ( !invWishart_sigma2 )
        {
          if ( random_effects )
	      {
#ifdef FIX1
            CVector tmpvec = ( fixed_predictors[ index ][j]->T() * ( response[ index ]->getColumn( j ) - ( (*(random_predictors[ index ][j])) * beta[ index ]->lastDraw().getVector() ) ) ) * precision;
            mTy[ index ]->add( tmpvec );
#else
            mTy[ index ]->add( ( fixed_predictors[ index ][j]->T() * ( response[ index ]->getColumn( j ) - ( (*(random_predictors[ index ][j])) * beta[ index ]->lastDraw().getVector() ) ) ) * precision );
#endif
//            mTy[ index ]->add( ( fixed_predictors[ index ][j]->T() * ( response[ index ]->getColumn( j ) - ( (*(random_predictors[ index ][j])) * beta[ index ]->lastDraw().getVector() ) ) ) * precision );
          }
          else
	      {
#ifdef FIX1
            CVector tmpvec = ( fixed_predictors[ index ][j]->T() * response[ index ]->getColumn( j ) ) * precision;
            mTy[ index ]->add( tmpvec );
#else
            mTy[ index ]->add( ( fixed_predictors[ index ][j]->T() * response[ index ]->getColumn( j ) ) * precision );
#endif
//            mTy[ index ]->add( ( fixed_predictors[ index ][j]->T() * response[ index ]->getColumn( j ) ) * precision );
          }
        }
        else if ( common_sigma )
	    {
          if ( random_effects )
	      {
#ifdef FIX1
            CVector tmpvec = ( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(random_predictors[ index ][j])) * beta[ index ]->lastDraw().getVector() ) ) ) ) * precision;
            mTy[ index ]->add( tmpvec );
#else
            mTy[ index ]->add( ( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(random_predictors[ index ][j])) * beta[ index ]->lastDraw().getVector() ) ) ) ) * precision );
#endif
//            mTy[ index ]->add( ( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(random_predictors[ index ][j])) * beta[ index ]->lastDraw().getVector() ) ) ) ) * precision );
          }
          else
	      {
#ifdef FIX1
            CVector tmpvec = ( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) ) ) * precision;
            mTy[ index ]->add( tmpvec );
#else
            mTy[ index ]->add( ( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) ) ) * precision );
#endif
//            mTy[ index ]->add( ( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) ) ) * precision );
          }
        }
        else
	    {
          if ( random_effects )
	      {
#ifdef FIX1
            CVector tmpvec = ( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(random_predictors[ index ][j])) * beta[ index ]->lastDraw().getVector() ) ) ) ) * precision;
            mTy[ index ]->add( tmpvec );
#else
            mTy[ index ]->add( ( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(random_predictors[ index ][j])) * beta[ index ]->lastDraw().getVector() ) ) ) ) * precision );
#endif
//            mTy[ index ]->add( ( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(random_predictors[ index ][j])) * beta[ index ]->lastDraw().getVector() ) ) ) ) * precision );
	      }
          else
	      {
#ifdef FIX1
            CVector tmpvec = ( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) ) ) * precision;
            mTy[ index ]->add( tmpvec );
#else
            mTy[ index ]->add( ( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) ) ) * precision );
#endif
//           mTy[ index ]->add( ( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) ) ) * precision );
	      }
        }
      }//end for j

      //update sum matrix: add updated value
      sMTm->add( ( *(mTm[ index ]) ) );
      sMTy->add( ( *(mTy[ index ]) ) );

    }//end if fixed effects
  }//end if t lklhd
  else if ( !strcmp( type, "random_effects" ) )
  {
    if ( invWishart_sigma2 )
    {
      xTx[ index ]->setToZero();
      for ( j = 0; j < response[ index ]->Col(); j++ )
      {
        if ( common_sigma )
        {
#ifdef FIX1
          CMatrix tmpmat1 = sigma2InvWishart[0]->inverseLastItemDrawn();
          CMatrix tmpmat2 = random_predictors[ index ][j]->xTransposedX( tmpmat1 );
          xTx[ index ]->add( tmpmat2 );
#else
          xTx[ index ]->add( random_predictors[ index ][j]->xTransposedX( sigma2InvWishart[0]->inverseLastItemDrawn() ) );
#endif
//          xTx[ index ]->add( random_predictors[ index ][j]->xTransposedX( sigma2InvWishart[0]->inverseLastItemDrawn() ) );
        }
        else
        {
#ifdef FIX1
          CMatrix tmpmat1 = sigma2InvWishart[ index ]->inverseLastItemDrawn();
          CMatrix tmpmat2 = random_predictors[ index ][j]->xTransposedX( tmpmat1 );
          xTx[ index ]->add( tmpmat2 );
#else
          xTx[ index ]->add( random_predictors[ index ][j]->xTransposedX( sigma2InvWishart[ index ]->inverseLastItemDrawn() ) );
#endif
//          xTx[ index ]->add( random_predictors[ index ][j]->xTransposedX( sigma2InvWishart[ index ]->inverseLastItemDrawn() ) );
        }
      }//end for j
    }//end if

    if ( fixed_effects )
    {
      xTy[ index ]->setToZero();
      for ( j = 0; j < response[ index ]->Col(); j++ )
      {
        if ( !invWishart_sigma2 )
        {
#ifdef FIX1
          CVector tmpvec = random_predictors[ index ][j]->T() * ( response[ index ]->getColumn( j ) - ( (*(fixed_predictors[ index ][j])) * gamma->lastDraw().getVector() ) );
          xTy[ index ]->add( tmpvec );
#else
          xTy[ index ]->add( random_predictors[ index ][j]->T() * ( response[ index ]->getColumn( j ) - ( (*(fixed_predictors[ index ][j])) * gamma->lastDraw().getVector() ) ) );
#endif
//          xTy[ index ]->add( random_predictors[ index ][j]->T() * ( response[ index ]->getColumn( j ) - ( (*(fixed_predictors[ index ][j])) * gamma->lastDraw().getVector() ) ) );
        }
        else if ( common_sigma )
        {
#ifdef FIX1
          CVector tmpvec = random_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(fixed_predictors[ index ][j])) * gamma->lastDraw().getVector() ) ) );
          xTy[ index ]->add( tmpvec );
#else
          xTy[ index ]->add( random_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(fixed_predictors[ index ][j])) * gamma->lastDraw().getVector() ) ) ) );
#endif
//          xTy[ index ]->add( random_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(fixed_predictors[ index ][j])) * gamma->lastDraw().getVector() ) ) ) );
        }
        else
        {
#ifdef FIX1
          CVector tmpvec = random_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(fixed_predictors[ index ][j])) * gamma->lastDraw().getVector() ) ) );
          xTy[ index ]->add( tmpvec );
#else
          xTy[ index ]->add( random_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(fixed_predictors[ index ][j])) * gamma->lastDraw().getVector() ) ) ) );
#endif
//          xTy[ index ]->add( random_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(fixed_predictors[ index ][j])) * gamma->lastDraw().getVector() ) ) ) );
        }
      }//end for j
    }
    else if ( invWishart_sigma2 )
    {
      xTy[ index ]->setToZero();
      for ( j = 0; j < response[ index ]->Col(); j++ )
      {
        if ( common_sigma )
        {
#ifdef FIX1
          CVector tmpvec = random_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) );
          xTy[ index ]->add( tmpvec );
#else
          xTy[ index ]->add( random_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) ) );
#endif
//          xTy[ index ]->add( random_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) ) );
        }
        else
        {
#ifdef FIX1
          CVector tmpvec = random_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) );
          xTy[ index ]->add( tmpvec );
#else
          xTy[ index ]->add( random_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) ) );
#endif
//          xTy[ index ]->add( random_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) ) );
        }
      }//end for j
    }
  }//end random effects
  else if ( !strcmp( type, "fixed_effects" ) )
  {
    if ( invWishart_sigma2 )
    {
#ifdef FIX1
      CMatrix tmpmat = ( -1.0 * ( *(mTm[ index ]) ) );
      sMTm->add( tmpmat );
#else
      sMTm->add( ( -1.0 * ( *(mTm[ index ]) ) ) );
#endif
//      sMTm->add( ( -1.0 * ( *(mTm[ index ]) ) ) );

      mTm[ index ]->setToZero();
      for ( j = 0; j < response[ index ]->Col(); j++ )
      {
        if ( common_sigma )
        {
#ifdef FIX1
          CMatrix tmpmat1 = sigma2InvWishart[0]->inverseLastItemDrawn();
          CMatrix tmpmat2 = fixed_predictors[ index ][j]->xTransposedX( tmpmat1 );
          mTm[ index ]->add( tmpmat2 );
#else
          mTm[ index ]->add( fixed_predictors[ index ][j]->xTransposedX( sigma2InvWishart[0]->inverseLastItemDrawn() ) );
#endif
//          mTm[ index ]->add( fixed_predictors[ index ][j]->xTransposedX( sigma2InvWishart[0]->inverseLastItemDrawn() ) );
        }
        else
        {
#ifdef FIX1
          CMatrix tmpmat1 = sigma2InvWishart[ index ]->inverseLastItemDrawn();
          CMatrix tmpmat2 = fixed_predictors[ index ][j]->xTransposedX( tmpmat1 );
          mTm[ index ]->add( tmpmat2 );
#else
          mTm[ index ]->add( fixed_predictors[ index ][j]->xTransposedX( sigma2InvWishart[ index ]->inverseLastItemDrawn() ) );
#endif
//          mTm[ index ]->add( fixed_predictors[ index ][j]->xTransposedX( sigma2InvWishart[ index ]->inverseLastItemDrawn() ) );
        }
      }//end for j

      sMTm->add( (*(mTm[ index ])) );
    }//end if


    if ( random_effects )
    {
      //update sum matrix: substract current value 
#ifdef FIX1
      CVector tmpvec = (*(mTy[ index ])) * ( -1.0 );
      sMTy->add( tmpvec );
#else
      sMTy->add( (*(mTy[ index ])) * ( -1.0 ) );
#endif
//      sMTy->add( (*(mTy[ index ])) * ( -1.0 ) );

      mTy[ index ]->setToZero();
      double precision = 1.0;
      for ( j = 0; j < response[ index ]->Col(); j++ )
      {
        if ( common_sigma || invWishart_sigma2 )
        {
          precision = 1.0;
        }
        else
	    {
          precision = 1.0 / sigma2->lastDrawFromComponent( index ).getScalar();
        }

        if ( !invWishart_sigma2 )
	    {
#ifdef FIX1
          CVector tmpvec = ( fixed_predictors[ index ][j]->T() * ( response[ index ]->getColumn( j ) - ( (*(random_predictors[ index ][j])) * beta[ index ]->lastDraw().getVector() ) ) ) * precision;
          mTy[ index ]->add( tmpvec );
#else
          mTy[ index ]->add( ( fixed_predictors[ index ][j]->T() * ( response[ index ]->getColumn( j ) - ( (*(random_predictors[ index ][j])) * beta[ index ]->lastDraw().getVector() ) ) ) * precision );
#endif
//          mTy[ index ]->add( ( fixed_predictors[ index ][j]->T() * ( response[ index ]->getColumn( j ) - ( (*(random_predictors[ index ][j])) * beta[ index ]->lastDraw().getVector() ) ) ) * precision );
    	}
        else if ( common_sigma )
	    {
#ifdef FIX1
          CVector tmpvec = ( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(random_predictors[ index ][j])) * beta[ index ]->lastDraw().getVector() ) ) ) ) * precision;
           mTy[ index ]->add( tmpvec );
#else
           mTy[ index ]->add( ( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(random_predictors[ index ][j])) * beta[ index ]->lastDraw().getVector() ) ) ) ) * precision );
#endif
//          mTy[ index ]->add( ( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(random_predictors[ index ][j])) * beta[ index ]->lastDraw().getVector() ) ) ) ) * precision );
        }
        else
	    {
#ifdef FIX1
          CVector tmpvec = ( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(random_predictors[ index ][j])) * beta[ index ]->lastDraw().getVector() ) ) ) ) * precision;
          mTy[ index ]->add( tmpvec );
#else
          mTy[ index ]->add( ( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(random_predictors[ index ][j])) * beta[ index ]->lastDraw().getVector() ) ) ) ) * precision );
#endif
//          mTy[ index ]->add( ( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * ( response[ index ]->getColumn( j ) - ( (*(random_predictors[ index ][j])) * beta[ index ]->lastDraw().getVector() ) ) ) ) * precision );
        }
      }//end for j

      //update sum matrix: add updated value
      sMTy->add( (*(mTy[ index ])) );

    }//end if random effects
    else if ( invWishart_sigma2 )
    {
      //update sum matrix: substract current value 
#ifdef FIX1
      CVector tmpvec = (*(mTy[ index ])) * (- 1.0 );
      sMTy->add( tmpvec );
#else
      sMTy->add( (*(mTy[ index ])) * (- 1.0 ) );
#endif
//      sMTy->add( (*(mTy[ index ])) * (- 1.0 ) );

      mTy[ index ]->setToZero();
      for ( j = 0; j < response[ index ]->Col(); j++ )
      {
        if ( common_sigma )
        {  
#ifdef FIX1
          CVector tmpvec = fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) );
          mTy[ index ]->add( tmpvec );
#else
          mTy[ index ]->add( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) ) );
#endif
//          mTy[ index ]->add( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[0]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) ) );
        }
        else
        {
#ifdef FIX1
          CVector tmpvec = fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) );
          mTy[ index ]->add( tmpvec );
#else
          mTy[ index ]->add( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) ) );
#endif
//          mTy[ index ]->add( fixed_predictors[ index ][j]->T() * ( sigma2InvWishart[ index ]->inverseLastItemDrawn() * response[ index ]->getColumn( j ) ) );
        } 
      }//end for j

      //update sum matrix: add updated value
     sMTy->add( ( *(mTy[ index ]) ) );

    }//end if invWishart sigma2
  }//end if fixed effects

}//end




void BayesianMVHLM::gibbsUpdateSecondStageWorkingMatrices()
{
  int index;

  zTvIb->setToZero();
  for ( index = 0; index < number_of_groups; index++ )
  {
    //Note: beta has changed initial covariance according to invWishart and t_beta if appropriate
    (*(zTb[ index ])) = beta_predictors[ index ]->T() * beta[ index ]->initialInverseCovariance() * beta[ index ]->lastDraw().getVector();
     zTvIb->add( (*(zTb[ index ])) );
  }//end for loop

  if ( invWishart_betaCov || t_beta )
  {
    zTvIz->setToZero();
    for ( index = 0; index < number_of_groups; index++ )
    {
#ifdef FIX1
      CMatrix tmpmat = beta[ index ]->initialInverseCovariance();
      (*(zTz[ index ])) = beta_predictors[ index ]->xTransposedX( tmpmat );
#else
      (*(zTz[ index ])) = beta_predictors[ index ]->xTransposedX( beta[ index ]->initialInverseCovariance() );
#endif
//      (*(zTz[ index ])) = beta_predictors[ index ]->xTransposedX( beta[ index ]->initialInverseCovariance() );

      zTvIz->add( (*(zTz[ index ])) );
    }
  }
  
}//end


void BayesianMVHLM::gibbsUpdateResiduals()
{
  int i;

  for ( i = 0; i < number_of_groups; i++ )
  {
    gibbsUpdateResiduals( i );
  }

}//end

void BayesianMVHLM::gibbsUpdateResiduals( int i )
{
  int j;

  if ( random_effects && fixed_effects )
  {
    for ( j = 0; j < response[i]->Col(); j++ )
    {
#ifdef FIX1
      CVector tmpvec = response[i]->getColumn( j ) - ( ( (*(random_predictors[i][j])) * beta[i]->lastDraw().getVector() ) + ( (*(fixed_predictors[i][j])) * gamma->lastDraw().getVector() ) );
      residuals[i]->setColumn( j, tmpvec );
#else
      residuals[i]->setColumn( j, response[i]->getColumn( j ) - ( ( (*(random_predictors[i][j])) * beta[i]->lastDraw().getVector() ) + ( (*(fixed_predictors[i][j])) * gamma->lastDraw().getVector() ) ) );
#endif
//      residuals[i]->setColumn( j, response[i]->getColumn( j ) - ( ( (*(random_predictors[i][j])) * beta[i]->lastDraw().getVector() ) + ( (*(fixed_predictors[i][j])) * gamma->lastDraw().getVector() ) ) );
    }

  }
  else if ( random_effects )
  {
    for ( j = 0; j < response[i]->Col(); j++ )
    {
#ifdef FIX1
      CVector tmpvec = response[i]->getColumn( j )  - ( (*(random_predictors[i][j])) * beta[i]->lastDraw().getVector() );
      residuals[i]->setColumn( j, tmpvec );
#else
      residuals[i]->setColumn( j, response[i]->getColumn( j )  - ( (*(random_predictors[i][j])) * beta[i]->lastDraw().getVector() ) );
#endif
//      residuals[i]->setColumn( j, response[i]->getColumn( j )  - ( (*(random_predictors[i][j])) * beta[i]->lastDraw().getVector() ) );
    }
  }
  else if ( fixed_effects )
  {
    for ( j = 0; j < response[i]->Col(); j++ )
    {
#ifdef FIX1
      CVector tmpvec = response[i]->getColumn( j ) - ( (*(fixed_predictors[i][j])) * gamma->lastDraw().getVector() );
      residuals[i]->setColumn( j, tmpvec );
#else
      residuals[i]->setColumn( j, response[i]->getColumn( j ) - ( (*(fixed_predictors[i][j])) * gamma->lastDraw().getVector() ) );
#endif
//      residuals[i]->setColumn( j, response[i]->getColumn( j ) - ( (*(fixed_predictors[i][j])) * gamma->lastDraw().getVector() ) );
    }
  }
}//end


void BayesianMVHLM::gibbsUpdateSecondStageResiduals()
{
  int i;

  for ( i = 0; i < number_of_groups; i++ )
  {
    gibbsUpdateSecondStageResiduals( i );
  }

}//end


void BayesianMVHLM::gibbsUpdateSecondStageResiduals( int index )
{
  if ( second_stage_effects )
  {
    (*(second_residuals[ index ])) = beta[ index ]->lastItemDrawn() - ( (*(beta_predictors[ index ])) * alpha->lastItemDrawn() );

    //printf("secondstage residuals[ %d ] = \n", index );
    //second_residuals[ index ]->Print();

  }
  else
  {
    (*(second_residuals[ index ])) = beta[ index ]->lastItemDrawn() - beta[ index ]->initialMean();
  }
}//end



void BayesianMVHLM::gibbsUpdateBeta( int index )
{
  CMatrix sigma_beta( beta[ index ]->dimension(), beta[ index ]->dimension() );
  CVector mu_beta( beta[ index ]->dimension() );

  //this update must be done first
  gibbsUpdateWorkingMatrices( index, "random_effects" );

  if ( invWishart_sigma2 )
  {
    mu_beta = (*(xTy[ index ]));
    sigma_beta = (*(xTx[ index ]));
  }
  else if ( common_sigma )
  {
    mu_beta = (*(xTy[ index ])) * ( 1 / sigma2->lastDraw().getScalar() );
    sigma_beta = ( 1 / sigma2->lastDraw().getScalar() ) * ( *(xTx[ index ]) );
  }
  else
  {
    mu_beta = (*(xTy[ index ])) * ( 1 / sigma2->lastDrawFromComponent( index ).getScalar() );
    sigma_beta = ( 1 / sigma2->lastDrawFromComponent( index ).getScalar() ) * ( *(xTx[ index ]) );
  }


  if ( invWishart_betaCov )
  {
#ifdef FIX1
    CMatrix tmpmat = tau2InvWishart->lastItemDrawn();
    beta[ index ]->updateVeryFirstCovariance( tmpmat );
#else
    beta[ index ]->updateVeryFirstCovariance( tau2InvWishart->lastItemDrawn() );
#endif
//    beta[ index ]->updateVeryFirstCovariance( tau2InvWishart->lastItemDrawn() );
  }


  if ( t_beta || invWishart_betaCov )
  {
    CMatrix * initCov = new CMatrix( beta[ index ]->veryFirstCovariance() );

    if ( t_beta )
    {
      initCov->multiplyByScalar( tau2_betas[ index ]->lastItemDrawn() );
    }

    beta[ index ]->updateInitialCovariance( initCov );
    delete initCov;
  }

  if ( second_stage_effects )
  {
    CVector * initMean = new CVector( (*(beta_predictors[ index ])) * alpha->lastItemDrawn() );
    beta[ index ]->updateInitialMean( initMean );

    //printf( "beta: initMean is \n");
    //initMean->Print();
    //printf( "beta: predictors[ %d ], and alpha are = \n", index );
    //beta_predictors[ index ]->Print();
    //alpha->lastItemDrawn().Print();

    delete initMean;
  }

  //printf("Bayes:updateBeta: new beta mean and cov are:\n" );
  //beta[index]->mean().Print();
  //beta[index]->covariance().Print();

  if ( invWishart_betaCov )
  {
    //printf( "beta to update: mu and sigma shifts are:\n");
    //mu_beta.Print();
    //sigma_beta.Print();

    beta[ index ]->update( mu_beta, sigma_beta );

    //printf("Bayes:updateBeta: after updating newest beta mean and cov are:\n" );
    //beta[index]->mean().Print();
    //beta[index]->covariance().Print();    

  }
  else
  {
    double precision_tau2 = 1 / tau2->lastDraw().getScalar();
    beta[ index ]->update( precision_tau2, mu_beta, sigma_beta );
  }

}//end


void BayesianMVHLM::gibbsUpdateGamma()
{
  int i;

  CMatrix sigma_gamma( gamma->dimension(), gamma->dimension() );
  CVector mu_gamma( gamma->dimension() );

  //this update must be done first
  for ( i = 0; i < number_of_groups; i++ )
  {
    gibbsUpdateWorkingMatrices( i, "fixed_effects" );
  }


  if ( invWishart_sigma2 )
  {
    mu_gamma = (*sMTy);
    sigma_gamma = (*sMTm);
  }
  else if ( !common_sigma )
  {
    mu_gamma = (*sMTy);
    sigma_gamma.setToZero();
    for ( i = 0; i < number_of_groups; i++ )
    {
      //done this in update working matrices each time a new sigma is drawn
      //mu_gamma.add( ( *(mTy[ i ]) ) * ( 1 / sigma2->lastDrawFromComponent( i ).getScalar() ) );

#ifdef FIX1
      CMatrix tmpmat = ( 1 / sigma2->lastDrawFromComponent( i ).getScalar() ) * ( *(mTm[ i ]) );
      sigma_gamma.add( tmpmat );
#else
      sigma_gamma.add( ( 1 / sigma2->lastDrawFromComponent( i ).getScalar() ) * ( *(mTm[ i ]) ) );
#endif
//      sigma_gamma.add( ( 1 / sigma2->lastDrawFromComponent( i ).getScalar() ) * ( *(mTm[ i ]) ) );
    }
  }
  else //if ( common_sigma )
  {
    mu_gamma = (*sMTy) * ( 1 / sigma2->lastDraw().getScalar() );
    sigma_gamma = ( 1 / sigma2->lastDraw().getScalar() ) * (*sMTm);
  }

  if ( t_gamma )
  {
    double precision_tau2 = 1.0 / tau2_gamma->lastItemDrawn();
    gamma->update( precision_tau2, mu_gamma, sigma_gamma );
  }
  else
  {

    //printf("updated mu_gamma and sigma_gamma are\n");
    //mu_gamma.Print();
    //sigma_gamma.Print();

    gamma->update( mu_gamma, sigma_gamma );
  }

}//end



void BayesianMVHLM::gibbsUpdateAlpha()
{
  CMatrix sigma_alpha( alpha->dimension(), alpha->dimension() );
  CVector mu_alpha( alpha->dimension() );

  gibbsUpdateSecondStageWorkingMatrices();

  if ( !invWishart_betaCov )
  {
    mu_alpha = (*zTvIb) * ( 1 / tau2->lastDraw().getScalar() );
    sigma_alpha = ( 1 / tau2->lastDraw().getScalar() ) * (*zTvIz);
  }
  else
  {
    mu_alpha = (*zTvIb);
    sigma_alpha = (*zTvIz);
  }

  if ( t_alpha )
  {
    double precision_tau2 = 1 / tau2_alpha->lastItemDrawn();
    alpha->update( precision_tau2, mu_alpha, sigma_alpha );
  }
  else
  {
    alpha->update( mu_alpha, sigma_alpha );
  }

  //printf( "alpha shift mu and cov are:\n" );
  //mu_alpha.Print();
  //sigma_alpha.Print();
  //printf("update alpha: new mean: \n" );
  //alpha->mean().Print();
  //printf("update alpha: new Cov: \n" );
  //alpha->covariance().Print();
  

}//end



void BayesianMVHLM::gibbsUpdateSigma2() throw( rtErr )
{
  int i, j;
  double scale, local_scale, add_df;

  if ( common_sigma )
  {
    if ( invChisq_sigma2 || nonInformativePower_sigma2 || properNonInfoPrior_sigma2 )
    {
      scale = 0;
      for ( i = 0; i < number_of_groups; i++ )
      {    
        for ( j = 0; j < response[i]->Col(); j++ )
	{
          local_scale = residuals[i]->getColumn( j ) * residuals[i]->getColumn( j );
          if ( t_lkhd || t_group_lkhd )
          {
            local_scale *= tau2_error_weights[i]->Val( j );
          }
          scale += local_scale;
        }
      }//end loop

      //printf( "BHL:gibbsUpdateSigma2: RSS = %f\n", scale );

      add_df = number_of_observations * response[0]->Row();

      if ( invChisq_sigma2 || nonInformativePower_sigma2 )
      {
        if ( nonInformativePower_sigma2 )
        {
          add_df += 2 * sigma2NIP[0];
          if ( scale == 0.0 )
	  {
            scale = sigma2ICS[0]->lastDraw().getScalar();
          }
        }

        sigma2ICS[0]->update( add_df, scale );
      }
      else 
      {
        sigma2PNIP[0]->setScale( scale );
      }
    }
    else if ( invWishart_sigma2 )
    {
      CMatrix scale_shift( response[0]->Row(), response[0]->Row() );
      scale_shift.setToZero();
      CMatrix local_shift( 1, response[0]->Row() );
      for ( i = 0; i < number_of_groups; i++ )
      {    
        for ( j = 0; j < response[i]->Col(); j++ )
	{
#ifdef FIX1
          CVector tmpvec = residuals[i]->getColumn( j );
          local_shift.setRow( 0, tmpvec );
#else
          local_shift.setRow( 0, residuals[i]->getColumn( j ) );
#endif
//          local_shift.setRow( 0, residuals[i]->getColumn( j ) );
          if ( t_lkhd || t_group_lkhd )
          {
#ifdef FIX1
            CMatrix tmpmat = tau2_error_weights[i]->Val( j ) * local_shift.xTransposedX();
            scale_shift.add( tmpmat );
#else
            scale_shift.add( tau2_error_weights[i]->Val( j ) * local_shift.xTransposedX() );
#endif
//            scale_shift.add( tau2_error_weights[i]->Val( j ) * local_shift.xTransposedX() );
          }
          else
          {
#ifdef FIX1
            CMatrix tmpmat = local_shift.xTransposedX();
            scale_shift.add( tmpmat );
#else
            scale_shift.add( local_shift.xTransposedX() );
#endif
//            scale_shift.add( local_shift.xTransposedX() );
          }   
        }
      }//end loop

      add_df = number_of_observations;
      sigma2InvWishart[0]->update( add_df, scale_shift );
    }
  }
  else
  {
    printf( "BayesianMVHLM::gibbsUpdateSigma2: Trying to update a common variance variable when there are number of groups variances in the model.\n" );
    char the_error[] = "BayesianMVHLM::gibbsUpdateSigma2: Trying to update a common variance variable when there are number of groups variances in the model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}//end



void BayesianMVHLM::gibbsUpdateSigma2( int index ) throw( rtErr )
{
  int j;
  double scale, add_df;

  if ( !common_sigma )
  {
    if ( invChisq_sigma2 || nonInformativePower_sigma2 || properNonInfoPrior_sigma2 )
    {
      scale = 0;
      for ( j = 0; j < response[ index ]->Col(); j++ )
      {
        if ( t_lkhd || t_group_lkhd )
        {
          scale += residuals[ index ]->getColumn( j ) * ( residuals[ index ]->getColumn( j ).weighted( ( *(tau2_error_weights[ index ]) ) ) );
        }
        else
        {
          scale += residuals[ index ]->getColumn( j ) * residuals[ index ]->getColumn( j );
        }
      }

      add_df = response[ index ]->Col() * response[ index ]->Row();

      if ( invChisq_sigma2 || nonInformativePower_sigma2 )
      {
        if ( nonInformativePower_sigma2 )
        {
          add_df += 2 * sigma2NIP[ index ];
          if ( scale == 0.0 )
	  {
            scale = sigma2ICS[ index ]->lastDraw().getScalar();
          }
        }

        sigma2ICS[ index ]->update( add_df, scale );
      }
      else 
      {
        sigma2PNIP[ index ]->setScale( scale );
      }
    }//end if
    else if ( invWishart_sigma2 )
    {
      CMatrix scale_shift( response[ index ]->Row(), response[ index ]->Row() );
      scale_shift.setToZero();
      CMatrix local_shift( 1, response[ index ]->Row() );
      for ( j = 0; j < response[ index ]->Col(); j++ )
      {
#ifdef FIX1
        CVector tmpvec = residuals[ index ]->getColumn( j );
        local_shift.setRow( 0, tmpvec );
#else
        local_shift.setRow( 0, residuals[ index ]->getColumn( j ) );
#endif
//        local_shift.setRow( 0, residuals[ index ]->getColumn( j ) );
        if ( t_lkhd || t_group_lkhd )
        {
#ifdef FIX1
          CMatrix tmpmat = tau2_error_weights[ index ]->Val( j ) * local_shift.xTransposedX();
          scale_shift.add( tmpmat );
#else
          scale_shift.add( tau2_error_weights[ index ]->Val( j ) * local_shift.xTransposedX()  );
#endif
//          scale_shift.add( tau2_error_weights[ index ]->Val( j ) * local_shift.xTransposedX()  );
        }
        else
        {
#ifdef FIX1
          CMatrix tmpmat = local_shift.xTransposedX();
          scale_shift.add( tmpmat );
#else
          scale_shift.add( local_shift.xTransposedX() );
#endif
//          scale_shift.add( local_shift.xTransposedX() );
        }   
      }//end loop

      add_df = response[ index ]->Col();
      sigma2InvWishart[ index ]->update( add_df, scale_shift );
    }
  }
  else
  {
    printf( "BayesianMVHLM::gibbsUpdateSigma2: Trying to update number of groups variances when there is a common variance in the model.\n" );
    char the_error[] = "BayesianMVHLM::gibbsUpdateSigma2: Trying to update number of groups variances when there is a common variance in the model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

}//end 


void BayesianMVHLM::gibbsUpdateTau2()
{
  int i;
  double scale, local_scale, add_df;

  if ( invChisq_tau2 || nonInformativePower_tau2 || properNonInfoPrior_tau2 )
  {
    scale = 0;
    for ( i = 0; i < number_of_groups; i++ )
    {    
      local_scale = M_D( beta[i]->initialInverseCovariance(), ( *(second_residuals[i]) ) );

      if ( t_beta )
      {
        local_scale /= tau2_betas[i]->lastItemDrawn();
      }

      scale += local_scale;
    }//end loop


    //printf("gibbsUpdateTau2: scale = %f\n", scale );

    add_df = number_of_groups * beta[0]->dimension();

    if ( invChisq_tau2 || nonInformativePower_tau2 )
    {
      if ( nonInformativePower_tau2 )
      {
        add_df += 2 * tau2NIP;
        if ( scale == 0.0 )
	{
          scale = tau2ICS->lastDraw().getScalar();
        }
      }

      tau2ICS->update( add_df, scale );
    }
    else
    {
      tau2PNIP->setScale( scale );
    }
  }
  else if ( invWishart_betaCov )
  {
    CMatrix scale_matrix( beta[0]->dimension(), beta[0]->dimension() );
    CMatrix resid( 1, beta[0]->dimension() );
    CVector diff( beta[0]->dimension() );

    scale_matrix.setToZero();
    for ( i = 0; i < number_of_groups; i++ )
    {
      diff = (*(second_residuals[i]));
      resid.setRow( 0, diff );

      if ( t_beta )
      {
#ifdef FIX1
        CMatrix tmpmat = ( 1.0 / tau2_betas[i]->lastItemDrawn() ) * resid.xTransposedX();
        scale_matrix.add( tmpmat );
#else
        scale_matrix.add( ( 1.0 / tau2_betas[i]->lastItemDrawn() ) * resid.xTransposedX() );
#endif
//        scale_matrix.add( ( 1.0 / tau2_betas[i]->lastItemDrawn() ) * resid.xTransposedX() );
      }
      else
      {
#ifdef FIX1
        CMatrix tmpmat = resid.xTransposedX();
        scale_matrix.add( tmpmat );
#else
        scale_matrix.add( resid.xTransposedX() );
#endif
//        scale_matrix.add( resid.xTransposedX() );
      }
    }//end loop

    add_df = number_of_groups;

    tau2InvWishart->update( add_df, scale_matrix );
  }

}//end



void BayesianMVHLM::gibbsUpdateTau2Gamma()
{
  double scale, add_df;
  CVector gamma_diff( gamma->dimension() );

  gamma_diff = gamma->lastItemDrawn() - gamma->initialMean();
  scale = M_D( gamma->initialInverseCovariance(), gamma_diff );

  add_df = gamma->dimension();

  tau2_gamma->update( add_df, scale );

}//end



void BayesianMVHLM::gibbsUpdateTau2Beta( int index )
{
  double scale, add_df;

  if ( !invWishart_betaCov )
  {
    scale = M_D( beta[ index ]->initialInverseCovariance(), ( *(second_residuals[ index ]) ) );
    scale /= tau2->lastDraw().getScalar();
  }
  else
  {
    scale = M_D( tau2->lastDraw().getMatrix().inverse(), (*(second_residuals[ index ])) );
  }


  add_df = beta[ index ]->dimension();

  tau2_betas[ index ]->update( add_df, scale );

}//end


void BayesianMVHLM::gibbsUpdateTau2Alpha()
{
  double scale, add_df;
  CVector alpha_diff( alpha->dimension() );

  alpha_diff = alpha->lastItemDrawn() - alpha->initialMean();
  scale = M_D( alpha->initialInverseCovariance(), alpha_diff );

  add_df = alpha->dimension();

  tau2_alpha->update( add_df, scale );

}//end


void BayesianMVHLM::gibbsUpdateTau2Group( int index )
{
  int j;
  double scale, add_df;

  scale = 0;
  for ( j = 0; j < response[ index ]->Col(); j++ )
  {
    if ( !invWishart_sigma2 )
    {
      scale += residuals[ index ]->getColumn( j ) * residuals[ index ]->getColumn( j );
    }
    else if ( common_sigma )
    {
      scale += M_D( sigma2InvWishart[0]->inverseLastItemDrawn(), residuals[ index ]->getColumn( j ) );
    }
    else
    {
      scale += M_D( sigma2InvWishart[ index ]->inverseLastItemDrawn(), residuals[ index ]->getColumn( j ) );
    }
  }

  if ( !invWishart_sigma2 )
  {
    if ( common_sigma )
    {
      scale /= sigma2->lastDraw().getScalar();
    }
    else
    {
      scale /= sigma2->lastDrawFromComponent( index ).getScalar();
    }
  }
 
  add_df = response[ index ]->Col() * response[ index ]->Row();

  tau2_groups[ index ]->update( add_df, scale );

}//end


void BayesianMVHLM::gibbsUpdateTau2Error( int index )
{
  int group, obs;
  double scale, add_df;

  DistributionParameter grp( groupIndex( index ) );

  group = (int) grp.getVector().Val(0);
  obs = (int) grp.getVector().Val(1);

  //printf(" group = %d,  obs = %d.  n groups = %d\n", group, obs, number_of_groups ); fflush(stdout);

  if ( !invWishart_sigma2 )
  {
    scale = residuals[ group ]->getColumn( obs ) * residuals[ group ]->getColumn( obs );
    if ( common_sigma )
    {
      scale /= sigma2->lastDraw().getScalar();
    }
    else
    {
      scale /= sigma2->lastDrawFromComponent( group ).getScalar();
    }
  }
  else if ( common_sigma )
  {
    scale = M_D( sigma2InvWishart[0]->inverseLastItemDrawn(), residuals[ group ]->getColumn( obs ) );
  }
  else
  {
    scale = M_D( sigma2InvWishart[ group ]->inverseLastItemDrawn(), residuals[ group ]->getColumn( obs ) );
  }

  add_df = response[ group ]->Row();

  tau2_errors[ index ]->update( add_df, scale );

}//end


void BayesianMVHLM::drawMissingResponse( int index )
{
  int group, obs, nmss, nobs;

  DistributionParameter grp( groupIndex( index ) );

  group = (int) grp.getVector().Val(0);
  obs = (int) grp.getVector().Val(1);

  nmss = (int) number_missing_response->Val( index );
  nobs = response[ group ]->Row() - nmss;

  //printf("drawing missing response [ %d ][ %d ]\n", group, obs ); fflush(stdout);
  //printf("missing components are:\n");
  //missing_response_components[ group ][ obs ]->Print();

  CVector mu_mss( nmss );
  CVector mu_obs( nobs );
  CVector mu_imputed( nmss );
  CMatrix Sigma_imputed( nmss, nmss );

  if ( nobs > 0 )
  {
    mu_mss = (response[ group ]->getColumn( obs ) - residuals[ group ]->getColumn( obs ) ).subVector( (*(missing_response_components[ group ][ obs ])) );
    mu_obs = residuals[ group ]->getColumn( obs ).subVectorComplement( (*(missing_response_components[ group ][ obs ])) );

    //printf(" beta[ %d ] is \n", group );
    //beta[group]->lastItemDrawn().Print();

    //printf("mu_mss is \n");
    //mu_mss.Print();

    //printf("residual is\n");
    //residuals[ group ]->getColumn( obs ).Print();

    //printf("response is \n");
    //response[group]->getColumn( obs ).Print();

    //printf("mu_obs is \n");
    //mu_obs.Print();

  }
  else
  {
    mu_mss = response[ group ]->getColumn( obs ) - residuals[ group ]->getColumn( obs );
  }

  if ( invWishart_sigma2 )
  {
    if ( nobs > 0 )
    {
      //form submatrix
      CMatrix Sigma_mss( nmss, nmss );
      CMatrix Sigma_obs( nobs, nobs );
      CMatrix Sigma_cross( nmss, nobs );

      if ( !common_sigma )
      {
        //printf( "Sigma[ %d ] = \n", group );
        //sigma2InvWishart[ group ]->lastItemDrawn().Print();

        

        Sigma_mss = sigma2InvWishart[ group ]->lastItemDrawn().subMatrix( (*(missing_response_components[ group ][ obs ])) );

        //printf( "SubSigma matrix is\n" );
        //Sigma_mss.Print();


        Sigma_obs = sigma2InvWishart[ group ]->lastItemDrawn().subMatrixComplement( (*(missing_response_components[ group ][ obs ])) );

        //printf( "SubSigma complement matrix is\n" );
        //Sigma_obs.Print();

        Sigma_cross = sigma2InvWishart[ group ]->lastItemDrawn().subMatrixCross( (*(missing_response_components[ group ][ obs ])) );

        //printf( "SubSigma crossed matrix is\n" );
        //Sigma_cross.Print();
      }
      else
      {
        Sigma_mss = sigma2InvWishart[0]->lastItemDrawn().subMatrix( (*(missing_response_components[ group ][ obs ])) );
        Sigma_obs = sigma2InvWishart[0]->lastItemDrawn().subMatrixComplement( (*(missing_response_components[ group ][ obs ])) );
        Sigma_cross = sigma2InvWishart[0]->lastItemDrawn().subMatrixCross( (*(missing_response_components[ group ][ obs ])) );
      }

      CMatrix inv_Sigma_obs( Sigma_obs.inverse() );
      mu_imputed = mu_mss + ( Sigma_cross * ( inv_Sigma_obs * mu_obs ) );
      Sigma_imputed = Sigma_mss - Sigma_cross.T().xTransposedX( inv_Sigma_obs );

      //printf("Sigma complement inv is \n" );
      //inv_Sigma_obs.Print();

      //printf("mu imputed is \n");
      //mu_imputed.Print();

      //printf("Sigma_imputed is\n");
      //Sigma_imputed.Print();

    }
    else
    {
      Sigma_imputed = sigma2InvWishart[0]->lastItemDrawn();
      mu_imputed = mu_mss;
    }
  }//end if invWishart
  else if ( common_sigma )
  {
    Sigma_imputed.setToZero();
    Sigma_imputed.setDiagonal( sigma2->lastDraw().getScalar() );
    mu_imputed = mu_mss;
  }
  else
  {
    Sigma_imputed.setToZero();
    Sigma_imputed.setDiagonal( sigma2->lastDrawFromComponent( group ).getScalar() );
    mu_imputed = mu_mss;
  }

  NormalDistribution y_imputed( nmss );
  y_imputed.setMean( &mu_imputed );
  y_imputed.setCovariance( &Sigma_imputed );
  y_imputed.draw();

  //printf("current response[ %d, %d ] is \n", group, obs );
  //response[ group ]->getColumn( obs ).Print();

#ifdef FIX1
  CVector tmpvec1 = (*(missing_response_components[ group ][ obs ])), tmpvec2 = y_imputed.lastItemDrawn();
  response[ group ]->setSubColumn( obs, tmpvec1, tmpvec2 );
#else
  response[ group ]->setSubColumn( obs, (*(missing_response_components[ group ][ obs ])), y_imputed.lastItemDrawn() );
#endif
//  response[ group ]->setSubColumn( obs, (*(missing_response_components[ group ][ obs ])), y_imputed.lastItemDrawn() );  


  //printf("y imputes is \n");
  //y_imputed.lastItemDrawn().Print();
  
  //printf("new response[ %d, %d ] is \n", group, obs );
  //response[ group ]->getColumn( obs ).Print();

  gibbsUpdateResiduals( group );

  //gibbsUpdateWorkingMatrices( group, "random_effects" );
  //gibbsUpdateWorkingMatrices( group, "fixed_effects" );

}//end



void BayesianMVHLM::drawMissingRandomPredictors( int index )
{
  int i, k, start_index, iVar, group, obs, nmss, p_dim;
  double mu1, mu2;

  DistributionParameter grp( groupIndex( index ) );

  group = (int) grp.getVector().Val(0);
  obs = (int) grp.getVector().Val(1);

  nmss = (int) number_missing_random_predictors->Val( index );
  //recall x predictors are original x repeated r times
  p_dim = (int) ( ((double) random_predictors[ group ][ obs ]->Col()) / ((double) response[ group ]->Row() ) );

  //printf("drawing missing random predictors [ %d ][ %d ]\n", group, obs ); fflush(stdout);
  //printf("missing components are:\n");
  //missing_random_predictors_components[ group ][ obs ]->Print();

  //printf(" p_dim = %d, nmss = %d\n", p_dim, nmss ); fflush(stdout);

  CVector mu_imputed( 1 );
  CMatrix Sigma_imputed( 1, 1 );

  //compute mean and variance of x (assuming normality)
  CMatrix beta_coef( p_dim, response[ group ]->Row() );
  CVector Sigma_x( p_dim );
  CMatrix Sigma_x_inv( p_dim, p_dim );
  CVector temp_x( p_dim );
  CVector temp_y( response[ group ]->Row() );
  CVector x_pred( p_dim );

  //get the covariates
  x_pred = random_predictors[ group ][ obs ]->getRow( 0 ).subVector( 0, p_dim - 1 );

  for ( i = 0; i < response[ group ]->Row(); i++ )
  {
#ifdef FIX1
    CVector tmpvec = beta[ group ]->lastItemDrawn().subVector( i * p_dim, (i+1) * p_dim - 1 );
    beta_coef.setColumn( i, tmpvec );
#else
    beta_coef.setColumn( i, beta[ group ]->lastItemDrawn().subVector( i * p_dim, (i+1) * p_dim - 1 ) );
#endif
//    beta_coef.setColumn( i, beta[ group ]->lastItemDrawn().subVector( i * p_dim, (i+1) * p_dim - 1 ) );
  }


  //printf( " x _pred is \n" );
  //x_pred.Print();
  //printf( "beta_coef is \n");
  //beta_coef.Print();
  //fflush(stdout);



  if ( !beta_coef.isZero( DELTA ) )
  {
    if ( fixed_effects )
    {
      temp_y = response[ group ]->getColumn( obs ) - ( (*(fixed_predictors[ group ][ obs ])) * gamma->lastItemDrawn() );
    }
    else
    {
      temp_y = response[ group ]->getColumn( obs );
    }


    if ( !invWishart_sigma2 )
    {
      if ( common_sigma )
      {
        Sigma_x_inv = ( 1.0 / sigma2->lastDraw().getScalar() ) * beta_coef.T().xTransposedX();
        temp_x = ( beta_coef * temp_y ) * ( 1.0 / sigma2->lastDraw().getScalar() ); 
      }
      else
      {
        Sigma_x_inv = ( 1.0 / sigma2->lastDrawFromComponent( group ).getScalar() ) * beta_coef.T().xTransposedX();
        temp_x = ( beta_coef * temp_y ) * ( 1.0 / sigma2->lastDrawFromComponent( group ).getScalar() ); 
      }
    }
    else if ( common_sigma )
    {
#ifdef FIX1
      CMatrix tmpmat = sigma2InvWishart[ 0 ]->inverseLastItemDrawn();
      Sigma_x_inv = beta_coef.T().xTransposedX( tmpmat );
#else
      Sigma_x_inv = beta_coef.T().xTransposedX( sigma2InvWishart[ 0 ]->inverseLastItemDrawn() );
#endif
//      Sigma_x_inv = beta_coef.T().xTransposedX( sigma2InvWishart[ 0 ]->inverseLastItemDrawn() );
      temp_x = beta_coef * ( sigma2InvWishart[ 0 ]->inverseLastItemDrawn() * temp_y );
    }
    else
    {
#ifdef FIX1
      CMatrix tmpmat = sigma2InvWishart[ group ]->inverseLastItemDrawn();
      Sigma_x_inv = beta_coef.T().xTransposedX( tmpmat );
#else
      Sigma_x_inv = beta_coef.T().xTransposedX( sigma2InvWishart[ group ]->inverseLastItemDrawn() );
#endif
//      Sigma_x_inv = beta_coef.T().xTransposedX( sigma2InvWishart[ group ]->inverseLastItemDrawn() );
      temp_x = beta_coef * ( sigma2InvWishart[ group ]->inverseLastItemDrawn() * temp_y );
    }


    //now compute sigma2_k and mu_k for each missing component
    Sigma_x.setToZero();
    for ( i = 0; i < p_dim; i++ )
    {
      if ( Sigma_x_inv.Val( i, i ) > 0.0 )
      {
        Sigma_x.Val( i ) = 1.0 / Sigma_x_inv.Val( i, i );
      }
    }

    //ready for imputation
    for ( i = 0; i < nmss; i++ )
    {
      //choose i at random
      start_index = randomStart( nmss );
      iVar = ( start_index + i ) %  nmss;
      k = (int) missing_random_predictors_components[ group ][ obs ]->Val( iVar );

      if ( Sigma_x.Val( k ) > 0 )
      {
        mu1 = ( temp_x.Val( k ) - ( Sigma_x_inv.getRow( k ) * x_pred ) + Sigma_x_inv.Val( k, k ) * x_pred.Val( k ) );

        mu2 = randomEffects_means[ group ]->Val( k ) / randomEffects_vars[ group ]->Val( k );

        Sigma_imputed.Val( 0, 0 ) =  1 / ( ( 1 / randomEffects_vars[ group ]->Val( k ) ) + ( 1 / Sigma_x.Val( k ) ) ); 

        mu_imputed.Val( 0 ) = Sigma_imputed.Val( 0, 0 ) * ( mu1 + mu2 );

        //mu_imputed.Val( 0 ) = Sigma_x.Val( k ) * ( temp_x.Val( k ) - ( Sigma_x_inv.getRow( k ) * x_pred ) + Sigma_x_inv.Val( k, k ) * x_pred.Val( k ) );
        //Sigma_imputed.Val( 0, 0 ) =  Sigma_x.Val( k );


        //ready for draw

        //printf(" mu imputed and Sigma_imputed are\n" );
        //mu_imputed.Print();
        //Sigma_imputed.Print();
        //fflush(stdout);

        NormalDistribution x_imputed( 1 );
        x_imputed.setMean( &mu_imputed );
        x_imputed.setCovariance( &Sigma_imputed );
        x_imputed.draw();

        //update missing values
        x_pred.Val( k ) = x_imputed.lastItemDrawn().Val( 0 );

        //printf( "x pred (%d) is \n", k );
        //x_pred.Print();
      }
    }//end for i

    //now set the predictors
    for ( i = 0; i < response[ group ]->Row(); i++ )
    {
      random_predictors[ group ][ obs ]->setSubRow( i, i * p_dim, (i+1) * p_dim, x_pred );
    }  

    //printf("new random_predictors[ %d ][ %d ] is \n", group, obs );
    //random_predictors[ group ][ obs ]->Print();

    gibbsUpdateResiduals( group );

    //gibbsUpdateWorkingMatrices( group, "random_effects" );
    //gibbsUpdateWorkingMatrices( group, "fixed_effects" );
  }//end if beta is not zero

}//end



void BayesianMVHLM::drawMissingFixedPredictors( int index )
{
  int i, k, start_index, iVar, group, obs, nmss, p_dim;
  double mu1, mu2;

  DistributionParameter grp( groupIndex( index ) );

  group = (int) grp.getVector().Val(0);
  obs = (int) grp.getVector().Val(1);

  nmss = (int) number_missing_fixed_predictors->Val( index );
  //recall x predictors are original x repeated r times
  p_dim = (int) ( ((double) fixed_predictors[ group ][ obs ]->Col()) / ((double) response[ group ]->Row() ) );

  //printf("drawing missing fixed predictors [ %d ][ %d ]\n", group, obs ); fflush(stdout);
  //printf("missing components are:\n");
  //missing_fixed_predictors_components[ group ][ obs ]->Print();

  //printf(" p_dim = %d, nmss = %d\n", p_dim, nmss ); fflush(stdout);

  CVector mu_imputed( 1 );
  CMatrix Sigma_imputed( 1, 1 );

  //compute mean and variance of x (assuming normality)
  CMatrix gamma_coef( p_dim, response[ group ]->Row() );
  CVector Sigma_x( p_dim );
  CMatrix Sigma_x_inv( p_dim, p_dim );
  CVector temp_x( p_dim );
  CVector temp_y( response[ group ]->Row() );
  CVector x_pred( p_dim );

  //get the covariates
  x_pred = fixed_predictors[ group ][ obs ]->getRow( 0 ).subVector( 0, p_dim - 1 );

  for ( i = 0; i < response[ group ]->Row(); i++ )
  {
#ifdef FIX1
    CVector tmpvec = gamma->lastItemDrawn().subVector( i * p_dim, (i+1) * p_dim - 1 );
    gamma_coef.setColumn( i, tmpvec );
#else
    gamma_coef.setColumn( i, gamma->lastItemDrawn().subVector( i * p_dim, (i+1) * p_dim - 1 ) );
#endif
//    gamma_coef.setColumn( i, gamma->lastItemDrawn().subVector( i * p_dim, (i+1) * p_dim - 1 ) );
  }


  //printf( " x _pred is \n" );
  //x_pred.Print();
  //printf( "gamma_coef is \n");
  //gamma_coef.Print();
  //fflush(stdout);



  if ( !gamma_coef.isZero( DELTA ) )
  {
    if ( random_effects )
    {
      temp_y = response[ group ]->getColumn( obs ) - ( (*(random_predictors[ group ][ obs ])) * beta[ group ]->lastItemDrawn() );
    }
    else
    {
      temp_y = response[ group ]->getColumn( obs );
    }


    if ( !invWishart_sigma2 )
    {
      if ( common_sigma )
      {
        Sigma_x_inv = ( 1.0 / sigma2->lastDraw().getScalar() ) * gamma_coef.T().xTransposedX();
        temp_x = ( gamma_coef * temp_y ) * ( 1.0 / sigma2->lastDraw().getScalar() ); 
      }
      else
      {
        Sigma_x_inv = ( 1.0 / sigma2->lastDrawFromComponent( group ).getScalar() ) * gamma_coef.T().xTransposedX();
        temp_x = ( gamma_coef * temp_y ) * ( 1.0 / sigma2->lastDrawFromComponent( group ).getScalar() ); 
      }
    }
    else if ( common_sigma )
    {
#ifdef FIX1
      CMatrix tmpmat = sigma2InvWishart[ 0 ]->inverseLastItemDrawn();
      Sigma_x_inv = gamma_coef.T().xTransposedX( tmpmat );
#else
      Sigma_x_inv = gamma_coef.T().xTransposedX( sigma2InvWishart[ 0 ]->inverseLastItemDrawn() );
#endif
//      Sigma_x_inv = gamma_coef.T().xTransposedX( sigma2InvWishart[ 0 ]->inverseLastItemDrawn() );
      temp_x = gamma_coef * ( sigma2InvWishart[ 0 ]->inverseLastItemDrawn() * temp_y );
    }
    else
    {
#ifdef FIX1
      CMatrix tmpmat = sigma2InvWishart[ group ]->inverseLastItemDrawn();
      Sigma_x_inv = gamma_coef.T().xTransposedX( tmpmat );
#else
      Sigma_x_inv = gamma_coef.T().xTransposedX( sigma2InvWishart[ group ]->inverseLastItemDrawn() );
#endif
//      Sigma_x_inv = gamma_coef.T().xTransposedX( sigma2InvWishart[ group ]->inverseLastItemDrawn() );
      temp_x = gamma_coef * ( sigma2InvWishart[ group ]->inverseLastItemDrawn() * temp_y );
    }


    //now compute sigma2_k and mu_k for each missing component
    Sigma_x.setToZero();
    for ( i = 0; i < p_dim; i++ )
    {
      if ( Sigma_x_inv.Val( i, i ) > 0.0 )
      {
        Sigma_x.Val( i ) = 1.0 / Sigma_x_inv.Val( i, i );
      }
    }

    //ready for imputation
    for ( i = 0; i < nmss; i++ )
    {
      //choose i at random
      start_index = randomStart( nmss );
      iVar = ( start_index + i ) %  nmss;
      k = (int) missing_fixed_predictors_components[ group ][ obs ]->Val( iVar );

      if ( Sigma_x.Val( k ) > 0 )
      {
        mu1 = ( temp_x.Val( k ) - ( Sigma_x_inv.getRow( k ) * x_pred ) + Sigma_x_inv.Val( k, k ) * x_pred.Val( k ) );

        mu2 = fixedEffects_means[ group ]->Val( k ) / fixedEffects_vars[ group ]->Val( k );

        Sigma_imputed.Val( 0, 0 ) =  1 / ( ( 1 / fixedEffects_vars[ group ]->Val( k ) ) + ( 1 / Sigma_x.Val( k ) ) ); 

        mu_imputed.Val( 0 ) = Sigma_imputed.Val( 0, 0 ) * ( mu1 + mu2 );

        //mu_imputed.Val( 0 ) = Sigma_x.Val( k ) * ( temp_x.Val( k ) - ( Sigma_x_inv.getRow( k ) * x_pred ) + Sigma_x_inv.Val( k, k ) * x_pred.Val( k ) );
        //Sigma_imputed.Val( 0, 0 ) =  Sigma_x.Val( k );
        //ready for draw

        //printf(" mu imputed and Sigma_imputed are\n" );
        //mu_imputed.Print();
        //Sigma_imputed.Print();
        //fflush(stdout);

        NormalDistribution x_imputed( 1 );
        x_imputed.setMean( &mu_imputed );
        x_imputed.setCovariance( &Sigma_imputed );
        x_imputed.draw();

        //update missing values
        x_pred.Val( k ) = x_imputed.lastItemDrawn().Val( 0 );

        //printf( "x pred (%d) is \n", k );
        //x_pred.Print();
      }
    }//end for i

    //now set the predictors
    for ( i = 0; i < response[ group ]->Row(); i++ )
    {
      fixed_predictors[ group ][ obs ]->setSubRow( i, i * p_dim, (i+1) * p_dim, x_pred );
    }  

    //printf("new fixed_predictors[ %d ][ %d ] is \n", group, obs );
    //fixed_predictors[ group ][ obs ]->Print();

    gibbsUpdateResiduals( group );

    //gibbsUpdateWorkingMatrices( group, "random_effects" );
    //gibbsUpdateWorkingMatrices( group, "fixed_effects" );
  }//end if gamma is not zero

}//end



void BayesianMVHLM::fullConditionalUpdateVariable( int index )
{
  //printf( "BHL:full conditional update: index is %d.\n  ", index ); fflush( stdout );

  if ( index < number_of_variables && index >= 0 )
  {
    //printf("BHL:full conditional update: var is [ %s ]\n", distr_map[ index ] ); fflush(stdout);

    if ( !strcmp( distr_map[ index ], "gamma" ) )
    {
      gibbsUpdateGamma();
    }
    else if ( !strcmp( distr_map[ index ], "alpha" ) )
    {
      gibbsUpdateAlpha();
    }
    else if ( !strcmp( distr_map[ index ], "sigma" ) )
    {
      gibbsUpdateSigma2();
    }
    else if ( !strcmp( distr_map[ index ], "tau2" ) )
    {
      gibbsUpdateTau2();
    }
    else if ( !strcmp( distr_map[ index ], "tau2a" ) )
    {
      gibbsUpdateTau2Alpha();
    }
    else if ( !strcmp( distr_map[ index ], "tau2gm" ) )
    {
      gibbsUpdateTau2Gamma();
    }
    else if ( !strncmp( distr_map[ index ], "beta", 4 ) 
              || !strncmp( distr_map[ index ], "sigma", 5 )
              || !strncmp( distr_map[ index ], "tau2b", 5 )
              || !strncmp( distr_map[ index ], "tau2g", 5 ) 
              || !strncmp( distr_map[ index ], "tau2e", 5 ) )
    {
      int group_index;
      char * ptr, * temp_distr;

      temp_distr = new char [ strlen( distr_map[ index ] ) + 1 ];
      sprintf( temp_distr, "%s", distr_map[ index ] );
      ptr = strtok( temp_distr, ":" );
      if ( ptr != NULL )
      {
        ptr = strtok( NULL, ":" );
        if ( ptr != NULL )
	{
          group_index = atoi( ptr );

          //printf( "BHL:full conditional: group index is %d\n", group_index );

          if ( !strncmp( distr_map[ index ], "tau2e", 5 ) )
          {
            //printf( "group index = %d, n obs is %d\n", group_index, number_of_observations ); fflush(stdout);

            if ( group_index >= 0 && group_index < number_of_observations )
	    {
              gibbsUpdateTau2Error( group_index );
            }
            else
	    {
              printf( " BayesianMVHLM::fullConditionalUpdateVariable: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
            }
          }
          else
	  {
            if ( group_index >= 0 && group_index < number_of_groups )
	    {
              if ( !strncmp( distr_map[ index ], "beta", 4 ) )
	      {
                gibbsUpdateBeta( group_index );
              }
              else if ( !strncmp( distr_map[ index ], "sigma", 5 ) )
	      {
                gibbsUpdateSigma2( group_index );
              }
              else if ( !strncmp( distr_map[ index ], "tau2b", 5 ) )
	      {
                gibbsUpdateTau2Beta( group_index );
              }          
              else if ( !strncmp( distr_map[ index ], "tau2g", 5 ) )
	      {
                gibbsUpdateTau2Group( group_index );
              }          
            }
            else
            {
              printf( " BayesianMVHLM::fullConditionalUpdateVariable: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
            }
          }
        }//end if not null
        else
        {
          printf( " BayesianMVHLM::fullConditionalUpdateVariable: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }
      else
      {
        printf( " BayesianMVHLM::fullConditionalUpdateVariable: Wrong argument in [%s]. Number expected.\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if 
  }//end if index is valid
  else
  {
    printf( "BayesianMVHLM::fullConditionalUpdateVariable: Variable index [%d] does not exist.\n", index );
  }

}//end



void BayesianMVHLM::drawVariable( int index )
{
  //printf( "BHL:draw variable: index is %d.\n  ", index ); fflush( stdout );

  if ( index < number_of_variables && index >= 0 )
  {
    if ( !strcmp( distr_map[ index ], "gamma" ) )
    {
      gamma->draw();

      //printf( "BHL: gamma draw = \n" );
      //gamma->lastItemDrawn().Print();

      gibbsUpdateResiduals();
    }
    else if ( !strcmp( distr_map[ index ], "alpha" ) )
    {
      alpha->draw();

      //printf("alpha draw is \n"); alpha->lastItemDrawn().Print();

      gibbsUpdateSecondStageResiduals();
    }
    else if ( !strcmp( distr_map[ index ], "sigma" ) )
    {
      sigma2->draw();

      //printf( "BHL: sigma draw = %f\n", sqrt( sigma2->lastDraw().getScalar() ) );
    }
    else if ( !strcmp( distr_map[ index ], "tau2" ) )
    {
      tau2->draw();
    }
    else if ( !strcmp( distr_map[ index ], "tau2a" ) )
    {
      tau2_alpha->draw();
    }
    else if ( !strcmp( distr_map[ index ], "tau2gm" ) )
    {
      tau2_gamma->draw();
    }
    else if ( !strncmp( distr_map[ index ], "beta", 4 ) 
              || !strncmp( distr_map[ index ], "sigma", 5 )
              || !strncmp( distr_map[ index ], "missing", 7 )
              || !strncmp( distr_map[ index ], "tau2b", 5 )
              || !strncmp( distr_map[ index ], "tau2g", 5 ) 
              || !strncmp( distr_map[ index ], "tau2e", 5 ) )
    {
      int group_index;
      char * ptr, * temp_distr;

      temp_distr = new char [ strlen( distr_map[ index ] ) + 1 ];
      sprintf( temp_distr, "%s", distr_map[ index ] );
      ptr = strtok( temp_distr, ":" );
      if ( ptr != NULL )
      {
        ptr = strtok( NULL, ":" );
        if ( ptr != NULL )
	{
          group_index = atoi( ptr );
          if ( !strncmp( distr_map[ index ], "tau2e", 5 ) )
          {
            if ( group_index >= 0 && group_index < number_of_observations )
	    {
              tau2_errors[ group_index ]->draw();
              updateRegressionWeight( group_index );
            }
            else
	    {
              printf( " BayesianMVHLM::drawVariable: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
            }
          }
          else if ( !strncmp( distr_map[ index ], "missingR:", 9 ) )
	  {
            if ( group_index >= 0 && group_index < number_of_observations )
	    {
              drawMissingResponse( group_index );
              //update of residuals takes place within drawMissingResponse
            }
            else
	    {
              printf( " BayesianMVHLM::drawVariable: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
            }
          }
          else if ( !strncmp( distr_map[ index ], "missingRP", 9 ) )
	  {
            if ( group_index >= 0 && group_index < number_of_observations )
	    {
              drawMissingRandomPredictors( group_index );
              //update of residuals takes place within drawMissingRandomPredictors
            }
            else
	    {
              printf( " BayesianMVHLM::drawVariable: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
            }
          }
          else if ( !strncmp( distr_map[ index ], "missingFP", 9 ) )
	  {
            if ( group_index >= 0 && group_index < number_of_observations )
	    {
              drawMissingFixedPredictors( group_index );
              //update of residuals takes place within drawMissingFixedPredictors
            }
            else
	    {
              printf( " BayesianMVHLM::drawVariable: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
            }
          }
          else
	  {
            if ( group_index >= 0 && group_index < number_of_groups )
	    {
              if ( !strncmp( distr_map[ index ], "beta", 4 ) )
	      {
                beta[ group_index ]->draw();

                //printf("beta[%d] draw is \n", group_index); 
		//beta[group_index]->lastItemDrawn().Print();

                gibbsUpdateSecondStageResiduals( group_index );
                gibbsUpdateResiduals( group_index );
              }
              else if ( !strncmp( distr_map[ index ], "sigma", 5 ) )
	      {
                sigma2->drawFromComponent( group_index );
              }
              else if ( !strncmp( distr_map[ index ], "tau2b", 5 ) )
	      {
                tau2_betas[ group_index ]->draw();
              }          
              else if ( !strncmp( distr_map[ index ], "tau2g", 5 ) )
	      {
                tau2_groups[ group_index ]->draw();

                //printf( "BHL: tau2 groups[%d] = %f\n", group_index, tau2_groups[ group_index ]->lastItemDrawn() );

                updateGroupRegressionWeight( group_index );
              }          
            }
            else
            {
              printf( " BayesianMVHLM::drawVariable: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
            }
          }
        }//end if not null
        else
        {
          printf( " BayesianMVHLM::drawVariable: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }
      else
      {
        printf( " BayesianMVHLM::drawVariable: Wrong argument in [%s]. Number expected.\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if 
  }//end if index is valid
  else
  {
    printf( "BayesianMVHLM::drawVariable: Variable index [%d] does not exist.\n", index );
  }

}//end



void BayesianMVHLM::dataAugmentationInitialDraws()
{
  int i;

  if ( fixed_effects && t_gamma )
  {
    tau2_gamma->draw();
  }

  if ( second_stage_effects && t_alpha )
  {
    tau2_alpha->draw();
  }

  if ( random_effects && t_beta )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      tau2_betas[ i ]->draw();
    }
  }

  if ( t_lkhd )
  {
    for ( i = 0; i < number_of_observations; i++ )
    {    
      tau2_errors[ i ]->draw();
      updateRegressionWeight( i );
    }
  }

  if ( t_group_lkhd )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      tau2_groups[ i ]->draw();
      updateGroupRegressionWeight( i );
    }
  }

}//end



void BayesianMVHLM::createOutput( int simulations_to_keep )
{
  int i;

  number_of_simulations = simulations_to_keep;

  if ( simulations_to_keep > 0 )
  {
    if ( random_effects )
    {
      simulated_beta = new CMatrix * [ number_of_groups ];
      for ( i = 0; i < number_of_groups; i++ )
      {
        simulated_beta[i] = new CMatrix( simulations_to_keep,  beta[i]->dimension() );
      }

      if ( t_beta && keep_tau2_error )
      {
        simulated_tau2_beta = new CMatrix( simulations_to_keep, number_of_groups );
      }

      simulated_tau2 = new CMatrix * [ simulations_to_keep ];

      if ( second_stage_effects )
      {
        simulated_alpha = new CMatrix( simulations_to_keep, alpha->dimension() );
        if ( t_alpha && keep_tau2_error )
        {
          simulated_tau2_alpha = new CVector( simulations_to_keep );
        }
      }
    }

    if ( fixed_effects )
    {
      simulated_gamma = new CMatrix( simulations_to_keep, gamma->dimension() );
      if ( t_gamma && keep_tau2_error )
      {
        simulated_tau2_gamma = new CVector( simulations_to_keep );
      }
    }

    if ( !known_sigma2 )
    {
      simulated_sigma2 = new CMatrix ** [ simulations_to_keep ];
    }

    if ( t_lkhd && keep_tau2_error )
    {
      simulated_tau2_errors = new CMatrix( simulations_to_keep, number_of_observations );
    }

    if ( t_group_lkhd && keep_tau2_error )
    {
      simulated_tau2_groups = new CMatrix( simulations_to_keep, number_of_groups );
    }

    if ( keep_missing_data )
    {
      if ( number_missing_response != NULL )
      {
        simulated_missingR = new CMatrix * [ number_of_observations ];
        for ( i = 0; i < number_of_observations; i++ )
        {
          if ( number_missing_response->Val( i ) > 0 )
	  {
            simulated_missingR[i] = new CMatrix( simulations_to_keep, (int) number_missing_response->Val( i ) );
          }
          else
	  {
            simulated_missingR[i] = new CMatrix( 1, 1);
          }
        }//end for i
      }//end if response

      if ( number_missing_random_predictors != NULL )
      {
        simulated_missingRP = new CMatrix * [ number_of_observations ];
        for ( i = 0; i < number_of_observations; i++ )
        {
          if ( number_missing_random_predictors->Val( i ) > 0 )
	  {
            simulated_missingRP[i] = new CMatrix( simulations_to_keep, (int) number_missing_random_predictors->Val( i ) );
          }
          else
	  {
            simulated_missingRP[i] = new CMatrix( 1, 1);
          }
        }//end for i
      }//end if random predictors


      if ( number_missing_fixed_predictors != NULL )
      {
        simulated_missingFP = new CMatrix * [ number_of_observations ];
        for ( i = 0; i < number_of_observations; i++ )
        {
          if ( number_missing_fixed_predictors->Val( i ) > 0 )
	  {
            simulated_missingFP[i] = new CMatrix( simulations_to_keep, (int) number_missing_fixed_predictors->Val( i ) );
          }
          else
	  {
            simulated_missingFP[i] = new CMatrix( 1, 1);
          }
        }//end for i
      }//end if random predictors

    }//end if missing


  }//end if keep any simulations

}//end



void BayesianMVHLM::keepSimulation( int simul_number )
{
  int i, group, obs;

  if ( simul_number < number_of_simulations )
  {
    if ( random_effects )
    {
      for ( i = 0; i < number_of_groups; i++ )
      {
#ifdef FIX1
        CVector tmpvec = beta[i]->lastItemDrawn();
        simulated_beta[i]->setRow( simul_number, tmpvec );
#else
        simulated_beta[i]->setRow( simul_number, beta[i]->lastItemDrawn() );
#endif
//        simulated_beta[i]->setRow( simul_number, beta[i]->lastItemDrawn() );
      }

      if ( t_beta && keep_tau2_error )
      {
        for ( i = 0; i < number_of_groups; i++ )
	{
          simulated_tau2_beta->Val( simul_number, i ) = tau2_betas[i]->lastItemDrawn();
	}
      }

      simulated_tau2[ simul_number ] = new CMatrix( tau2->lastDraw().getMatrix() );

      if ( second_stage_effects )
      {
#ifdef FIX1
        CVector tmpvec = alpha->lastItemDrawn();
        simulated_alpha->setRow( simul_number, tmpvec );
#else
        simulated_alpha->setRow( simul_number, alpha->lastItemDrawn() );
#endif
//        simulated_alpha->setRow( simul_number, alpha->lastItemDrawn() );

        if ( t_alpha && keep_tau2_error )
        {
          simulated_tau2_alpha->Val( simul_number ) = tau2_alpha->lastItemDrawn();
        }
      }
    }

    if ( fixed_effects )
    {
#ifdef FIX1
      CVector tmpvec = gamma->lastItemDrawn();
      simulated_gamma->setRow( simul_number, tmpvec );
#else
      simulated_gamma->setRow( simul_number, gamma->lastItemDrawn() );
#endif
//      simulated_gamma->setRow( simul_number, gamma->lastItemDrawn() );

      if ( t_gamma && keep_tau2_error )
      {
        simulated_tau2_gamma->Val( simul_number ) = tau2_gamma->lastItemDrawn();
      }
    }


    if ( !known_sigma2 )
    {
      if ( common_sigma )
      {
        simulated_sigma2[ simul_number ] = new CMatrix * [1];
        simulated_sigma2[ simul_number ][ 0 ] =  new CMatrix( sigma2->lastDraw().getMatrix() );
      }
      else
      {
        simulated_sigma2[ simul_number ] = new CMatrix * [ number_of_groups ];
        for ( i = 0; i < number_of_groups; i++ )
        {
          simulated_sigma2[ simul_number ][ i ] = new CMatrix( sigma2->lastDrawFromComponent( i ).getMatrix() );
        }
      }
    }

    if ( t_lkhd && keep_tau2_error )
    {
      for ( i = 0; i < number_of_observations; i++ )
      {
        simulated_tau2_errors->Val( simul_number, i ) = tau2_errors[i]->lastItemDrawn();
      }
    }

    if ( t_group_lkhd && keep_tau2_error )
    {
      for ( i = 0; i < number_of_groups; i++ )
      {
        simulated_tau2_groups->Val( simul_number, i ) = tau2_groups[i]->lastItemDrawn();
      }
    }

    if ( keep_missing_data )
    {

      if ( number_missing_response != NULL )
      {
        for ( i = 0; i< number_of_observations; i++ )
        {
          if ( number_missing_response->Val(i) > 0 )
	  {
            DistributionParameter grp( groupIndex( i ) );
            group = (int) grp.getVector().Val(0);
            obs = (int) grp.getVector().Val(1);          

            //printf("keep missing data: [%d][%d] missing are:\n", group, obs ); fflush(stdout);
            //missing_response_components[group][obs]->Print(); fflush(stdout);

#ifdef FIX1
            CVector tmpvec = response[ group ]->getColumn( obs ).subVector( (*(missing_response_components[ group ][ obs ])) );
            simulated_missingR[i]->setRow( simul_number, tmpvec );
#else
            simulated_missingR[i]->setRow( simul_number, response[ group ]->getColumn( obs ).subVector( (*(missing_response_components[ group ][ obs ])) ) );
#endif
//            simulated_missingR[i]->setRow( simul_number, response[ group ]->getColumn( obs ).subVector( (*(missing_response_components[ group ][ obs ])) ) );
          }
        }//end for i
      }//end response

      if ( number_missing_random_predictors != NULL )
      {
        for ( i = 0; i< number_of_observations; i++ )
        {
          if ( number_missing_random_predictors->Val(i) > 0 )
	  {
            DistributionParameter grp( groupIndex( i ) );
            group = (int) grp.getVector().Val(0);
            obs = (int) grp.getVector().Val(1);          

            //printf("keep missing random predictors data: [%d][%d] missing are:\n", group, obs ); fflush(stdout);
            //missing_random_predictors_components[group][obs]->Print(); fflush(stdout);

#ifdef FIX1
            CVector tmpvec = random_predictors[ group ][ obs ]->getRow( 0 ).subVector( (*(missing_random_predictors_components[ group ][ obs ])) );
            simulated_missingRP[i]->setRow( simul_number, tmpvec );
#else
            simulated_missingRP[i]->setRow( simul_number, random_predictors[ group ][ obs ]->getRow( 0 ).subVector( (*(missing_random_predictors_components[ group ][ obs ])) ) );
#endif
//            simulated_missingRP[i]->setRow( simul_number, random_predictors[ group ][ obs ]->getRow( 0 ).subVector( (*(missing_random_predictors_components[ group ][ obs ])) ) );
          }
        }//end for i
      }//end random predictors

      if ( number_missing_fixed_predictors != NULL )
      {
        for ( i = 0; i< number_of_observations; i++ )
        {
          if ( number_missing_fixed_predictors->Val(i) > 0 )
	  {
            DistributionParameter grp( groupIndex( i ) );
            group = (int) grp.getVector().Val(0);
            obs = (int) grp.getVector().Val(1);          

            //printf("keep missing fixed predictors data: [%d][%d] missing are:\n", group, obs ); fflush(stdout);
            //missing_fixed_predictors_components[group][obs]->Print(); fflush(stdout);

#ifdef FIX1
            CVector tmpvec = fixed_predictors[ group ][ obs ]->getRow( 0 ).subVector( (*(missing_fixed_predictors_components[ group ][ obs ])) );
            simulated_missingFP[i]->setRow( simul_number, tmpvec );
#else
            simulated_missingFP[i]->setRow( simul_number, fixed_predictors[ group ][ obs ]->getRow( 0 ).subVector( (*(missing_fixed_predictors_components[ group ][ obs ])) ) );
#endif
//            simulated_missingFP[i]->setRow( simul_number, fixed_predictors[ group ][ obs ]->getRow( 0 ).subVector( (*(missing_fixed_predictors_components[ group ][ obs ])) ) );
          }
        }//end for i
      }//end fixed predictors

    }//end missing data

  }//end if keep any simulations

}//end
    


void BayesianMVHLM::simulationsToArray( double * simul_output, int simulations_to_keep )
{
  int i, j, k, q, total_dim, start_dim, g, v, s;

  start_dim = 0;
  if ( random_effects )
  {
    total_dim = simulations_to_keep * beta[0]->dimension();
    for ( g = 0; g < number_of_groups; g++ )
    {
      for ( v = 0; v < beta[g]->dimension(); v++ )
      {
        for ( s = 0; s < simulations_to_keep; s++ )
	{
          simul_output[ g * total_dim + v * simulations_to_keep + s ] = simulated_beta[g]->Val(s,v);
        }
      }
    }
    start_dim = number_of_groups * total_dim;
  }

  if ( fixed_effects )
  {
    total_dim = simulations_to_keep;
    for ( v = 0; v < gamma->dimension(); v++ )
    {
      for ( s = 0; s < simulations_to_keep; s++ )
      {
        simul_output[ start_dim + v * total_dim + s ] = simulated_gamma->Val(s,v);
      }
    }
    start_dim += gamma->dimension() * total_dim;
  }

  if ( second_stage_effects )
  {
    total_dim = simulations_to_keep;
    for ( v = 0; v < alpha->dimension(); v++ )
    {
      for ( s = 0; s < simulations_to_keep; s++ )
      {
        simul_output[ start_dim + v * total_dim + s ] = simulated_alpha->Val(s,v);
      }
    }
    start_dim += alpha->dimension() * total_dim;
  }


  if ( !known_sigma2 )
  {
    if ( common_sigma )
    {
      q = simulated_sigma2[0][0]->Col();
      total_dim = q * simulations_to_keep;
      for ( i = 0; i < q; i++ )
      {
        for ( j = 0; j < q; j++ )
        {
          for ( s = 0; s < simulations_to_keep; s++ )
          {
            if ( invWishart_sigma2 )
	    {
              simul_output[ start_dim + i * total_dim + j * simulations_to_keep + s ] = simulated_sigma2[s][0]->Val( i, j );
            }
            else
	    {
              simul_output[ start_dim + i * total_dim + j * simulations_to_keep + s ] = sqrt( simulated_sigma2[s][0]->Val( i, j ) );
            }
          }//end for s
        }
      }//end for i
      start_dim += q * total_dim;
    }
    else
    {
      q = simulated_sigma2[0][0]->Col();
      total_dim = simulations_to_keep * q*q;
      int partial_dim = simulations_to_keep * q;
      for ( i = 0; i < number_of_groups; i++ )
      {
        for ( k = 0; k < q; k++ )
        {
          for ( j = 0; j < q; j++ )
	  {
            for ( s = 0; s < simulations_to_keep; s++ )
            {
              if ( invWishart_sigma2 )
	      {
                simul_output[ start_dim + i * total_dim + k * partial_dim + j * simulations_to_keep + s ] = simulated_sigma2[s][i]->Val( k, j );
              }
              else
	      {
                simul_output[ start_dim + i * total_dim + k * partial_dim + j * simulations_to_keep + s ] = sqrt( simulated_sigma2[s][i]->Val( k, j ) );
              }
            }
          }
        }
      }//end for i
      start_dim += number_of_groups * total_dim;
    }
  }  

  if ( random_effects )
  {
    int q = simulated_tau2[0]->Col();
    total_dim = simulations_to_keep * q;
    for ( i = 0; i < q; i++ )
    {
      for ( j = 0; j < q; j++ )
      {   
        for ( s = 0; s < simulations_to_keep; s++ )
        {
          if ( invWishart_betaCov )
	  {
            simul_output[ start_dim + i * total_dim + j * simulations_to_keep + s ] = simulated_tau2[s]->Val(i,j);
          }
          else
	  {
            simul_output[ start_dim + i * total_dim + j * simulations_to_keep + s ] = sqrt( simulated_tau2[s]->Val(i,j) );
          }
        }
      }
    }//end for i
    start_dim += q * total_dim;
  }//end if random effects

  if ( random_effects && t_beta && keep_tau2_error )
  {
    total_dim = simulations_to_keep;
    for ( g = 0; g < number_of_groups; g++ )
    {    
      for ( s = 0; s < simulations_to_keep; s++ )
      {
        simul_output[ start_dim + g * total_dim + s ] =  simulated_tau2_beta->Val(s,g);
      }
    }    
    start_dim += number_of_groups * total_dim;
  }

  if ( fixed_effects && t_gamma && keep_tau2_error )
  {
    for ( s = 0; s < simulations_to_keep; s++ )
    {
      simul_output[ start_dim +  s ] = simulated_tau2_gamma->Val(s);
    }
    start_dim += simulations_to_keep;
  }

  if ( second_stage_effects && t_alpha && keep_tau2_error )
  {
    for ( s = 0; s < simulations_to_keep; s++ )
    {
      simul_output[ start_dim +  s ] = simulated_tau2_alpha->Val(s);
    }
    start_dim += simulations_to_keep;
  }

  if ( t_group_lkhd && keep_tau2_error )
  {
    total_dim = simulations_to_keep;
    for ( g = 0; g < number_of_groups; g++ )
    {    
      for ( s = 0; s < simulations_to_keep; s++ )
      {
        simul_output[ start_dim + g * total_dim + s ] =  simulated_tau2_groups->Val(s,g);
      }
      
    }    
    start_dim += number_of_groups * total_dim;
  }

  if ( t_lkhd && keep_tau2_error )
  {
    total_dim = simulations_to_keep;
    for ( i = 0; i < number_of_observations; i++ )
    {
      for ( s = 0; s < simulations_to_keep; s++ )
      { 
        simul_output[ start_dim + i * total_dim + s ] = simulated_tau2_errors->Val(s,i);
      }
    }
    start_dim += number_of_observations * total_dim;
  }

  if ( keep_missing_data )
  {
    k = 0;    
    if ( number_missing_response != NULL )
    {
      for ( i = 0; i < number_of_observations; i++ )
      {
        //printf( "Bayes: number missing[%d] = %f\n", i, number_missing_response->Val(i) );
        if ( number_missing_response->Val(i) > 0 )
        {
          for ( j = 0; j < number_missing_response->Val(i); j++ )
	  {
            for ( s = 0; s < simulations_to_keep; s++ )
	    {
              //printf( "missingR[%d].Val( %d, %d ) = %f\n", i,s,j, simulated_missingR[i]->Val( s, j ) );
              //fflush(stdout);

              simul_output[ start_dim + k ] = simulated_missingR[i]->Val( s, j );
              k++;
            }//end for s
          }//end for j
        }//end if
      }//end for i
      start_dim += k;
    }//end if response

    k = 0;
    if ( number_missing_random_predictors != NULL )
    {
      for ( i = 0; i < number_of_observations; i++ )
      {
        //printf( "Bayes: number random_predictors missing[%d] = %f\n", i, number_missing_random_predictors->Val(i) );
        if ( number_missing_random_predictors->Val(i) > 0 )
        {
          for ( j = 0; j < number_missing_random_predictors->Val(i); j++ )
	  {
            for ( s = 0; s < simulations_to_keep; s++ )
	    {
              //printf( "missingRP[%d].Val( %d, %d ) = %f\n", i,s,j, simulated_missingRP[i]->Val( s, j ) );
              //fflush(stdout);

              simul_output[ start_dim + k ] = simulated_missingRP[i]->Val( s, j );
              k++;
            }//end for s
          }//end for j
        }//end if
      }//end for i
      start_dim += k;
    }//end if random predictors

    k = 0;
    if ( number_missing_fixed_predictors != NULL )
    {
      for ( i = 0; i < number_of_observations; i++ )
      {
        //printf( "Bayes: number fixed_predictors missing[%d] = %f\n", i, number_missing_fixed_predictors->Val(i) );
        if ( number_missing_fixed_predictors->Val(i) > 0 )
        {
          for ( j = 0; j < number_missing_fixed_predictors->Val(i); j++ )
	  {
            for ( s = 0; s < simulations_to_keep; s++ )
	    {
              //printf( "missingFP[%d].Val( %d, %d ) = %f\n", i,s,j, simulated_missingFP[i]->Val( s, j ) );
              //fflush(stdout);

              simul_output[ start_dim + k ] = simulated_missingFP[i]->Val( s, j );
              k++;
            }//end for s
          }//end for j
        }//end if
      }//end for i
      start_dim += k;
    }//end if fixed predictors

  }//end if missing

  //printf( "Bayes: end: start dim = %d\n", start_dim); fflush(stdout);
}//end

