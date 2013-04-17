#include <math.h>
#include <string.h>
#include "R.h"
#include "rtErr.h"

#include "Const.h"
#include "DistributionParameter.h"
#include "BayesianHierarchicalLinearModel.h"

#define SQR(x) ((x)*(x))

extern int randomStart( int n );

/* constructor.
   returns: an empty Distribution of type BayesianLinearModel
*/
BayesianHierarchicalLinearModel::BayesianHierarchicalLinearModel()
{

  emptyModel();

}//end


void BayesianHierarchicalLinearModel::emptyModel()
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

  //second stage
  alpha = NULL;
  tau2 = NULL;
  tau2ICS = NULL;
  tau2PNIP = NULL;
  tau2NIP = NULL;
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
  dim_beta = 0;

  t_lkhd = false;
  t_beta = false;
  t_gamma = false;
  t_alpha = false;
  t_group_lkhd = false;
  non_informative_beta = false;
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



BayesianHierarchicalLinearModel::~BayesianHierarchicalLinearModel()
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
      for ( j = 0; j < response[i]->Len(); j++ )
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
      for ( j = 0; j < response[i]->Len(); j++ )
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
      for ( j = 0; j < response[i]->Len(); j++ )
      {
        delete missing_fixed_predictors_components[i][j];
      }
      delete [] missing_fixed_predictors_components[i];
    }
    delete [] missing_fixed_predictors_components;
  }


  if ( response != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete response[i];
    }
    delete [] response;
  }

  if ( random_effects )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete random_predictors[i];
    }
    delete [] random_predictors;
  }

  if ( fixed_effects )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete fixed_predictors[i];
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


  //second stage
  if ( alpha != NULL )
  {
    delete alpha;
  }

  if ( tau2 != NULL )
  {
  	if( invWishart_betaCov )
  		delete tau2[0];
  	else 
  		for( i = 0; i < dim_beta; i++ )
        delete tau2[i];
    delete [] tau2;
  }

  if ( tau2ICS != NULL )
  {
  	for( i = 0; i < dim_beta; i++ )
      delete tau2ICS[i];
    delete [] tau2ICS;
  }

  if ( tau2PNIP != NULL )
  {
  	for( i = 0; i < dim_beta; i++ )
  	  delete tau2PNIP[i];
    delete [] tau2PNIP;
  }

  if ( tau2InvWishart != NULL )
  {
    delete tau2InvWishart;
  }

  if ( tau2NIP != NULL )
  {
    delete [] tau2NIP;
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
  	if( invWishart_betaCov ){
  		delete tau2_first_draw[0];
    } else {
    	for( i = 0; i < dim_beta; i++ )
        delete tau2_first_draw[i];
    }
    delete [] tau2_first_draw;
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
    delete simulated_sigma2;
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



void BayesianHierarchicalLinearModel::initialize( CVector ** y, int n_groups )
{
  int i;

  number_of_groups = n_groups;
  response = new CVector * [ number_of_groups ];
  number_of_observations = 0;

  final_group_index = new CVector( number_of_groups );
  for ( i = 0; i < number_of_groups; i++ )
  {
    response[i] = new CVector( (*y[i]) );

    if ( i > 0 )
    {
      final_group_index->Val( i ) = final_group_index->Val( i - 1 ) + response[i]->Len();
    }
    else
    {
      final_group_index->Val( i ) = response[i]->Len() - 1;
    }

    number_of_observations += response[i]->Len();

    //printf("BHL: response[%d] = \n", i );
    //response[i]->Print();

  }
}//end



void BayesianHierarchicalLinearModel::initializeResponseMissingData( bool keep_them )
{
  int i, j;

  keep_missing_data = keep_them;

  number_missing_response = new CVector( number_of_observations ); 
  number_missing_response->setToZero();
  missing_response_components = new CVector ** [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    missing_response_components[i] = new CVector * [ response[i]->Len() ];
    for ( j = 0; j < response[i]->Len(); j++ )
    {
      missing_response_components[i][j] = new CVector( 1 );
      missing_response_components[i][j]->setToZero();
    }
  }
}//end


void BayesianHierarchicalLinearModel::responseMissingData( int index, CVector & m_comp )
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


void BayesianHierarchicalLinearModel::initializeRandomPredictorsMissingData( bool keep_them )
{
  int i, j;

  keep_missing_data = keep_them;

  number_missing_random_predictors = new CVector( number_of_observations ); 
  number_missing_random_predictors->setToZero();
  missing_random_predictors_components = new CVector ** [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    missing_random_predictors_components[i] = new CVector * [ response[i]->Len() ];
    for ( j = 0; j < response[i]->Len(); j++ )
    {
      missing_random_predictors_components[i][j] = new CVector( 1 );
      missing_random_predictors_components[i][j]->setToZero();
    }
  }
}//end


void BayesianHierarchicalLinearModel::randomPredictorsMissingData( int index, CVector & m_comp )
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


void BayesianHierarchicalLinearModel::initializeFixedPredictorsMissingData( bool keep_them )
{
  int i, j;

  keep_missing_data = keep_them;

  number_missing_fixed_predictors = new CVector( number_of_observations ); 
  number_missing_fixed_predictors->setToZero();
  missing_fixed_predictors_components = new CVector ** [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    missing_fixed_predictors_components[i] = new CVector * [ response[i]->Len() ];
    for ( j = 0; j < response[i]->Len(); j++ )
    {
      missing_fixed_predictors_components[i][j] = new CVector( 1 );
      missing_fixed_predictors_components[i][j]->setToZero();
    }
  }
}//end


void BayesianHierarchicalLinearModel::fixedPredictorsMissingData( int index, CVector & m_comp )
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



void BayesianHierarchicalLinearModel::randomEffects( CMatrix ** x )
{
  int i;

  random_predictors = new CMatrix * [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    random_predictors[i] = new CMatrix( (*x[i]) );

    //printf("BHL: random predictors[%d] = \n", i );
    //random_predictors[i]->Print();
  }

  number_of_variables += number_of_groups;

  random_effects = true;  
  dim_beta = random_predictors[0]->Col();
}//end


void BayesianHierarchicalLinearModel::fixedEffects( CMatrix ** m )
{
  int i;

  fixed_predictors = new CMatrix * [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    fixed_predictors[i] = new CMatrix( (*m[i]) );

    //printf("BHL: fixed predictors[%d] = \n", i );
    //fixed_predictors[i]->Print();
  }
  //count gamma
  number_of_variables++;

  fixed_effects = true;  
}//end


void BayesianHierarchicalLinearModel::secondStageRandomEffects( CMatrix ** z )
{
  int i;

  beta_predictors = new CMatrix * [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    beta_predictors[i] = new CMatrix( (*z[i]) );

    //printf("BHL: level2 predictors[%d] = \n", i );
    //beta_predictors[i]->Print();
  }
  //count alpha
  number_of_variables++;

  second_stage_effects = true;
}//end


void BayesianHierarchicalLinearModel::betaPrior( CMatrix * p_betaCov )
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


void BayesianHierarchicalLinearModel::betaPrior( CVector * p_beta, CMatrix * p_betaCov )
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


void BayesianHierarchicalLinearModel::betaPrior( CMatrix * p_beta, CMatrix * p_betaCov )
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


void BayesianHierarchicalLinearModel::betaPriorNonInformative( int dim )
{
  int i;

  beta = new NormalDistribution * [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    beta[i] = new NormalDistribution( dim );
    beta[i]->setNonInformative();
  }

  non_informative_beta = true;

}//end


void BayesianHierarchicalLinearModel::gammaPriorNonInformative( int dim )
{
  gamma = new NormalDistribution( dim );
  gamma->setNonInformative();

  non_informative_gamma = true;

  //printf( "BHL: noninfo gamma prior\n" );
  //gamma->initialMean().Print();
  //gamma->initialCovariance().Print();

}//end


void BayesianHierarchicalLinearModel::gammaPrior( CVector * p_gamma, CMatrix * p_gammaCov )    
{
  gamma = new NormalDistribution( p_gamma->Len() );
  gamma->setMean( p_gamma );
  gamma->setCovariance( p_gammaCov );

  //printf( "BHL: gamma prior\n" );
  //gamma->initialMean().Print();
  //gamma->initialCovariance().Print();

}//end



void BayesianHierarchicalLinearModel::alphaPrior( CVector * p_alpha, CMatrix * p_alphaCov ) throw( rtErr )   
{
  if ( second_stage_effects )
  {
    alpha = new NormalDistribution( p_alpha->Len() );
    alpha->setMean( p_alpha );
    alpha->setCovariance( p_alphaCov );
  }
  else
  {
    printf( "BayesianHierarchicalLinearModel::alphaPrior: alpha is not a parameter of this model.\n" );
    char the_error[] = "BayesianHierarchicalLinearModel::alphaPrior: alpha is not a parameter of this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

  //printf( "BHL: alpha prior\n" );
  //alpha->initialMean().Print();
  //alpha->initialCovariance().Print();

}//end


void BayesianHierarchicalLinearModel::alphaPriorNonInformative( int dim ) throw( rtErr )
{
  if ( second_stage_effects )
  {
    alpha = new NormalDistribution( dim );
    alpha->setNonInformative();
    non_informative_alpha = true;
  }
  else
  {
    printf( "BayesianHierarchicalLinearModel::alphaPriorNonInformative: alpha is not a parameter of this model.\n" );
    char the_error[] = "BayesianHierarchicalLinearModel::alphaPriorNonInformative: alpha is not a parameter of this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

  //printf( "BHL: alpha prior: non informative\n" );
  //alpha->initialMean().Print();
  //alpha->initialCovariance().Print();

}//end



void BayesianHierarchicalLinearModel::sigma2CommonPrior( bool val )
{
  common_sigma = val;
}//end



void BayesianHierarchicalLinearModel::sigma2PriorInvChisq( double p_nuSigma2, double p_sigma2 )
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



void BayesianHierarchicalLinearModel::sigma2PriorInvChisq( double * p_nuSigma2, double * p_sigma2 )
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


void BayesianHierarchicalLinearModel::sigma2PriorDuMouchel( double p_sigma2 )
{
  double df;

  df = number_of_observations + 1;
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


void BayesianHierarchicalLinearModel::sigma2PriorDuMouchel( double * p_sigma2 )
{
  double df;

  df = number_of_observations + 1;
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



void BayesianHierarchicalLinearModel::sigma2PriorUniformShrinkage( double p_sigma2 )
{
  double df;

  df = number_of_observations;
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


void BayesianHierarchicalLinearModel::sigma2PriorUniformShrinkage( double * p_sigma2 )
{
  double df;

  df = number_of_observations;
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



void BayesianHierarchicalLinearModel::sigma2PriorNonInformative( double p_power ) throw( rtErr )
{

  nonInformativePower_sigma2 = true;

  if ( common_sigma )
  {
    if ( p_power <= -0.5 * number_of_observations )
    {
      printf( "BayesianHierarchicalLinearModel::sigma2PriorNonInformative: Power [%f] in noninformative prior is not valid.\n", p_power );  
      char the_error[] = "BayesianHierarchicalLinearModel::sigma2PriorNonInformative: Power in noninformative prior is not valid.";
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
      if ( p_power <= -0.5 * response[i]->Len() )
      {
        printf( "BayesianHierarchicalLinearModel::sigma2PriorNonInformative: Power [%f] in noninformative prior is not valid.\n", p_power );  
	char the_error[] = "BayesianHierarchicalLinearModel::sigma2PriorNonInformative: Power in noninformative prior is not valid.";
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


void BayesianHierarchicalLinearModel::sigma2PriorNonInformative( double * p_power ) throw( rtErr )
{

  nonInformativePower_sigma2 = true;

  if ( common_sigma )
  {
    if ( p_power[0] <= -0.5 * number_of_observations )
    {
      printf( "BayesianHierarchicalLinearModel::sigma2PriorNonInformative: Power [%f] in noninformative prior is not valid.\n", p_power );  
      char the_error[] = "BayesianHierarchicalLinearModel::sigma2PriorNonInformative: Power in noninformative prior is not valid.";
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
      if ( p_power[i] <= -0.5 * response[i]->Len() )
      {
        printf( "BayesianHierarchicalLinearModel::sigma2PriorNonInformative: Power [%f] in noninformative prior is not valid.\n", p_power );  
        char the_error[] = "BayesianHierarchicalLinearModel::sigma2PriorNonInformative: Power in noninformative prior is not valid.";
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


void BayesianHierarchicalLinearModel::sigma2Known( double * p_sigma2 )
{
  known_sigma2 = true;

  if ( common_sigma )
  {
    sigma2 = new DistributionMixture( 1 );
    sigma2ICS = new InvChisqDistribution * [ 1 ];
    sigma2ICS[0] = new InvChisqDistribution( 1 );
    sigma2ICS[0]->setScale( p_sigma2[0] );
    sigma2_first_draw = new DistributionParameter * [ 1 ];
    sigma2_first_draw[0] = new DistributionParameter( p_sigma2[0] );
    sigma2->set( 0, sigma2ICS[0], 1.0 );
  }
  else
  {
    int i;

    sigma2 = new DistributionMixture( number_of_groups );
    sigma2ICS = new InvChisqDistribution * [ number_of_groups ];
    sigma2_first_draw = new DistributionParameter * [ number_of_groups ];
    for ( i = 0; i < number_of_groups; i++ )
    {
      sigma2ICS[i] = new InvChisqDistribution( 1 );
      sigma2ICS[i]->setScale( p_sigma2[i] );
      sigma2_first_draw[i] = new DistributionParameter( p_sigma2[i] );
      sigma2->set( i, sigma2ICS[i], 1.0 );
    }
  }
}//end


void BayesianHierarchicalLinearModel::sigma2Known( double p_sigma2 )
{
  known_sigma2 = true;

  if ( common_sigma )
  {
    sigma2 = new DistributionMixture( 1 );
    sigma2ICS = new InvChisqDistribution * [ 1 ];
    sigma2ICS[0] = new InvChisqDistribution( 1 );
    sigma2ICS[0]->setScale( p_sigma2 );
    sigma2_first_draw = new DistributionParameter * [ 1 ];
    sigma2_first_draw[0] = new DistributionParameter( p_sigma2 );
    sigma2->set( 0, sigma2ICS[0], 1.0 );
  }
  else
  {
    int i;

    sigma2 = new DistributionMixture( number_of_groups );
    sigma2ICS = new InvChisqDistribution * [ number_of_groups ];
    sigma2_first_draw = new DistributionParameter * [ number_of_groups ];
    for ( i = 0; i < number_of_groups; i++ )
    {
      sigma2ICS[i] = new InvChisqDistribution( 1 );
      sigma2ICS[i]->setScale( p_sigma2 );
      sigma2_first_draw[i] = new DistributionParameter( p_sigma2 );
      sigma2->set( i, sigma2ICS[i], 1.0 );
    }
  }
}//end




void BayesianHierarchicalLinearModel::tau2PriorInvChisq( double *p_nuTau2, double *p_tau2 ) throw( rtErr )
{
  if ( random_effects )
  {
  	int i;
  	tau2 = new DistributionMixture * [ dim_beta ];
  	tau2ICS = new InvChisqDistribution * [ dim_beta ];
  	tau2_first_draw = new DistributionParameter * [ dim_beta ];
    for( i = 0; i < dim_beta; i++ ){
    	// for each random effect, specify the prior for the variance
      tau2[i] = new DistributionMixture( 1 );    	
      tau2ICS[i] = new InvChisqDistribution( p_nuTau2[i] );
      tau2ICS[i]->setScale( p_tau2[i] );
      tau2[i]->set( 0, tau2ICS[i], 1.0 );
      tau2_first_draw[i] = new DistributionParameter( p_tau2[i] );
    }
    invChisq_tau2 = true;
    // Only increase the number of variables by 1, since the tau2s will 
    // be updated together.
    number_of_variables++;
  }
  else
  {
    printf( "BayesianHierarchicalLinearModel::tau2PriorInvChisq: There are no random effects in the model. Tau2 does not make sense in this model.\n" );
    char the_error[] = "BayesianHierarchicalLinearModel::tau2PriorInvChisq: There are no random effects in the model. Tau2 does not make sense in this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

}//end


void BayesianHierarchicalLinearModel::tau2PriorDuMouchel( double *p_tau2 ) throw( rtErr )
{
  double df;

  if ( random_effects )
  {
  	int i;
  	tau2 = new DistributionMixture * [ dim_beta ];
  	tau2PNIP = new ProperNonInfoPosteriorHLM * [ dim_beta ];
  	tau2_first_draw = new DistributionParameter * [ dim_beta ];
    df = number_of_groups * dim_beta + 1;
    for( i = 0; i < dim_beta; i++ ){
    	// for each random effect, specify the prior for the variance
      tau2[i] = new DistributionMixture( 1 );    	
      tau2PNIP[i] = new ProperNonInfoPosteriorHLM( p_tau2[i], df );
      tau2PNIP[i]->setPriorDistribution( "du_mouchel" );
      tau2_first_draw[i] = new DistributionParameter( p_tau2[i] );
      tau2[i]->set( 0, tau2PNIP[i], 1.0 );
  	}
    duMouchel_tau2 = true;
    properNonInfoPrior_tau2 = true;
    // Only increase the number of variables by 1, since the tau2s will 
    // be updated together.
    number_of_variables++;
  }
  else
  {
    printf( "BayesianHierarchicalLinearModel::tau2PriorDuMouchel: There are no random effects in the model. Tau2 does not make sense in this model.\n" );
    char the_error[] = "BayesianHierarchicalLinearModel::tau2PriorDuMouchel: There are no random effects in the model. Tau2 does not make sense in this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}//end


void BayesianHierarchicalLinearModel::tau2PriorUniformShrinkage( double *p_tau2 ) throw( rtErr )
{
  double df;

  if ( random_effects )
  {
  	int i;
  	tau2 = new DistributionMixture * [ dim_beta ];
  	tau2PNIP = new ProperNonInfoPosteriorHLM * [ dim_beta ];
  	tau2_first_draw = new DistributionParameter * [ dim_beta ];
    df = number_of_groups * dim_beta;
    for( i = 0; i < dim_beta; i++ ){
    	// for each random effect, specify the prior for the variance
      tau2[i] = new DistributionMixture( 1 );    	
      tau2PNIP[i] = new ProperNonInfoPosteriorHLM( p_tau2[i], df );
      tau2PNIP[i]->setPriorDistribution( "uniform_shrinkage" );
      tau2_first_draw[i] = new DistributionParameter( p_tau2[i] );
      tau2[i]->set( 0, tau2PNIP[i], 1.0 );
    }
   	  
    uniformShrinkage_tau2 = true;
    properNonInfoPrior_tau2 = true;
    // Only increase the number of variables by 1, since the tau2s will 
    // be updated together.
    number_of_variables++;
  }
  else
  {
    printf( "BayesianHierarchicalLinearModel::tau2PriorUniformShrinkage: There are no random effects in the model. Tau2 does not make sense in this model.\n" );
    char the_error[] = "BayesianHierarchicalLinearModel::tau2PriorUniformShrinkage: There are no random effects in the model. Tau2 does not make sense in this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}//end


void BayesianHierarchicalLinearModel::tau2PriorNonInformative( double *p_power ) throw( rtErr )
{

  if ( random_effects )
  {
    //printf("BHLM: p_power = %f, n groups = %d,  p = %d\n", p_power, number_of_groups, random_predictors[0]->Col()); fflush(stdout);
  	int i;
  	tau2 = new DistributionMixture * [ dim_beta ];
    tau2ICS = new InvChisqDistribution * [ dim_beta ];
  	tau2NIP = new double [ dim_beta ];
  	tau2_first_draw = new DistributionParameter * [ dim_beta ];
    int df = number_of_groups * dim_beta + 1;

    for( i = 0; i < dim_beta; i++ ){
    	// for each random effect, specify the prior for the variance
  	
      if ( p_power[i] <= -0.5 * number_of_groups * dim_beta )
      {
        printf( "BayesianHierarchicalLinearModel::tau2PriorNonInformative: Power [%f] in noninformative prior is not valid.\n", p_power[i] );  
        char the_error[] = "BayesianHierarchicalLinearModel::tau2PriorNonInformative: Power in noninformative prior is not valid.";
        rtErr runtime_error( the_error );
        throw runtime_error;
      }

      tau2NIP[i] = p_power[i];
      tau2[i] = new DistributionMixture( 1 );
      tau2ICS[i] = new InvChisqDistribution( 0 );
      tau2ICS[i]->setScale( 1.0 );
      tau2_first_draw[i] = new DistributionParameter( 1.0 );
      tau2[i]->set( 0, tau2ICS[i], 1.0 );
    }
    nonInformativePower_tau2 = true;
    // Only increase the number of variables by 1, since the tau2s will 
    // be updated together.
    number_of_variables++;
  }
  else
  {
    printf( "BayesianHierarchicalLinearModel::tau2PriorNonInformative: There are no random effects in the model. Tau2 does not make sense in this model.\n" );
    char the_error[] = "BayesianHierarchicalLinearModel::tau2PriorNonInformative: There are no random effects in the model. Tau2 does not make sense in this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

}//end



void BayesianHierarchicalLinearModel::betaCovPriorInvWishart( double p_nuV, CMatrix * p_VCov ) throw( rtErr )
{
  if ( random_effects )
  {
  	tau2 = new DistributionMixture *[ 1 ];
    tau2[0] = new DistributionMixture( 1 );
    tau2InvWishart = new InvWishartDistribution( p_nuV );
    tau2InvWishart->setScaleMatrix( p_VCov );
    tau2_first_draw = new DistributionParameter *[ 1 ];
    tau2_first_draw[0] = new DistributionParameter( (*p_VCov) );
    tau2[0]->set( 0, tau2InvWishart, 1.0 );
    invWishart_betaCov = true;
  
    number_of_variables++;
  }
  else
  {
    printf( "BayesianHierarchicalLinearModel::betaCovPriorInvWishart: There are no random effects in the model. Tau2 does not make sense in this model.\n" );
    char the_error[] = "BayesianHierarchicalLinearModel::betaCovPriorInvWishart: There are no random effects in the model. Tau2 does not make sense in this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

  //printf( "bhlm: set prior Inv Wishart. Initial Cov is \n" ); fflush(stdout);
  //tau2InvWishart->scaleMatrix().Print();
}//end



void BayesianHierarchicalLinearModel::tLikelihood( double p_nuError )
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



void BayesianHierarchicalLinearModel::groupTLikelihood( double p_nuGroup )
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



void BayesianHierarchicalLinearModel::gammaTPrior( double p_nuGamma ) throw( rtErr )
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
    printf( "BayesianHierarchicalLinearModel::gammaTPrior: gamma is not a parameter of this model.\n" );
    char the_error[] = "BayesianHierarchicalLinearModel::gammaTPrior: gamma is not a parameter of this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

}//end



void BayesianHierarchicalLinearModel::betaTPrior( double p_nuBeta ) throw( rtErr )
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
    printf( "BayesianHierarchicalLinearModel::betaTPrior: beta is not a parameter of this model.\n" );
    char the_error[] = "BayesianHierarchicalLinearModel::betaTPrior: beta is not a parameter of this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}//end


void BayesianHierarchicalLinearModel::alphaTPrior( double p_nuAlpha ) throw( rtErr )
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
    printf( "BayesianHierarchicalLinearModel::alphaTPrior: alpha is not a parameter of this model.\n" );
    char the_error[] = "BayesianHierarchicalLinearModel::alphaTPrior: alpha is not a parameter of this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}//end



void BayesianHierarchicalLinearModel::samplerDefaultInitialPoint()
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

    if ( !non_informative_beta )
    {
    	if( invWishart_betaCov ){
    		tau2[0]->setLastDraw( *(tau2_first_draw[0]) );
      } else {
      	for( i=0; i<dim_beta; i++ )
          tau2[i]->setLastDraw( *(tau2_first_draw[i]) );
      }
    }

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


void BayesianHierarchicalLinearModel::samplerSigma2InitialPoint( double init_sigma2 )
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


void BayesianHierarchicalLinearModel::samplerSigma2InitialPoint( CVector & init_sigma2 )
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



void BayesianHierarchicalLinearModel::samplerTau2InitialPoint( double *init_tau2 )
{
  if ( !invWishart_betaCov )
  {
  	int i;
  	tau2_first_draw = new DistributionParameter * [ dim_beta ];
  	for( i = 0; i < dim_beta; i++ ){
      tau2_first_draw[i] = new DistributionParameter( init_tau2[i] );
      tau2[i]->setLastDraw( *(tau2_first_draw[i]) );
    }
  }
  else
  {
    printf( "BayesianHierarchicalLinearModel::samplerTau2InitialPoint: Initial point should be a matrix, not a vector.\n" );
  }
}//end


void BayesianHierarchicalLinearModel::samplerTau2InitialPoint( CMatrix & init_tau2 )
{
  if ( invWishart_betaCov )
  {
  	tau2_first_draw = new DistributionParameter * [ 1 ];
    tau2_first_draw[0] = new DistributionParameter( init_tau2 );
    tau2[0]->setLastDraw( *(tau2_first_draw[0]) );

    //printf( "bhlm: tau2 set last draw is \n" ); fflush(stdout);
    //tau2[0]->lastDraw().getMatrix().Print();
  }
  else
  {
    printf( "BayesianHierarchicalLinearModel::samplerTau2InitialPoint: Initial point should not be a matrix, but a scalar.\n" );
  }
}//end


void BayesianHierarchicalLinearModel::samplerAlphaInitialPoint( CVector & init_alpha )
{
  alpha->setLastDraw( init_alpha );
}//end


void BayesianHierarchicalLinearModel::samplerBetaInitialPoint( CVector & init_beta )
{
  int i;
  
  for ( i = 0; i < number_of_groups; i++ )
  {
    beta[ i ]->setLastDraw( init_beta );
  }
}//end


void BayesianHierarchicalLinearModel::samplerBetaInitialPoint( CMatrix & init_beta )
{
  int i;
  
  for ( i = 0; i < number_of_groups; i++ )
  {
    beta[ i ]->setLastDraw( init_beta.getColumn( i ) );
  }
}//end



void BayesianHierarchicalLinearModel::samplerGammaInitialPoint( CVector & init_gamma )
{
  gamma->setLastDraw( init_gamma );

}//end



void BayesianHierarchicalLinearModel::samplerMissingVariablesInitialPoint()
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



void BayesianHierarchicalLinearModel::samplerMissingResponseInitialPoint()
{
  int index, group, obs, nmss;

  //compute means for each group ignoring missing data

  DistributionParameter * grp;
  CVector response_means( number_of_groups );;
  CVector n_items( number_of_groups );
  CVector missing_in_group( number_of_groups );
  
  missing_in_group.setToZero();
  response_means.setToZero();
  n_items.setToZero();

  index = 0;
  group = 0;
  obs = 0;
  nmss = (int) number_missing_response->Val( 0 );

  while ( index < number_of_observations )
  {
    if ( nmss > 0 )
    {
      missing_in_group.Val( group ) += 1;
    }
    else
    {
      response_means.Val( group ) += response[ group ]->Val( obs );
      n_items.Val( group ) += 1;
    }

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
      if ( n_items.Val( group ) > 0 )
      {
        response_means.Val( group ) = response_means.Val( group ) / n_items.Val( group );
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
      response[ group ]->Val( obs ) = response_means.Val( group );
   
      //printf("initial response[%d][%d] is %f \n", group, obs, response[ group ]->Val( obs )  );
    }

    index++;
    grp = new DistributionParameter( groupIndex( index ) );
    group = (int) grp->getVector().Val(0);
    obs = (int) grp->getVector().Val(1);
    nmss = (int) number_missing_response->Val( index );

    delete grp;
  }//end while index


}//end



void BayesianHierarchicalLinearModel::samplerMissingRandomPredictorsInitialPoint()
{
  int i, k, index, group, obs, nmss, p_dim;
  double diff;

  //compute means for each group ignoring missing data
  p_dim = random_predictors[ 0 ]->Col();

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
      x_pred = random_predictors[ group ]->getRow( obs );

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
      x_pred = random_predictors[ group ]->getRow( obs );

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
        random_predictors[ group ]->Val( obs, k ) = randomEffects_means[ group ]->Val( k );

        //printf("initial random predictors[%d][%d] is \n", group, obs );
        //random_predictors[ group ]->getRow( obs ).Print();  

        //printf("initial random predictors variances[%d][%d] is \n", group, obs );
        //randomEffects_vars[ group ]->Print();  

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



void BayesianHierarchicalLinearModel::samplerMissingFixedPredictorsInitialPoint()
{
  int i, k, index, group, obs, nmss, p_dim;
  double diff;

  //compute means for each group ignoring missing data
  p_dim = fixed_predictors[ 0 ]->Col();

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
      x_pred = fixed_predictors[ group ]->getRow( obs );

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
      x_pred = fixed_predictors[ group ]->getRow( obs );

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
        fixed_predictors[ group ]->Val( obs, k ) = fixedEffects_means[ group ]->Val( k );

        //printf("initial fixed predcitors[%d][%d] is \n", group, obs );
        //fixed_predictors[ group ]->getRow( obs ).Print();  
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





DistributionParameter BayesianHierarchicalLinearModel::groupIndex( int index )
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



void BayesianHierarchicalLinearModel::updateRegressionWeight( int index ) throw( rtErr )
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
    printf( "BayesianHierarchicalLinearModel::updateRegressionWeight: Index [%d] out of range.\n", index );
    char the_error[] = "BayesianHierarchicalLinearModel::updateRegressionWeight: Index out of range.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}//end


void BayesianHierarchicalLinearModel::updateGroupRegressionWeight( int grp )
{
  int i;

  for ( i = 0; i < response[ grp ]->Len(); i++ )
  {
    tau2_error_weights[ grp ]->Val( i ) = 1 / tau2_groups[ grp ]->lastItemDrawn();
  }
}//end


void BayesianHierarchicalLinearModel::initializeTemporaryStructures()
{
  int i, k, len_digits_vars, len_digits_obs;

  //create working matrices
  if ( random_effects )
  {
    xTx = new CMatrix * [ number_of_groups ];
    xTy = new CVector * [ number_of_groups ];
    for ( i = 0; i < number_of_groups; i++ )
    {
      xTx[i] = new CMatrix( random_predictors[i]->Col(), random_predictors[i]->Col() );
      xTy[i] = new CVector( random_predictors[i]->Col() ); 

     ( *(xTx[i]) ) = random_predictors[i]->xTransposedX();
     ( *(xTy[i]) ) = random_predictors[i]->T() * ( *(response[i]) );

    }//end loop
  }

  if ( fixed_effects )
  {
    sMTm = new CMatrix( fixed_predictors[0]->Col(), fixed_predictors[0]->Col() );
    sMTy = new CVector( fixed_predictors[0]->Col() ); 
    sMTm->setToZero();
    sMTy->setToZero();

    mTm = new CMatrix * [ number_of_groups ];
    mTy = new CVector * [ number_of_groups ];
    for ( i = 0; i < number_of_groups; i++ )
    {
      mTm[i] = new CMatrix( fixed_predictors[i]->Col(), fixed_predictors[i]->Col() );
      mTy[i] = new CVector( fixed_predictors[i]->Col() ); 

     ( *(mTm[i]) ) = fixed_predictors[i]->xTransposedX();
     ( *(mTy[i]) ) = fixed_predictors[i]->T() * ( *(response[i]) );

     sMTm->add( ( *(mTm[i]) ) );
     sMTy->add( ( *(mTy[i]) ) );

    }//end loop
  }

  //second stage
  if ( second_stage_effects )
  {
    zTvIz = new CMatrix( beta_predictors[0]->Col(), beta_predictors[0]->Col() );
    zTvIb = new CVector( beta_predictors[0]->Col() );
    zTvIz->setToZero();
    zTvIb->setToZero();

    zTz = new CMatrix * [ number_of_groups ];
    zTb = new CVector * [ number_of_groups ];
    for ( i = 0; i < number_of_groups; i++ )
    {
      zTz[i] = new CMatrix( beta_predictors[i]->Col(), beta_predictors[i]->Col() );
      zTb[i] = new CVector( beta_predictors[i]->Col() ); 

      ( *(zTz[i]) ) = beta_predictors[i]->T() * ( beta[i]->initialInverseCovariance() * ( (*beta_predictors[i]) ) );
      ( *(zTb[i]) ) = beta_predictors[i]->T() * ( beta[i]->initialInverseCovariance() * beta[i]->lastDraw().getVector() );

      zTvIz->add( ( *(zTz[i]) ) );
      zTvIb->add( ( *(zTb[i]) ) );
    }//end loop
  }


  //initialize residuals
  residuals = new CVector * [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    residuals[i] = new CVector( response[i]->Len() );
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

  //printf( "initializeTemporaryStructures: number of variables is %d\n", number_of_variables ); fflush(stdout);

  len_digits_vars = ( (int) ( log( (double) number_of_variables ) / log( 10.0 ) ) ) + 1;
  len_digits_obs = ( (int) ( log( (double) number_of_observations ) / log( 10.0 ) ) ) + 1;

  distr_map = new char * [ number_of_variables ];
  k = 0;

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


  if ( random_effects && !non_informative_beta )
  {
    distr_map[k] = new char [ 5 ];	
    sprintf( distr_map[k], "tau2" );
    k++;
  }

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
      tau2_error_weights[i] = new CVector( response[i]->Len() );
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
      tau2_error_weights[i] = new CVector( response[i]->Len() );
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


  //printf( "BHL: total vars = %d.  Number vars = %d.x The map is:\n", k, number_of_variables  );
  //for ( i = 0; i < k; i++ )
  //{
  //  printf("map[%d] = %s\n", i, distr_map[ i ] );
  //}
  //fflush( stdout );


}//end




void BayesianHierarchicalLinearModel::gibbsUpdateWorkingMatrices()
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


void BayesianHierarchicalLinearModel::gibbsUpdateWorkingMatrices( int index, char * type )
{
  //printf( "bhlm: gibbsUpdateWorkingMatrices\n" ); fflush(stdout);

  if ( t_lkhd || t_group_lkhd )
  {
    //weighted (or transformed) predictors. 
    //update must be made before updating beta
    if ( !strcmp( type, "random_effects" ) )
    {
      ( *(xTx[ index ]) ) = random_predictors[ index ]->xTransposedX( ( *(tau2_error_weights[ index ]) ) );

      if ( fixed_effects )
      {
        ( *(xTy[ index ]) ) = ( random_predictors[ index ]->T() ) * ( (*(response[ index ])) - ( (*(fixed_predictors[ index ])) * gamma->lastDraw().getVector() ) ).weighted( ( *(tau2_error_weights[ index ]) ) );
      }
      else
      {
        ( *(xTy[ index ]) ) = ( random_predictors[ index ]->T() ) * ( response[ index ]->weighted( ( *(tau2_error_weights[ index ]) ) ) );
      }
    }

    if ( !strcmp( type, "fixed_effects" ) )
    {
      //update sum matrix: substract current value 
#ifdef FIX1
      CVector tmpvec = ( ( *(mTy[ index ]) ) * -1.0 );
      CMatrix tmpmat = ( -1.0 * ( *(mTm[ index ]) ) );
      sMTy->add( tmpvec );
      sMTm->add( tmpmat );
#else
      sMTy->add( ( ( *(mTy[ index ]) ) * -1.0 ) );
      sMTm->add( ( -1.0 * ( *(mTm[ index ]) ) ) );
#endif
//      sMTy->add( ( ( *(mTy[ index ]) ) * -1.0 ) );
//      sMTm->add( ( -1.0 * ( *(mTm[ index ]) ) ) );

      if ( common_sigma )
      {
        ( *(mTm[ index ]) ) = fixed_predictors[ index ]->xTransposedX( (*tau2_error_weights[ index ]) );
        if ( random_effects )
        {
          ( *(mTy[ index ]) ) = ( fixed_predictors[ index ]->T() ) * ( (*(response[ index ])) - ( (*(random_predictors[ index ])) * beta[ index ]->lastDraw().getVector() ) ).weighted( (*(tau2_error_weights[ index ])) );
        }
        else
        {
          ( *(mTy[ index ]) ) = ( fixed_predictors[ index ]->T() ) * ( response[ index ]->weighted( (*tau2_error_weights[ index ]) ) );
        }
      }
      else
      {
        //printf( "working matrix: sigma2->lastDrawFromComponent( %d ) = %f\n", index,  sigma2->lastDrawFromComponent( index ).getScalar() ); fflush( stdout );

        ( *(mTm[ index ]) ) = ( 1.0 / sigma2->lastDrawFromComponent( index ).getScalar() ) * fixed_predictors[ index ]->xTransposedX( (*tau2_error_weights[ index ]) );
        if ( random_effects )
        {
          ( *(mTy[ index ]) ) = ( fixed_predictors[ index ]->T() * ( (*(response[ index ])) - ( (*(random_predictors[ index ])) * beta[ index ]->lastDraw().getVector() ) ).weighted( (*(tau2_error_weights[ index ])) ) ) * ( 1.0 / sigma2->lastDrawFromComponent( index ).getScalar() );
        }
        else
        {
          ( *(mTy[ index ]) ) = ( fixed_predictors[ index ]->T() * ( response[ index ]->weighted( (*(tau2_error_weights[ index ])) ) ) ) * ( 1.0 / sigma2->lastDrawFromComponent( index ).getScalar() );
        }
      }

      //update sum matrix: add updated value
      sMTm->add( ( *(mTm[ index ]) ) );
      sMTy->add( ( *(mTy[ index ]) ) );

    }
  }//end if t lklhd
  else if ( !strcmp( type, "random_effects" ) && fixed_effects )
  {
    ( *(xTy[ index ]) ) = random_predictors[ index ]->T() * ( (*(response[ index ])) - ( (*(fixed_predictors[ index ])) * gamma->lastDraw().getVector() ) );
  }
  else if ( !strcmp( type, "fixed_effects" ) && random_effects )
  {
    //update sum matrix: substract current value 
#ifdef FIX1
      CVector tmpvec = ( ( *(mTy[ index ]) ) * (-1.0) );
      sMTy->add( tmpvec );
#else
      sMTy->add( ( ( *(mTy[ index ]) ) * (-1.0) ) );
#endif
//    sMTy->add( ( ( *(mTy[ index ]) ) * (-1.0) ) );

    if ( common_sigma )
    {
      ( *(mTy[ index ]) ) = ( fixed_predictors[ index ]->T() ) * ( (*(response[ index ])) - ( (*(random_predictors[ index ])) * beta[ index ]->lastDraw().getVector() ) );
    }
    else
    {
      ( *(mTy[ index ]) ) = ( ( fixed_predictors[ index ]->T() ) * ( (*(response[ index ])) - ( (*(random_predictors[ index ])) * beta[ index ]->lastDraw().getVector() ) ) ) * ( 1.0 / sigma2->lastDrawFromComponent( index ).getScalar() );
    }

    //update sum matrix: add updated value
    sMTy->add( ( *(mTy[ index ]) ) );
  }

}//end




void BayesianHierarchicalLinearModel::gibbsUpdateSecondStageWorkingMatrices()
{
  int index, j;

  zTvIb->setToZero();
  CVector tmpVec ( dim_beta );
  CMatrix tmpMat ( dim_beta, dim_beta );
  
  for ( index = 0; index < number_of_groups; index++ )
  {
  	tmpVec = beta[ index ]->lastDraw().getVector();
  	#ifdef DEBUGALPHA
  	  printf( "Last beta[%d] sample = \n", index ); tmpVec.Print();
  	#endif
  	if ( !invWishart_betaCov ) {
      for( j = 0; j < dim_beta; j++ )
        tmpVec.Val(j) = tmpVec.Val(j) / tau2[ j ]->lastDraw().getScalar();
    }
    #ifdef DEBUGALPHA
  	  printf( "After dividing each element by tau2[%d] = \n", index ); tmpVec.Print();
  	#endif
  	// In the following line, for invWishart_betaCov = F it simplifies to: beta_predictors[index]->T() * tmpVec
    ( *(zTb[ index ]) ) = beta_predictors[ index ]->T() * ( beta[ index ]->initialInverseCovariance() * tmpVec );

    zTvIb->add( ( *(zTb[ index ]) ) );
  }//end for loop

//  if ( invWishart_betaCov || t_beta )
//  {
  zTvIz->setToZero();
  for ( index = 0; index < number_of_groups; index++ )
  {
    	// In the following line, for invWishart_betaCov = F tmpMat is just the identity matrix.
    	tmpMat = beta[ index ]->initialInverseCovariance();
    	if ( !invWishart_betaCov ) {
    		// Divide the diagonal elements of tmpMat by the corresponding tau2 parameter.
        for( j = 0; j < dim_beta; j++ )
          tmpMat.Val(j, j) = tmpMat.Val(j, j) / tau2[ j ]->lastDraw().getScalar();
      }
  	  #ifdef DEBUGALPHA
  	    printf( "The diagonal matrix with elements 1/tau2[j], or the matrix tau2^-1 = \n" ); tmpMat.Print();
  	    printf( "beta predictors = \n" ); beta_predictors[index]->Print();
  	  #endif
      ( *(zTz[ index ]) ) = beta_predictors[ index ]->T() * ( tmpMat * (*beta_predictors[ index ] ) );

      // Sum over the groups
      zTvIz->add( ( *(zTz[ index ]) ) );
  }
//  }
  
}//end


void BayesianHierarchicalLinearModel::gibbsUpdateResiduals()
{
  int i;

  for ( i = 0; i < number_of_groups; i++ )
  {
    gibbsUpdateResiduals( i );
  }

}//end

void BayesianHierarchicalLinearModel::gibbsUpdateResiduals( int i )
{

  if ( random_effects && fixed_effects )
  {
    ( *(residuals[i]) ) = ( *(response[i]) ) - ( ( ( *(random_predictors[i]) ) * beta[i]->lastDraw().getVector() ) + ( ( *(fixed_predictors[i]) ) * gamma->lastDraw().getVector() ) );

    //printf("BHLM: residuals:[%d] =\n" );
    //residuals[i]->Print();

    //printf( "response[%d] is\n");
    //response[i]->Print();
    //printf( "random predictors[%d] =\n",i);
    //random_predictors[i]->Print();
    //printf(" beta is \n");
    //beta[i]->lastDraw().getVector().Print();
    //printf(" fixed predictors[%d] = \n", i);
    //fixed_predictors[i]->Print();
    //printf("gamma is \n");
    //gamma->lastDraw().getVector().Print();
  }
  else if ( random_effects )
  {
    ( *(residuals[i]) ) = ( *(response[i]) ) - ( ( *(random_predictors[i]) ) * beta[i]->lastDraw().getVector() );
  }
  else if ( fixed_effects )
  {
    ( *(residuals[i]) ) = ( *(response[i]) ) - ( ( *(fixed_predictors[i]) ) * gamma->lastDraw().getVector() );
  }
}//end


void BayesianHierarchicalLinearModel::gibbsUpdateSecondStageResiduals()
{
  int i;

  for ( i = 0; i < number_of_groups; i++ )
  {
    gibbsUpdateSecondStageResiduals( i );
  }

}//end


void BayesianHierarchicalLinearModel::gibbsUpdateSecondStageResiduals( int index )
{
  if ( second_stage_effects )
  {
    ( *(second_residuals[ index ]) ) = beta[ index ]->lastItemDrawn() - ( ( *(beta_predictors[ index ]) ) * alpha->lastItemDrawn() );

    //printf("secondstage residuals[ %d ] = \n", index );
    //second_residuals[ index ]->Print();

  }
  else
  {
    ( *(second_residuals[ index ]) ) = beta[ index ]->lastItemDrawn() - beta[ index ]->initialMean();
  }
}//end



void BayesianHierarchicalLinearModel::gibbsUpdateBeta( int index )
{
  CMatrix sigma_beta( beta[ index ]->dimension(), beta[ index ]->dimension() );
  CVector mu_beta( beta[ index ]->dimension() );

  //this update must be done first
  gibbsUpdateWorkingMatrices( index, "random_effects" );

  if ( common_sigma )
  {
    mu_beta = ( *(xTy[ index ]) ) * ( 1 / sigma2->lastDraw().getScalar() );
    sigma_beta = ( 1 / sigma2->lastDraw().getScalar() ) * ( *(xTx[ index ]) );

    //printf("gibbsUpdateBeta: sigma2 is %f,  xTy[ %d ] is \n", sigma2->lastDraw().getScalar(), index );
    //xTy[index]->Print();
    //printf("xTx[ %d ] is \n", index );
    //xTx[index]->Print();

  }
  else
  {
    mu_beta = ( *(xTy[ index ]) ) * ( 1 / sigma2->lastDrawFromComponent( index ).getScalar() );
    sigma_beta = ( 1 / sigma2->lastDrawFromComponent( index ).getScalar() ) * ( *(xTx[ index ]) );
  }

  if ( invWishart_betaCov )
  {
      CMatrix tmpmat = tau2[0]->lastDraw().getMatrix();
      beta[ index ]->updateVeryFirstCovariance( tmpmat );
      #ifdef DEBUGBETA
        printf("tau2 covar matrix:\n");
        tmpmat.Print(); 
      #endif
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
    #ifdef DEBUGBETA
      printf("beta predictors * alpha vec:\n");
      initMean->Print(); 
    #endif
    delete initMean;
  }

  if ( invWishart_betaCov )
  {
    #ifdef DEBUGBETA
      printf("beta[%d] prior precision matrix:\n", index );
      beta[ index ]->initialInverseCovariance().Print();
      printf("and prec mat * prior mean vector:\n" );
      beta[ index ]->initialInvCovBeta().Print();  
      printf("and updating precision matrix:\n" );
      sigma_beta.Print();
      printf("and updating mean vector:\n" );
      mu_beta.Print();  fflush(stdout);
    #endif
    beta[ index ]->update( mu_beta, sigma_beta );
  }
  else if ( non_informative_beta )
  {
    beta[ index ]->update( mu_beta, sigma_beta );
  }
  else
  {
    CVector *precision_tau2 = new CVector( dim_beta );
    for( int i=0; i < dim_beta; i++ )
      precision_tau2->Val(i) = 1 / tau2[i]->lastDraw().getScalar();
    beta[ index ]->update( (*precision_tau2), mu_beta, sigma_beta );
    delete( precision_tau2 );
  }

  //printf( "bhlm: updated beta. Initial Cov is \n" ); fflush( stdout);
  //beta[index]->initialCovariance().Print();
}//end


void BayesianHierarchicalLinearModel::gibbsUpdateGamma()
{
  int i;

  CMatrix sigma_gamma( gamma->dimension(), gamma->dimension() );
  CVector mu_gamma( gamma->dimension() );

  //this update must be done first
  for ( i = 0; i < number_of_groups; i++ )
  {
    gibbsUpdateWorkingMatrices( i, "fixed_effects" );
  }

  if ( common_sigma )
  {
    mu_gamma = (*sMTy) * ( 1 / sigma2->lastDraw().getScalar() );
    sigma_gamma = ( 1 / sigma2->lastDraw().getScalar() ) * (*sMTm);
  }
  else
  {
    mu_gamma = (*sMTy);

    //sigma_gamma = (*sMTm);

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

  if ( t_gamma )
  {
    double precision_tau2 = 1.0 / tau2_gamma->lastItemDrawn();
    gamma->update( precision_tau2, mu_gamma, sigma_gamma );
  }
  else
  {
  	#ifdef DEBUGGAMMA
  	  printf( "Updating gamma with mu and sigma = \n" );
  	  mu_gamma.Print();
  	  sigma_gamma.Print();
  	  fflush( stdout );
  	#endif
    gamma->update( mu_gamma, sigma_gamma );
  }

}//end



void BayesianHierarchicalLinearModel::gibbsUpdateAlpha()
{
  CMatrix sigma_alpha( alpha->dimension(), alpha->dimension() );
  CVector mu_alpha( alpha->dimension() );

  gibbsUpdateSecondStageWorkingMatrices();

/*  if ( !invWishart_betaCov )
  {
    mu_alpha = (*zTvIb) * ( 1 / tau2->lastDraw().getScalar() );
    sigma_alpha = ( 1 / tau2->lastDraw().getScalar() ) * (*zTvIz);
  }
  else
  {*/
    mu_alpha = (*zTvIb);
    sigma_alpha = (*zTvIz);
  //}

  //printf( "gibbsUpdateAlpha: alpha current mean and variance are:\n" );
  //alpha->mean().Print();
  //alpha->covariance().Print();
  //printf( "Alpha: shifts means and covariance are:\n" );
  //mu_alpha.Print();
  //sigma_alpha.Print();

  if ( t_alpha )
  {
    double precision_tau2 = 1 / tau2_alpha->lastItemDrawn();
    alpha->update( precision_tau2, mu_alpha, sigma_alpha );
  }
  else
  {
    alpha->update( mu_alpha, sigma_alpha );
  }

  #ifdef DEBUGALPHA
    printf("update alpha: new mean: \n" );
    alpha->mean().Print();
    printf("update alpha: new Cov: \n" );
    alpha->covariance().Print();
  #endif
  

}//end



void BayesianHierarchicalLinearModel::gibbsUpdateSigma2() throw( rtErr )
{
  int i;
  double scale, local_scale, add_df;

  if ( common_sigma )
  {
    if ( invChisq_sigma2 || nonInformativePower_sigma2 || properNonInfoPrior_sigma2 )
    {
      scale = 0;
      for ( i = 0; i < number_of_groups; i++ )
      {    
        if ( t_lkhd || t_group_lkhd )
        {
          local_scale = ( *(residuals[i]) ) * residuals[i]->weighted( ( *(tau2_error_weights[i]) ) );
        }
        else
        {
          local_scale = ( *(residuals[i]) * ( *(residuals[i]) ) );
        }

        scale += local_scale;
      }//end loop

      //printf( "BHL:gibbsUpdateSigma2: RSS = %f\n", scale );

      add_df = number_of_observations;

      if ( invChisq_sigma2 || nonInformativePower_sigma2 )
      {
        if ( nonInformativePower_sigma2 )
        {
        	// Bug fix by Dawn: this was add_df += 2 * sigma2NIP[0];
          add_df += 2 * ( sigma2NIP[0] + 1 );
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
  }
  else
  {
    printf( "BayesianHierarchicalLinearModel::gibbsUpdateSigma2: Trying to update a common variance variable when there are number of groups variances in the model.\n" );
    char the_error[] = "BayesianHierarchicalLinearModel::gibbsUpdateSigma2: Trying to update a common variance variable when there are number of groups variances in the model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}//end



void BayesianHierarchicalLinearModel::gibbsUpdateSigma2( int index ) throw( rtErr )
{
  double scale, add_df;

  if ( !common_sigma )
  {
    if ( invChisq_sigma2 || nonInformativePower_sigma2 || properNonInfoPrior_sigma2 )
    {
      if ( t_lkhd || t_group_lkhd )
      {
        scale = ( *(residuals[ index ]) ) * residuals[ index ]->weighted( ( *(tau2_error_weights[ index ]) ) );
      }
      else
      {
        scale = ( *(residuals[ index ]) * ( *(residuals[ index ]) ) );
      }

      add_df = response[ index ]->Len();

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

        //printf( "Bayes: update sigma2: df = %f, scale = %f\n", add_df, scale );
        //printf("residual[%d] is \n", index );
        //residuals[index]->Print();

        sigma2ICS[ index ]->update( add_df, scale );
      }
      else 
      {
        sigma2PNIP[ index ]->setScale( scale );
      }
    }//end if
  }
  else
  {
    printf( "BayesianHierarchicalLinearModel::gibbsUpdateSigma2: Trying to update number of groups variances when there is a common variance in the model.\n" );
    char the_error[] = "BayesianHierarchicalLinearModel::gibbsUpdateSigma2: Trying to update number of groups variances when there is a common variance in the model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}//end 


void BayesianHierarchicalLinearModel::gibbsUpdateTau2()
{
  int i,j;
  double scale, local_scale, add_df;
  
  if ( invChisq_tau2 || nonInformativePower_tau2 || properNonInfoPrior_tau2 )
  {
  	for ( j = 0; j < dim_beta; j++ ){
  		scale = 0;
      for ( i = 0; i < number_of_groups; i++ )
      {    
        local_scale = SQR( second_residuals[i]->Val(j) );

        if ( t_beta )
        {
          local_scale /= tau2_betas[i]->lastItemDrawn();
        }

        scale += local_scale;
      }//end i loop

      //printf("gibbsUpdateTau2: scale = %f\n", scale );

      add_df = number_of_groups;
      if ( invChisq_tau2 || nonInformativePower_tau2 )
      {
        if ( nonInformativePower_tau2 )
        {
          add_df += 2 * tau2NIP[j];
          if ( scale == 0.0 )
          {
            scale = tau2ICS[j]->lastDraw().getScalar();
          }
        }

        tau2ICS[j]->update( add_df, scale );
      }
      else
      {
        tau2PNIP[j]->setScale( scale );
      }
    }// end j loop 
  } else if ( invWishart_betaCov ) {
    CMatrix scale_matrix( beta[0]->dimension(), beta[0]->dimension() );
    CMatrix resid( 1, beta[0]->dimension() );
    CVector diff( beta[0]->dimension() );

    scale_matrix.setToZero();
    for ( i = 0; i < number_of_groups; i++ )
    {
      diff = (*(second_residuals[i]));

      resid.setRow( 0, diff );

      #ifdef DEBUGTAU2
        printf( "bhlm: resid is \n");
        resid.Print();
        printf("bhlm: update: resid.xTx is \n" );
        resid.xTransposedX().Print();
      #endif

      if ( t_beta )
      {
        #ifdef FIX1
          CMatrix tmpmat = ( 1.0 / tau2_betas[i]->lastItemDrawn() ) * resid.xTransposedX();
          scale_matrix.add( tmpmat );
        #else
          scale_matrix.add( ( 1.0 / tau2_betas[i]->lastItemDrawn() ) * resid.xTransposedX() );
        #endif
      } else {
        CMatrix tmpmat = resid.xTransposedX();
        scale_matrix.add( tmpmat );
      }

      //printf( "bhlm:update tau2: in loop[%d]: scale matrix is \n", i ); fflush(stdout);
      //scale_matrix.Print(); fflush(stdout);
    }//end loop

    add_df = number_of_groups;
    #ifdef DEBUGTAU2
      printf( "bhlm:update tau2: add_df = %f, scale matrix is \n", add_df ); fflush(stdout);
      scale_matrix.Print(); fflush(stdout);
      printf( "bhlm:update tau2: prior_df = %f, prior precision matrix is \n", tau2InvWishart->initialDegreesOfFreedom() ); fflush(stdout);
      tau2InvWishart->inverseInitialScaleMatrix().Print();
    #endif

    tau2InvWishart->update( add_df, scale_matrix );

    //printf("bhlm: back from update\n"); fflush(stdout);
    //printf( "bhlm:update tau2: now scale matrix is \n" ); fflush(stdout);
    //scale_matrix.Print(); fflush(stdout);   
  }

  //printf("bhlm: out of update Tau2\n"); fflush(stdout);
}//end



void BayesianHierarchicalLinearModel::gibbsUpdateTau2Gamma()
{
  double scale, add_df;
  CVector gamma_diff( gamma->dimension() );

  gamma_diff = gamma->lastItemDrawn() - gamma->initialMean();
  scale = M_D( gamma->initialInverseCovariance(), gamma_diff );

  add_df = gamma->dimension();

  tau2_gamma->update( add_df, scale );

}//end



void BayesianHierarchicalLinearModel::gibbsUpdateTau2Beta( int index )
{
  double scale, add_df;

  if ( !invWishart_betaCov )
  {
    scale = M_D( beta[ index ]->initialInverseCovariance(), ( *(second_residuals[ index ]) ) );
    scale /= tau2[0]->lastDraw().getScalar();
  }
  else
  {
    scale = M_D( tau2[0]->lastDraw().getMatrix().inverse(), (*(second_residuals[ index ])) );
  }

  add_df = beta[ index ]->dimension();

  tau2_betas[ index ]->update( add_df, scale );

}//end


void BayesianHierarchicalLinearModel::gibbsUpdateTau2Alpha()
{
  double scale, add_df;
  CVector alpha_diff( alpha->dimension() );

  alpha_diff = alpha->lastItemDrawn() - alpha->initialMean();
  scale = M_D( alpha->initialInverseCovariance(), alpha_diff );

  add_df = alpha->dimension();

  tau2_alpha->update( add_df, scale );

}//end


void BayesianHierarchicalLinearModel::gibbsUpdateTau2Group( int index )
{
  double scale, add_df;

  scale = ( *(residuals[ index ]) ) * ( *(residuals[ index ]) );
  if ( common_sigma )
  {
    scale /= sigma2->lastDraw().getScalar();
  }
  else
  {
    scale /= sigma2->lastDrawFromComponent( index ).getScalar();
  }
 
  add_df = response[ index ]->Len();

  tau2_groups[ index ]->update( add_df, scale );

}//end


void BayesianHierarchicalLinearModel::gibbsUpdateTau2Error( int index )
{
  int group, obs;
  double scale, add_df;

  DistributionParameter grp( groupIndex( index ) );

  group = (int) grp.getVector().Val(0);
  obs = (int) grp.getVector().Val(1);

  scale = residuals[ group ]->Val( obs ) * residuals[ group ]->Val( obs );
  if ( common_sigma )
  {
    scale /= sigma2->lastDraw().getScalar();
  }
  else
  {
    scale /= sigma2->lastDrawFromComponent( group ).getScalar();
  }
 
  add_df = 1.0;

  tau2_errors[ index ]->update( add_df, scale );

}//end



void BayesianHierarchicalLinearModel::drawMissingResponse( int index )
{
  int group, obs;

  DistributionParameter grp( groupIndex( index ) );

  group = (int) grp.getVector().Val(0);
  obs = (int) grp.getVector().Val(1);

  //printf("drawing missing response [ %d ][ %d ]\n", group, obs ); fflush(stdout);
  //printf("missing components are:\n");
  //missing_response_components[ group ][ obs ]->Print();

  CVector mu_imputed( 1 );
  CMatrix Sigma_imputed( 1, 1 );

  //just the mean
  mu_imputed.Val( 0 ) = response[ group ]->Val( obs ) - residuals[ group ]->Val( obs );

  if ( common_sigma )
  {
    Sigma_imputed.Val( 0, 0 ) = sigma2->lastDraw().getScalar();
  }
  else
  {
    Sigma_imputed.Val( 0, 0 ) = sigma2->lastDrawFromComponent( group ).getScalar();
  }

  NormalDistribution y_imputed( 1 );
  y_imputed.setMean( &mu_imputed );
  y_imputed.setCovariance( &Sigma_imputed );
  y_imputed.draw();

  //printf(" y mean and variance are %f and %f\n", mu_imputed.Val(0), Sigma_imputed.Val(0,0) );
  //printf("current response[ %d, %d ] is %f\n", group, obs, response[ group ]->Val( obs ) );

  response[ group ]->Val( obs ) = y_imputed.lastItemDrawn().Val( 0 );  


  //printf("y imputed is \n");
  //y_imputed.lastItemDrawn().Print();
  
  //printf("new response[ %d, %d ] is %f\n", group, obs, response[ group ]->Val( obs ) );

  gibbsUpdateResiduals( group );

}//end



void BayesianHierarchicalLinearModel::drawMissingRandomPredictors( int index )
{
  int i, k, start_index, iVar, group, obs, nmss, p_dim;
  double Sigma_x, mu1, mu2;

  DistributionParameter grp( groupIndex( index ) );

  group = (int) grp.getVector().Val(0);
  obs = (int) grp.getVector().Val(1);

  nmss = (int) number_missing_random_predictors->Val( index );
  p_dim = random_predictors[ group ]->Col();

  //printf("drawing missing random predictors [ %d ][ %d ]\n", group, obs ); fflush(stdout);
  //printf("missing components are:\n");
  //missing_random_predictors_components[ group ][ obs ]->Print();

  //printf(" p_dim = %d, nmss = %d\n", p_dim, nmss ); fflush(stdout);

  CVector mu_imputed( 1 );
  CMatrix Sigma_imputed( 1, 1 );

  //compute mean and variance of x (assuming normality)
  CVector beta_coef( beta[ group ]->lastItemDrawn() );
  CVector x_pred( p_dim );

  //get the covariates
  x_pred = random_predictors[ group ]->getRow( obs );

  //printf( " x _pred is \n" );
  //x_pred.Print();
  //printf( "beta_coef is \n");
  //beta_coef.Print();
  //fflush(stdout);


  if ( !beta_coef.isZero( DELTA ) )
  {
    if ( common_sigma )
    {
      Sigma_x = sigma2->lastDraw().getScalar();
    }
    else
    {
      Sigma_x = sigma2->lastDrawFromComponent( group ).getScalar();
    }

    //ready for imputation
    for ( i = 0; i < nmss; i++ )
    {
      //choose i at random
      start_index = randomStart( nmss );
      iVar = ( start_index + i ) %  nmss;
      k = (int) missing_random_predictors_components[ group ][ obs ]->Val( iVar );

      if ( fabs( beta_coef.Val( k ) ) > DELTA &&  randomEffects_vars[ group ]->Val( k ) > DELTA )
      {
        mu1 = ( residuals[ group ]->Val( obs ) + x_pred.Val( k ) * beta_coef.Val( k ) ) * beta_coef.Val( k ) / Sigma_x;
        mu2 = randomEffects_means[ group ]->Val( k ) / randomEffects_vars[ group ]->Val( k );

        Sigma_imputed.Val( 0, 0 ) =  1 / ( ( 1 / randomEffects_vars[ group ]->Val( k ) ) + ( ( beta_coef.Val( k ) * beta_coef.Val( k ) ) / Sigma_x ) ); 

        mu_imputed.Val( 0 ) = Sigma_imputed.Val( 0, 0 ) * ( mu1 + mu2 );

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
    random_predictors[ group ]->setRow( obs , x_pred );

    //printf("new random_predictors[ %d ][ %d ] is \n", group, obs );
    //random_predictors[ group ]->getRow( obs ).Print();

    gibbsUpdateResiduals( group );

    //gibbsUpdateWorkingMatrices( group, "random_effects" );
    //gibbsUpdateWorkingMatrices( group, "fixed_effects" );
  }//end if beta is not zero

}//end




void BayesianHierarchicalLinearModel::drawMissingFixedPredictors( int index )
{
  int i, k, start_index, iVar, group, obs, nmss, p_dim;
  double Sigma_x, mu1, mu2;

  DistributionParameter grp( groupIndex( index ) );

  group = (int) grp.getVector().Val(0);
  obs = (int) grp.getVector().Val(1);

  nmss = (int) number_missing_fixed_predictors->Val( index );
  p_dim = fixed_predictors[ group ]->Col();

  //printf("drawing missing fixed predictors [ %d ][ %d ]\n", group, obs ); fflush(stdout);
  //printf("missing components are:\n");
  //missing_fixed_predictors_components[ group ][ obs ]->Print();

  //printf(" p_dim = %d, nmss = %d\n", p_dim, nmss ); fflush(stdout);

  CVector mu_imputed( 1 );
  CMatrix Sigma_imputed( 1, 1 );

  //compute mean and variance of x (assuming normality)
  CVector gamma_coef( gamma->lastItemDrawn() );
  CVector x_pred( p_dim );

  //get the covariates
  x_pred = fixed_predictors[ group ]->getRow( obs );

  //printf( " x _pred is \n" );
  //x_pred.Print();
  //printf( "gamma_coef is \n");
  //gamma_coef.Print();
  //fflush(stdout);


  if ( !gamma_coef.isZero( DELTA ) )
  {
    if ( common_sigma )
    {
      Sigma_x = sigma2->lastDraw().getScalar();
    }
    else
    {
      Sigma_x = sigma2->lastDrawFromComponent( group ).getScalar();
    }

    //ready for imputation
    for ( i = 0; i < nmss; i++ )
    {
      //choose i at random
      start_index = randomStart( nmss );
      iVar = ( start_index + i ) %  nmss;
      k = (int) missing_fixed_predictors_components[ group ][ obs ]->Val( iVar );

      if ( fabs( gamma_coef.Val( k ) ) > DELTA &&  fixedEffects_vars[ group ]->Val( k ) > DELTA )
      {
        mu1 = ( residuals[ group ]->Val( obs ) + x_pred.Val( k ) * gamma_coef.Val( k ) ) * gamma_coef.Val( k ) / Sigma_x;
        mu2 = fixedEffects_means[ group ]->Val( k ) / fixedEffects_vars[ group ]->Val( k );

        Sigma_imputed.Val( 0, 0 ) =  1 / ( ( 1 / fixedEffects_vars[ group ]->Val( k ) ) + ( ( gamma_coef.Val( k ) * gamma_coef.Val( k ) ) / Sigma_x ) ); 

        mu_imputed.Val( 0 ) = Sigma_imputed.Val( 0, 0 ) * ( mu1 + mu2 );

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
    fixed_predictors[ group ]->setRow( obs , x_pred );

    //printf("new fixed_predictors[ %d ][ %d ] is \n", group, obs );
    //fixed_predictors[ group ]->getRow( obs ).Print();

    gibbsUpdateResiduals( group );

    //gibbsUpdateWorkingMatrices( group, "random_effects" );
    //gibbsUpdateWorkingMatrices( group, "fixed_effects" );
  }//end if gamma is not zero

}//end



void BayesianHierarchicalLinearModel::fullConditionalUpdateVariable( int index )
{
	#ifdef DEBUG1
    printf( "BHL:full conditional update for variable %s.\n  ", distr_map[ index ] ); fflush( stdout );
  #endif

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
            if ( group_index >= 0 && group_index < number_of_observations )
	    {
              gibbsUpdateTau2Error( group_index );
            }
            else
	    {
              printf( " BayesianHierarchicalLinearModel::fullConditionalUpdateVariable: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
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
              printf( " BayesianHierarchicalLinearModel::fullConditionalUpdateVariable: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
            }
          }
        }//end if not null
        else
        {
          printf( " BayesianHierarchicalLinearModel::fullConditionalUpdateVariable: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }
      else
      {
        printf( " BayesianHierarchicalLinearModel::fullConditionalUpdateVariable: Wrong argument in [%s]. Number expected.\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if 
  }//end if index is valid
  else
  {
    printf( "BayesianHierarchicalLinearModel::fullConditionalUpdateVariable: Variable index [%d] does not exist.\n", index );
  }

  //printf("bhlm: update full conditinal: out\n\n"); fflush( stdout );
}//end



void BayesianHierarchicalLinearModel::drawVariable( int index )
{
	int i;
  #ifdef DEBUG2
    printf( "bhlm: will draw variable %d: [ %s ]\n", index, distr_map[ index ] ); fflush(stdout);
  #endif

  if ( index < number_of_variables && index >= 0 )
  {
    if ( !strcmp( distr_map[ index ], "gamma" ) )
    {
      gamma->draw();

    	#ifdef DEBUGGAMMA
    	  printf("drew gamma with mean and covar = \n"); 
    	  gamma->mean().Print();  gamma->covariance().Print();  
    	  printf("new gamma = \n");  gamma->lastItemDrawn().Print();
    	  fflush( stdout );
    	#endif

      //printf( "BHL: gamma draw = \n" );
      //gamma->lastItemDrawn().Print();

      gibbsUpdateResiduals();
    }
    else if ( !strcmp( distr_map[ index ], "alpha" ) )
    {
      alpha->draw();
      gibbsUpdateSecondStageResiduals();
      #ifdef DEBUGALPHA
        printf("drew alpha with mean \n"); alpha->mean().getVector().Print();
        printf("alpha = \n"); alpha->lastItemDrawn().Print();
      #endif
      #ifdef DEBUGTAU2
        printf("drew alpha with mean \n"); alpha->mean().getVector().Print();
        printf("alpha = \n"); alpha->lastItemDrawn().Print();
      #endif
    }
    else if ( !strcmp( distr_map[ index ], "sigma" ) )
    {
      sigma2->draw();

      #ifdef DEBUGSIGMA
        printf("drew sigma with df = %f and scale = %f\n", sigma2ICS[0]->degreesOfFreedom(), sigma2ICS[0]->scale() );
        printf( "BHL: sigma draw = %f\n", sqrt( sigma2->lastDraw().getScalar() ) );
      #endif
    }
    else if ( !strcmp( distr_map[ index ], "tau2" ) )
    {
    	if( invWishart_betaCov ){
    		tau2[0]->draw();
    		#ifdef DEBUGTAU2
        	printf("drew tau2 matrix with df = %f\n", tau2InvWishart->degreesOfFreedom() );
        	printf("sample[0][0]: %f \n", tau2[0]->lastDraw().getScalar() );
        #endif
    		#ifdef DEBUGALPHA
        	printf("drew tau2 matrix = \n" );
          tau2[0]->lastDraw().getMatrix().Print();
        #endif
    	} else {
        for( i=0; i<dim_beta; i++ ){
          tau2[i]->draw();
          #ifdef DEBUGTAU2
            if( invChisq_tau2 ){
          	  printf("drew tau2[%d] with df = %f and scale = %f\n", i, tau2ICS[i]->degreesOfFreedom(), tau2ICS[i]->scale() );
        	    printf("sample: %f \n", sqrt( tau2[i]->lastDraw().getScalar() ) );
        	  }
        	#endif
        }
      }
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
              printf( " BayesianHierarchicalLinearModel::drawVariable: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
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
              printf( " BayesianHierarchicalLinearModel::drawVariable: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
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
              printf( " BayesianHierarchicalLinearModel::drawVariable: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
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
              printf( " BayesianHierarchicalLinearModel::drawVariable: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
            }
          }
          else
	  {
            if ( group_index >= 0 && group_index < number_of_groups )
	    {
              if ( !strncmp( distr_map[ index ], "beta", 4 ) )
  	          {
                beta[ group_index ]->draw();
                gibbsUpdateSecondStageResiduals( group_index );
                gibbsUpdateResiduals( group_index );
                #ifdef DEBUGTAU2
                  printf("drew beta[%d] with mean \n", group_index); beta[ group_index ]->mean().getVector().Print();
                  printf("beta = \n"); beta[ group_index ]->lastItemDrawn().Print();
                #endif
              }
              else if ( !strncmp( distr_map[ index ], "sigma", 5 ) )
	      {
                sigma2->drawFromComponent( group_index );

                //printf( "Bayes: sigma2 draw[%d] = %f\n", group_index, sigma2->lastDraw().getScalar() );
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
              printf( " BayesianHierarchicalLinearModel::drawVariable: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
            }
          }
        }//end if not null
        else
        {
          printf( " BayesianHierarchicalLinearModel::drawVariable: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }
      else
      {
        printf( " BayesianHierarchicalLinearModel::drawVariable: Wrong argument in [%s]. Number expected.\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if 
  }//end if index is valid
  else
  {
    printf( "BayesianHierarchicalLinearModel::drawVariable: Variable index [%d] does not exist.\n", index );
  }

}//end



void BayesianHierarchicalLinearModel::dataAugmentationInitialDraws()
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



void BayesianHierarchicalLinearModel::createOutput( int simulations_to_keep )
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
        simulated_beta[i] = new CMatrix( simulations_to_keep, random_predictors[0]->Col() );
      }

      if ( t_beta && keep_tau2_error )
      {
        simulated_tau2_beta = new CMatrix( simulations_to_keep, number_of_groups );
      }


      if ( !non_informative_beta )
      {
        simulated_tau2 = new CMatrix * [ simulations_to_keep ];
      }

      if ( second_stage_effects )
      {
        simulated_alpha = new CMatrix( simulations_to_keep, beta_predictors[0]->Col() );
        if ( t_alpha && keep_tau2_error )
        {
          simulated_tau2_alpha = new CVector( simulations_to_keep );
        }
      }
    }

    if ( fixed_effects )
    {
      simulated_gamma = new CMatrix( simulations_to_keep, fixed_predictors[0]->Col() );
      if ( t_gamma && keep_tau2_error )
      {
        simulated_tau2_gamma = new CVector( simulations_to_keep );
      }
    }

    if ( !known_sigma2 )
    {
      if ( common_sigma )
      {
        simulated_sigma2 = new CMatrix( simulations_to_keep, 1 );
      }
      else
      {
        simulated_sigma2 = new CMatrix( simulations_to_keep, number_of_groups );
      }
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
            //printf( "number_missing_response->Val( %d ) = %d\n", i, (int) number_missing_response->Val( i ) ); fflush(stdout);

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



void BayesianHierarchicalLinearModel::keepSimulation( int simul_number )
{
  int i, group, obs;

  if ( simul_number < number_of_simulations )
  {
    //printf( "** simul number is %d **\n", simul_number );

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


      if ( !non_informative_beta )
      {
        //printf("tau2 is %f\n", tau2->lastDraw().getScalar() ); fflush(stdout);

        if( invWishart_betaCov ){
          simulated_tau2[ simul_number ] = new CMatrix( tau2[0]->lastDraw().getMatrix() );
        } else {
        	simulated_tau2[ simul_number ] = new CMatrix( 1, dim_beta );
        	for( i=0; i<dim_beta; i++ ){
        		simulated_tau2[ simul_number ]->Val( 0, i ) = tau2[i]->lastDraw().getScalar();
        	}
        }
      }


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
        simulated_sigma2->Val( simul_number, 0 ) = sigma2->lastDraw().getScalar();
      }
      else
      {
        for ( i = 0; i < number_of_groups; i++ )
        {
          simulated_sigma2->Val( simul_number, i ) = sigma2->lastDrawFromComponent( i ).getScalar();
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

            simulated_missingR[i]->Val( simul_number, 0 ) = response[ group ]->Val( obs );
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
            CVector tmpvec = random_predictors[ group ]->getRow( obs ).subVector( (*(missing_random_predictors_components[ group ][ obs ])) );
            simulated_missingRP[i]->setRow( simul_number, tmpvec );
#else
            simulated_missingRP[i]->setRow( simul_number, random_predictors[ group ]->getRow( obs ).subVector( (*(missing_random_predictors_components[ group ][ obs ])) ) );
#endif
//            simulated_missingRP[i]->setRow( simul_number, random_predictors[ group ]->getRow( obs ).subVector( (*(missing_random_predictors_components[ group ][ obs ])) ) );
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
            CVector tmpvec = fixed_predictors[ group ]->getRow( obs ).subVector( (*(missing_fixed_predictors_components[ group ][ obs ])) );
            simulated_missingFP[i]->setRow( simul_number, tmpvec );
#else
            simulated_missingFP[i]->setRow( simul_number, fixed_predictors[ group ]->getRow( obs ).subVector( (*(missing_fixed_predictors_components[ group ][ obs ])) ) );
#endif
//            simulated_missingFP[i]->setRow( simul_number, fixed_predictors[ group ]->getRow( obs ).subVector( (*(missing_fixed_predictors_components[ group ][ obs ])) ) );
          }
        }//end for i
      }//end fixed predictors

    }//end missing data

  }//end if keep any simulations

}//end
    


void BayesianHierarchicalLinearModel::simulationsToArray( double * simul_output, int simulations_to_keep )
{
  int i, j, k, total_dim, start_dim, g, v, s;

  start_dim = 0;
  if ( random_effects )
  {
    total_dim = simulations_to_keep * random_predictors[0]->Col();
    for ( g = 0; g < number_of_groups; g++ )
    {
      for ( v = 0; v < random_predictors[0]->Col(); v++ )
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
    for ( v = 0; v < fixed_predictors[0]->Col(); v++ )
    {
      for ( s = 0; s < simulations_to_keep; s++ )
      {
        simul_output[ start_dim + v * total_dim + s ] = simulated_gamma->Val(s,v);
      }
    }
    start_dim += fixed_predictors[0]->Col() * total_dim;
  }

  if ( second_stage_effects )
  {
    total_dim = simulations_to_keep;
    for ( v = 0; v < beta_predictors[0]->Col(); v++ )
    {
      for ( s = 0; s < simulations_to_keep; s++ )
      {
        simul_output[ start_dim + v * total_dim + s ] = simulated_alpha->Val(s,v);
      }
    }
    start_dim += beta_predictors[0]->Col() * total_dim;
  }


  if ( !known_sigma2 )
  {
    if ( common_sigma )
    {
      for ( s = 0; s < simulations_to_keep; s++ )
      {
        simul_output[ start_dim + s ] = sqrt( simulated_sigma2->Val( s, 0 ) );
      }
      start_dim += simulations_to_keep;
    }
    else
    {
      total_dim = simulations_to_keep;
      for ( i = 0; i < number_of_groups; i++ )
      {
        for ( s = 0; s < simulations_to_keep; s++ )
        {
          simul_output[ start_dim + i * total_dim + s ] =  sqrt( simulated_sigma2->Val( s, i ) );
        }
      }
      start_dim += number_of_groups * total_dim;
    }
  }  


  if ( random_effects && !non_informative_beta )
  {
    if( invWishart_betaCov ){
      int q = simulated_tau2[0]->Col();
      total_dim = simulations_to_keep * q;
      for ( i = 0; i < q; i++ )
      {
        for ( j = 0; j < q; j++ )
        {   
          for ( s = 0; s < simulations_to_keep; s++ )
          {
            simul_output[ start_dim + i * total_dim + j * simulations_to_keep + s ] = simulated_tau2[s]->Val(i,j);
          }
        }
      }
      start_dim += q * total_dim;
    } else {
    	int q = simulated_tau2[0]->Col();
      for ( i = 0; i < q; i++ )
      {
        for ( s = 0; s < simulations_to_keep; s++ )
        {
          simul_output[ start_dim + i * simulations_to_keep + s ] = sqrt( simulated_tau2[s]->Val(0,i) );
        }
      }
      start_dim += simulations_to_keep * q;
    }
    
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
        simul_output[ start_dim + g * total_dim + s ] =  sqrt( simulated_tau2_groups->Val(s,g) );
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
        simul_output[ start_dim + i * total_dim + s ] = sqrt( simulated_tau2_errors->Val(s,i) );
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




  //printf(" simul to array: start_index is %d\n", start_dim);

}//end
