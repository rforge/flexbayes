#include <math.h>
#include <string.h>
#include "R.h"
#include "rtErr.h"

#include "Const.h"
#include "BayesianGlmModel.h"
#include "DistributionParameter.h"
#include "BayesianHierarchicalGlmModel.h"

#define SQR(x) ((x)*(x))


/* constructor.
   returns: an empty Distribution of type BayesianHierarchicalGlmModel
*/
BayesianHierarchicalGlmModel::BayesianHierarchicalGlmModel()
{
  emptyModel();
}//end


void BayesianHierarchicalGlmModel::emptyModel()
{
  random_predictors = NULL;
  fixed_predictors = NULL;
  beta_predictors = NULL;

  //distributions
  //first stage
  gamma = NULL;
  tau2_gamma = NULL;

  beta = NULL;
  dim_beta = 0 ;
  //t prior
  tau2_betas = NULL;    

  //second stage
  alpha = NULL;
  tau2 = NULL;
  tau2ICS = NULL;
  tau2PNIP = NULL;
  tau2NIP = NULL;
  tau2InvWishart = NULL;
  tau2_alpha = NULL;

  tau2_first_draw = NULL;

  number_of_variables = 0;
  number_of_groups = 0;
  number_of_observations = 0;

  t_beta = false;
  t_gamma = false;
  t_alpha = false;
  non_informative_alpha = false;
  non_informative_beta = false;
  non_informative_gamma = false;
  invWishart_betaCov = false;
  keep_tau2_error = true;
  random_effects = false;
  fixed_effects = false;
  second_stage_effects = false;

  invChisq_tau2 = false;
  duMouchel_tau2 = false;
  uniformShrinkage_tau2 = false;
  properNonInfoPrior_tau2 = false;
  nonInformativePower_tau2 = false;

  glm_model = NULL;
  update_for_hessians = true;

  //simulation samples
  number_of_simulations = 0;
  simulated_beta = NULL;
  simulated_gamma = NULL;
  simulated_alpha = NULL;
  simulated_tau2 = NULL;
  simulated_tau2_beta = NULL;
  simulated_tau2_gamma = NULL;
  simulated_tau2_alpha = NULL;

  //working matrices for Gibbs sampler
  //second stage
  zTz = NULL;
  zTb = NULL;
  zTvIz = NULL;
  zTvIb = NULL;

  second_residuals = NULL;

  //array of distributions map (maps distributions names to actual objects)
  distr_map = NULL;

}//end



BayesianHierarchicalGlmModel::~BayesianHierarchicalGlmModel()
{
  int i;

  if ( random_effects )
  {
    for ( i = 0; i < number_of_groups + 2; i++ )
    {
      delete random_predictors[i];
    }
    delete [] random_predictors;
  }

  if ( fixed_effects )
  {

    for ( i = 0; i < number_of_groups + 2; i++ )
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

  if ( tau2NIP != NULL )
  {
    delete [] tau2NIP;
  }

  if ( tau2InvWishart != NULL )
  {
    delete tau2InvWishart;
  }

  if ( tau2_alpha != NULL )
  {
    delete tau2_alpha;
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

  if ( glm_model != NULL )
  {
    delete [] glm_model;
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

  //working matrices for Gibbs sampler
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

}//end



void BayesianHierarchicalGlmModel::initialize( CVector ** y, int n_groups )
{
  int i;

  number_of_groups = n_groups;
  number_of_observations = 0;

  for ( i = 0; i < number_of_groups; i++ )
  {
    number_of_observations += y[i]->Len();
  }
}//end


void BayesianHierarchicalGlmModel::randomEffects( CMatrix ** x )
{
  int i;

  random_predictors = new DistributionParameter * [ number_of_groups + 2 ];
  //indicate random effects
  random_predictors[0] = new DistributionParameter( 0 );
  //points to first group
  random_predictors[1] = new DistributionParameter( 0 );

  for ( i = 0; i < number_of_groups; i++ )
  {
    random_predictors[ i + 2 ] = new DistributionParameter( (*x[i]) );
  }

  number_of_variables += number_of_groups;

  random_effects = true;  
  dim_beta = random_predictors[2]->getMatrix().Col();

}//end


void BayesianHierarchicalGlmModel::fixedEffects( CMatrix ** m )
{
  int i;

  fixed_predictors = new DistributionParameter * [ number_of_groups + 2 ];
  //indicate fixed effects
  fixed_predictors[0] = new DistributionParameter( 1 );
  //points to first group
  fixed_predictors[1] = new DistributionParameter( 0 );

  for ( i = 0; i < number_of_groups; i++ )
  {
    fixed_predictors[ i + 2 ] = new DistributionParameter( (*m[i]) );

    //printf("fixed predictors[%d] = \n", i );
    //fixed_predictors[ i + 2 ]->Print();

  }
  //count gamma
  number_of_variables++;

  fixed_effects = true;  

}//end


void BayesianHierarchicalGlmModel::secondStageRandomEffects( CMatrix ** z )
{
  int i;

  beta_predictors = new CMatrix * [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    beta_predictors[i] = new CMatrix( (*z[i]) );
  }
  //count alpha
  number_of_variables++;

  second_stage_effects = true;
}//end



//attach the glm variables and model
void BayesianHierarchicalGlmModel::glmModel( BayesianGlmModel * glm_vars )
{
  glm_model = new BayesianGlmModel * [1];
  glm_model[0] = glm_vars;
  number_of_variables += glm_model[0]->numberOfVariables();

  //initialize counts based on response
  int gamma_dim;

  gamma_dim = 0;

  if ( fixed_effects )
  {
    gamma_dim = fixed_predictors[2]->getMatrix().Col();
  }

  glm_model[0]->setCoefficientDimension( dim_beta, gamma_dim );
  glm_model[0]->considerCounts();
  glm_model[0]->setUpdateHessians( update_for_hessians );

}//end 


void BayesianHierarchicalGlmModel::setUpdateProposalCovariance( bool value )
{
  update_for_hessians = value;
  if ( glm_model != NULL )
  {
    glm_model[0]->setUpdateHessians( update_for_hessians );
  }

}//end


void BayesianHierarchicalGlmModel::betaPrior( CMatrix * p_betaCov )
{
  int i;

  beta = new NormalDistribution * [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    beta[i] = new NormalDistribution( p_betaCov->Col() );
    beta[i]->setCovariance( p_betaCov );
  }

}//end


void BayesianHierarchicalGlmModel::betaPrior( CVector * p_beta, CMatrix * p_betaCov )
{
  int i;
  beta = new NormalDistribution * [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    beta[i] = new NormalDistribution( p_betaCov->Col() );
    beta[i]->setCovariance( p_betaCov );
    beta[i]->setMean( p_beta );
  }

}//end


void BayesianHierarchicalGlmModel::betaPrior( CMatrix * p_beta, CMatrix * p_betaCov )
{
  int i;

  beta = new NormalDistribution * [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    beta[i] = new NormalDistribution( p_betaCov->Col() );
    beta[i]->setCovariance( p_betaCov );
    CVector mean_val( p_beta->getColumn( i ) );
    beta[i]->setMean( &mean_val );
  }

}//end


void BayesianHierarchicalGlmModel::betaPriorNonInformative( int dim ) throw( rtErr )
{
  if ( !second_stage_effects )
  {
    int i;

    beta = new NormalDistribution * [ number_of_groups ];
    for ( i = 0; i < number_of_groups; i++ )
    {
      beta[i] = new NormalDistribution( dim );
      beta[i]->setNonInformative();
    }
    non_informative_beta = true;
  }
  else
  {
    printf( "BayesianHierarchicalGlmModel::betaPriorNonInformative: prior for beta cannot be non-informative when second stage effects are present in the model.\n" );
    char the_error[] = "BayesianHierarchicalGlmModel::betaPriorNonInformative: prior for beta cannot be non-informative when second stage effects are present in the model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

}//end


void BayesianHierarchicalGlmModel::gammaPriorNonInformative( int dim )
{
  gamma = new NormalDistribution( dim );
  gamma->setNonInformative();
  non_informative_gamma = true;

}//end


void BayesianHierarchicalGlmModel::gammaPrior( CVector * p_gamma, CMatrix * p_gammaCov )    
{
  gamma = new NormalDistribution( p_gamma->Len() );
  gamma->setMean( p_gamma );
  gamma->setCovariance( p_gammaCov );

}//end



void BayesianHierarchicalGlmModel::alphaPrior( CVector * p_alpha, CMatrix * p_alphaCov ) throw( rtErr )   
{
  if ( second_stage_effects )
  {
    alpha = new NormalDistribution( p_alpha->Len() );
    alpha->setMean( p_alpha );
    alpha->setCovariance( p_alphaCov );
  }
  else
  {
    printf( "BayesianHierarchicalGlmModel::alphaPrior: alpha is not a parameter of this model.\n" );
    char the_error[] = "BayesianHierarchicalGlmModel::alphaPrior: alpha is not a parameter of this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

}//end


void BayesianHierarchicalGlmModel::alphaPriorNonInformative( int dim ) throw( rtErr )
{
  if ( second_stage_effects )
  {
    alpha = new NormalDistribution( dim );
    alpha->setNonInformative();
    non_informative_alpha = true;
  }
  else
  {
    printf( "BayesianHierarchicalGlmModel::alphaPriorNonInformative: alpha is not a parameter of this model.\n" );
    char the_error[] = "BayesianHierarchicalGlmModel::alphaPriorNonInformative: alpha is not a parameter of this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

}//end


void BayesianHierarchicalGlmModel::tau2PriorInvChisq( double *p_nuTau2, double *p_tau2 ) throw( rtErr )
{
  if ( random_effects )
  {
  	int i;
  	tau2 = new DistributionMixture * [ dim_beta ];
  	tau2ICS = new InvChisqDistribution * [ dim_beta ];
  	tau2_first_draw = new DistributionParameter * [ dim_beta ];
  	#ifdef DEBUG1
      printf("HGLM: Specifying tau2 prior; dim_beta = %d\n", dim_beta ); fflush(stdout);
    #endif
    for( i = 0; i < dim_beta; i++ ){
    	// for each random effect, specify the prior for the variance
      tau2[i] = new DistributionMixture( 1 );    	
      tau2ICS[i] = new InvChisqDistribution( p_nuTau2[i] );
      tau2ICS[i]->setScale( p_tau2[i] );
      tau2[i]->set( 0, tau2ICS[i], 1.0 );
      tau2_first_draw[i] = new DistributionParameter( p_tau2[i] );
      #ifdef DEBUG1
        printf("Specifying inverse-chisquared prior with nu=%f, scale=%f, for tau2[%d]\n", 
          p_nuTau2[i], p_tau2[i], i ); fflush(stdout);
      #endif
    }
    invChisq_tau2 = true;
    // Only increase the number of variables by 1, since the tau2s will 
    // be updated together.
    number_of_variables++;
  }
  else
  {
    printf( "BayesianHierarchicalGlmModel::tau2PriorInvChisq: There are no random effects in the model. Tau2 does not make sense in this model.\n" );
    char the_error[] = "BayesianHierarchicalGlmModel::tau2PriorInvChisq: There are no random effects in the model. Tau2 does not make sense in this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

}//end


void BayesianHierarchicalGlmModel::tau2PriorDuMouchel( double *p_tau2 ) throw( rtErr )
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
    printf( "BayesianHierarchicalGlmModel::tau2PriorDuMouchel: There are no random effects in the model. Tau2 does not make sense in this model.\n" );
    char the_error[] = "BayesianHierarchicalGlmModel::tau2PriorDuMouchel: There are no random effects in the model. Tau2 does not make sense in this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}//end


void BayesianHierarchicalGlmModel::tau2PriorUniformShrinkage( double *p_tau2 ) throw( rtErr )
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
    printf( "BayesianHierarchicalGlmModel::tau2PriorUniformShrinkage: There are no random effects in the model. Tau2 does not make sense in this model.\n" );
    char the_error[] = "BayesianHierarchicalGlmModel::tau2PriorUniformShrinkage: There are no random effects in the model. Tau2 does not make sense in this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}//end


void BayesianHierarchicalGlmModel::tau2PriorNonInformative( double *p_power ) throw( rtErr )
{

  if ( random_effects )
  {
    //printf("BHLM: p_power = %f, n groups = %d,  p = %d\n", p_power, number_of_groups, random_predictors[2]->getMatrix().Col()); fflush(stdout);
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
    printf( "BayesianHierarchicalGlmModel::tau2PriorNonInformative: There are no random effects in the model. Tau2 does not make sense in this model.\n" );
    char the_error[] = "BayesianHierarchicalGlmModel::tau2PriorNonInformative: There are no random effects in the model. Tau2 does not make sense in this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

}//end



void BayesianHierarchicalGlmModel::betaCovPriorInvWishart( double p_nuV, CMatrix * p_VCov ) throw( rtErr )
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
    printf( "BayesianHierarchicalGlmModel::betaCovPriorInvWishart: There are no random effects in the model. Tau2 does not make sense in this model.\n" );
    char the_error[] = "BayesianHierarchicalGlmModel::betaCovPriorInvWishart: There are no random effects in the model. Tau2 does not make sense in this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

}//end



void BayesianHierarchicalGlmModel::gammaTPrior( double p_nuGamma ) throw( rtErr )
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
    printf( "BayesianHierarchicalGlmModel::gammaTPrior: gamma is not a parameter of this model.\n" );
    char the_error[] = "BayesianHierarchicalGlmModel::gammaTPrior: gamma is not a parameter of this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

}//end



void BayesianHierarchicalGlmModel::betaTPrior( double p_nuBeta ) throw( rtErr )
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
    printf( "BayesianHierarchicalGlmModel::betaTPrior: beta is not a parameter of this model.\n" );
    char the_error[] = "BayesianHierarchicalGlmModel::betaTPrior: beta is not a parameter of this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}//end


void BayesianHierarchicalGlmModel::alphaTPrior( double p_nuAlpha ) throw( rtErr )
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
    printf( "BayesianHierarchicalGlmModel::alphaTPrior: alpha is not a parameter of this model.\n" );
    char the_error[] = "BayesianHierarchicalGlmModel::alphaTPrior: alpha is not a parameter of this model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}//end



void BayesianHierarchicalGlmModel::samplerDefaultInitialPoint()
{
  int i;

  if ( random_effects )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
#ifdef FIX1
      DistributionParameter tmpmean = beta[i]->mean();
      beta[i]->setLastDraw( tmpmean );
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
      DistributionParameter tmpmean = alpha->mean();
      alpha->setLastDraw( tmpmean );
#else
      alpha->setLastDraw( alpha->mean() );
#endif
//      alpha->setLastDraw( alpha->mean() );
    }
  }

  if ( fixed_effects )
  {
#ifdef FIX1
    DistributionParameter tmpmean = gamma->mean();
    gamma->setLastDraw( tmpmean );
#else
    gamma->setLastDraw( gamma->mean() );
#endif
//    gamma->setLastDraw( gamma->mean() );
  }

}//end


void BayesianHierarchicalGlmModel::samplerTau2InitialPoint( double * init_tau2 )
{
  if ( !invWishart_betaCov )
  {
  	int i;
  	tau2_first_draw = new DistributionParameter * [ dim_beta ];
  	for( i = 0; i < dim_beta; i++ ){
      tau2_first_draw[i] = new DistributionParameter( init_tau2[i] );
      tau2[i]->setLastDraw( *(tau2_first_draw[i]) );
    }
    #ifdef DEBUG1
      printf( "HGLM: tau2 initialized to: " );
      for( i = 0; i < dim_beta; i++ ){
        printf( " %f ", init_tau2[i] );
      }
      printf( "\n" );
    #endif
  }
  else
  {
    printf( "BayesianHierarchicalGlmModel::samplerTau2InitialPoint: Initial point should be a matrix, not a scalar.\n" );
  }
}//end


void BayesianHierarchicalGlmModel::samplerTau2InitialPoint( CMatrix & init_tau2 )
{
  if ( invWishart_betaCov )
  {
  	tau2_first_draw = new DistributionParameter * [ 1 ];
    tau2_first_draw[0] = new DistributionParameter( init_tau2 );
    tau2[0]->setLastDraw( *(tau2_first_draw[0]) );
  }
  else
  {
    printf( "BayesianHierarchicalGlmModel::samplerTau2InitialPoint: Initial point should not be a matrix, but a scalar.\n" );
  }
}//end


void BayesianHierarchicalGlmModel::samplerAlphaInitialPoint( CVector & init_alpha )
{
  alpha->setLastDraw( init_alpha );
}//end


void BayesianHierarchicalGlmModel::samplerBetaInitialPoint( CVector & init_beta )
{
  int i;
  
  for ( i = 0; i < number_of_groups; i++ )
  {
    beta[ i ]->setLastDraw( init_beta );
  }
}//end


void BayesianHierarchicalGlmModel::samplerBetaInitialPoint( CMatrix & init_beta )
{
  int i;
  
  for ( i = 0; i < number_of_groups; i++ )
  {
    beta[ i ]->setLastDraw( init_beta.getColumn( i ) );
  }
}//end



void BayesianHierarchicalGlmModel::samplerGammaInitialPoint( CVector & init_gamma )
{
  gamma->setLastDraw( init_gamma );

}//end

// "index" is the group index
void BayesianHierarchicalGlmModel::updateLinearResponse( int index )
{
  int i, dim;

  dim = 1;
  if ( random_effects )
  {
    dim = random_predictors[ index + 2 ]->getMatrix().Row();
  }
  else if ( fixed_effects )
  {
    dim = fixed_predictors[ index + 2 ]->getMatrix().Row();
  }

  CVector mu( dim );

  if ( random_effects && fixed_effects )
  {
  	// These random / fixed predictors matrices have # rows = # obs in this group, # cols = # predictors.
    mu = ( random_predictors[ index + 2 ]->getMatrix() * beta[ index ]->lastDraw().getVector() ) + ( fixed_predictors[ index + 2 ]->getMatrix() * gamma->lastDraw().getVector() );
  }
  else if ( random_effects )
  {
    mu = random_predictors[ index + 2 ]->getMatrix() * beta[ index ]->lastDraw().getVector();
  }
  else
  {
    mu = fixed_predictors[ index + 2 ]->getMatrix() * gamma->lastDraw().getVector();
  }
  // The vector mu is of length = # obs in this group

  //printf("gamma = \n");
  //gamma->lastDraw().Print();

  // The vector mu_plus is equal to the vector mu, except the zero-th element is equal to the group index
  // and the rest of the elements are moved one index up.
  CVector mu_plus( mu.Len() + 1 );
  mu_plus.Val(0) = index;
  for ( i = 0; i < mu.Len(); i++ )
  {
    mu_plus.Val( i + 1 ) = mu.Val( i );
  }
  DistributionParameter mu_par( mu_plus );

  glm_model[0]->metropolisHastingsUpdateModel( &mu_par );

}//end


void BayesianHierarchicalGlmModel::updateLinearResponse()
{
  int i;
  
  for ( i = 0; i < number_of_groups; i++ )
  {
    updateLinearResponse( i );
  }

}//end


void BayesianHierarchicalGlmModel::gibbsUpdateSecondStageResiduals()
{
  int i;

  for ( i = 0; i < number_of_groups; i++ )
  {
    gibbsUpdateSecondStageResiduals( i );
  }

}//end


void BayesianHierarchicalGlmModel::gibbsUpdateSecondStageResiduals( int index )
{
  if ( second_stage_effects )
  {
    ( *(second_residuals[ index ]) ) = beta[ index ]->lastItemDrawn() - ( ( *(beta_predictors[ index ]) ) * alpha->lastItemDrawn() );

  }
  else
  {
    ( *(second_residuals[ index ]) ) = beta[ index ]->lastItemDrawn() - beta[ index ]->veryFirstMean();
  }
}//end


// The argument "index" here is the group index, although the result
// should not actually depend on this index.
CMatrix BayesianHierarchicalGlmModel::betaInverseCovariance( int index )
{
  CMatrix beta_inv_cov( beta[ index ]->dimension(), beta[ index ]->dimension() );

  if ( !non_informative_beta )
  {
    if ( invWishart_betaCov )
    {
      beta_inv_cov = tau2[0]->lastDraw().getMatrix().inverse();
    }
    else
    {
      beta_inv_cov = beta[ index ]->veryFirstInverseCovariance();
      // The result of the above call is the identity matrix.  Divide the 
      // diagonal elements by the associated tau2 last draw values.
      for( int i=0; i<dim_beta; i++ ){
      	// the following code may cause a problem if the Val function does not 
      	// successfully update beta_inv_cov
      	beta_inv_cov.Val( i, i ) = 1.0 / tau2[i]->lastDraw().getScalar();
      }
    }

    if ( t_beta )
    {
      beta_inv_cov.multiplyByScalar( ( 1.0 / tau2_betas[ index ]->lastItemDrawn() ) );
    }
  }
  else
  {
    beta_inv_cov = beta[ index ]->veryFirstInverseCovariance();
  }

  return (beta_inv_cov);
}//end



void BayesianHierarchicalGlmModel::gibbsUpdateBeta( int index )
{
  //Gibbs proposal

  //get shift mean and cov from current values
  //random_predictors[1] indicates the index of the group that is currently being considered.
  //Fix by Dawn: the following was changed from   random_predictors[1]->parameter().Val( 0, 0 ) = index;
  //that code fails to change the value of the parameter, instead changing the value of a copy of the
  //parameter.
  random_predictors[1]->updateParameter( (double) index );
  CMatrix beta_shift_cov( (-1.0) * glm_model[0]->metropolisHastingsHessian( 3, random_predictors ) );

  DistributionParameter ** beta_vectors;
  beta_vectors = new DistributionParameter * [1];
#ifdef FIX1
  CVector tmpvec = beta[ index ]->lastItemDrawn();
  beta_vectors[0] = new DistributionParameter( tmpvec );
#else
  beta_vectors[0] = new DistributionParameter( beta[ index ]->lastItemDrawn() );
#endif
//  beta_vectors[0] = new DistributionParameter( beta[ index ]->lastItemDrawn() );

  CVector beta_shift_mean( glm_model[0]->metropolisHastingsMean( 3, random_predictors, beta_vectors ) );

  if ( !non_informative_beta )
  {
    if ( invWishart_betaCov || t_beta )
    {
      if ( invWishart_betaCov )
      {
#ifdef FIX1
        CMatrix tmpmat = tau2[0]->lastDraw().getMatrix();
        beta[ index ]->updateVeryFirstCovariance( tmpmat );
#else
        beta[ index ]->updateVeryFirstCovariance( tau2[0]->lastDraw().getMatrix() );
#endif

#ifdef DEBUGBETA
        printf("Updating beta[%d]->very_first_cov and beta[%d]->very_first_inv_cov to \n", index, index);
        tmpmat.Print();  fflush(stdout);
#endif
      }

      CMatrix init_cov( beta[ index ]->veryFirstCovariance() );

      if ( t_beta )
      {
        init_cov.multiplyByScalar( tau2_betas[ index ]->lastItemDrawn() );
      }
 
      beta[ index ]->updateInitialCovariance( &init_cov );
#ifdef DEBUGBETA
      printf("Calling beta[%d]->updateInitialCovariance with \n", index);
      init_cov.Print();  fflush(stdout);
#endif
    }

    if ( second_stage_effects )
    {

#ifdef FIX1
      CMatrix tmpmat = *(beta_predictors[ index ]);
      CVector init_mean = tmpmat * alpha->lastItemDrawn();
#else
      CVector init_mean( (*(beta_predictors[ index ])) * alpha->lastItemDrawn() );
#endif
      beta[ index ]->updateInitialMean( &init_mean );
#ifdef DEBUGBETA
      printf("Calling beta[%d]->updateInitialMean with \n", index);
      init_mean.Print();  fflush(stdout);
#endif
    }
  }

  //ready to update mean and covariance
  // Bug fix by Dawn: handling when tau2 is a scalar.  
  // The following code is inspired by that in 
  // BayesianHierarchicalLinearModel::gibbsUpdateBeta
  if ( invWishart_betaCov || t_beta ){
    beta[ index ]->update( beta_shift_mean, beta_shift_cov );    
#ifdef DEBUGBETA
    printf("Calling beta[%d]->update with \n", index );
    beta_shift_mean.Print(); beta_shift_cov.Print(); fflush(stdout);
#endif
  } else {
    CVector *precision_tau2 = new CVector( dim_beta );
    for( int i=0; i < dim_beta; i++ )
      precision_tau2->Val(i) = 1 / tau2[i]->lastDraw().getScalar();

  	beta[ index ]->update( (*precision_tau2), beta_shift_mean, beta_shift_cov );    
  	delete( precision_tau2 );
#ifdef DEBUGBETA
    printf("Calling beta[%d]->update( %f, \n", index, precision_tau2);
    beta_shift_mean.Print(); beta_shift_cov.Print(); fflush(stdout);
#endif
  }

  delete beta_vectors[0];
  delete [] beta_vectors;

}//end




void BayesianHierarchicalGlmModel::gibbsUpdateGamma()
{
  //Gibbs proposal

  //get shift mean and cov from current values
  CMatrix gamma_shift_cov( (-1.0) * glm_model[0]->metropolisHastingsHessian( number_of_groups + 1, fixed_predictors ) );


  DistributionParameter ** gamma_vectors;
  gamma_vectors = new DistributionParameter * [1];
#ifdef FIX1
  CVector tmpvec = gamma->lastItemDrawn();
  gamma_vectors[0] = new DistributionParameter( tmpvec );
#else
  gamma_vectors[0] = new DistributionParameter( gamma->lastItemDrawn() );
#endif
//  gamma_vectors[0] = new DistributionParameter( gamma->lastItemDrawn() );

  CVector gamma_shift_mean( glm_model[0]->metropolisHastingsMean( number_of_groups + 1, fixed_predictors, gamma_vectors ) );

  //printf("gibbs update gamma: shift mean and cov are \n" );
  //gamma_shift_mean.Print();
  //gamma_shift_cov.Print();


  if ( !non_informative_gamma && t_gamma )
  {
    CMatrix init_cov( gamma->veryFirstCovariance() );
    init_cov.multiplyByScalar( tau2_gamma->lastItemDrawn() );
    gamma->updateInitialCovariance( &init_cov );
  }

  //ready to update mean and covariance
  gamma->update( gamma_shift_mean, gamma_shift_cov );    

  delete gamma_vectors[0];
  delete [] gamma_vectors;

}//end


void BayesianHierarchicalGlmModel::metropolisHastingsUpdateBeta( int index )
{
  int dimB;
  CVector * cholesky_dec;

  if ( glm_model[0]->gibbsDrawingsForCoefficients() )
  {
    gibbsUpdateBeta( index );
  }
  else
  {
    //get Hessian from Glm model 
    // Bug fix by DBW: the following was previously:
    // random_predictors[1]->parameter().Val(0,0) = index;
    // which fails to change the value of random_predictors[1].  The following is 
    // the corrected version.
    random_predictors[1]->updateParameter( (double) index );
    CMatrix beta_hessian( (-1.0) * glm_model[0]->metropolisHastingsHessian( 3, random_predictors ) );

    if ( !non_informative_beta )
    {
#ifdef FIX1
      CMatrix tmpmat = betaInverseCovariance( index );
      beta_hessian.add( tmpmat );
#else
      beta_hessian.add( betaInverseCovariance( index ) );
#endif
//      beta_hessian.add( betaInverseCovariance( index ) );
    }

    //CMatrix beta_cov( beta_hessian.inverse() );

    CVector beta_mean( beta[ index ]->lastItemDrawn() );

    //keep this in case we need it
    CMatrix prev_beta_cov( beta[ index ]->covariance() );

    dimB = beta[ index ]->dimension();
    cholesky_dec = new CVector( ( dimB * ( dimB + 1 ) ) / 2 );
    // Bug fix by Dawn: added " && update_for_hessians" in the following line
    if ( ( beta_hessian.exploreCholeskyDecomposition( cholesky_dec ) ) && update_for_hessians )
    {
      //CMatrix beta_cov( beta_hessian.inverse() );

      CMatrix inverse_cholLT( dimB, dimB );
      inverse_cholLT.assignInverseOfLowerTriangular( (*cholesky_dec) );
      CMatrix beta_cov( inverse_cholLT.T() * inverse_cholLT );

      CMatrix beta_ID( dimB, dimB );
      beta_ID.setDiagonal( 1.0 );
      CMatrix is_ID_matrix( beta_cov * beta_hessian );
      if ( (is_ID_matrix - beta_ID).isZero( 0.0001 ) )
      {
        //accept beta_cov
        beta[ index ]->updateInitialMean( &beta_mean );
        beta[ index ]->updateInitialCovariance( &beta_cov );
      }
      else
      {
        //try another covariance
        beta[ index ]->updateInitialMean( &beta_mean );
        beta[ index ]->updateInitialCovariance( &prev_beta_cov );
      }
    }
    else
    {
      //try another covariance
      beta[ index ]->updateInitialMean( &beta_mean );
      beta[ index ]->updateInitialCovariance( &prev_beta_cov );

#ifdef DEBUGBETA
      printf( "BayesianHierarchicalGlmModel::metropolisHastingsUpdateBeta: Using previous covariance\n" );
      fflush( stdout );
#endif
    }
    delete cholesky_dec;
  }

}//end


void BayesianHierarchicalGlmModel::metropolisHastingsUpdateGamma()
{
  int dimG;
  CVector * cholesky_dec;

  if ( glm_model[0]->gibbsDrawingsForCoefficients() )
  {
    gibbsUpdateGamma();
  }
  else
  {
    CMatrix gamma_hessian_correction( gamma->dimension(), gamma->dimension() );

    //get Hessian from Glm model 
    CMatrix gamma_hessian( (-1.0) * glm_model[0]->metropolisHastingsHessian( number_of_groups + 1, fixed_predictors ) );
 
    if ( !non_informative_gamma )
    {  
      gamma_hessian_correction = gamma->veryFirstInverseCovariance();

      if ( t_gamma )
      {
        gamma_hessian_correction.multiplyByScalar( ( 1.0 / tau2_gamma->lastItemDrawn() ) );
      }

      gamma_hessian.add( gamma_hessian_correction );
    }

    //CMatrix gamma_cov( gamma_hessian.inverse() );

    CVector gamma_mean( gamma->lastItemDrawn() );

    //keep this in case we need it
    CMatrix prev_gamma_cov( gamma->covariance() );

    dimG = gamma->dimension();
    cholesky_dec = new CVector( ( dimG * ( dimG + 1 ) ) / 2 );
    if ( gamma_hessian.exploreCholeskyDecomposition( cholesky_dec ) )
    {
      //CMatrix gamma_cov( gamma_hessian.inverse() );

      CMatrix inverse_cholLT( dimG, dimG );
      inverse_cholLT.assignInverseOfLowerTriangular( (*cholesky_dec) );
      CMatrix gamma_cov( inverse_cholLT.T() * inverse_cholLT );

      CMatrix gamma_ID( dimG, dimG );
      gamma_ID.setDiagonal( 1.0 );
      CMatrix is_ID_matrix( gamma_cov * gamma_hessian );
      if ( (is_ID_matrix - gamma_ID).isZero( 0.0001 ) )
      {
        gamma->updateInitialMean( &gamma_mean );
        gamma->updateInitialCovariance( &gamma_cov );
      }
      else
      {
        //try another covariance
        gamma->updateInitialMean( &gamma_mean );
        gamma->updateInitialCovariance( &prev_gamma_cov );

        //printf( "BayesianHierarchicalGlmModel::metropolisHastingsUpdateGamma: Using previous covariance\n" );
        //fflush( stdout );
      }
    }
    else
    {
      //try another covariance
      gamma->updateInitialMean( &gamma_mean );
      gamma->updateInitialCovariance( &prev_gamma_cov );

      //printf( "BayesianHierarchicalGlmModel::metropolisHastingsUpdateGamma: Using previous covariance\n" );
      //fflush( stdout );
    }
    delete cholesky_dec;
  }

}//end




void BayesianHierarchicalGlmModel::gibbsUpdateSecondStageWorkingMatrices()
{
  int index, j;

  CVector tmpVec ( dim_beta );
  CMatrix beta_inv_cov( beta[ 0 ]->dimension(), beta[ 0 ]->dimension() );

  zTvIb->setToZero();
  
  
#ifdef DEBUGALPHA
  printf("zTb = ");
#endif
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
  	// Bug fix by Dawn: changed to match BayesianHierarchicalLinearModel.  The following line was changed from:
    //beta_inv_cov = betaInverseCovariance( index );
    beta_inv_cov = beta[ index ]->veryFirstInverseCovariance();

    ( *(zTb[ index ]) ) = beta_predictors[ index ]->T() * ( beta_inv_cov * tmpVec );

    zTvIb->add( ( *(zTb[ index ]) ) );

#ifdef DEBUGALPHA
    zTb[index]->Print(); 
    printf("beta_predictors[index] = \n");  beta_predictors[ index ]->Print();
#endif
  }
  zTvIz->setToZero();
  /*if ( invWishart_betaCov || t_beta )
  {*/
  for ( index = 0; index < number_of_groups; index++ )
  {
  	// In the following line, for invWishart_betaCov = F, beta_inv_cov is supposed to be the identity matrix.
  	beta_inv_cov = beta[ index ]->veryFirstInverseCovariance();
    
    if ( !invWishart_betaCov ) {
    	// Divide the diagonal elements of beta_inv_cov by the corresponding tau2 parameter.
      for( j = 0; j < dim_beta; j++ )
        beta_inv_cov.Val(j, j) = beta_inv_cov.Val(j, j) / tau2[ j ]->lastDraw().getScalar();
    }

    ( *(zTz[ index ]) ) = beta_predictors[ index ]->T() * ( beta_inv_cov * (*beta_predictors[ index ] ) );

    zTvIz->add( ( *(zTz[ index ]) ) );
  }//end for loop
#ifdef DEBUGALPHA
  printf("last beta_inv_cov = \n");  beta[ number_of_groups - 1 ]->veryFirstInverseCovariance().Print();
#endif
  
}//end


void BayesianHierarchicalGlmModel::gibbsUpdateAlpha()
{
  CMatrix sigma_alpha( alpha->dimension(), alpha->dimension() );
  CVector mu_alpha( alpha->dimension() );

  gibbsUpdateSecondStageWorkingMatrices();

  // Bug fix by Dawn: previously this was
  // mu_alpha = (*zTvIb);
  // sigma_alpha = (*zTvIz);
  // which omits a factor of tau^2.  The following code was
  // copied from BayesianHierarchicalLinearModel::gibbsUpdateAlpha()
  /*if ( !invWishart_betaCov )
  {
    mu_alpha = (*zTvIb) * ( 1 / tau2->lastDraw().getScalar() );
    sigma_alpha = ( 1 / tau2->lastDraw().getScalar() ) * (*zTvIz);
  }
  else
  {*/
    mu_alpha = (*zTvIb);
    sigma_alpha = (*zTvIz);
  //}

#ifdef DEBUGALPHA
  printf("mu_alpha = \n");  mu_alpha.Print(); 
  printf("sigma_alpha = \n");  sigma_alpha.Print();
  printf("zTvIz = \n"); zTvIz->Print();
  fflush(stdout);
#endif

  if ( t_alpha )
  {
    double precision_tau2 = 1 / tau2_alpha->lastItemDrawn();
    alpha->update( precision_tau2, mu_alpha, sigma_alpha );
  }
  else
  {
    alpha->update( mu_alpha, sigma_alpha );
  }

}//end



void BayesianHierarchicalGlmModel::gibbsUpdateTau2()
{
  int i, j;
  double scale, local_scale, add_df;

  if ( invChisq_tau2 || nonInformativePower_tau2 || properNonInfoPrior_tau2 )
  {
  	for ( j = 0; j < dim_beta; j++ ){
      scale = 0;
      #ifdef DEBUGTAU2
        printf("beta[0]->veryFirstInvCov() = "); beta[0]->veryFirstInverseCovariance().Print(); 
        printf("second_residuals = \n");
      #endif
      for ( i = 0; i < number_of_groups; i++ ){    
        local_scale = SQR( second_residuals[i]->Val(j) );
        #ifdef DEBUGTAU2
          second_residuals[i]->Print(); 
          printf("local_scale = %f\n", local_scale);
        #endif

        if ( t_beta )
        {
          local_scale /= tau2_betas[i]->lastItemDrawn();
        }

        scale += local_scale;
      }//end i loop


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
        #ifdef DEBUGTAU2
          printf("Calling tau2ICS[j]->update( %f, %f )\n", add_df, scale );
          fflush(stdout);
        #endif
      }
      else
      {
        tau2PNIP[j]->setScale( scale );
        #ifdef DEBUGTAU2
          printf("Calling tau2PNIP[j]->setScale( %f )\n", scale ); fflush(stdout);
        #endif
      }
    }// end j loop 
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
#ifdef DEBUGTAU2
    printf("Calling tau2InvWishart->update( %f, \n", add_df );
    scale_matrix.Print(); fflush(stdout);
#endif

  }

}//end


void BayesianHierarchicalGlmModel::gibbsUpdateTau2Gamma()
{
  double scale, add_df;
  CVector gamma_diff( gamma->dimension() );

  gamma_diff = gamma->lastItemDrawn() - gamma->initialMean();
  scale = M_D( gamma->veryFirstInverseCovariance(), gamma_diff );

  add_df = gamma->dimension();

  tau2_gamma->update( add_df, scale );

}//end



void BayesianHierarchicalGlmModel::gibbsUpdateTau2Beta( int index )
{
  double scale, add_df;

  if ( !invWishart_betaCov )
  {
    scale = M_D( beta[ index ]->veryFirstInverseCovariance(), ( *(second_residuals[ index ]) ) );
    scale /= tau2[0]->lastDraw().getScalar();
  }
  else
  {
    scale = M_D( tau2[0]->lastDraw().getMatrix().inverse(), (*(second_residuals[ index ])) );
  }

  add_df = beta[ index ]->dimension();

  tau2_betas[ index ]->update( add_df, scale );

}//end


void BayesianHierarchicalGlmModel::gibbsUpdateTau2Alpha()
{
  double scale, add_df;
  CVector alpha_diff( alpha->dimension() );

  alpha_diff = alpha->lastItemDrawn() - alpha->initialMean();
  scale = M_D( alpha->initialInverseCovariance(), alpha_diff );

  add_df = alpha->dimension();

  tau2_alpha->update( add_df, scale );

}//end


void BayesianHierarchicalGlmModel::updateVariableForProposal( int index )
{
	
  #ifdef DEBUG2
    printf( "HGLM:update for proposal: index is %d.\n  ", index ); fflush( stdout );
  #endif
  if ( index < number_of_variables && index >= 0 )
  {
    
    if ( !strncmp( distr_map[ index ], "glm:", 4 ) )
    {
      int glm_index;
      char * ptr, * temp_distr;

      temp_distr = new char [ strlen( distr_map[ index ] ) + 1 ];
      sprintf( temp_distr, "%s", distr_map[ index ] );
      ptr = strtok( temp_distr, ":" );
      if ( ptr != NULL )
      {
        ptr = strtok( NULL, ":" );
        if ( ptr != NULL )
        {
          glm_index = atoi( ptr );
          glm_model[0]->updateVariableForProposal( glm_index );
        }
        else
        {
          printf( " BayesianHierarchicalGlmModel::updateVariableForProposal: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }//end if null
      else
      {
        printf( " BayesianHierarchicalGlmModel::updateVariableForProposal: Wrong argument in [%s].\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if glm
    else if ( !strcmp( distr_map[ index ], "gamma" ) )
    {
      metropolisHastingsUpdateGamma();
    }
    else if ( !strcmp( distr_map[ index ], "alpha" ) )
    {
      gibbsUpdateAlpha();
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
              || !strncmp( distr_map[ index ], "tau2b", 5 ) )
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

          //printf( "BGlm:full conditional: group index is %d\n", group_index );

          if ( group_index >= 0 && group_index < number_of_groups )
          {
            if ( !strncmp( distr_map[ index ], "beta", 4 ) )
            {
              metropolisHastingsUpdateBeta( group_index );
            }
            else if ( !strncmp( distr_map[ index ], "tau2b", 5 ) )
            {
              gibbsUpdateTau2Beta( group_index );
            }          
            else
            {
              printf( " BayesianHierarchicalGlmModel::updateVariableForProposal: Wrong argument in [%s].\n", distr_map[ index ] );
            }
          }//end if
          else
          {
            printf( " BayesianHierarchicalGlmModel::updateVariableForProposal: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
          } 
        }//end if not null
        else
        {
          printf( " BayesianHierarchicalGlmModel::updateVariableForProposal: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }
      else
      {
        printf( " BayesianHierarchicalGlmModel::updateVariableForProposal: Wrong argument in [%s]. Number expected.\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if 
  }//end if index is valid
  else
  {
    printf( "BayesianHierarchicalGlmModel::updateVariableForProposal: Variable index [%d] does not exist.\n", index );
  }

}//end



void BayesianHierarchicalGlmModel::drawVariableFromProposal( int index )
{
#ifdef DEBUG1
  printf("Updating variable %s (index %d)\n", distr_map[ index ], index); fflush(stdout);
#endif

  if ( index < number_of_variables && index >= 0 )
  {
    if ( !strncmp( distr_map[ index ], "glm:", 4 ) )
    {
      int glm_index;
      char * ptr, * temp_distr;

      temp_distr = new char [ strlen( distr_map[ index ] ) + 1 ];
      sprintf( temp_distr, "%s", distr_map[ index ] );
      ptr = strtok( temp_distr, ":" );
      if ( ptr != NULL )
      {
        ptr = strtok( NULL, ":" );
        if ( ptr != NULL )
	{
          glm_index = atoi( ptr );
          glm_model[0]->drawVariableFromProposal( glm_index );
        }
        else
        {
          printf( " BayesianHierarchicalGlmModel::drawVariableFromProposal: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }//end if null
      else
      {
        printf( " BayesianHierarchicalGlmModel::drawVariableFromProposal: Wrong argument in [%s].\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if glm
    else if ( !strcmp( distr_map[ index ], "gamma" ) )
    {
#ifdef DEBUGGAMMA
      printf("Calling gamma->draw() where mean and Cholesky inv cov = \n");
      gamma->meanVec().Print(); 
      gamma->choleskyInvCov().Print();  fflush(stdout); 
#endif
      gamma->draw();
    }
    else if ( !strcmp( distr_map[ index ], "alpha" ) )
    {
      alpha->draw();
#ifdef DEBUGALPHA
      printf("Calling alpha->draw()\n"); 
      printf("mean vector = \n"); alpha->mean().Print();
      printf("Cholesky decomp of inv cov = \n"); alpha->choleskyInvCov().Print();
      fflush(stdout); 
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
        for( int i=0; i<dim_beta; i++ ){
          tau2[i]->draw();
          #ifdef DEBUG1
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
#ifdef DEBUGTAU2A
      printf("Calling tau2_alpha->draw()\n"); fflush(stdout); 
#endif
    }
    else if ( !strcmp( distr_map[ index ], "tau2gm" ) )
    {
      tau2_gamma->draw();
#ifdef DEBUGTAU2GM
      printf("Calling tau2_gamma->draw()\n"); fflush(stdout); 
#endif
    }
    else if ( !strncmp( distr_map[ index ], "beta", 4 ) 
              || !strncmp( distr_map[ index ], "tau2b", 5 ) )
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

          if ( group_index >= 0 && group_index < number_of_groups )
          {
            if ( !strncmp( distr_map[ index ], "beta", 4 ) )
            {
              beta[ group_index ]->draw();
#ifdef DEBUGBETA
              printf("Calling beta[%d]->draw()\n", group_index ); 
              printf("mean vector = \n"); beta[ group_index ]->mean().Print();
              printf("Cholesky decomp of inv cov = \n"); beta[ group_index ]->choleskyInvCov().Print();
              fflush(stdout); 
#endif
            }
            else if ( !strncmp( distr_map[ index ], "tau2b", 5 ) )
            {
              tau2_betas[ group_index ]->draw();
#ifdef DEBUGTAU2B
              printf("Calling tau2_betas[%d]->draw()\n", group_index); fflush(stdout); 
#endif
            }          
          }
          else
          {
            printf( " BayesianHierarchicalGlmModel::drawVariableFromProposal: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
          }          
        }//end if not null
        else
        {
          printf( " BayesianHierarchicalGlmModel::drawVariableFromProposal: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }
      else
      {
        printf( " BayesianHierarchicalGlmModel::drawVariableFromProposal: Wrong argument in [%s]. Number expected.\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if 
  }//end if index is valid
  else
  {
    printf( "BayesianHierarchicalGlmModel::drawVariableFromProposal: Variable index [%d] does not exist.\n", index );
  }

}//end



void BayesianHierarchicalGlmModel::initializeTemporaryStructures()
{
  int i, k, len_digits_vars, len_digits_obs;

  //initialize the linear part of the link function (the mu vector). 
  if ( random_effects || fixed_effects )
  {
    updateLinearResponse();
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

  if ( random_effects && !non_informative_beta )
  {
    distr_map[k] = new char [ 5 ];
    sprintf( distr_map[k], "tau2" );
    k++;
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


  //add glm variables
  glm_model[0]->initializeTemporaryStructures();

  for ( i = 0; i < glm_model[0]->numberOfVariables(); i++ )
  {
    distr_map[ k ] = new char [ 7 + len_digits_vars ];
    sprintf( distr_map[ k ], "glm:%d", i );
    k++;
  }

  #ifdef DEBUG1
    printf( "HGLM: Finished initializing temporary structures.\n" );
    printf( "Total vars = %d.  Number vars = %d. The map is:\n", k, number_of_variables  );
    for ( i = 0; i < k; i++ ){
      printf("map[%d] = %s\n", i, distr_map[ i ] );
    }
    fflush( stdout );
  #endif


}//end



void BayesianHierarchicalGlmModel::dataAugmentationInitialDraws()
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

  glm_model[0]->dataAugmentationInitialDraws();

  #ifdef DEBUG1
    printf( "Finished data augmentation initial draws\n" ); fflush(stdout);
  #endif
}//end



void BayesianHierarchicalGlmModel::createOutput( int simulations_to_keep )
{
  int i;

  number_of_simulations = simulations_to_keep;

  #ifdef DEBUG1
    printf("HGLM: Creating objects to hold simulation output\n");fflush(stdout);
  #endif

  if ( simulations_to_keep > 0 )
  {
    if ( random_effects )
    {
      simulated_beta = new CMatrix * [ number_of_groups ];
      for ( i = 0; i < number_of_groups; i++ )
      {
        simulated_beta[i] = new CMatrix( simulations_to_keep, random_predictors[2]->getMatrix().Col() );
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
      simulated_gamma = new CMatrix( simulations_to_keep, fixed_predictors[2]->getMatrix().Col() );
      if ( t_gamma && keep_tau2_error )
      {
        simulated_tau2_gamma = new CVector( simulations_to_keep );
      }
    }

    glm_model[0]->createOutput( simulations_to_keep );

  }//end if keep any simulations
  
  #ifdef DEBUG1
    printf("HGLM: Finished creating objects to hold simulation output\n");  fflush(stdout);
  #endif

}//end


void BayesianHierarchicalGlmModel::keepSimulation( int simul_number )
{
  int i;

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


    glm_model[0]->keepSimulation( simul_number );

  }//end if keep any simulations

}//end
    


void BayesianHierarchicalGlmModel::simulationsToArray( double * simul_output, int simulations_to_keep )
{
  int i, j, total_dim, start_dim, g, v, s;

  start_dim = 0;
  if ( random_effects )
  {
    total_dim = simulations_to_keep * random_predictors[2]->getMatrix().Col();
    for ( g = 0; g < number_of_groups; g++ )
    {
      for ( v = 0; v < random_predictors[2]->getMatrix().Col(); v++ )
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
    for ( v = 0; v < fixed_predictors[2]->getMatrix().Col(); v++ )
    {
      for ( s = 0; s < simulations_to_keep; s++ )
      {
        simul_output[ start_dim + v * total_dim + s ] = simulated_gamma->Val(s,v);
      }
    }
    start_dim += fixed_predictors[2]->getMatrix().Col() * total_dim;
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


  glm_model[0]->simulationsToArray( simul_output, simulations_to_keep, start_dim );

}//end



void BayesianHierarchicalGlmModel::simulationsToArrayWithMeans( double * simul_output, int simulations_to_keep )
{
  int i, j, total_dim, start_dim, g, v, s, start_i, end_i, g_length;

  #ifdef DEBUG1
    printf( "HGLM: Converting samples into an array with means \n" );
  #endif
  
  start_dim = 0;
  if ( random_effects )
  {
    total_dim = simulations_to_keep * random_predictors[2]->getMatrix().Col();
    for ( g = 0; g < number_of_groups; g++ )
    {
      for ( v = 0; v < random_predictors[2]->getMatrix().Col(); v++ )
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
    for ( v = 0; v < fixed_predictors[2]->getMatrix().Col(); v++ )
    {
      for ( s = 0; s < simulations_to_keep; s++ )
      {
        simul_output[ start_dim + v * total_dim + s ] = simulated_gamma->Val(s,v);
      }
    }
    start_dim += fixed_predictors[2]->getMatrix().Col() * total_dim;
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


  
  CMatrix simulated_mu( simulations_to_keep, number_of_observations );  
  for ( s = 0; s < simulations_to_keep; s++ )
  {
    start_i = 0;
    end_i = 0;
    for ( g = 0; g < number_of_groups; g++ )
    {
      if ( random_effects )
      {     
        g_length = random_predictors[ g + 2 ]->getMatrix().Row();
      }
      else if ( fixed_effects )
      {
        g_length = fixed_predictors[ g + 2 ]->getMatrix().Row();
      }

      end_i += g_length;
      CVector mu( g_length );

      if ( random_effects && fixed_effects )
      {
        mu = random_predictors[ g + 2 ]->getMatrix() * simulated_beta[ g ]->getRow( s ) + fixed_predictors[ g + 2 ]->getMatrix() * simulated_gamma->getRow( s );
      }
      else if ( random_effects )
      {
        mu = random_predictors[ g + 2 ]->getMatrix() * simulated_beta[ g ]->getRow( s );
      }
      else if ( fixed_effects )
      {
        mu = fixed_predictors[ g + 2 ]->getMatrix() * simulated_gamma->getRow( s );
      }

      simulated_mu.setSubRow( s, start_i, end_i, mu );

      //printf( "HGlm: sim to array:  s = %d,  start_i = %d, end_i = %d, mu = \n", s, start_i, end_i );
      //simulated_mu.Print();

      start_i = end_i;
    }//end for g
  }//end for s

  glm_model[0]->simulationsToArray( simul_output, start_dim, &simulated_mu );

  #ifdef DEBUG1
    printf( "HGLM: Finished converting samples into an array with means \n" );
  #endif

}//end




double BayesianHierarchicalGlmModel::logRatioProposal( int index )
{
  double ratio;
  int dimB;
  CVector * cholesky_dec;

  //printf( "BGlm:logRatio for proposal: index is %d.\n  ", index ); fflush( stdout );

  ratio = 0;

  if ( index < number_of_variables && index >= 0 )
  {
    if ( !strncmp( distr_map[ index ], "glm:", 4 ) )
    {
      int glm_index;
      char * ptr, * temp_distr;

      temp_distr = new char [ strlen( distr_map[ index ] ) + 1 ];
      sprintf( temp_distr, "%s", distr_map[ index ] );
      ptr = strtok( temp_distr, ":" );
      if ( ptr != NULL )
      {
        ptr = strtok( NULL, ":" );
        if ( ptr != NULL )
        {
          glm_index = atoi( ptr );
          ratio = glm_model[0]->logRatioProposal( glm_index );
        }
        else
        {
          printf( " BayesianHierarchicalGlmModel::logRatioProposal: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }//end if null
      else
      {
        printf( " BayesianHierarchicalGlmModel::logRatioProposal: Wrong argument in [%s].\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if glm
    else if ( !strcmp( distr_map[ index ], "gamma" ) && update_for_hessians && !glm_model[0]->gibbsDrawingsForCoefficients() )
    {
      DistributionParameter ** gamma_vectors;
      gamma_vectors = new DistributionParameter * [2];
#ifdef FIX1
      CVector tmpvec1 = gamma->lastItemDrawn();
      CVector tmpvec2 = gamma->initialMean();
      gamma_vectors[0] = new DistributionParameter( tmpvec1 );
      gamma_vectors[1] = new DistributionParameter( tmpvec2 );
#else
      gamma_vectors[0] = new DistributionParameter( gamma->lastItemDrawn() );
      gamma_vectors[1] = new DistributionParameter( gamma->initialMean() );
#endif
//      gamma_vectors[0] = new DistributionParameter( gamma->lastItemDrawn() );
//      gamma_vectors[1] = new DistributionParameter( gamma->initialMean() );

      CMatrix glm_proposal_hessian( (-1.0) * glm_model[0]->metropolisHastingsProposalHessian( number_of_groups + 1, fixed_predictors, gamma_vectors ) );

      //CMatrix gamma_cov( glm_proposal_hessian.inverse() );
      CVector proposal_mean( gamma->lastItemDrawn() );
      NormalDistribution proposal( gamma->dimension() );

      //keep this in case we need it
      //CMatrix prev_gamma_cov( gamma->covariance() );

      dimB = gamma->dimension();
      cholesky_dec = new CVector( ( dimB * ( dimB + 1 ) ) / 2 );
      if ( glm_proposal_hessian.exploreCholeskyDecomposition( cholesky_dec ) )
      {
        //CMatrix gamma_cov( glm_proposal_hessian.inverse() );

        CMatrix inverse_cholLT( dimB, dimB );
        inverse_cholLT.assignInverseOfLowerTriangular( (*cholesky_dec) );
        CMatrix gamma_cov( inverse_cholLT.T() * inverse_cholLT );

        CMatrix gamma_ID( dimB, dimB );
        gamma_ID.setDiagonal( 1.0 );
        CMatrix is_ID_matrix( gamma_cov * glm_proposal_hessian );
        if ( (is_ID_matrix - gamma_ID).isZero( 0.0001 ) )
        {
          proposal.setMean( &proposal_mean );
          proposal.setCovariance( &gamma_cov );
          ratio = gamma->logDensity( (*(gamma_vectors[0])) ) - proposal.logDensity( (*(gamma_vectors[1])) );
        }
        else
        {
          //printf( "BayesianHierarchicalGlmModel::logRatioProposal: Using covariance from gamma for proposal.\n" );
          //fflush( stdout );

          ratio = 0;
        }            
      }
      else
      {
        //try another covariance
        //proposal.setMean( &proposal_mean );
        //proposal.setCovariance( &prev_gamma_cov );


        //printf( "BayesianHierarchicalGlmModel::logRatioProposal: Using covariance from gamma for proposal.\n" );
        //fflush( stdout );

        ratio = 0;
      }    

      delete cholesky_dec;

      //ratio = gamma->logDensity( (*(gamma_vectors[0])) ) - proposal.logDensity( (*(gamma_vectors[1])) );

      delete gamma_vectors[0];
      delete gamma_vectors[1];  
      delete [] gamma_vectors;
    }// end gamma
    else if ( !strcmp( distr_map[ index ], "alpha" ) )
    {
      ratio = 0;
    }
    else if ( !strcmp( distr_map[ index ], "tau2" ) )
    {
      ratio = 0;
    }
    else if ( !strcmp( distr_map[ index ], "tau2a" ) )
    {
      ratio = 0;
    }
    else if ( !strcmp( distr_map[ index ], "tau2gm" ) )
    {
      ratio = 0;
    }
    else if ( !strncmp( distr_map[ index ], "beta", 4 ) 
              || !strncmp( distr_map[ index ], "tau2b", 5 ) )
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

          if ( group_index >= 0 && group_index < number_of_groups )
	  {
            if ( !strncmp( distr_map[ index ], "beta", 4 ) && update_for_hessians && !glm_model[0]->gibbsDrawingsForCoefficients() )
            {
              DistributionParameter ** beta_vectors;
              beta_vectors = new DistributionParameter * [2];
#ifdef FIX1
              CVector tmpvec1 = beta[ group_index ]->lastItemDrawn();
              CVector tmpvec2 = beta[ group_index ]->initialMean();
              beta_vectors[0] = new DistributionParameter( tmpvec1 );
              beta_vectors[1] = new DistributionParameter( tmpvec2 );
#else
              beta_vectors[0] = new DistributionParameter( beta[ group_index ]->lastItemDrawn() );
              beta_vectors[1] = new DistributionParameter( beta[ group_index ]->initialMean() );
#endif
//              beta_vectors[0] = new DistributionParameter( beta[ group_index ]->lastItemDrawn() );
//              beta_vectors[1] = new DistributionParameter( beta[ group_index ]->initialMean() );

              // Bug fix by Dawn: the following was previously:
              // random_predictors[1]->parameter().Val(0,0) = group_index;
              // which fails to change the value of random_predictors[1].  The following is 
              // the corrected version.
              random_predictors[1]->updateParameter( (double) group_index );

              CMatrix glm_proposal_hessian( (-1.0) * glm_model[0]->metropolisHastingsProposalHessian( 2, random_predictors, beta_vectors ) );

              //CMatrix beta_cov( glm_proposal_hessian.inverse() );

#ifdef DEBUGBETA
              printf("logRatioProposal: glm proposal beta matrix = \n" );
              glm_proposal_hessian.Print();
#endif

              CVector proposal_mean( beta[ group_index ]->lastItemDrawn() );
              NormalDistribution proposal( beta[ group_index ]->dimension() );

              //keep this in case we need it
              //CMatrix prev_beta_cov( beta[ group_index ]->covariance() );

              dimB = beta[ group_index ]->dimension();
              cholesky_dec = new CVector( ( dimB * ( dimB + 1 ) ) / 2 );
              if ( glm_proposal_hessian.exploreCholeskyDecomposition( cholesky_dec ) )
              {
                CMatrix inverse_cholLT( dimB, dimB );
                inverse_cholLT.assignInverseOfLowerTriangular( (*cholesky_dec) );
                CMatrix beta_cov( inverse_cholLT.T() * inverse_cholLT );

                CMatrix beta_ID( dimB, dimB );
                beta_ID.setDiagonal( 1.0 );
                CMatrix is_ID_matrix( beta_cov * glm_proposal_hessian );
                if ( (is_ID_matrix - beta_ID).isZero( 0.0001 ) )
                {
                  //accept beta_cov
                  proposal.setMean( &proposal_mean );
                  proposal.setCovariance( &beta_cov );
                  ratio = beta[ group_index ]->logDensity( (*(beta_vectors[0])) ) - proposal.logDensity( (*(beta_vectors[1])) );
                } else {
#ifdef DEBUGBETA
                  printf( "BayesianHierarchicalGlmModel::logRatioProposal: Using covariance from beta for proposal.\n" );
                  fflush( stdout );
#endif

                  ratio = 0;                  
                }
              }
              else
              {
                //try another covariance
                //proposal.setMean( &proposal_mean );
                //proposal.setCovariance( &prev_beta_cov );

#ifdef DEBUGBETA
                printf( "BayesianHierarchicalGlmModel::logRatioProposal: Using covariance from beta for proposal.\n" );
                fflush( stdout );
#endif

                ratio = 0;
              }    

              delete cholesky_dec;

              //ratio = beta[ group_index ]->logDensity( (*(beta_vectors[0])) ) - proposal.logDensity( (*(beta_vectors[1])) );

              delete beta_vectors[0];
              delete beta_vectors[1];  
              delete [] beta_vectors;
            }
            else if ( !strncmp( distr_map[ index ], "tau2b", 5 ) )
            {
              ratio = 0;
            }          
          }
          else
          {
            printf( " BayesianHierarchicalGlmModel::logRatioProposal: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
          }          
        }//end if not null
        else
        {
          printf( " BayesianHierarchicalGlmModel::logRatioProposal: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }
      else
      {
        printf( " BayesianHierarchicalGlmModel::logRatioProposal: Wrong argument in [%s]. Number expected.\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if 
  }//end if index is valid
  else
  {
    printf( "BayesianHierarchicalGlmModel::logRatioProposal: Variable index [%d] does not exist.\n", index );
  }

#ifdef DEBUG1
  printf("proposal: ratio[%d] = %f\n", index, ratio ); fflush(stdout);
#endif

  return ( ratio );


}//end



double BayesianHierarchicalGlmModel::logRatioTargetDensity( int index )
{
  double ratio, glm_ratio, scale;

  //printf( "BGlm:logRatio for Target Density: index is %d.\n  ", index ); fflush( stdout );  

  ratio = 0;
  if ( index < number_of_variables && index >= 0 )
  {
    if ( !strncmp( distr_map[ index ], "glm:", 4 ) )
    {
      int glm_index;
      char * ptr, * temp_distr;

      temp_distr = new char [ strlen( distr_map[ index ] ) + 1 ];
      sprintf( temp_distr, "%s", distr_map[ index ] );
      ptr = strtok( temp_distr, ":" );
      if ( ptr != NULL )
      {
        ptr = strtok( NULL, ":" );
        if ( ptr != NULL )
	{
          glm_index = atoi( ptr );

          ratio = glm_model[0]->logRatioTargetDensity( glm_index );
        }
        else
        {
          printf( " BayesianHierarchicalGlmModel::logRatioTargetDensity: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }//end if null
      else
      {
        printf( " BayesianHierarchicalGlmModel::logRatioTargetDensity: Wrong argument in [%s].\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if glm
    else if ( !strcmp( distr_map[ index ], "gamma" ) && !glm_model[0]->gibbsDrawingsForCoefficients() )
    {
      DistributionParameter ** gamma_vectors;
      gamma_vectors = new DistributionParameter * [2];
#ifdef FIX1
      CVector tmpvec1 = gamma->lastItemDrawn();
      CVector tmpvec2 = gamma->initialMean();
      gamma_vectors[0] = new DistributionParameter( tmpvec1 );
      gamma_vectors[1] = new DistributionParameter( tmpvec2 );
#else
      gamma_vectors[0] = new DistributionParameter( gamma->lastItemDrawn() );
      gamma_vectors[1] = new DistributionParameter( gamma->initialMean() );
#endif
//      gamma_vectors[0] = new DistributionParameter( gamma->lastItemDrawn() );
//      gamma_vectors[1] = new DistributionParameter( gamma->initialMean() );

      glm_ratio = glm_model[0]->logRatioTargetDensity( 2, fixed_predictors, gamma_vectors );

      delete gamma_vectors[0];
      delete gamma_vectors[1];  
      delete [] gamma_vectors;


      CVector gamma_diff( gamma->dimension() );
      if ( non_informative_gamma )
      {
        ratio = glm_ratio;
      }
      else
      {
        gamma_diff = gamma->lastItemDrawn() - gamma->veryFirstMean();
        scale = M_D( gamma->veryFirstInverseCovariance(), gamma_diff );
        gamma_diff = gamma->initialMean() - gamma->veryFirstMean();
        ratio = 0.5 * ( M_D( gamma->veryFirstInverseCovariance(), gamma_diff ) - scale );

        if ( t_gamma )
        {
          ratio /= tau2_gamma->lastItemDrawn();
        }

        ratio += glm_ratio;
      }
    }
    else if ( !strcmp( distr_map[ index ], "alpha" ) )
    {
      ratio = 0;
    }
    else if ( !strcmp( distr_map[ index ], "tau2" ) )
    {
      ratio = 0;
    }
    else if ( !strcmp( distr_map[ index ], "tau2a" ) )
    {
      ratio = 0;
    }
    else if ( !strcmp( distr_map[ index ], "tau2gm" ) )
    {
      ratio = 0;
    }
    else if ( !strncmp( distr_map[ index ], "beta", 4 ) 
              || !strncmp( distr_map[ index ], "tau2b", 5 ) )
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

          if ( group_index >= 0 && group_index < number_of_groups )
	  {
            if ( !strncmp( distr_map[ index ], "beta", 4 ) && !glm_model[0]->gibbsDrawingsForCoefficients() )
            {
              DistributionParameter ** beta_vectors;
              beta_vectors = new DistributionParameter * [2];
#ifdef FIX1
              CVector tmpvec1 = beta[ group_index ]->lastItemDrawn();
              CVector tmpvec2 = beta[ group_index ]->initialMean();
              beta_vectors[0] = new DistributionParameter( tmpvec1 );
              beta_vectors[1] = new DistributionParameter( tmpvec2 );
#else
              beta_vectors[0] = new DistributionParameter( beta[ group_index ]->lastItemDrawn() );
              beta_vectors[1] = new DistributionParameter( beta[ group_index ]->initialMean() );
#endif
//              beta_vectors[0] = new DistributionParameter( beta[ group_index ]->lastItemDrawn() );
//              beta_vectors[1] = new DistributionParameter( beta[ group_index ]->initialMean() );

              // Bug fix by Dawn: the following was previously:
              // random_predictors[1]->parameter().Val(0,0) = group_index;
              // which fails to change the value of random_predictors[1].  The following is 
              // the corrected version.
              random_predictors[1]->updateParameter( (double) group_index );
              

              glm_ratio = glm_model[0]->logRatioTargetDensity( 2, random_predictors, beta_vectors );

              //printf( "glm_ratio[ %d ] = %f\n", index, glm_ratio ); fflush(stdout);

              delete beta_vectors[0];
              delete beta_vectors[1];  
              delete [] beta_vectors;

              if ( non_informative_beta ){
                ratio = glm_ratio;
              }
              else {
                CMatrix beta_inv_cov( beta[ group_index ]->dimension(), beta[ group_index ]->dimension() );
                CVector beta_diff( beta[ group_index ]->dimension() );
                beta_inv_cov = betaInverseCovariance( group_index );
                //Bug fix by Dawn: was:  
                // beta_diff = beta[ group_index ]->lastItemDrawn() - beta[ group_index ]->veryFirstMean();
                beta_diff = *(second_residuals[ group_index ]);
                scale = M_D( beta_inv_cov, beta_diff );
                //Bug fix by Dawn: was:  
                //beta_diff = beta[ group_index ]->initialMean() - beta[ group_index ]->veryFirstMean();
                beta_diff = beta_diff + beta[ group_index ]->lastItemDrawn();
                beta_diff = beta_diff - beta[ group_index ]->initialMean();
                ratio = 0.5 * ( scale - M_D( beta_inv_cov, beta_diff ) );
                ratio += glm_ratio;
                #ifdef DEBUGBETA
                  printf("beta_inv_cov = \n");  beta_inv_cov.Print(); 
                  printf("second_residuals = \n");  second_residuals[ group_index ]->Print();
                  printf("beta_diff = \n");  beta_diff.Print();
                  printf("beta = \n");  beta[ group_index ]->initialMean().Print();  
                  printf("beta star = \n");  beta[ group_index ]->lastItemDrawn().Print();  fflush(stdout);
                #endif
              }
            }
            else if ( !strncmp( distr_map[ index ], "tau2b", 5 ) )
	    {
              ratio = 0;
            }          
          }
          else
          {
            printf( " BayesianHierarchicalGlmModel::logRatioTargetDensity: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
          }          
        }//end if not null
        else
        {
          printf( " BayesianHierarchicalGlmModel::logRatioTargetDensity: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }
      else
      {
        printf( " BayesianHierarchicalGlmModel::logRatioTargetDensity: Wrong argument in [%s]. Number expected.\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if 
  }//end if index is valid
  else
  {
    printf( "BayesianHierarchicalGlmModel::logRatioTargetDensity: Variable index [%d] does not exist.\n", index );
  }

#ifdef DEBUG1
  printf("target density: ratio[%d] = %f\n", index, ratio ); fflush(stdout);
#endif

  return ( ratio );

}//end



void BayesianHierarchicalGlmModel::setCurrentVariableFromProposal( int index )
{
  //do nothing: keep the last draws as current draws
  if ( index < number_of_variables && index >= 0 )
  {
    if ( !strncmp( distr_map[ index ], "glm:", 4 ) ) 
    	// For the Poisson log-normal and binomial logit-normal models, this updates 
    	// either sigma^2, xi, or one of the lambdas
    {
      int glm_index;
      char * ptr, * temp_distr;

      temp_distr = new char [ strlen( distr_map[ index ] ) + 1 ];
      sprintf( temp_distr, "%s", distr_map[ index ] );
      ptr = strtok( temp_distr, ":" );
      if ( ptr != NULL )
      {
        ptr = strtok( NULL, ":" );
        if ( ptr != NULL )
	      {
          glm_index = atoi( ptr );
          glm_model[0]->setCurrentVariableFromProposal( glm_index );
        }
        else
        {
          printf( " BayesianHierarchicalGlmModel::setCurrentVariableFromProposal: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }//end if null
      else
      {
        printf( " BayesianHierarchicalGlmModel::setCurrentVariableFromProposal: Wrong argument in [%s].\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if glm
    else if ( !strcmp( distr_map[ index ], "gamma" ) )
    {
      updateLinearResponse();
#ifdef DEBUG1
      printf("Gamma draw accepted.  New value:\n");
      gamma->lastItemDrawn().Print(); fflush(stdout);
#endif
    }
    else if ( !strcmp( distr_map[ index ], "alpha" ) )
    {
      gibbsUpdateSecondStageResiduals();
#ifdef DEBUG1
      printf("New value of alpha is:\n" );  alpha->lastItemDrawn().Print();  fflush(stdout);
#endif
    }
    else if ( !strcmp( distr_map[ index ], "tau2" ) )
    {
      //nothing to do: this is is a Gibbs draw.
    }
    else if ( !strcmp( distr_map[ index ], "tau2a" ) )
    {
      //nothing to do: this is is a Gibbs draw.
    }
    else if ( !strcmp( distr_map[ index ], "tau2gm" ) )
    {
      //nothing to do: this is is a Gibbs draw.
    }
    else if ( !strncmp( distr_map[ index ], "beta", 4 ) 
              || !strncmp( distr_map[ index ], "tau2b", 5 ) )
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

          if ( group_index >= 0 && group_index < number_of_groups )
	  {
            if ( !strncmp( distr_map[ index ], "beta", 4 ) )
            {
              updateLinearResponse( group_index );
              gibbsUpdateSecondStageResiduals( group_index );
#ifdef DEBUG1
              printf("beta[%d] accepted.  New value: \n", group_index);
              beta[ group_index ]->lastItemDrawn().Print();
#endif
            }
            else if ( !strncmp( distr_map[ index ], "tau2b", 5 ) )
	    {
              //nothing to do: this is is a Gibbs draw.
            }          
          }
          else
          {
            printf( " BayesianHierarchicalGlmModel::setCurrentVariableFromProposal: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
          }          
        }//end if not null
        else
        {
          printf( " BayesianHierarchicalGlmModel::setCurrentVariableFromProposal: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }
      else
      {
        printf( " BayesianHierarchicalGlmModel::setCurrentVariableFromProposal: Wrong argument in [%s]. Number expected.\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if 
  }//end if index is valid
  else
  {
    printf( "BayesianHierarchicalGlmModel::setCurrentVariableFromProposal: Variable index [%d] does not exist.\n", index );
  }

}//end


void BayesianHierarchicalGlmModel::keepCurrentVariable( int index )
{
  if ( index < number_of_variables && index >= 0 )
  {
    if ( !strncmp( distr_map[ index ], "glm:", 4 ) )
    {
      int glm_index;
      char * ptr, * temp_distr;

      temp_distr = new char [ strlen( distr_map[ index ] ) + 1 ];
      sprintf( temp_distr, "%s", distr_map[ index ] );
      ptr = strtok( temp_distr, ":" );
      if ( ptr != NULL )
      {
        ptr = strtok( NULL, ":" );
        if ( ptr != NULL )
        {
          glm_index = atoi( ptr );
          glm_model[0]->keepCurrentVariable( glm_index );
        }
        else
        {
          printf( " BayesianHierarchicalGlmModel::keepCurrentVariable: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }//end if null
      else
      {
        printf( " BayesianHierarchicalGlmModel::keepCurrentVariable: Wrong argument in [%s].\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if glm
    else if ( !strcmp( distr_map[ index ], "gamma" ) )
    {
      gamma->setLastDraw( gamma->initialMean() );
#ifdef DEBUG1
      printf("Gamma draw rejected.  Keeping old value:\n");
      gamma->lastItemDrawn().Print(); fflush(stdout);
#endif
    }
    else if ( !strcmp( distr_map[ index ], "alpha" ) )
    {
      //this will not occur: this is is a Gibbs draw.
      printf( " BayesianHierarchicalGlmModel::keepCurrentVariable: Wrong call to [%s].\n", distr_map[ index ] );
    }
    else if ( !strcmp( distr_map[ index ], "tau2" ) )
    {
      //this will not occur: this is is a Gibbs draw.    
      printf( " BayesianHierarchicalGlmModel::keepCurrentVariable: Wrong call to [%s].\n", distr_map[ index ] );
    }
    else if ( !strcmp( distr_map[ index ], "tau2a" ) )
    {
      //this will not occur: this is is a Gibbs draw.
      printf( " BayesianHierarchicalGlmModel::keepCurrentVariable: Wrong call to [%s].\n", distr_map[ index ] );
    }
    else if ( !strcmp( distr_map[ index ], "tau2gm" ) )
    {
      //this will not occur: this is is a Gibbs draw.
      printf( " BayesianHierarchicalGlmModel::keepCurrentVariable: Wrong call to [%s].\n", distr_map[ index ] ); 
    }
    else if ( !strncmp( distr_map[ index ], "beta", 4 ) 
              || !strncmp( distr_map[ index ], "tau2b", 5 ) )
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

          if ( group_index >= 0 && group_index < number_of_groups )
	  {
            if ( !strncmp( distr_map[ index ], "beta", 4 ) )
            {
              beta[ group_index ]->setLastDraw( beta[ group_index ]->initialMean() );
#ifdef DEBUG1
              printf("Setting beta[%d] back to its previous value: \n", group_index);
              beta[ group_index ]->lastItemDrawn().Print();
#endif
            }
            else if ( !strncmp( distr_map[ index ], "tau2b", 5 ) )
	    {
              //this will not occur: this is is a Gibbs draw.
              printf( " BayesianHierarchicalGlmModel::keepCurrentVariable: Wrong call to [%s].\n", distr_map[ index ] ); 
            }          
          }
          else
          {
            printf( " BayesianHierarchicalGlmModel::keepCurrentVariable: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
          }          
        }//end if not null
        else
        {
          printf( " BayesianHierarchicalGlmModel::keepCurrentVariable: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }
      else
      {
        printf( " BayesianHierarchicalGlmModel::keepCurrentVariable: Wrong argument in [%s]. Number expected.\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if 
  }//end if index is valid
  else
  {
    printf( "BayesianHierarchicalGlmModel::keepCurrentVariable: Variable index [%d] does not exist.\n", index );
  }

}//end

