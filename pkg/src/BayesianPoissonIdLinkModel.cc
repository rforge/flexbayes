#include <math.h>
#include <string.h>
#include "S.h"

#include "Const.h"
#include "DistributionParameter.h"
#include "BayesianPoissonIdLinkModel.h"


/* constructor.
   returns: an empty Distribution of type BayesianPoissonIdLinkModel
*/
BayesianPoissonIdLinkModel::BayesianPoissonIdLinkModel()
{
  emptyModel();
}//end


void BayesianPoissonIdLinkModel::emptyModel()
{
  exposures = NULL;
  counts = NULL;
  mu_linear = NULL;
  residuals = NULL;
  sigma2_first_draw = NULL;
  log_lambda_previous = NULL;

  log_lambda = NULL;

  beta_dim = 0; 
  gamma_dim = 0;

  sigma2 = NULL;
  sigma2ICS = NULL;
  sigma2PNIP = NULL;
  sigma2NIP = NULL;

  common_sigma = false;
  invChisq_sigma2 = false;
  duMouchel_sigma2 = false;
  uniformShrinkage_sigma2 = false;
  properNonInfoPrior_sigma2 = false;
  known_sigma2 = false;
  nonInformativePower_sigma2 = false;
  
  number_of_observations = 0;
  number_of_variables = 0;
  number_of_groups = 0;

  number_of_simulations = 0;
  simulated_lambda = NULL;
  simulated_sigma2 = NULL;

  update_hessians = true;
  gibbs_for_coefs = true;

  //array of distributions map (maps distributions names to actual objects)
  distr_map = NULL;
  model_type = NULL;
}//end



BayesianPoissonIdLinkModel::~BayesianPoissonIdLinkModel()
{
  int i, j;

  if ( exposures != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete exposures[i];
    }
    delete [] exposures;
  }


  if ( mu_linear != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete mu_linear[i];
    }
    delete [] mu_linear;
  }

  if ( residuals != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete residuals[i];
    }
    delete [] residuals;
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


  if ( log_lambda_previous != NULL )
  {
    for ( j = 0; j < number_of_groups; j++ )
    {
      for ( i = 0; i < counts[j]->Len(); i++ )
      {
        delete log_lambda_previous[j][i];
      }
      delete [] log_lambda_previous[j];
    }
    delete [] log_lambda_previous;
  }

  if ( log_lambda != NULL )
  {
    for ( j = 0; j < number_of_groups; j++ )
    {
      for ( i = 0; i < counts[j]->Len(); i++ )
      {
        delete log_lambda[j][i];
      }
      delete [] log_lambda[j];
    }
    delete [] log_lambda;
  }


  if ( simulated_sigma2 != NULL )
  {
    delete simulated_sigma2;
  }

  if ( simulated_lambda != NULL )
  {
    for ( j = 0; j < number_of_groups; j++ )
    {
      for ( i = 0; i < counts[j]->Len(); i++ )
      {
        delete simulated_lambda[j][i];
      }
      delete [] simulated_lambda[j];
    }
    delete [] simulated_lambda;  
  }


  if ( counts != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete counts[i];
    }
    delete [] counts;
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

  if ( model_type != NULL )
  {
    delete [] model_type;
  }

}//end




void BayesianPoissonIdLinkModel::initialize( CVector ** count_v, CVector ** expos, int n_groups  )
{
	// count_v contains the outcome counts.  count_v[i] contains the counts for group i, for each
	// group (0, ..., #groups-1).  count_v[i]->Len() is the number of observations in group i.
	// Similar for expos.
  int i;

  number_of_groups = n_groups;
  number_of_observations = 0;

  counts = new CVector * [ number_of_groups ];
  exposures = new CVector * [ number_of_groups ];
  mu_linear = new CVector * [ number_of_groups ];
  residuals = new CVector * [ number_of_groups ];

  for ( i = 0; i < number_of_groups; i++ )
  {
    counts[i] = new CVector( (*count_v[i]) );
    exposures[i] = new CVector( (*expos[i]) );

    mu_linear[i] = new CVector( count_v[i]->Len() );
    residuals[i] = new CVector( count_v[i]->Len() );
    
    number_of_observations += exposures[i]->Len();
  }

  setLambdaPrior();
  setModelType( "Poisson" );

}//end



void BayesianPoissonIdLinkModel::initialize( CVector ** count_v, int n_groups )
{
  int i;

  number_of_groups = n_groups;
  number_of_observations = 0;  

  counts = new CVector * [ number_of_groups ];
  exposures = new CVector * [ number_of_groups ];
  mu_linear = new CVector * [ number_of_groups ];
  residuals = new CVector * [ number_of_groups ];

  for ( i = 0; i < number_of_groups; i++ )
  {
    counts[i] = new CVector( (*count_v[i]) );
    exposures[i] = new CVector( count_v[i]->Len() );
    exposures[i]->setTo( 1.0 );

    mu_linear[i] = new CVector( count_v[i]->Len() );
    residuals[i] = new CVector( count_v[i]->Len() );

    number_of_observations += count_v[i]->Len();
  }

  setLambdaPrior();
  setModelType( "Poisson" );

}//end



void BayesianPoissonIdLinkModel::setLambdaPrior()
{
  int i, j;

  log_lambda = new NormalDistribution ** [ number_of_groups ];
  log_lambda_previous = new DistributionParameter ** [ number_of_groups ];

  for ( j = 0; j < number_of_groups; j++ )
  {
    log_lambda[j] = new NormalDistribution * [ counts[j]->Len() ];
    log_lambda_previous[j] = new DistributionParameter * [ counts[j]->Len() ];
    for ( i = 0; i < counts[j]->Len(); i++ )
    {
      log_lambda[j][i] = new NormalDistribution( 1 );
      log_lambda_previous[j][i] = new DistributionParameter( 1, 1 );
    }
  }

  number_of_variables += number_of_observations;
}//end



void BayesianPoissonIdLinkModel::setModelType( char * type )
{
  if ( model_type != NULL )
  {
    if ( strcmp( model_type, type ) )
    {
      delete [] model_type;
      if ( !strcmp( type, "Poisson" ) )
      {
        model_type = new char [ 7 + 1 ];
        sprintf( model_type, "%s", "Poisson" );
      }
      else 
      {
        model_type = new char [ 8 + 1 ];
        sprintf( model_type, "%s", "Binomial" );
      }
    }
  }
  else 
  {
    if ( !strcmp( type, "Poisson" ) )
    {
      model_type = new char [ 7 + 1 ];
      sprintf( model_type, "%s", "Poisson" );
    }
    else 
    {
      model_type = new char [ 8 + 1 ];
      sprintf( model_type, "%s", "Binomial" );
    }
  }    

}//end


void BayesianPoissonIdLinkModel::sigma2PriorInvChisq( double p_nuSigma2, double p_sigma2 )
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


void BayesianPoissonIdLinkModel::sigma2PriorInvChisq( double * p_nuSigma2, double * p_sigma2 )
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
}//end


void BayesianPoissonIdLinkModel::sigma2PriorDuMouchel( double p_sigma2 )
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


void BayesianPoissonIdLinkModel::sigma2PriorDuMouchel( double * p_sigma2 )
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



void BayesianPoissonIdLinkModel::sigma2PriorUniformShrinkage( double p_sigma2 )
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


void BayesianPoissonIdLinkModel::sigma2PriorUniformShrinkage( double * p_sigma2 )
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



void BayesianPoissonIdLinkModel::sigma2PriorNonInformative( double p_power ) throw( rtErr )
{

  nonInformativePower_sigma2 = true;

  if ( common_sigma )
  {
    if ( p_power <= -0.5 * number_of_observations )
    {
      printf( "BayesianPoissonIdLinkModel::sigma2PriorNonInformative: Power [%f] in noninformative prior is not valid.\n", p_power );  
      char the_error[] = "BayesianPoissonIdLinkModel::sigma2PriorNonInformative: Power in noninformative prior is not valid.";
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
      if ( p_power <= -0.5 * counts[i]->Len() )
      {
        printf( "BayesianPoissonIdLinkModel::sigma2PriorNonInformative: Power [%f] in noninformative prior is not valid.\n", p_power );  
	char the_error[] = "BayesianPoissonIdLinkModel::sigma2PriorNonInformative: Power in noninformative prior is not valid.";
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


void BayesianPoissonIdLinkModel::sigma2PriorNonInformative( double * p_power ) throw( rtErr )
{

  nonInformativePower_sigma2 = true;

  if ( common_sigma )
  {
    if ( p_power[0] <= -0.5 * number_of_observations )
    {
      printf( "BayesianPoissonIdLinkModel::sigma2PriorNonInformative: Power [%f] in noninformative prior is not valid.\n", p_power );  
      char the_error[] = "BayesianPoissonIdLinkModel::sigma2PriorNonInformative: Power in noninformative prior is not valid.";
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
      if ( p_power[i] <= -0.5 * counts[i]->Len() )
      {
        printf( "BayesianPoissonIdLinkModel::sigma2PriorNonInformative: Power [%f] in noninformative prior is not valid.\n", p_power );  
        char the_error[] = "BayesianPoissonIdLinkModel::sigma2PriorNonInformative: Power in noninformative prior is not valid.";
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


void BayesianPoissonIdLinkModel::sigma2Known( double * p_sigma2 )
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


void BayesianPoissonIdLinkModel::sigma2Known( double p_sigma2 )
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



void BayesianPoissonIdLinkModel::samplerDefaultInitialPoint()
{
  int i;

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


void BayesianPoissonIdLinkModel::samplerSigma2InitialPoint( double init_sigma2 )
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


void BayesianPoissonIdLinkModel::samplerSigma2InitialPoint( CVector & init_sigma2 )
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



void BayesianPoissonIdLinkModel::setCoefficientDimension( int p_random, int p_fixed )
{
  beta_dim = p_random;
  gamma_dim = p_fixed;

}//end


void BayesianPoissonIdLinkModel::metropolisHastingsUpdateModel( DistributionParameter * par_vals )
{

  if ( par_vals->length() == 1 )
  {
    metropolisHastingsUpdateModel( ((int) par_vals->getScalar()) );
  }
  else
  {
#ifdef FIX1
    // The mu vector is passed to this function, and the first element simply indexes the group.  The
    // rest of the vector contains the mu values for each observation in the group
    CVector tmpvec = par_vals->getVector().subVector( 1, par_vals->length() - 1 );
    metropolisHastingsUpdateModel( ((int) par_vals->getVector().Val( 0 )), tmpvec );
#else
    metropolisHastingsUpdateModel( ((int) par_vals->getVector().Val( 0 )), par_vals->getVector().subVector( 1, par_vals->length() - 1 ) );
#endif
//    metropolisHastingsUpdateModel( ((int) par_vals->getVector().Val( 0 )), par_vals->getVector().subVector( 1, par_vals->length() - 1 ) );
  }
  // The mu vectors should have been updated in this function.  
}//end



void BayesianPoissonIdLinkModel::updateResiduals( int j, int i )
{
  residuals[j]->Val( i ) = ( log_lambda[j][i]->lastDraw().getScalar() - mu_linear[j]->Val(i) );

#ifdef DEBUGBETA
  printf("residuals[%d][%d] changed to %f\n", j, i, residuals[j]->Val( i ) );
#endif
}//end


void BayesianPoissonIdLinkModel::updateResiduals( int j )
{
  int i;

  for ( i = 0; i < mu_linear[ j ]->Len(); i++ )
  {
    residuals[j]->Val( i ) = ( log_lambda[j][i]->lastDraw().getScalar() - mu_linear[j]->Val(i) );
  }

#ifdef DEBUGBETA
  printf("residuals[%d] changed to:\n", j );  residuals[j]->Print(); fflush(stdout);
#endif
}//end


void BayesianPoissonIdLinkModel::updateResiduals()
{
  int j;

  for ( j = 0; j < number_of_groups; j++ )
  {
    updateResiduals( j );
  }

}//end



void BayesianPoissonIdLinkModel::gibbsUpdateSigma2() throw( rtErr )
{
  int i;
  double scale, local_scale, add_df;

  if ( common_sigma )
  {
    if ( invChisq_sigma2 || nonInformativePower_sigma2 || properNonInfoPrior_sigma2 )
    {
      scale = 0;
#ifdef DEBUGGLM
      printf("residuals = \n");
#endif
      for ( i = 0; i < number_of_groups; i++ )
      {    
        local_scale = ( *(residuals[i]) * ( *(residuals[i]) ) );
        scale += local_scale;
#ifdef DEBUGGLM
        residuals[i]->Print(); printf("local_scale = %f\n", local_scale); fflush(stdout);
#endif
      }//end loop

      add_df = number_of_observations;

      //printf("update sigma2: scale = %f, df = %f\n", scale, add_df ); fflush(stdout);

      if ( invChisq_sigma2 || nonInformativePower_sigma2 )
      {
        if ( nonInformativePower_sigma2 )
        {
          add_df += 2 * sigma2NIP[0];
        }

        sigma2ICS[0]->update( add_df, scale );
#ifdef DEBUGGLM
        printf("Calling sigma2ICS[0]->update( %f, %f )\n", add_df, scale ); fflush(stdout);
#endif
      }
      else 
      {
        sigma2PNIP[0]->setScale( scale );
#ifdef DEBUGGLM
        printf("Calling sigma2PNIP[0]->setScale( %f )\n", scale ); fflush(stdout);
#endif
      }
    }
  }
  else
  {
    printf( "BayesianPoissonIdLinkModel::gibbsUpdateSigma2: Trying to update a common variance variable when there are number of groups variances in the model.\n" );
    char the_error[] = "BayesianPoissonIdLinkModel::gibbsUpdateSigma2: Trying to update a common variance variable when there are number of groups variances in the model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}//end



void BayesianPoissonIdLinkModel::gibbsUpdateSigma2( int index ) throw( rtErr )
{
  double scale, add_df;

  if ( !common_sigma )
  {
    if ( invChisq_sigma2 || nonInformativePower_sigma2 || properNonInfoPrior_sigma2 )
    {
      scale = ( *(residuals[ index ]) * ( *(residuals[ index ]) ) );
      add_df = counts[ index ]->Len();

      if ( invChisq_sigma2 || nonInformativePower_sigma2 )
      {
        if ( nonInformativePower_sigma2 )
        {
          add_df += 2 * sigma2NIP[ index ];
        }

        sigma2ICS[ index ]->update( add_df, scale );
#ifdef DEBUGSIGMA2
        printf("Calling sigma2ICS[%d]->update( %f, %f )\n", index, add_df, scale ); fflush(stdout);
#endif
      }
      else 
      {
        sigma2PNIP[ index ]->setScale( scale );
#ifdef DEBUGSIGMA2
        printf("Calling sigma2PNIP[%d]->setScale( %f )\n", index, scale ); fflush(stdout);
#endif
      }
    }//end if
  }
  else
  {
    printf( "BayesianPoissonIdLinkModel::gibbsUpdateSigma2: Trying to update number of groups variances when there is a common variance in the model.\n" );
    char the_error[] = "BayesianPoissonIdLinkModel::gibbsUpdateSigma2: Trying to update number of groups variances when there is a common variance in the model.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}//end 

// update the residuals vector and the mu_linear (mean) vector
// j is the group index.  mu is the vector of mu values for each observation in the group.
void BayesianPoissonIdLinkModel::metropolisHastingsUpdateModel( int j, CVector & mu )
{
  //update beta or gamma
  int i;

  for ( i = 0; i < mu_linear[ j ]->Len(); i++ )
  {
    mu_linear[ j ]->Val( i ) = mu.Val( i );
    residuals[j]->Val( i ) = ( log_lambda[j][i]->lastDraw().getScalar() - mu.Val(i) );


    CVector lmean( 1 );
    CMatrix lcov( 1, 1 );

    lmean.Val(0) = mu.Val(i);

    if ( common_sigma )
    {
      lcov.Val( 0, 0 ) = sigma2->lastDraw().getScalar();
      //lmean.Val(0) = mu.Val(i) + sigma2->lastDraw().getScalar() * counts[j]->Val( i ); 
    }
    else
    {
      lcov.Val( 0, 0 ) = sigma2->lastDrawFromComponent( j ).getScalar();
      //lmean.Val(0) = mu.Val(i) + sigma2->lastDrawFromComponent( j ).getScalar() * counts[j]->Val( i ); 
    }

    // Following is a bug fix by Dawn:  log_lambda mean and covariance do not need to get set here.  They
    // are set in BayesianPoissonIdLinkModel::MetropolisHastingsUpdateLogLambda(...), which is called 
    // from BayesianPoissonIdLinkModel::updateVariableForProposal(index) when index is the variable index
    // for this lambda.
    // The problem with updating this here is that it resets log_lambda[j][i]->last.item.drawn to be 
    // equal to &lmean.
    //log_lambda[j][i]->setMean( &lmean );
    //log_lambda[j][i]->setCovariance( &lcov );
  }

#ifdef DEBUGBETA
  printf("\n\nresiduals[%d] changed to \n", j );
  residuals[j]->Print();
#endif
#ifdef DEBUGGAMMA
  printf("\n\nresiduals[%d] changed to \n", j );
  residuals[j]->Print();

  printf( "mu linear[%d] changed to \n", j );
  mu_linear[j]->Print();
#endif

}//end


// updates log_lambda proposal mean and covariance 
void BayesianPoissonIdLinkModel::metropolisHastingsUpdateModel( int j )
{
  //update sigma j
  int i;

  for ( i = 0; i < mu_linear[ j ]->Len(); i++ )
  {
    CVector lmean( 1 );
    CMatrix lcov( 1, 1 );

    lmean.Val( 0 ) = mu_linear[ j ]->Val( i );

    if ( common_sigma )
    {
      lcov.Val( 0, 0 ) = sigma2->lastDraw().getScalar();
      //lmean.Val( 0 ) = mu_linear[ j ]->Val( i ) + sigma2->lastDraw().getScalar() * counts[ j ]->Val( i ); 
    }
    else
    {
      lcov.Val( 0, 0 ) = sigma2->lastDrawFromComponent( j ).getScalar();
      //lmean.Val( 0 ) = mu_linear[ j ]->Val( i ) + sigma2->lastDrawFromComponent( j ).getScalar() * counts[ j ]->Val( i ); 
    }

    log_lambda[j][i]->setMean( &lmean );
    log_lambda[j][i]->setCovariance( &lcov );
  }      

}//end


void BayesianPoissonIdLinkModel::metropolisHastingsUpdateLogLambda( int g, int obs )
{
  // prepare lambda for drawing.  The proposal for each log lambda parameter is its 
  // prior distribution conditional on the other parameters.

  CVector lmean( 1 );
  CMatrix lcov( 1, 1 );

  lmean.Val( 0 ) = mu_linear[ g ]->Val( obs );

  if ( common_sigma )
  {
    lcov.Val( 0, 0 ) = sigma2->lastDraw().getScalar();
    //lmean.Val( 0 ) = mu_linear[ g ]->Val( obs ) + sigma2->lastDraw().getScalar() * counts[ g ]->Val( obs ); 
  }
  else
  {
    lcov.Val( 0, 0 ) = sigma2->lastDrawFromComponent( g ).getScalar();
    //lmean.Val( 0 ) = mu_linear[ g ]->Val( obs ) + sigma2->lastDrawFromComponent( g ).getScalar() * counts[ g ]->Val( obs ); 
  }


  //printf("update log lambda[%d, %d]: mean = and lcov = \n", g, obs );
  //lmean.Print();
  //lcov.Print();

#ifdef FIX1
  // FIX1 breaks up the computation into two steps.  Here there 
  // is an additional bug fix by Dawn: previously this
  // code used "setParameter" instead of "updateParameter", which
  // caused a memory leak
  CVector tmpvec = log_lambda[ g ][ obs ]->lastItemDrawn();
  log_lambda_previous[ g ][ obs ]->updateParameter( tmpvec );
#else
  log_lambda_previous[ g ][ obs ]->updateParameter( log_lambda[ g ][ obs ]->lastItemDrawn() );
#endif
//  log_lambda_previous[ g ][ obs ]->setParameter( log_lambda[ g ][ obs ]->lastItemDrawn() );
  log_lambda[ g ][ obs ]->setMean( &lmean );
  log_lambda[g ][ obs ]->setCovariance( &lcov );
  
#ifdef DEBUGGLM
  printf("Calling log_lambda_previous[%d][%d]->updateParameter with \n", g, obs);
  tmpvec.Print();
  printf("Calling log_lambda[%d][%d]->setMean with \n", g, obs);
  lmean.Print();
  printf("Calling log_lambda[%d][%d]->setCovariance with \n", g, obs);
  lcov.Print();
  printf("mu_linear[ g ]->Val( obs ) = %f\n", mu_linear[ g ]->Val( obs ));
  printf("log lambda old = %f\n", log_lambda_previous[ g ][ obs ]->getScalar());
#endif

}//end



void BayesianPoissonIdLinkModel::metropolisHastingsUpdateModel( int g, int obs )
{
  updateResiduals( g, obs );
}//end



CMatrix BayesianPoissonIdLinkModel::metropolisHastingsProposalHessian( int n_pars, DistributionParameter ** par_vals, DistributionParameter ** par_vector )
{
  CMatrix null_hessian( 1, 1 );
  null_hessian.setToZero();
  return ( null_hessian );
}//end



CVector BayesianPoissonIdLinkModel::metropolisHastingsMean( int n_pars, DistributionParameter ** par_vals, DistributionParameter ** par_vector )
{
  int j;
  double scale;

  if ( n_pars > 0 )
  {
    if ( par_vector != NULL )
    {
      CVector b_star( par_vector[0]->getVector() );

      if ( par_vals != NULL )
      {
        if ( par_vals[0]->getScalar() == 0 )
        {
          //beta
          j = ((int) par_vals[1]->getScalar());
          if ( common_sigma )
          {
            scale = sigma2->lastDraw().getScalar();
          }
          else
          {
            scale = sigma2->lastDrawFromComponent( j ).getScalar();
          }

          CVector theta( (*(residuals[j])) + (par_vals[ j + 2 ]->getMatrix() * b_star ) );
          CVector proposal_mean( par_vals[ j + 2 ]->getMatrix().T() * theta );
          proposal_mean.multiplyByScalar( 1.0 / scale );
#ifdef DEBUGBETA
          printf("residuals[%d] = \n", j);
          residuals[j]->Print();
          printf("current value of beta = \n");
          b_star.Print();
          printf("sigma2= %f\n", sigma2->lastDraw().getScalar());
          printf("random predictors matrix = \n"); par_vals[ j + 2 ]->getMatrix().Print();
          printf("proposal_mean = \n"); proposal_mean.Print(); fflush(stdout);
#endif       
          return (proposal_mean);
        }
        else
        {
          //gamma
          if ( n_pars == number_of_groups + 1 )
	  {
            if ( common_sigma )
            {
              scale = sigma2->lastDraw().getScalar();
            }
            else
            {
              scale = sigma2->lastDrawFromComponent( 0 ).getScalar();
            }
#ifdef DEBUGGAMMA
            printf("b_star = \n");
            b_star.Print();
            printf("residuals[0] = \n" );
            (*residuals[0]).Print();
            printf("par_vals[ 2 ]->getMatrix() = \n");
            par_vals[ 2 ]->getMatrix().Print();
#endif
            // Here par_vals[ 2 ]->getMatrix() is the fixed effect predictor matrix for the 
            // first group.  Should have # rows = # obs in group, # cols = # predictors
            CVector theta( (*(residuals[0])) + ( par_vals[ 2 ]->getMatrix() * b_star ) );
            CVector proposal_mean( par_vals[ 2 ]->getMatrix().T() * theta );
            proposal_mean.multiplyByScalar( 1.0 / scale );

            if ( number_of_groups > 1 )
            {
              for ( j = 1; j < number_of_groups; j++ )
              {   
                if ( !common_sigma )
                {
                  scale = sigma2->lastDrawFromComponent( j ).getScalar();
                }
                theta = (*(residuals[j])) + ( par_vals[ j + 2 ]->getMatrix() * b_star );
                theta.multiplyByScalar( 1.0 / scale );
#ifdef FIX1
                // FIX1 breaks up following computations into two steps.
                // Additional bugfix by Dawn: was par_vals[ j + 2 ]->getMatrix() * theta, which yields
                // an ill-defined matrix multiplication.  Change to
                // match previous computations in this function.
                CVector tmpvec = par_vals[ j + 2 ]->getMatrix().T() * theta;
#ifdef DEBUGGAMMA
                printf("proposal_mean = ");
                proposal_mean.Print();
                printf("tmpvec = ");
                tmpvec.Print();
#endif
                proposal_mean.add( tmpvec );
#else
                proposal_mean.add( par_vals[ j + 2 ]->getMatrix().T() * theta );
#endif
//                proposal_mean.add( par_vals[ j + 2 ]->getMatrix() * theta );
              }
            } 
  
            return ( proposal_mean );
          }
          else
	  {
            printf( " BayesianPoissonIdLinkModel::metropolisHastingsMean: Wrong number of parameters for Gamma [%d].\n", n_pars );
          }
        }//end gamma
      }//end if not null
      else
      {
        printf( " BayesianPoissonIdLinkModel::metropolisHastingsMean: Null array.\n" );
      }
    }//end if par_vector
    else
    {
      printf( " BayesianPoissonIdLinkModel::metropolisHastingsMean: Null array.\n" );
    }
  }//end if > 0
  else
  {
    printf( " BayesianPoissonIdLinkModel::metropolisHastingsMean: No parameters passed.\n" );
  }

  CVector null_mean( 1 );
  null_mean.setToZero();
  return ( null_mean );
}//end
  


CMatrix BayesianPoissonIdLinkModel::metropolisHastingsHessian( int n_pars, DistributionParameter ** par_vals )
{
  int j;
  double scale;

  if ( n_pars > 0 )
  {
    if ( par_vals != NULL )
    {
      if ( par_vals[0]->getScalar() == 0 )
      {
        //beta
        j = ((int) par_vals[1]->getScalar());
        if ( common_sigma )
        {
          scale = sigma2->lastDraw().getScalar();
        }
        else
        {
          scale = sigma2->lastDrawFromComponent( j ).getScalar();
        }

        CMatrix proposal_hessian( (-1.0 / scale ) * par_vals[ j + 2 ]->getMatrix().xTransposedX() );
     
        return (proposal_hessian);
      }
      else
      {
        //gamma
        if ( n_pars == number_of_groups + 1 )
	  {
          if ( common_sigma )
          {
            scale = sigma2->lastDraw().getScalar();
          }
          else
          {
            scale = sigma2->lastDrawFromComponent( 0 ).getScalar();
          }

          CMatrix proposal_hessian( ( 1.0 / scale ) * par_vals[ 2 ]->getMatrix().xTransposedX() );
          if ( number_of_groups > 1 )
          {
            for ( j = 1; j < number_of_groups; j++ )
            {   
              if ( !common_sigma )
              {
                scale = sigma2->lastDrawFromComponent( j ).getScalar();
              }
#ifdef FIX1
              CMatrix tmpmat = ( 1.0 / scale ) * par_vals[ j + 2 ]->getMatrix().xTransposedX();
              proposal_hessian.add( tmpmat );
#else
              proposal_hessian.add( ( 1.0 / scale ) * par_vals[ j + 2 ]->getMatrix().xTransposedX() );
#endif
//              proposal_hessian.add( ( 1.0 / scale ) * par_vals[ j + 2 ]->getMatrix().xTransposedX() );
            }
          } 
          proposal_hessian.multiplyByScalar( -1.0 );   

          return ( proposal_hessian );
        }
        else
	  {
          printf( " BayesianPoissonIdLinkModel::metropolisHastingsProposalHessian: Wrong number of parameters for Gamma [%d].\n", n_pars );
        }
      }//end gamma
    }//end if not null
    else
    {
      printf( " BayesianPoissonIdLinkModel::metropolisHastingsProposalHessian: Null array.\n" );
    }
  }//end if > 0
  else
  {
    printf( " BayesianPoissonIdLinkModel::metropolisHastingsProposalHessian: No parameters passed.\n" );
  }

  CMatrix null_hessian( 1, 1 );
  null_hessian.setToZero();
  return ( null_hessian );
}//end



void BayesianPoissonIdLinkModel::simulationsToArray( double * simul_output, int simulations_to_keep, int start_index )
{
  int i, j, k, s, start_dim, total_dim;

  start_dim = start_index;

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


  //now keep lambda's
  k = 0;
  for ( j = 0; j < number_of_groups; j++ )
  {
    for ( i = 0; i < counts[j]->Len(); i++ )
    {
      for ( s = 0; s < simulations_to_keep; s++ )
      {
        simul_output[ start_dim + k ] = simulated_lambda[j][i]->Val( s );
        k++;
      }
    }
  }

}//end


void BayesianPoissonIdLinkModel::keepSimulation( int simul_number )
{
  int i, j;
  double w;

  if ( simul_number < number_of_simulations )
  {
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

    //now keep lambda's
    for ( j = 0; j < number_of_groups; j++ )
    {
      for ( i = 0; i < counts[j]->Len(); i++ )
      {
        w = exp( log_lambda[j][i]->lastDraw().getScalar() );
        if ( !strcmp( model_type, "Poisson" ) )
          simulated_lambda[j][i]->Val( simul_number ) = w;
        else
          simulated_lambda[j][i]->Val( simul_number ) = w / ( 1 + w );
      }
    }
  }

}//end



void BayesianPoissonIdLinkModel::createOutput( int sim_to_keep )
{
  int i,j;

  number_of_simulations = sim_to_keep;

  if ( !known_sigma2 )
  {
    if ( common_sigma )
    {
      simulated_sigma2 = new CMatrix( sim_to_keep, 1 );
    }
    else
    {
      simulated_sigma2 = new CMatrix( sim_to_keep, number_of_groups );
    }
  }

  simulated_lambda = new CVector ** [ number_of_groups ];
  for ( j = 0; j < number_of_groups; j++ )
  {
    simulated_lambda[j] = new CVector * [ counts[j]->Len() ];
    for ( i = 0; i < counts[j]->Len(); i++ )
    {
      simulated_lambda[j][i] = new CVector( sim_to_keep );
    }
  }
}//end


void BayesianPoissonIdLinkModel::initializeTemporaryStructures()
{
  int i,j, k, len_digits_vars, len_digits_obs;

  len_digits_vars = ( (int) ( log( (double) number_of_variables ) / log( 10.0 ) ) ) + 1;
  len_digits_obs = ( (int) ( log( (double) number_of_observations ) / log( 10.0 ) ) ) + 1;

  distr_map = new char * [ number_of_variables ];
  k = 0;

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

  //now lambda's
  for ( j = 0; j < number_of_groups; j++ )
  {
    for ( i = 0; i < counts[j]->Len(); i++ )
    {
      distr_map[k] = new char [ 7 + 2 + 2 * len_digits_obs ];
      if ( !strcmp( model_type, "Poisson" ) )
      {
        sprintf( distr_map[k], "lambda:%d:%d", j, i );
      }
      else
      {
        //Binomial
        sprintf( distr_map[k], "theta:%d:%d", j, i );
      }
      k++;
    }
  }

  //printf( "BHPM: total vars = %d.  Number vars = %d. The map is:\n", k, number_of_variables  );
  //for ( i = 0; i < k; i++ )
  //{
  //printf("map[%d] = %s\n", i, distr_map[ i ] );
  //}
  //fflush( stdout );


}//end



void BayesianPoissonIdLinkModel::dataAugmentationInitialDraws()
{
  int i, j;

  #ifdef DEBUG1
    printf( "PoissonIDModel: Sampling initial lambdas.\n" ); fflush( stdout );
  #endif

  //draw initial lambda's
  for ( j = 0; j < number_of_groups; j++ )
  {
    for ( i = 0; i < counts[j]->Len(); i++ )
    {
    	//Bug fix by Dawn: The following code was moved from 
    	// the function metropolisHastingsUpdateModel because
    	// it only needs to be called at the beginning of the 
    	// MCMC, not every iteration after updating the regression
    	// coefficients.
    	CVector lmean( 1 );
      CMatrix lcov( 1, 1 );

      lmean.Val(0) = mu_linear[ j ]->Val( i );

      if ( common_sigma )
      {
        lcov.Val( 0, 0 ) = sigma2->lastDraw().getScalar();
      }
      else
      {
        lcov.Val( 0, 0 ) = sigma2->lastDrawFromComponent( j ).getScalar();
      }
      #ifdef DEBUG_LAMBDA
        printf( "PoissonIDModel: log_lambda variance = %f \n", lcov.Val(0,0) ); fflush( stdout );
      #endif

      log_lambda[j][i]->setMean( &lmean );
      log_lambda[j][i]->setCovariance( &lcov );
      // end bug fix by Dawn
      
      log_lambda[j][i]->draw();
#ifdef FIX1
      CVector tmpvec = log_lambda[j][i]->lastItemDrawn();
      log_lambda_previous[j][i]->setParameter( tmpvec );
#else
      log_lambda_previous[j][i]->setParameter( log_lambda[j][i]->lastItemDrawn() );
#endif
    }
  }
  #ifdef DEBUG1
    printf( "PoissonIDModel: Finished sampling initial lambdas\n" ); fflush( stdout );
  #endif

  updateResiduals();

}//end



void BayesianPoissonIdLinkModel::drawVariableFromProposal( int index )
{
  if ( index < number_of_variables && index >= 0 )
  {
    if ( !strcmp( distr_map[ index ], "sigma" ) )
    {
      sigma2->draw();
#ifdef DEBUGGLM
      printf("Calling sigma2->draw()\n"); 
      printf("sigma2->mean() = %f\n", sigma2->mean().getScalar());
      printf("sigma2->variance() = %f\n", sigma2->variance().getScalar());
      fflush(stdout); 
#endif
    }
    else if ( !strncmp( distr_map[ index ], "sigma:", 6 )  ||
              !strncmp( distr_map[ index ], "lambda", 6 )  ||
              !strncmp( distr_map[ index ], "theta", 5 ) )
    {
      int group_index, obs;
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
	          if ( !strncmp(  distr_map[ index ], "sigma", 5 ) )
	          {
              sigma2->drawFromComponent( group_index );
#ifdef DEBUGGLM
              printf("Calling sigma2->drawFromComponent( %d )\n", group_index); fflush(stdout); 
#endif
            }
            else if ( !strncmp( distr_map[ index ], "lambda", 6 ) || !strncmp( distr_map[ index ], "theta", 5 ) )
	    {
              //get observation number  
              ptr = strtok( NULL, ":" );

              if ( ptr != NULL )
	      {
                obs = atoi( ptr );
                if ( obs >= 0 && obs < number_of_observations )
		{
                  log_lambda[ group_index ][ obs ]->draw();
#ifdef DEBUGGLM
                  printf("Calling log_lambda[%d][%d]->draw()\n", group_index, obs); fflush(stdout); 
#endif
                }
                else
		{
                  printf( "BayesianPoissonIdLinkModel::drawVariableFromProposal: Wrong argument in [%s]. Observation number [%d] out of range.\n", distr_map[ index ], obs );
                }
              }
              else
              {
                printf( "BayesianPoissonIdLinkModel::drawVariableFromProposal: Wrong argument in [%s].\n", distr_map[ index ] );
              }
            }//end if lambda
            else
	    {
              printf( "BayesianPoissonIdLinkModel::drawVariableFromProposal: Unknown variable in [%s].\n", distr_map[ index ] );
            }
          }//end if group index
          else
	  {
            printf( "BayesianPoissonIdLinkModel::drawVariableFromProposal: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
          }
        }//end if null
        else
        {
          printf( "BayesianPoissonIdLinkModel::drawVariableFromProposal: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }//end if not null
      else
      {
        printf( "BayesianPoissonIdLinkModel::drawVariableFromProposal: Wrong argument in [%s].\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if valid variable
    else
    {
      printf( " BayesianPoissonIdLinkModel::drawVariableFromProposal: Wrong argument index [%d].\n", index );
    }
  }
  else
  {
    printf( "BayesianPoissonIdLinkModel::drawVariableFromProposal: Index [%d] out of range.\n", index );
  }
}//end



void BayesianPoissonIdLinkModel::updateVariableForProposal( int index )
{
  if ( index < number_of_variables && index >= 0 )
  {
    if ( !strcmp( distr_map[ index ], "sigma" ) )
    {
      gibbsUpdateSigma2();
    }
    else if ( !strncmp( distr_map[ index ], "sigma:", 6 )  ||
              !strncmp( distr_map[ index ], "lambda", 6 )  ||
              !strncmp( distr_map[ index ], "theta", 5 ) )
    {
      int group_index, obs;
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
	    if ( !strncmp(  distr_map[ index ], "sigma", 5 ) )
	    {
              gibbsUpdateSigma2( group_index );
            }
            else if ( !strncmp( distr_map[ index ], "lambda", 6 ) || !strncmp( distr_map[ index ], "theta", 5 ))
	    {
              //get observation number  
              ptr = strtok( NULL, ":" );
              if ( ptr != NULL )
	      {
                obs = atoi( ptr );
                if ( obs >= 0 && obs < number_of_observations )
		{
                  metropolisHastingsUpdateLogLambda( group_index, obs );
                }
                else
		{
                  printf( "BayesianPoissonIdLinkModel::updateVariableForProposal: Wrong argument in [%s]. Observation number [%d] out of range.\n", distr_map[ index ], obs );
                }
              }
              else
              {
                printf( "BayesianPoissonIdLinkModel::updateVariableForProposal: Wrong argument in [%s].\n", distr_map[ index ] );
              }
            }//end if lambda
            else
	    {
              printf( "BayesianPoissonIdLinkModel::updateVariableForProposal: Unknown variable in [%s].\n", distr_map[ index ] );
            }
          }//end if group index
          else
	  {
            printf( "BayesianPoissonIdLinkModel::updateVariableForProposal: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
          }
        }//end if null
        else
        {
          printf( "BayesianPoissonIdLinkModel::updateVariableForProposal: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }//end if not null
      else
      {
        printf( "BayesianPoissonIdLinkModel::updateVariableForProposal: Wrong argument in [%s].\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if valid variable
    else
    {
      printf( " BayesianPoissonIdLinkModel::updateVariableForProposal: Wrong argument index [%d].\n", index );
    }
  }
  else
  {
    printf( "BayesianPoissonIdLinkModel::updateVariableForProposal: Index [%d] out of range.\n", index );
  }
}//end



void BayesianPoissonIdLinkModel::setCurrentVariableFromProposal( int index )
{

  if ( index < number_of_variables && index >= 0 )
  {
    if ( !strcmp( distr_map[ index ], "sigma" ) )
    {
      //do nothing
#ifdef DEBUG1
      printf("New value of sigma^2 is %f\n ", sigma2->lastDraw().getScalar() ); fflush(stdout);
#endif
    }
    else if ( !strncmp( distr_map[ index ], "sigma:", 6 )  ||
              !strncmp( distr_map[ index ], "lambda", 6 )  ||
              !strncmp( distr_map[ index ], "theta", 5 ) )
    {
      int group_index, obs;
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
            if ( !strncmp(  distr_map[ index ], "sigma", 5 ) )
            {
              //do nothing
            }
            else if ( !strncmp( distr_map[ index ], "lambda", 6 ) || !strncmp( distr_map[ index ], "theta", 5 ) )
            {
              //get observation number  
              ptr = strtok( NULL, ":" );
              if ( ptr != NULL )
              {
                obs = atoi( ptr );
                if ( obs >= 0 && obs < number_of_observations )
                {
                  //just update the residuals
                  metropolisHastingsUpdateModel( group_index, obs );
#ifdef DEBUG1
                  printf("Updating lambda [%d][%d]\n", group_index, obs);
                  printf("New log lambda = %f \n", log_lambda[ group_index ][ obs ]->lastDraw().getScalar());
#endif
                }
                else
                {
                  printf( "BayesianPoissonIdLinkModel::setCurrentVariableFromProposal: Wrong argument in [%s]. Observation number [%d] out of range.\n", distr_map[ index ], obs );
                }
              }
              else
              {
                printf( "BayesianPoissonIdLinkModel::setCurrentVariableFromProposal: Wrong argument in [%s].\n", distr_map[ index ] );
              }
            }//end if lambda
            else
	    {
              printf( "BayesianPoissonIdLinkModel::setCurrentVariableFromProposal: Unknown variable in [%s].\n", distr_map[ index ] );
            }
          }//end if group index
          else
	  {
            printf( "BayesianPoissonIdLinkModel::setCurrentVariableFromProposal: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
          }
        }//end if null
        else
        {
          printf( "BayesianPoissonIdLinkModel::setCurrentVariableFromProposal: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }//end if not null
      else
      {
        printf( "BayesianPoissonIdLinkModel::setCurrentVariableFromProposal: Wrong argument in [%s].\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if valid variable
    else
    {
      printf( " BayesianPoissonIdLinkModel::setCurrentVariableFromProposal: Wrong argument index [%d].\n", index );
    }
  }
  else
  {
    printf( "BayesianPoissonIdLinkModel::setCurrentVariableFromProposal: Index [%d] out of range.\n", index );
  }
}//end



void BayesianPoissonIdLinkModel::keepCurrentVariable( int index )
{
  if ( index < number_of_variables && index >= 0 )
  {
    if ( !strcmp( distr_map[ index ], "sigma" ) )
    {
      //do nothing; this is a Gibbs sampling
#ifdef DEBUG1 
      printf("Rejecting proposed value of sigma^2, which is wrong since it was a Gibbs step.\n"); 
      fflush(stdout);
#endif
    }
    else if ( !strncmp( distr_map[ index ], "sigma:", 6 )  ||
              !strncmp( distr_map[ index ], "lambda", 6 )  ||
              !strncmp( distr_map[ index ], "theta", 5 ) )
    {
      int group_index, obs;
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
	    if ( !strncmp(  distr_map[ index ], "sigma", 5 ) )
	    {
              //do nothing; this is a Gibbs sampling
            }
            else if ( !strncmp( distr_map[ index ], "lambda", 6 ) || !strncmp( distr_map[ index ], "theta", 5 ) )
	    {
              //get observation number  
              ptr = strtok( NULL, ":" );
              if ( ptr != NULL )
	      {
                obs = atoi( ptr );
                if ( obs >= 0 && obs < number_of_observations )
		{
			            // set log_lambda back to the previous value
                  log_lambda[ group_index ][ obs ]->setLastDraw( (*(log_lambda_previous[ group_index ][ obs ])) );
#ifdef DEBUG1
                  printf("Setting log_lambda[%d][%d] back to previous value, %f\n", group_index, obs, 
                    log_lambda_previous[ group_index ][ obs ]->getScalar() );
#endif
                }
                else
		{
                  printf( "BayesianPoissonIdLinkModel::keepCurrentVariable: Wrong argument in [%s]. Observation number [%d] out of range.\n", distr_map[ index ], obs );
                }
              }
              else
              {
                printf( "BayesianPoissonIdLinkModel::keepCurrentVariable: Wrong argument in [%s].\n", distr_map[ index ] );
              }
            }//end if lambda
            else
	    {
              printf( "BayesianPoissonIdLinkModel::keepCurrentVariable: Unknown variable in [%s].\n", distr_map[ index ] );
            }
          }//end if group index
          else
	  {
            printf( "BayesianPoissonIdLinkModel::keepCurrentVariable: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
          }
        }//end if null
        else
        {
          printf( "BayesianPoissonIdLinkModel::keepCurrentVariable: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }//end if not null
      else
      {
        printf( "BayesianPoissonIdLinkModel::keepCurrentVariable: Wrong argument in [%s].\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if valid variable
    else
    {
      printf( " BayesianPoissonIdLinkModel::keepCurrentVariable: Wrong argument index [%d].\n", index );
    }
  }
  else
  {
    printf( "BayesianPoissonIdLinkModel::keepCurrentVariable: Index [%d] out of range.\n", index );
  }
}//end



double BayesianPoissonIdLinkModel::logRatioProposal( int index )
{
  double ratio, diff, sum, scale;

  ratio = 0;

  if ( index < number_of_variables && index >= 0 )
  {
    if ( !strcmp( distr_map[ index ], "sigma" ) )
    {
      //do nothing; this is a Gibbs sampling
    }
    else if ( !strncmp( distr_map[ index ], "sigma:", 6 )  ||
              !strncmp( distr_map[ index ], "lambda", 6 )  ||
              !strncmp( distr_map[ index ], "theta", 5 ) )
    {
      int group_index, obs;
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
	          if ( !strncmp(  distr_map[ index ], "sigma", 5 ) )
            {
              //do nothing; this is a Gibbs sampling
            }
            else if ( !strncmp( distr_map[ index ], "lambda", 6 ) || !strncmp( distr_map[ index ], "theta", 5 ) )
            {
              //get observation number  
              ptr = strtok( NULL, ":" );
              if ( ptr != NULL )
              {
                obs = atoi( ptr );
                if ( obs >= 0 && obs < number_of_observations )
                {
                 
                  if ( common_sigma )
                  {
                    scale = sigma2->lastDraw().getScalar();
                  }
                  else
                  {
                    scale = sigma2->lastDrawFromComponent( group_index ).getScalar();
                  }

                  diff = log_lambda[ group_index ][ obs ]->lastDraw().getScalar() - log_lambda_previous[ group_index ][ obs ]->getScalar();
                  sum = log_lambda[ group_index ][ obs ]->lastDraw().getScalar() + log_lambda_previous[ group_index ][ obs ]->getScalar();
                  //ratio = (-0.5) * diff * ( sum - 2 * log_lambda[ group_index ][ obs ]->mean().getScalar() );
                  ratio = (-0.5) * diff * ( sum - 2 * mu_linear[ group_index ]->Val( obs ) );
                  ratio /= scale;
#ifdef DEBUGGLM
  printf("group=%d, obs = %d\n", group_index, obs);
  printf("mu = %f, sigma^2 = %f\n", mu_linear[ group_index ]->Val( obs ), scale);
  printf("log lambda new %f and old %f \n", 
    log_lambda[ group_index ][ obs ]->lastDraw().getScalar(), 
    log_lambda_previous[ group_index ][ obs ]->getScalar());
  printf("log proposal ratio = %f\n", ratio);
#endif
                }
                else
		{
                  printf( "BayesianPoissonIdLinkModel::logRatioProposal: Wrong argument in [%s]. Observation number [%d] out of range.\n", distr_map[ index ], obs );
                }
              }
              else
              {
                printf( "BayesianPoissonIdLinkModel::logRatioProposal: Wrong argument in [%s].\n", distr_map[ index ] );
              }
            }//end if lambda
            else
	    {
              printf( "BayesianPoissonIdLinkModel::logRatioProposal: Unknown variable in [%s].\n", distr_map[ index ] );
            }
          }//end if group index
          else
	  {
            printf( "BayesianPoissonIdLinkModel::logRatioProposal: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
          }
        }//end if null
        else
        {
          printf( "BayesianPoissonIdLinkModel::logRatioProposal: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }//end if not null
      else
      {
        printf( "BayesianPoissonIdLinkModel::logRatioProposal: Wrong argument in [%s].\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if valid variable
    else
    {
      printf( " BayesianPoissonIdLinkModel::logRatioProposal: Wrong argument index [%d].\n", index );
    }
  }
  else
  {
    printf( "BayesianPoissonIdLinkModel::logRatioProposal: Index [%d] out of range.\n", index );
  }

  return ( ratio );
}//end


double BayesianPoissonIdLinkModel::logRatioTargetDensity( int index )
{
  double ratio, diff, sum, scale, mean_shift;

  ratio = 0;

  if ( index < number_of_variables && index >= 0 )
  {
    if ( !strcmp( distr_map[ index ], "sigma" ) )
    {
      //do nothing; this is a Gibbs sampling
    }
    else if ( !strncmp( distr_map[ index ], "sigma:", 6 )  ||
              !strncmp( distr_map[ index ], "lambda", 6 )  ||
              !strncmp( distr_map[ index ], "theta", 5 ) )
    {
      int group_index, obs;
      double wx, wy;
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
	    if ( !strncmp(  distr_map[ index ], "sigma", 5 ) )
	    {
              //do nothing; this is a Gibbs sampling
            }
            else if ( !strncmp( distr_map[ index ], "lambda", 6 ) || !strncmp( distr_map[ index ], "theta", 5 ) )
	    {
              //get observation number  
              ptr = strtok( NULL, ":" );
              if ( ptr != NULL )
	      {
                obs = atoi( ptr );
                if ( obs >= 0 && obs < number_of_observations )
		{
                 
                  if ( common_sigma )
                  {
                    scale = sigma2->lastDraw().getScalar();
                  }
                  else
                  {
                    scale = sigma2->lastDrawFromComponent( group_index ).getScalar();
                  }

                  diff = log_lambda[ group_index ][ obs ]->lastDraw().getScalar() - log_lambda_previous[ group_index ][ obs ]->getScalar();
                  sum = log_lambda[ group_index ][ obs ]->lastDraw().getScalar() + log_lambda_previous[ group_index ][ obs ]->getScalar();
                  //ratio = (-0.5) * diff * ( sum - 2 * log_lambda[ group_index ][ obs ]->mean().getScalar() );
                  mean_shift = mu_linear[ group_index ]->Val( obs ) + scale * counts[ group_index ]->Val( obs );
                  ratio = (-0.5) * diff * ( sum - 2 * mean_shift );
                  ratio /= scale;
 
                  if ( !strcmp( model_type, "Poisson" ) )
		              {
                    ratio -= ( exposures[ group_index ]->Val( obs ) * ( exp( log_lambda[ group_index ][ obs ]->lastDraw().getScalar() ) - exp( log_lambda_previous[ group_index ][ obs ]->getScalar() ) ) );
                  }
                  else
		  {
                    //Binomial
                    wy = 1 + exp( log_lambda[ group_index ][ obs ]->lastDraw().getScalar() );

                    wx = 1 + exp( log_lambda_previous[ group_index ][ obs ]->getScalar() );

                    ratio -= ( exposures[ group_index ]->Val( obs ) * log( wy / wx ) );
                  }
#ifdef DEBUGGLM
  printf("log target ratio = %f\n", ratio);
#endif
                }
                else
		{
                  printf( "BayesianPoissonIdLinkModel::logRatioTargetDensity: Wrong argument in [%s]. Observation number [%d] out of range.\n", distr_map[ index ], obs );
                }
              }
              else
              {
                printf( "BayesianPoissonIdLinkModel::logRatioTargetDensity: Wrong argument in [%s].\n", distr_map[ index ] );
              }
            }//end if lambda
            else
	    {
              printf( "BayesianPoissonIdLinkModel::logRatioTargetDensity: Unknown variable in [%s].\n", distr_map[ index ] );
            }
          }//end if group index
          else
	  {
            printf( "BayesianPoissonIdLinkModel::logRatioTargetDensity: Wrong argument in [%s]. Index [%d]out of range.\n", distr_map[ index ], group_index );
          }
        }//end if null
        else
        {
          printf( "BayesianPoissonIdLinkModel::logRatioTargetDensity: Wrong argument in [%s].\n", distr_map[ index ] );
        }
      }//end if not null
      else
      {
        printf( "BayesianPoissonIdLinkModel::logRatioTargetDensity: Wrong argument in [%s].\n", distr_map[ index ] );
      }
      delete [] temp_distr;
    }//end if valid variable
    else
    {
      printf( " BayesianPoissonIdLinkModel::logRatioTargetDensity: Wrong argument index [%d].\n", index );
    }
  }
  else
  {
    printf( "BayesianPoissonIdLinkModel::logRatioTargetDensity: Index [%d] out of range.\n", index );
  }

  return ( ratio );
}//end



