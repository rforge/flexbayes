#include "R.h"
#include "Rmath.h"

#include "Const.h"
#include "DistributionParameter.h"
#include "BayesianPoissonLogLinkModel.h"


/* constructor.
   returns: an empty Distribution of type BayesianPoissonLogLinkModel
*/
BayesianPoissonLogLinkModel::BayesianPoissonLogLinkModel()
{
  emptyModel();
}//end


void BayesianPoissonLogLinkModel::emptyModel()
{
  exposures = NULL;
  counts = NULL;
  mu_linear = NULL;
  xi_first_draw = NULL;

  log_xi = NULL;
  xi_z0 = NULL;
  xi_type = 2;   // default is no overdispersion parameter

  beta_dim = 0; 
  gamma_dim = 0;

  number_of_observations = 0;
  number_of_variables = 0;
  number_of_groups = 0;

  max_counts = NULL;
  mean_counts = NULL;
  n_equal = NULL;
  n_larger_equal = NULL;
  log_xi_hessian = NULL;
  m0 = NULL;
  b = NULL;
  b_denom = NULL;
  d = NULL;
 
  number_of_simulations = 0;
  simulated_xi = NULL;
  simulated_lambda = NULL;

  computed_starting_hessian_logXi = NULL;
  computed_starting_hessian_beta = NULL;
  computed_starting_hessian_gamma = 0;

  hessian_beta = NULL;
  hessian_gamma = NULL;
  update_hessians = true;

  gibbs_for_coefs = false;

  computed_starting_logXi_variance = NULL;
  logXi_variance = NULL;

}//end


BayesianPoissonLogLinkModel::~BayesianPoissonLogLinkModel()
{
  int i;

  if ( exposures != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete exposures[i];
    }
    delete [] exposures;
  }

  if ( counts != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete counts[i];
    }
    delete [] counts;
  }

  if ( mu_linear != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete mu_linear[i];
    }
    delete [] mu_linear;
  }

  if ( xi_first_draw != NULL )
  {
    if ( xi_type == 1 )   // single overdispersion param.
    {
      delete xi_first_draw[0];
    }
    else if( xi_type == 0 )  // one overdispersion param. per group
    {
      for ( i = 0; i < number_of_groups; i++ )
      {
        delete xi_first_draw[i];
      }
    }
    delete [] xi_first_draw;
  }

  if ( log_xi != NULL )
  {
    if ( xi_type == 1 )   // single overdispersion param.
    {
      delete log_xi[0];
    }
    else if( xi_type == 0 )   // one overdispersion param. per group
    {
      for ( i = 0; i < number_of_groups; i++ )
      {
        delete log_xi[i];
      }
    }
    delete [] log_xi;
  }

  if ( max_counts != NULL )
  {
    delete max_counts;
  }

  if ( mean_counts != NULL )
  {
    delete mean_counts;
  }

  if ( n_equal != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete n_equal[i];
    }
    delete [] n_equal;
  }

  if ( n_larger_equal != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete n_larger_equal[i];
    }
    delete [] n_larger_equal;
  }

  if ( log_xi_hessian != NULL )
  {
    delete log_xi_hessian;
  }

  if ( m0 != NULL )
  {
    delete m0;
  }

  if ( b != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete b[i];
    }
    delete [] b;
  }

  if ( b_denom != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete b_denom[i];
    }
    delete [] b_denom;
  }

  if ( d != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete d[i];
    }
    delete [] d;
  }

  if ( simulated_xi != NULL )
  {
    delete simulated_xi;
  }
     
  if ( simulated_lambda != NULL )
  {
    delete simulated_lambda;
  } 

  if ( computed_starting_hessian_logXi != NULL )
  {
    delete computed_starting_hessian_logXi;
  }

  if ( computed_starting_hessian_beta != NULL )
  {
    delete computed_starting_hessian_beta;
  }

  if ( hessian_beta != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete hessian_beta[i];
    }
    delete [] hessian_beta;
  }

  if ( hessian_gamma != NULL )
  {
    delete hessian_gamma;
  }


  if ( computed_starting_logXi_variance != NULL )
  {
    delete computed_starting_logXi_variance;
  }

  if ( logXi_variance != NULL )
  {
    delete logXi_variance;
  }
  if ( xi_z0 != NULL)
  {
    // should be from new double[] so we don't need to delete its elements
    delete [] xi_z0 ;
  }

}//end


void BayesianPoissonLogLinkModel::initialize( CVector ** count_v, CVector ** expos, int n_groups  )
{
  int i;

  number_of_groups = n_groups;
  number_of_observations = 0;

  counts = new CVector * [ number_of_groups ];
  exposures = new CVector * [ number_of_groups ];
  mu_linear = new CVector * [ number_of_groups ];

  for ( i = 0; i < number_of_groups; i++ )
  {
    counts[i] = new CVector( (*count_v[i]) );
    exposures[i] = new CVector( (*expos[i]) );

    mu_linear[i] = new CVector( count_v[i]->Len() );

    number_of_observations += exposures[i]->Len();
  }

	// The following creates objects that are large enough for the case of 
	// separate xi parameters for each group.  In the case that there is a
	// single xi parameter, only the first element of this vector will ever
	// be used, and if there is no overdispersion parameter then these
	// vectors will never be used.
  computed_starting_hessian_logXi = new CVector( number_of_groups );
  computed_starting_hessian_logXi->setToZero();

  computed_starting_logXi_variance = new CVector( number_of_groups );
  computed_starting_logXi_variance->setToZero();

  logXi_variance = new CVector( number_of_groups );
  logXi_variance->setTo( 1.0 );

}//end


void BayesianPoissonLogLinkModel::initialize( CVector ** count_v, int n_groups )
{
  int i;

  number_of_groups = n_groups;
  number_of_observations = 0;  

  counts = new CVector * [ number_of_groups ];
  exposures = new CVector * [ number_of_groups ];
  mu_linear = new CVector * [ number_of_groups ];

  for ( i = 0; i < number_of_groups; i++ )
  {
    counts[i] = new CVector( (*count_v[i]) );
    exposures[i] = new CVector( count_v[i]->Len() );
    exposures[i]->setTo( 1.0 );

    mu_linear[i] = new CVector( count_v[i]->Len() );

    number_of_observations += count_v[i]->Len();
  }

	// The following creates objects that are large enough for the case of 
	// separate xi parameters for each group.  In the case that there is a
	// single xi parameter, only the first element of these vectors may ever
	// be used, and if there is no overdispersion parameter then these
	// vectors will never be used.
  computed_starting_hessian_logXi = new CVector( number_of_groups );
  computed_starting_hessian_logXi->setToZero();

  computed_starting_logXi_variance = new CVector( number_of_groups );
  computed_starting_logXi_variance->setToZero();

  logXi_variance = new CVector( number_of_groups );
  logXi_variance->setTo( 1.0 );

}//end


void BayesianPoissonLogLinkModel::initializeCountStorage()
{
  int i;

  max_counts = new CVector( number_of_groups );
  mean_counts = new CVector( number_of_groups );
  n_equal = new CVector * [ number_of_groups ];
  n_larger_equal = new CVector * [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    n_equal[i] = NULL;
    n_larger_equal[i] = NULL;
  }


  //initialize working structures
  log_xi_hessian = new CVector( number_of_groups );
  log_xi_hessian->setToZero();

  m0 = new CVector( number_of_groups );

  d = new CVector * [ number_of_groups ];
  b = new CVector * [ number_of_groups ];
  b_denom = new CVector * [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    d[i] = new CVector( exposures[i]->Len() );
    b[i] = new CVector( exposures[i]->Len() );
    b_denom[i] = new CVector( exposures[i]->Len() );
  }  

}//end



void BayesianPoissonLogLinkModel::logXiPrior( double p_xi )
{
  int i;
#ifdef DEBUG1
  printf( "p_xi = %f\n", p_xi );  fflush(stdout);
#endif

  if( xi_type == 1 ){  // single overdispersion param
    xi_z0 = new double[ 1 ];
    xi_z0[ 0 ] = p_xi;
  	log_xi = new NormalDistribution * [ 1 ];
  	log_xi[ 0 ] = new NormalDistribution( 1 );
  	xi_first_draw = new DistributionParameter *[ 1 ];
  	xi_first_draw[ 0 ] = new DistributionParameter( p_xi );
  	number_of_variables = 1;
  	
  } else if( xi_type == 0 ){   // 1 overdispersion param / group
    xi_z0 = new double[ number_of_groups ];
    log_xi = new NormalDistribution * [ number_of_groups ];
    xi_first_draw = new DistributionParameter *[ number_of_groups ];

    for ( i = 0; i < number_of_groups; i++ )
    {
    	xi_z0[ i ] = p_xi;
      xi_first_draw[i] = new DistributionParameter( p_xi );
      log_xi[i] = new NormalDistribution( 1 );
    }
    number_of_variables = number_of_groups;
  }

}//end


void BayesianPoissonLogLinkModel::logXiPrior( double * p_xi )
{
	int i;
	
	if ( xi_type == 1 )  // single overdispersion param
  {
  	log_xi = new NormalDistribution * [ 1 ];
  	log_xi[0] = new NormalDistribution( 1 );
  	xi_z0 = new double[ 1 ];
  	xi_z0[0] = p_xi[0];
    xi_first_draw = new DistributionParameter *[ 1 ];
  	xi_first_draw[0] = new DistributionParameter( xi_z0[0] );
  	number_of_variables = 1;

  } else if( xi_type == 0 ){   // 1 overdispersion param / group
  	log_xi = new NormalDistribution * [ number_of_groups ];
  	xi_first_draw = new DistributionParameter *[ number_of_groups ];
  	xi_z0 = new double[ number_of_groups ];
  	
  	for ( i = 0; i < number_of_groups; i++ )
    {
      xi_first_draw[i] = new DistributionParameter( p_xi[i] );
      log_xi[i] = new NormalDistribution( 1 );
      xi_z0[i] = p_xi[i];
    }
    number_of_variables = number_of_groups;
  }
}

void BayesianPoissonLogLinkModel::samplerDefaultInitialPoint()
{
  int i;

  if ( xi_type == 1 ){  // single overdispersion param.
  	DistributionParameter log_xi_mean( log( xi_first_draw[0]->getScalar() ) );
    log_xi[0]->setLastDraw( log_xi_mean );

  } else if( xi_type == 0 ) {   // one overdispersion param / group
    for ( i = 0; i < number_of_groups; i++ )
    {
      DistributionParameter log_xi_mean( log( xi_first_draw[i]->getScalar() ) );
      log_xi[i]->setLastDraw( log_xi_mean );
    }
  }  // otherwise there is no overdispersion parameter
  
}//end


void BayesianPoissonLogLinkModel::samplerXiInitialPoint( double init_xi )
{
  int i;
#ifdef DEBUG1
  printf( "init_xi = %f\n", init_xi );  fflush(stdout);
#endif
  if ( xi_type == 1 ){    // single overdispersion param.
  	DistributionParameter init_log_xi( log( init_xi ) );
    log_xi[ 0 ]->setLastDraw( init_log_xi );
    
  } else if( xi_type == 0 ){     // one overdispersion param / group
    for ( i = 0; i < number_of_groups; i++ )
    {
      DistributionParameter init_log_xi( log( init_xi ) );
      log_xi[ i ]->setLastDraw( init_log_xi );
    }
  }
}//end


void BayesianPoissonLogLinkModel::samplerXiInitialPoint( CVector & init_xi )
{
  int i;
  
  if ( xi_type == 1 ){    // single overdispersion param.
  	DistributionParameter init_log_xi( log( init_xi.Val( 0 ) ) );
    log_xi[ 0 ]->setLastDraw( init_log_xi );
  } else if( xi_type == 0 ){     // one overdispersion param / group
    for ( i = 0; i < number_of_groups; i++ )
    {
#ifdef DEBUG1
      printf( "init_xi[%d] = %f\n", i, init_xi.Val(i) );  fflush(stdout);
#endif
      DistributionParameter init_log_xi( log( init_xi.Val( i ) ) );
      log_xi[ i ]->setLastDraw( init_log_xi );
    }
  }
}//end



void BayesianPoissonLogLinkModel::setCoefficientDimension( int p_random, int p_fixed )
{
  int j;

  beta_dim = p_random;
  gamma_dim = p_fixed;

  if ( beta_dim > 0 )
  {
    computed_starting_hessian_beta = new CVector( number_of_groups );
    computed_starting_hessian_beta->setToZero();
    hessian_beta = new CMatrix * [ number_of_groups ];
    for ( j = 0; j < number_of_groups; j++ )
    {  
      hessian_beta[j] = new CMatrix( beta_dim, beta_dim );
    }
  }

  if ( gamma_dim > 0 )
  {
    computed_starting_hessian_gamma = 0;
    hessian_gamma = new CMatrix( gamma_dim, gamma_dim );
  }
}//end



void BayesianPoissonLogLinkModel::considerCounts()
{
  int j, k;

  initializeCountStorage();
  
  for ( j = 0; j < number_of_groups; j++ )
  {
    max_counts->Val( j ) = counts[ j ]->max();
    mean_counts->Val( j ) = counts[ j ]->Mean();

    n_equal[j] = new CVector( ((int) max_counts->Val( j )) );
    n_larger_equal[j] = new CVector( ((int) max_counts->Val( j )) );
    
    for ( k = ((int) max_counts->Val( j )); k >= 1; k-- )
    {
      n_equal[ j ]->Val( k - 1 ) = counts[ j ]->howManyEqual( ((double) k) );
      if ( k < ((int) max_counts->Val( j )) )
      {
        n_larger_equal[ j ]->Val( k - 1 ) = n_equal[ j ]->Val( k - 1 ) + n_larger_equal[ j ]->Val( k );
      }
      else
      {
        n_larger_equal[ j ]->Val( k - 1 ) = n_equal[ j ]->Val( k - 1 );
      }
    }//end for k
  } // end for j

  if( xi_type == 1 ){  // common xi parameter for all groups
  	// Pool the counts into a single max_counts value, single n_larger_equal vector,
  	// and single mean_counts vector.
  	// First find the overall max outcome (count) measurement
  	double mxCount = max_counts->max();
  	// Create a vector that contains the total numbers of each outcome (count) value
  	CVector *nLargerEqual = new CVector( (int) mxCount );
  	// Create a variable to calculate the overall outcome (count) mean
  	double meanCount = 0.0;
  	nLargerEqual->setToZero();
  	for( j = 0; j < number_of_groups; j++ ){
  	  for( k = 0; k < n_larger_equal[ j ]->Len(); k++ ){
  		  nLargerEqual->Val( k ) = nLargerEqual->Val( k ) + n_larger_equal[ j ]->Val( k );
  	  }
  	  meanCount += mean_counts->Val( j ) * exposures[ j ]->Len();
  	}
  	max_counts->Val( 0 ) = mxCount;
        if (n_larger_equal[0] != NULL) {
          delete n_larger_equal[0] ;
        }
  	n_larger_equal[ 0 ] = nLargerEqual;
  	meanCount = meanCount / (double) number_of_observations;
  	mean_counts->Val( 0 ) = meanCount;
  } 

  //compute moj for correction of logXi Hessian
  for ( j = 0; j < number_of_groups; j++ )
  {
    m0->Val(j) = counts[j]->Mean() / exposures[j]->Mean();
  }

}//end



void BayesianPoissonLogLinkModel::metropolisHastingsUpdateModel( DistributionParameter * par_vals )
{

  if ( par_vals->length() == 1 )
  {
    metropolisHastingsUpdateModel( ((int) par_vals->getScalar()) );
  }
  else
  {
#ifdef FIX1
    CVector tmpvec = par_vals->getVector().subVector( 1, par_vals->length() - 1 );
    metropolisHastingsUpdateModel( ((int) par_vals->getVector().Val( 0 )), tmpvec );
#else
    metropolisHastingsUpdateModel( ((int) par_vals->getVector().Val( 0 )), 
      par_vals->getVector().subVector( 1, par_vals->length() - 1 ) );
#endif
//    metropolisHastingsUpdateModel( ((int) par_vals->getVector().Val( 0 )), par_vals->getVector().subVector( 1, par_vals->length() - 1 ) );
  }
}//end



void BayesianPoissonLogLinkModel::metropolisHastingsUpdateModel( int j, CVector & mu )
{
  int i;
  double xi, lambda;

  if( xi_type == 1 )  // common overdispersion parameter case
    xi = exp( log_xi[ 0 ]->lastDraw().getScalar() );
  else if( xi_type == 0 )   // group-specific overdispersion parameter case
    xi = exp( log_xi[ j ]->lastDraw().getScalar() );
  else          // no overdispersion case
  	xi = 0.0;

  #ifdef DEBUGGLM
    printf("mhUpdateModel: enter j = %d, xi = %f, mu = \n", j, xi  );  
    mu.Print();  fflush(stdout);
  #endif

  for ( i = 0; i < mu_linear[ j ]->Len(); i++ )
  {
    mu_linear[ j ]->Val( i ) = exp( mu.Val( i ) );
    b_denom[ j ]->Val( i ) = exposures[ j ]->Val( i ) * mu_linear[ j ]->Val( i ) + xi;
    b[ j ]->Val( i ) = xi / b_denom[ j ]->Val( i );

    //update diagonal covariance for hessians of beta and gamma
    lambda = ( 1 - b[ j ]->Val( i ) ) * ( counts[ j ]->Val( i ) / exposures[ j ]->Val( i ) ) + b[ j ]->Val( i ) * mu_linear[ j ]->Val( i );
    d[ j ]->Val( i )  = exposures[ j ]->Val( i ) * b[ j ]->Val( i ) * lambda;
  }

  //printf("mhUpdateModel: mu_linear = and d = \n");
  //mu_linear[j]->Print();
  //d[j]->Print();

  //update (gradient and) hessian of logXi
  if ( ( computed_starting_hessian_logXi->Val( j ) == 0 ) && 
  	( xi_type != 2 ) )   // check for overdispersion parameters
    computeHessianForLogXi( j, xi );

}//end


void BayesianPoissonLogLinkModel::metropolisHastingsUpdateModel( int j )
{
  int i;
  double xi, lambda;
  
  if( xi_type == 1 )  // common overdispersion parameter case
    xi = exp( log_xi[ 0 ]->lastDraw().getScalar() );
  else if( xi_type == 0 )   // group-specific overdispersion parameter case
    xi = exp( log_xi[ j ]->lastDraw().getScalar() );
  else          // no overdispersion case
  	xi = 0.0;
  	
  #ifdef DEBUGGLM
    printf("mhUpdateModel: enter j = %d, xi = %f, mu_linear = \n", j, xi  );  
    mu_linear[j]->Print();  fflush(stdout);
  #endif
  for ( i = 0; i < mu_linear[ j ]->Len(); i++ )
  {
    b_denom[ j ]->Val( i ) = exposures[ j ]->Val( i ) * mu_linear[ j ]->Val( i ) + xi;
    b[ j ]->Val( i ) = xi / b_denom[ j ]->Val( i );

    //update diagonal covariance for hessians of beta and gamma
    lambda = ( 1 - b[ j ]->Val( i ) ) * ( counts[ j ]->Val( i ) / exposures[ j ]->Val( i ) ) + b[ j ]->Val( i ) * mu_linear[ j ]->Val( i );
    d[ j ]->Val( i )  = exposures[ j ]->Val( i ) * b[ j ]->Val( i ) * lambda;
  }

  //printf("mhUpdateModel: d = \n");
  //d[j]->Print();

  //update (gradient and) hessian of logXi
  if ( ( computed_starting_hessian_logXi->Val( j ) == 0 ) && 
  	( xi_type != 2 ) )   // check for overdispersion parameters
    computeHessianForLogXi( j, xi );

}//end


CVector BayesianPoissonLogLinkModel::computeWeights( int var_index, DistributionParameter ** predictors, CVector & b_star, CVector & b_0 )
{
  int i;
  double xi, zxb, proposal_mu_linear, proposal_b_denom, proposal_b, proposal_lambda;

  CVector proposal_d( exposures[ var_index ]->Len() );

  CVector bdiff( b_star - b_0 );
  CVector xb( predictors[ var_index + 2 ]->getMatrix() * bdiff );
  zxb = (*counts[ var_index ]) * xb;
   
  if( xi_type == 1 )  // common overdispersion param. case
    xi = exp( log_xi[ 0 ]->lastDraw().getScalar() );
  else if( xi_type == 0 )   // group-specific overdispersion param. case
    xi = exp( log_xi[ var_index ]->lastDraw().getScalar() );
  else     // no overdispersion param. case.  Crease a pretend overdispersion parameter for this calculation.
  	xi = 1;     // Need a less arbitrary choice here 
    
  for ( i = 0; i < exposures[ var_index ]->Len(); i++ )
  {
    proposal_mu_linear = mu_linear[ var_index ]->Val( i ) * exp( xb.Val( i ) );
    proposal_b_denom = xi + ( b_denom[ var_index ]->Val( i ) - xi ) * exp( xb.Val( i ) );
    proposal_b = xi / proposal_b_denom;

    //update diagonal covariance for hessians of beta and gamma
    proposal_lambda = ( 1 - proposal_b ) * ( counts[ var_index ]->Val( i ) / exposures[ var_index ]->Val( i ) ) + proposal_b * proposal_mu_linear;

    proposal_d.Val( i ) = exposures[ var_index ]->Val( i ) * proposal_b * proposal_lambda;
  }

  return (proposal_d);

}//end




double BayesianPoissonLogLinkModel::logXiGradient( int j, double xi )
{
  int i, k;
  double sum1, sum2, sum3, gradient;

  //update gradient of logXi
  sum1 = 0.0;
  for ( k = 0; k < ((int) max_counts->Val(j)); k++ )
  {
    sum1 += n_larger_equal[j]->Val(k) * ( xi / ( xi + (k + 1) - 1 ) );
  }

  sum2 = 0.0;
  sum3 = 0.0;
  for ( i = 0; i < exposures[j]->Len(); i++ )
  {
    sum2 += log( b[j]->Val(i) );
    sum3 += ( counts[j]->Val(i) - exposures[j]->Val(i) * mu_linear[ j ]->Val(i) ) * b[j]->Val(i);
  }
  sum2 *= xi;

  //printf("in logXi gradient: xi = %f, sum1 = %f, sum2 = %f, sum3 = %f\n", xi, sum1, sum2, sum3 ); fflush(stdout);
  //printf( "in logXi gradient: b[%d] = \n", j);
  //b[j]->Print();

  gradient = sum1 + sum2 - sum3;

  return (gradient);

}//end



void BayesianPoissonLogLinkModel::computeHessianForLogXi(  int j, double xi )
{
  if( xi_type == 2 ){
    printf( "Attempt to compute Hessian for nonexistent parameter xi[%d]\n", j );  fflush(stdout);
  	MESSAGE "" ERROR;
  } 

  int i, k;
  double sum1, sum2, sum3;

  //printf( "in compute Hessian for logXi: j = %d, xi = %f\n", j, xi ); fflush(stdout);

  sum1 = 0.0;
  for ( k = 0; k < ((int) max_counts->Val(j)); k++ )
  {
    sum1 += n_larger_equal[j]->Val(k) * ( xi / ( xi + (k + 1) - 1 ) ) * ( xi / ( xi + (k + 1) - 1 ) );
  }
  
  sum2 = 0.0;
  sum3 = 0.0;
  for ( i = 0; i< exposures[j]->Len(); i++ )
  {
    sum2 += counts[j]->Val(i) *  b[j]->Val(i);
    sum3 += ( counts[j]->Val(i) - exposures[j]->Val(i) * mu_linear[ j ]->Val(i) ) * b[j]->Val(i) * ( 1 - b[j]->Val(i) );
  }

  //printf("in compute Hessian logXi: sum1 = %f, sum2 = %f, sum3 = %f \n", sum1, sum2, sum3 ); fflush(stdout);
  //printf( "in compute Hessian logXi: b[%d] = \n", j);
  //b[j]->Print();
  double z0;
  if( xi_type == 1 )   // common overdispersion parameter case
  	z0 = xi_z0[ 0 ];
  else                 // group-specific overdispersion parameter case
  	z0 = xi_z0[ j ];
  log_xi_hessian->Val(j) = -sum1 + sum2 - sum3 - 2 * ( ( xi * z0 ) / ( ( xi + z0 ) * ( xi + z0 ) ) );
  log_xi_hessian->Val(j) += logXiGradient( j, xi );

  computed_starting_hessian_logXi->Val( j ) = 1;

}//end



CMatrix BayesianPoissonLogLinkModel::metropolisHastingsProposalHessian( int n_pars, 
	DistributionParameter ** par_vals, DistributionParameter ** par_vector )
{
  int j;

  if ( n_pars > 0 )
  {
    if ( par_vector != NULL )
    {
      CVector b_star( par_vector[0]->getVector() );
      CVector b_0( par_vector[1]->getVector() );

      if ( par_vals != NULL )
      {
        if ( par_vals[0]->getScalar() == 0 )
        {
          //beta
          j = ((int) par_vals[1]->getScalar());

          CVector proposal_d( computeWeights( j, par_vals, b_star, b_0 ) );
          CMatrix proposal_hessian( (-1.0) * par_vals[ j + 2 ]->getMatrix().xTransposedX( proposal_d ) );
       
          return (proposal_hessian);
        }
        else
        {
          //gamma
          if ( n_pars == number_of_groups + 1 )
	  {
            CVector proposal_d( computeWeights( 0, par_vals, b_star, b_0 ) );
            CMatrix proposal_hessian( par_vals[ 2 ]->getMatrix().xTransposedX( proposal_d ) );
            if ( number_of_groups > 1 )
            {
              for ( j = 1; j < number_of_groups; j++ )
              {   
                proposal_d = computeWeights( j, par_vals, b_star, b_0 );
#ifdef FIX1
                CMatrix tmpmat = par_vals[ j + 2 ]->getMatrix().xTransposedX( proposal_d );
                proposal_hessian.add( tmpmat );
#else
                proposal_hessian.add( par_vals[ j + 2 ]->getMatrix().xTransposedX( proposal_d ) );
#endif
//                proposal_hessian.add( par_vals[ j + 2 ]->getMatrix().xTransposedX( proposal_d ) );
              }
            } 
            proposal_hessian.multiplyByScalar( -1.0 );   
  
            return ( proposal_hessian );
          }
          else
	  {
            printf( " BayesianPoissonLogLinkModel::metropolisHastingsProposalHessian: Wrong number of parameters for Gamma [%d].\n", n_pars );
          }
        }//end gamma
      }//end if not null
      else
      {
        printf( " BayesianPoissonLogLinkModel::metropolisHastingsProposalHessian: Null array.\n" );
      }
    }//end if par_vector
    else
    {
      printf( " BayesianPoissonLogLinkModel::metropolisHastingsProposalHessian: Null array.\n" );
    }
  }//end if > 0
  else
  {
    printf( " BayesianPoissonLogLinkModel::metropolisHastingsProposalHessian: No parameters passed.\n" );
  }

  CMatrix null_hessian( 1, 1 );
  null_hessian.setToZero();
  return ( null_hessian );
}//end



CVector BayesianPoissonLogLinkModel::metropolisHastingsMean( int n_pars, 
	DistributionParameter ** par_vals, DistributionParameter ** par_vector )
{
  //not used

  CVector null_vector( 1 );
  null_vector.setToZero();
  return ( null_vector );

}//end




CMatrix BayesianPoissonLogLinkModel::metropolisHastingsHessian( int n_pars, DistributionParameter ** par_vals )
{
  int j;

  if ( n_pars > 0 )
  {
    if ( par_vals != NULL )
    {
      if ( par_vals[0]->getScalar() == 0 )
      {
        //beta
        j = ((int) par_vals[1]->getScalar());

        //printf("in Hessian: j = %d, predictors = \n", j );
        //par_vals[ j + 2]->Print();
        //printf("in Hessian d = \n");
        //d[j]->Print();

        if ( !update_hessians && computed_starting_hessian_beta->Val( j ) == 0 )
	{
          (*(hessian_beta[ j ])) = ( (-1.0) * par_vals[ j + 2 ]->getMatrix().xTransposedX( (*d[ j ]) ) );
          computed_starting_hessian_beta->Val( j ) = 1;
        }
        else if ( update_hessians )
	{
          (*(hessian_beta[ j ])) = ( (-1.0) * par_vals[ j + 2 ]->getMatrix().xTransposedX( (*d[ j ]) ) );
        }
       
        return ( (*(hessian_beta[j])) );
      }
      else
      {
        //gamma
        if ( n_pars == number_of_groups + 1 )
	{
          if ( !update_hessians && computed_starting_hessian_gamma == 0 )
	  {
            (*hessian_gamma) = ( par_vals[ 2 ]->getMatrix().xTransposedX( (*d[ 0 ]) ) );
            if ( number_of_groups > 1 )
            {
              for ( j = 1; j < number_of_groups; j++ )
              {  
#ifdef FIX1
                CMatrix tmpmat = par_vals[ j + 2 ]->getMatrix().xTransposedX( (*d[ j ]) );
                hessian_gamma->add( tmpmat );
#else
                hessian_gamma->add( par_vals[ j + 2 ]->getMatrix().xTransposedX( (*d[ j ]) ) );
#endif
              }
            } 
            hessian_gamma->multiplyByScalar( -1.0 );   
            computed_starting_hessian_gamma = 1;
          }
          else if ( update_hessians )
	  {
            (*hessian_gamma) = ( par_vals[ 2 ]->getMatrix().xTransposedX( (*d[ 0 ]) ) );
            if ( number_of_groups > 1 )
            {
              for ( j = 1; j < number_of_groups; j++ )
              {  
#ifdef FIX1
                CMatrix tmpmat = par_vals[ j + 2 ]->getMatrix().xTransposedX( (*d[ j ]) );
                hessian_gamma->add( tmpmat );
#else
                hessian_gamma->add( par_vals[ j + 2 ]->getMatrix().xTransposedX( (*d[ j ]) ) );
#endif
//                hessian_gamma->add( par_vals[ j + 2 ]->getMatrix().xTransposedX( (*d[ j ]) ) );
              }
            } 
            hessian_gamma->multiplyByScalar( -1.0 );   
          }
  
          return ( (*hessian_gamma) );
        }
        else
	{
          printf( " BayesianPoissonLogLinkModel::metropolisHastingsHessian: Wrong number of parameters for Gamma [%d].\n", n_pars );
        }
      }//end gamma
    }//end if not null
    else
    {
      printf( " BayesianPoissonLogLinkModel::metropolisHastingsHessian: Null array.\n" );
    }
  }//end if > 0
  else
  {
    printf( " BayesianPoissonLogLinkModel::metropolisHastingsHessian: No parameters passed.\n" );
  }

  CMatrix null_hessian( 1, 1 );
  null_hessian.setToZero();
  return ( null_hessian );
}//end


// This function never gets used
double BayesianPoissonLogLinkModel::hessianForlogXi( int i )
{
  if( xi_type == 2 ){
    printf( "Attempt to compute Hessian for nonexistent parameter xi[%d]\n", i );  fflush(stdout);
  	MESSAGE "" ERROR;
  } 
  if( xi_type == 1 )  // common overdispersion param. case
    return ( log_xi_hessian->Val(0) );
  else    //  group-specific overdispersion param. case
    return ( log_xi_hessian->Val(i) );
}//end


// Log xi[j] has a proposal distribution the variance of which is set in the first iteration of 
// the sampler and never changed.  The mean of the proposal distribution is the last draw of log xi
// (apparently).  In the case where there is a common xi parameter, the log xi proposal variance 
// is calculated just using data from the first group; this should be fixed.  The Markov chain is
// still correct but may not be as efficient as possible.
// If (xi_type == 2), so that there is no overdispersion parameter, this function never gets called.
void BayesianPoissonLogLinkModel::metropolisHastingsUpdateLogXi( int j )
{
  int i;
  double xi, correction_term, denom, var_log_xi;
  
  if( ( xi_type == 2 ) ||
  	( (j > 0) && (xi_type == 1) ) ){
  	MESSAGE "Attempt to update nonexistent parameter xi[j]" ERROR;
  } 

  xi = exp( log_xi[ j ]->lastDraw().getScalar() );

  #ifdef DEBUGGLM
    printf( "Update Log Xi: xi = %f, m0 = %f, exposure = \n", xi, m0->Val(j) ); fflush(stdout);
    exposures[j]->Print();
  #endif

  var_log_xi = 1;
  if ( computed_starting_logXi_variance->Val( j ) == 0 )  // If you have not done this computation before for this group
  {
      correction_term = 0.0;
      
      for ( i = 0; i < exposures[j]->Len(); i++ )
      {
        denom = ( xi + exposures[j]->Val(i) * m0->Val(j) ) * ( xi + exposures[j]->Val(i) * m0->Val(j) );
        correction_term += ( exposures[j]->Val(i) / denom );
      }
      correction_term *= m0->Val(j) * xi;
      correction_term *= ( ((double) beta_dim ) / ((double) 2 * exposures[j]->Len() ) );

      denom = log_xi_hessian->Val(j) + correction_term;
      if ( denom < 0.0 && ( denom <= - 0.1 ) )   // the upper bound for denom has been changed from -DELTA
      {
        var_log_xi = - 1.0 / denom;
      }
      else 
      {
        double w;

        //too small or negative variance: use gradient estimate of variance
        correction_term = 0.0;
        for ( i = 0; i < exposures[j]->Len(); i++ )
        {
          denom = ( xi + exposures[j]->Val(i) * m0->Val(j) );
          correction_term += ( exposures[j]->Val(i) / denom );
        }
        correction_term *= m0->Val(j);
        correction_term *= ( ((double) beta_dim ) / ((double) 2 * exposures[j]->Len() ) );

        double z0;
        if( xi_type == 1 )   // common overdispersion param. case
        {
        	z0 = xi_z0[ 0 ];
        }
        else if( xi_type == 0 )
        {
        	z0 = xi_z0[ j ];   // group-specific overdispersion param. case
        }
        else
        {
                error("internal error: unexpected xi_type=%d, z0 not set", xi_type) ;
        }
        w =  logXiGradient( j, xi )  - ( ( xi - z0 ) / ( xi + z0 ) ) - correction_term;
        w = w * w;

        if ( w  >=  DELTA )
        {
          var_log_xi = 1.0 / w;
        }
        else
        {
          var_log_xi = 1.0 / DELTA;
        }
      }
    
      logXi_variance->Val( j ) = var_log_xi;    // Store for all future calls to this function
      computed_starting_logXi_variance->Val( j ) = 1;   

      //printf( "correction_term = %f\n", correction_term ); fflush(stdout);
  }
  else
  {
  	if( xi_type == 1 )  // common overdispersion param. case
  		var_log_xi = logXi_variance->Val( 0 );
    else if( xi_type == 0 )   // group-specific overdispersion param. case
      var_log_xi = logXi_variance->Val( j );
  }

  CVector xi_mean( 1 );
  CMatrix xi_var( 1, 1 );
  
  if( xi_type == 1 )   //  common overdispersion param. case
    xi_mean.Val( 0 ) = log_xi[ 0 ]->lastDraw().getScalar();
  else if( xi_type == 0 )   // group-specific overdispersion param. case
    xi_mean.Val( 0 ) = log_xi[ j ]->lastDraw().getScalar();
  
  xi_var.Val( 0, 0 ) = var_log_xi;

  
  //printf( "log_xi_hessian = %f\n", log_xi_hessian->Val(j) ); fflush(stdout);
  //printf( "xi_var = %f, mean = %f\n", xi_var.Val( 0, 0 ), xi_mean.Val(0) ); fflush(stdout);
  if( xi_type == 1 )   //  common overdispersion param. case
    log_xi[0]->updateInitialMean( &xi_mean );
  else if( xi_type == 0 )   // group-specific overdispersion param. case
  	log_xi[j]->updateInitialMean( &xi_mean );
  	
  if ( xi_var.Val(0,0) > 0 )
  {
    if( xi_type == 1 )   //  common overdispersion param. case
      log_xi[0]->updateInitialCovariance( &xi_var );
    else if( xi_type == 0 )   // group-specific overdispersion param. case
    	log_xi[j]->updateInitialCovariance( &xi_var );
  }
  else
  {
    //keep the same covariance
    printf( "BayesianPoissonLogLinkModel::metropolisHastingsUpdateLogXi: Keeping previous covariance. Current one is %f\n", 
      xi_var.Val(0,0) );
    fflush( stdout );
  }
  
  #ifdef DEBUGGLM
    printf( "log xi proposal mean = %f, var = %f\n", xi_mean.Val(0), xi_var.Val(0,0) );  fflush(stdout);
  #endif
  
}//end


// This function is never used
void BayesianPoissonLogLinkModel::simulationsToArray( double * simul_output, 
	int simulations_to_keep, int start_index )
{
  int j, s, start_dim;

  start_dim = start_index;
  
  if( xi_type == 1 ){  // common overdispersion param. case
  	for ( s = 0; s < simulations_to_keep; s++ ){
  	  simul_output[ start_dim + s ] = simulated_xi->Val( s, 0 );
  	}  
  } else if( xi_type == 0){   // group-specific overdispersion param. case
    for ( j = 0; j < number_of_groups; j++ )
    {
      for ( s = 0; s < simulations_to_keep; s++ )
      {
        simul_output[ start_dim + j * simulations_to_keep + s ] = simulated_xi->Val( s, j );
        //printf( "simul xi[%d] = %f\n", s, simulated_xi->Val( s, j ) ); fflush(stdout);
      }
    }
  }

}//end


// if( xi_type == 2 ), so that there is no overdispersion, this function should not
// alter simul_output
void BayesianPoissonLogLinkModel::simulationsToArray( double * simul_output, 
	int start_index, CMatrix * simulated_mu )
{
  int j, s, start_dim, simulations_to_keep;

  simulations_to_keep = simulated_mu->Row();
  start_dim = start_index;
  if( xi_type == 1 ){   // common overdispersion param. case
    for ( s = 0; s < simulations_to_keep; s++ )
    {
      simul_output[ start_dim + s ] = simulated_xi->Val( s, 0 );
    }
    start_dim += simulations_to_keep;
    
  } else if( xi_type == 0 ){    // group-specific overdispersion param. case
    for ( j = 0; j < number_of_groups; j++ )
    {
      for ( s = 0; s < simulations_to_keep; s++ )
      {
        simul_output[ start_dim + j * simulations_to_keep + s ] = simulated_xi->Val( s, j );
        //printf( "simul xi[%d] = %f\n", s, simulated_xi->Val( s, j ) ); fflush(stdout);
      }
    }
    start_dim += number_of_groups * simulations_to_keep;
    
  } else {  // no overdispersion param. case
  	// nothing to do
  }
  
  if( xi_type != 2 ){   // if there is overdispersion
    generateLambdas( simulated_mu );
    //keep lambdas
    for ( j = 0; j < number_of_observations; j++ )
    {
      for ( s = 0; s < simulations_to_keep; s++ )
      {
        simul_output[ start_dim + j * simulations_to_keep + s ] = simulated_lambda->Val( s, j );
      }
    }
  }
  
}//end


// If (xi_type == 2) then this function should never be called since there is no 
// overdispersion.
void BayesianPoissonLogLinkModel::generateLambdas( CMatrix * simulated_mu )
{
	if( xi_type == 2 ){
		printf("BayesianPoissonLogLinkModel::drawVariableFromProposal: There is no parameter to draw\n");
		return;
	}
  int i, j, k, s, simulations_to_keep;
  double xi, wmu, w1, w2, wgamma;

  simulations_to_keep = simulated_mu->Row();
  k = 0;
  for ( j = 0; j < number_of_groups; j++ )
  {
    //get alpha and beta
    for ( i = 0; i < exposures[j]->Len(); i++ )
    {
      for ( s = 0; s < simulations_to_keep; s++ )
      {
      	if( xi_type == 1 )  // common overdispersion param. case
          xi = simulated_xi->Val( s, 0 );
        else               // group-specific overdispersion case
        	xi = simulated_xi->Val( s, j );
        wmu = exp( simulated_mu->Val( s, k ) );
        w1 = counts[ j ]->Val( i ) + xi;
        w2 = exposures[ j ]->Val( i ) + ( xi / wmu );

        //wgamma/w2 ~ Gamma( w1, w2 )
        wgamma = rgamma( w1, 1.0 );
        simulated_lambda->Val( s, k ) = wgamma / w2;

        //printf( "generate lambdas: k = %d, s = %d, xi = %f, wx = %f, wmu = %f, w1 = %f, w2 = %f, lambda = %f\n", k, s, xi, wx, wmu, w1, w2, simulated_lambda->Val( s, k ) ); fflush(stdout);

      }//end for s
      k++;
    }//end for i
  }//end for j

  //printf("Lambdas are\n" );
  //simulated_lambda->Print();

}//end


void BayesianPoissonLogLinkModel::keepSimulation( int simul_number )
{
  int j;

  if ( simul_number < number_of_simulations )
  {
  	if( xi_type == 1 ){  // common overdispersion param. case
  		simulated_xi->Val( simul_number, 0 ) = exp( log_xi[0]->lastDraw().getScalar() );
  	} else if ( xi_type == 0 ) {   // separate overdispersion param. / group
      for ( j = 0; j < number_of_groups; j++ )
      {
        simulated_xi->Val( simul_number, j ) = exp( log_xi[j]->lastDraw().getScalar() );
      }
    } else {    // no overdispersion param. case
    	//  do nothing
    }
  }    

}//end



void BayesianPoissonLogLinkModel::createOutput( int sim_to_keep )
{
  number_of_simulations = sim_to_keep;
  
  if ( xi_type == 1 )   // common overdispersion param. case
  {
    simulated_xi = new CMatrix( sim_to_keep, 1 );
  }
  else if( xi_type == 0 )  // group-specific overdispersion param. case
  {
    simulated_xi = new CMatrix( sim_to_keep, number_of_groups );
  } 

  if( xi_type != 2 )    // if there is overdispersion
    simulated_lambda = new CMatrix( number_of_simulations, number_of_observations );
    
}//end


void BayesianPoissonLogLinkModel::initializeTemporaryStructures()
{
  //nothing to do here
  //except computing starting (and perhaps fixed) hessian for xi
  
  

}//end


void BayesianPoissonLogLinkModel::dataAugmentationInitialDraws()
{
  //nothing to do here because the lambda parameters are integrated 
  // out in the sampling procedure; no data augmentation needs to be used.
}//end


// If (xi_type == 2), so that there is no overdispersion parameter, this function never gets called
void BayesianPoissonLogLinkModel::drawVariableFromProposal( int var_index )
{
	if( xi_type == 2 ){
		printf("BayesianPoissonLogLinkModel::drawVariableFromProposal: There is no parameter to draw\n");
		return;
	}
	if( (var_index > 0) && (xi_type == 1) ){
		MESSAGE "Attempt to draw nonexistant parameter xi[j]" ERROR;
	}
	
  if ( var_index < number_of_variables )
  {
    log_xi[ var_index ]->draw();
  }
  else
  {
    printf( " BayesianPoissonLogLinkModel::drawVariableFromProposal: Wrong argument index [%d].\n", var_index );
  }
}//end


// If (xi_type == 2), so that there is no overdispersion parameter, this function never gets called
void BayesianPoissonLogLinkModel::updateVariableForProposal( int var_index )
{
	if( xi_type == 2 ){
		printf("BayesianPoissonLogLinkModel::updateVariableForProposal: There is no parameter to update\n");
		return;
	}
  if ( var_index < number_of_variables )
  {
    metropolisHastingsUpdateLogXi( var_index );
  }
  else
  {
    printf( " BayesianPoissonLogLinkModel::updateVariableForProposal: Wrong argument index [%d].\n", var_index );
  }
}//end



// If (xi_type == 2), so that there is no overdispersion parameter, this function never 
// gets called.
void BayesianPoissonLogLinkModel::setCurrentVariableFromProposal( int var_index )
{
  if( xi_type == 2 ){
		printf("BayesianPoissonLogLinkModel::updateVariableForProposal: There is no parameter to update\n");
		return;
	}
  if ( var_index < number_of_variables )
  {
  	// We have changed the value of xi so we need to update b_denom, b, etc.
  	// If we have updated the common xi parameter, then we need to update b_denom,
  	// b, etc. for all groups, and otherwise we just need to update them for the
  	// group whose xi parameter we have changed.
  	if( xi_type == 1 ){  // common overdispersion param. case
  		for( int i=0; i<number_of_groups; i++ )
  		  metropolisHastingsUpdateModel( i );
  	} else if( xi_type == 0 ) {
      metropolisHastingsUpdateModel( var_index );
    }
    #ifdef DEBUG1
      printf(" Accepted proposed change in xi[%d].  New value is %f.\n", var_index, 
        exp( log_xi[ var_index ]->lastDraw().getScalar() ) );
    #endif
  }
  else
  {
    printf( " BayesianPoissonLogLinkModel::setCurrentVariableFromProposal: Wrong argument index [%d].\n", var_index );
  }
}//end


// If (xi_type == 2), so that there is no overdispersion parameter, this function never 
// gets called.
void BayesianPoissonLogLinkModel::keepCurrentVariable( int var_index )
{
	if( xi_type == 2 ){
		printf("BayesianPoissonLogLinkModel::keepCurrentVariable: There is no parameter to keep\n");
		return;
	}
  if ( var_index < number_of_variables )
  {
  	#ifdef DEBUG1
  	  printf("rejected proposed change in xi.  restoring xi[var_index] = %f\n", 
  	    exp( log_xi[ var_index ]->mean().getScalar() ) );  fflush(stdout);
  	#endif
  	
#ifdef FIX1
    DistributionParameter tmpdist = log_xi[ var_index ]->mean();
    log_xi[ var_index ]->setLastDraw( tmpdist );
#else
    log_xi[ var_index ]->setLastDraw( log_xi[ var_index ]->mean() );
#endif
//    log_xi[ var_index ]->setLastDraw( log_xi[ var_index ]->mean() );
  }
  else
  {
    printf( " BayesianPoissonLogLinkModel::keepCurrentVariable: Wrong argument index [%d].\n", var_index );
  }
}//end


// Calculates the log ratio target density for a proposed change in xi.  The ratio is
// P(new) / P(old).
// The calculation is done as described in Christiansen and Morris (1997), p. 621.
// If (xi_type == 2), so that there is no overdispersion parameter, this function never 
// gets called.
double BayesianPoissonLogLinkModel::logRatioTargetDensity( int var_index )
{
	if( xi_type == 2 ){
		MESSAGE "BayesianPoissonLogLinkModel::logRatioTargetDensity: There is no parameter to calculate\n" ERROR;
	}
	if( (xi_type == 1) && (var_index > 0) )  // common overdispersion parameter case
		MESSAGE "logRatioTargetDensity: var_index out of bounds" ERROR;

  int i, j, k;
  double ratio, sum1, sum2, sum3, b_star, xi, xi_star;

  ratio = 0;
  if ( var_index < number_of_variables )
  {
    xi_star = exp( log_xi[ var_index ]->lastDraw().getScalar() );
    xi = exp( log_xi[ var_index ]->mean().getScalar() );
#ifdef DEBUGGLM
  printf("xi = %f \n", xi );
  printf("xi_star = %f \n", xi_star );
  printf("var_index = %d, xi_z0[var_index] = %f, max_counts[var_index] = %f \n", var_index, xi_z0[var_index], max_counts->Val(var_index) );
  printf("n_larger_equal[0] = \n");  n_larger_equal[0]->Print();
  printf("number_of_observations = %d, mean_counts->Val(var_index) = %f\n", number_of_observations, mean_counts->Val(var_index));
  printf("b_denom[<first group>] = \n");  b_denom[0]->Print();
  printf("b[<first group>] = \n");  b[0]->Print();
  if( number_of_groups > 1 ){
    printf("b_denom[<2nd group>] = \n");  b_denom[1]->Print();
    printf("b[<2nd group>] = \n");  b[1]->Print();
  }
  fflush(stdout);
#endif
    sum1 = 0.0;
    for ( k = 0; k < ((int) max_counts->Val( var_index )); k++ )
    {
      sum1 += n_larger_equal[ var_index ]->Val( k ) * log( ( xi_star + (k + 1) - 1 ) / ( xi + (k + 1) - 1 ) );
    }

    if( xi_type == 1 )   // common overdispersion param. case
    	sum1 -= number_of_observations * mean_counts->Val( var_index ) * log( xi_star / xi );
    else                 // group-specific overdispersion case
      sum1 -= exposures[ var_index ]->Len() * mean_counts->Val( var_index ) * log( xi_star / xi );
    
    sum2 = 0.0;
    sum3 = 0.0;
    if( xi_type == 1 ){  // common overdispersion param. case
    	// have to take into account all observations
    	for ( j = 0; j < number_of_groups; j++ ){
    		for ( i = 0; i < exposures[ j ]->Len(); i++ ){
          b_star = xi_star / ( b_denom[ j ]->Val( i ) - xi + xi_star );
          sum2 += xi_star * log( 1.0 / b_star );
          sum2 -= xi * log( 1.0 / b[ j ]->Val( i ) );
          sum3 += counts[ j ]->Val( i ) * log( b[ j ]->Val( i ) / b_star ); 	
    		}
    	}
    } else {
    	// just take into account the observations for the group with index var_index
      for ( i = 0; i < exposures[ var_index ]->Len(); i++ )
      {
        b_star = xi_star / ( b_denom[ var_index ]->Val( i ) - xi + xi_star );
        sum2 += xi_star * log( 1.0 / b_star );
        sum2 -= xi * log( 1.0 / b[ var_index ]->Val( i ) );
        sum3 += counts[ var_index ]->Val( i ) * log( b[ var_index ]->Val( i ) / b_star ); 
      }
    } 
    #ifdef DEBUGGLM
      printf("sum3 = %f \n", sum3);
    #endif
    
    sum3 += 2 * log( ( xi_star + xi_z0[ var_index ] ) / ( xi + xi_z0[ var_index ] ) );

    ratio = sum1 - sum2 - sum3;
    #ifdef DEBUGGLM
      printf("sum1 = %f, sum2 = %f, sum3 = %f \n", sum1, sum2, sum3);
      printf("ratio = %f \n", ratio);  fflush(stdout);
    #endif
  }
  else
  {
    printf( " BayesianPoissonLogLinkModel::logRatioTargetDensity: Wrong argument index [%d].\n", var_index );
  }

  return ( ratio );

}//end


// Calculates the log likelihood ratio for a proposed change in beta (the random effects)
// or gamma (the fixed effects), for the case where there is overdispersion.
// The calculation is done as described in Christiansen and Morris (1997), p. 621
double BayesianPoissonLogLinkModel::logRatioTargetDensity( DistributionParameter ** predictors, 
	CVector & b_star, CVector & b_0 )
{
  int i, j, var_index;
  double sum, numer, ratio1, ratio, xi, zxb;

  ratio = 0.0;
  if ( ((int) predictors[0]->getScalar()) == 0 )
  {
    //beta
    var_index = ((int) predictors[1]->getScalar());

    CVector bdiff( beta_dim );
    CVector xb( counts[ var_index ]->Len() );

    bdiff = b_star - b_0;
    xb = predictors[ var_index + 2 ]->getMatrix() * bdiff;
    #ifdef DEBUGBETA
      printf("xb = \n");  xb.Print();  
      printf("counts[ var_index ] \n");  counts[ var_index ]->Print();  fflush(stdout);
    #endif
    zxb = (*counts[ var_index ]) * xb;
    
    if( xi_type == 1 )     // common overdispersion param. case
    	xi = exp( log_xi[ 0 ]->lastDraw().getScalar() );
    else if( xi_type == 0)    // group-specific overdispersion case
      xi = exp( log_xi[ var_index ]->lastDraw().getScalar() );

    sum = 0.0;
    for ( i = 0; i < exposures[ var_index ]->Len(); i++ )
    {
      numer = xi + ( b_denom[ var_index ]->Val( i ) - xi ) * exp( xb.Val( i ) );
      ratio1 = numer /  b_denom[ var_index ]->Val( i );

      //printf(" [%d: e:%f r:%f bd:%f] ", i, exposures[ var_index ]->Val( i ), ratio1, b_denom[ var_index ]->Val( i ) );

      sum += ( xi + counts[ var_index ]->Val( i ) ) * log( ratio1 );
    }

    ratio = zxb - sum;

#ifdef DEBUGBETA
    printf("\nlogRatio: var_index = %d, xi = %f, zxb = %f, sum = %f, likelihood ratio = %f\n", var_index, xi, zxb, sum, ratio );
    printf("b_star and b_0 = \n" );
    b_star.Print();
    b_0.Print();  
    printf("b_denom[%d]: \n", var_index);  b_denom[ var_index ]->Print();  
    fflush(stdout);
#endif
    //bdiff.Print();
    //xb.Print();

  }//end if
  else if ( ((int) predictors[0]->getScalar()) == 1 )
  {
    //gamma
    CVector bdiff( gamma_dim );

    bdiff = b_star - b_0;
    for ( j = 0; j < number_of_groups; j++ )
    {
      CVector xb( counts[ j ]->Len() );

      xb = predictors[ j + 2 ]->getMatrix() * bdiff;
      zxb = (*counts[ j ]) * xb;

      if( xi_type == 1 )     // common overdispersion param. case
        xi = exp( log_xi[ 0 ]->lastDraw().getScalar() );
      else if( xi_type == 0)    // group-specific overdispersion case
        xi = exp( log_xi[ j ]->lastDraw().getScalar() );
      sum = 0.0;
      for ( i = 0; i < exposures[ j ]->Len(); i++ )
      {
        numer = xi + ( b_denom[ j ]->Val( i ) - xi ) * exp( xb.Val( i ) );
        ratio1 = numer /  b_denom[ j ]->Val( i );
        sum += ( xi + counts[ j ]->Val( i ) ) * log( ratio1 );
      }//end for i

      ratio += (zxb - sum);
    }//end for j
#ifdef DEBUGGAMMA
    printf("\nlogRatio: xi = %f, zxb = %f, sum = %f, likelihood ratio = %f\n", xi, zxb, sum, ratio );
    printf("b_star and b_0 = \n" );
    b_star.Print();
    b_0.Print();  fflush(stdout);
#endif

  }//end if
  else
  {
    printf( " BayesianPoissonLogLinkModel::logRatioTargetDensity: Wrong argument type [%d].\n", ((int) predictors[0]->getScalar())  );
  }
  return ( ratio );
}//end


// Calculates the log likelihood ratio for a proposed change in beta (the random effects)
// or gamma (the fixed effects), for the case where there is no overdispersion.
double BayesianPoissonLogLinkModel::logRatioTargetNoOverdisperse( 
	DistributionParameter ** predictors, 
	CVector & b_star, CVector & b_0 )
{
  int i, j, var_index;
  double sum, ratio, zxb, mu_0, mu_star; // removed unused numer, ratio1, xi. -wwd.

  ratio = 0.0;
  if ( ((int) predictors[0]->getScalar()) == 0 )
  {
    //beta
    var_index = ((int) predictors[1]->getScalar());

    CVector bdiff( beta_dim );
    CVector xb( counts[ var_index ]->Len() );

    bdiff = b_star - b_0;
    xb = predictors[ var_index + 2 ]->getMatrix() * bdiff;
    #ifdef DEBUGBETA
      printf("xb = \n");  xb.Print();  
      printf("counts[ var_index ] \n");  counts[ var_index ]->Print();  fflush(stdout);
    #endif
    zxb = (*counts[ var_index ]) * xb;
    
    sum = 0.0;
    for ( i = 0; i < exposures[ var_index ]->Len(); i++ )
    {
      mu_0 = b_denom[ var_index ]->Val( i ) / exposures[ var_index ]->Val( i );
      mu_star = mu_0 * exp( xb.Val( i ) );

      sum += exposures[ var_index ]->Val( i ) * ( mu_star - mu_0 );
    }

    ratio = zxb - sum;

#ifdef DEBUGBETA
    printf("\nlogRatio: var_index = %d, zxb = %f, sum = %f, likelihood ratio = %f\n", var_index, zxb, sum, ratio );
    printf("b_star and b_0 = \n" );
    b_star.Print();
    b_0.Print();  
    printf("b_denom[%d]: \n", var_index);  b_denom[ var_index ]->Print();  
    fflush(stdout);
#endif
    //bdiff.Print();
    //xb.Print();

  }//end if
  else if ( ((int) predictors[0]->getScalar()) == 1 )
  {
    //gamma
    CVector bdiff( gamma_dim );

    bdiff = b_star - b_0;
    
    for ( j = 0; j < number_of_groups; j++ )
    {
      CVector xb( counts[ j ]->Len() );

      xb = predictors[ j + 2 ]->getMatrix() * bdiff;
      zxb = (*counts[ j ]) * xb;
      #ifdef DEBUGGAMMA
        printf("counts[ j ] \n");  counts[ j ]->Print();  
        printf("xb = \n");  xb.Print();  printf("zxb = %f\n", zxb );  
        fflush(stdout);
      #endif

      sum = 0.0;
      for ( i = 0; i < exposures[ j ]->Len(); i++ )
      {
        if( exposures[ j ]->Val( i ) > 0 ){
        	mu_0 = b_denom[ j ]->Val(i) / exposures[ j ]->Val( i );
          mu_star = mu_0 * exp( xb.Val( i ) );
          sum += exposures[ j ]->Val( i ) * ( mu_star - mu_0 );
          #ifdef DEBUGGAMMA
            printf("mu_0 = %f, mu_star = %f, exposure[j,i] = %f, \n", mu_0, mu_star, 
              exposures[ j ]->Val( i ) );
          #endif
        }
      }//end for i

      ratio += (zxb - sum);
    }//end for j
    #ifdef DEBUGGAMMA
      printf("\nlogRatio: zxb = %f, sum = %f, likelihood ratio = %f\n", zxb, sum, ratio );
      printf("b_star and b_0 = \n" );
      b_star.Print();
      b_0.Print();  fflush(stdout);
    #endif

  }//end if
  else
  {
    printf( " BayesianPoissonLogLinkModel::logRatioTargetDensity: Wrong argument type [%d].\n", ((int) predictors[0]->getScalar())  );
  }
  return ( ratio );
}//end


// Calculates the log ratio target density for a proposed change in beta (or gamma?)
double BayesianPoissonLogLinkModel::logRatioTargetDensity( int n_pars, 
	DistributionParameter ** pars_pred, DistributionParameter ** pars_vector )
{
  double ratio;

  ratio = 0;

  if ( n_pars > 0 )
  {
    if ( n_pars == 1 )
    {
      // Get the ratio for a proposed change in xi
      ratio = logRatioTargetDensity( ((int) pars_pred[0]->getScalar() ) );
    }
    else 
    {
    	// Get the ratio for a proposed change in gamma / beta
      CVector tmpvec1 = pars_vector[0]->getVector(), tmpvec2 = pars_vector[1]->getVector();      
      if( xi_type == 2 ){    // no overdispersion case
        ratio = logRatioTargetNoOverdisperse( pars_pred, tmpvec1, tmpvec2 );
      } else {
        ratio = logRatioTargetDensity( pars_pred, tmpvec1, tmpvec2 );
      }
    }//end
  }
  else
  {
    printf( " BayesianPoissonLogLinkModel::logRatioTargetDensity: No arguments passed.\n" );
  }

  return ( ratio );

}//end


// Calculate the ratio of the proposal density q(xi_star | xi) to 
// q(xi | xi_star).  The proposal is normal on the log-xi space, so 
// the proposal on the xi space is not normal.  We must account for
// the asymmetry, so the log ratio is not zero.
// If (xi_type == 2), so that there is no overdispersion parameter, this function never 
// gets called.
double BayesianPoissonLogLinkModel::logRatioProposal( int var_index )
{	
	if( xi_type == 2 ){
		MESSAGE "BayesianPoissonLogLinkModel::logRatioProposal: There is no parameter to calculate\n" ERROR;
	}
	if( (xi_type == 1) && (var_index > 0) )
		MESSAGE "logRatioProposal: var_index out of bounds" ERROR;
		
  double xi_star = exp( log_xi[ var_index ]->lastDraw().getScalar() ),
    xi = exp( log_xi[ var_index ]->mean().getScalar() );

	double ratio = log( xi / xi_star );
	
	return( ratio );
}//end
