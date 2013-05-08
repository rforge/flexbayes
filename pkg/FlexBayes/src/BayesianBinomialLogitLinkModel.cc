#include "R.h"
#include "Rmath.h"

#include "Const.h"
#include "DistributionParameter.h"
#include "BayesianBinomialLogitLinkModel.h"


/* constructor.
   returns: an empty Distribution of type BayesianBinomialLogitLinkModel
*/
BayesianBinomialLogitLinkModel::BayesianBinomialLogitLinkModel()
{
  emptyModel();
}//end


void BayesianBinomialLogitLinkModel::emptyModel()
{
  trials = NULL;
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

  max_trials = NULL;
  n_equal = NULL;
  n_larger_equal = NULL;
  log_xi_hessian = NULL;
  b1 = NULL;
  b2 = NULL;
  d = NULL;
 
  number_of_simulations = 0;
  simulated_xi = NULL;
  simulated_theta = NULL;

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


BayesianBinomialLogitLinkModel::~BayesianBinomialLogitLinkModel()
{
  int i;

  if ( trials != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete trials[i];
    }
    delete [] trials;
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

  if ( max_trials != NULL )
  {
    delete max_trials;
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

  if ( b1 != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete b1[i];
    }
    delete [] b1;
  }

  if ( b2 != NULL )
  {
    for ( i = 0; i < number_of_groups; i++ )
    {
      delete b2[i];
    }
    delete [] b2;
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
      
  if ( simulated_theta != NULL )
  {
    delete simulated_theta;
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

}//end


void BayesianBinomialLogitLinkModel::initialize( CVector ** count_v, CVector ** n_trials, int n_groups  )
{
  int i;

  number_of_groups = n_groups;
  number_of_observations = 0;

  counts = new CVector * [ number_of_groups ];
  trials = new CVector * [ number_of_groups ];
  mu_linear = new CVector * [ number_of_groups ];

  for ( i = 0; i < number_of_groups; i++ )
  {
    counts[i] = new CVector( (*count_v[i]) );
    trials[i] = new CVector( (*n_trials[i]) );

    mu_linear[i] = new CVector( count_v[i]->Len() );

    number_of_observations += trials[i]->Len();
  }

	// The following creates objects that are large enough for the case of 
	// separate xi parameters for each group.  In the case that there is a
	// single xi parameter, only the first element of this vector will ever
	// be used.
  computed_starting_hessian_logXi = new CVector( number_of_groups );
  computed_starting_hessian_logXi->setToZero();

  computed_starting_logXi_variance = new CVector( number_of_groups );
  computed_starting_logXi_variance->setToZero();

  logXi_variance = new CVector( number_of_groups );
  logXi_variance->setTo( 1.0 );

}//end


void BayesianBinomialLogitLinkModel::initializeCountStorage()
{
  int i;

  max_trials = new CVector( number_of_groups );
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

  d = new CVector * [ number_of_groups ];
  b1 = new CVector * [ number_of_groups ];
  b2 = new CVector * [ number_of_groups ];
  for ( i = 0; i < number_of_groups; i++ )
  {
    d[i] = new CVector( trials[i]->Len() );
    b1[i] = new CVector( trials[i]->Len() );
    b2[i] = new CVector( trials[i]->Len() );
  }  

}//end



void BayesianBinomialLogitLinkModel::logXiPrior( double p_xi )
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


void BayesianBinomialLogitLinkModel::logXiPrior( double * p_xi )
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

void BayesianBinomialLogitLinkModel::samplerDefaultInitialPoint()
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


void BayesianBinomialLogitLinkModel::samplerXiInitialPoint( double init_xi )
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


void BayesianBinomialLogitLinkModel::samplerXiInitialPoint( CVector & init_xi )
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



void BayesianBinomialLogitLinkModel::setCoefficientDimension( int p_random, int p_fixed )
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



void BayesianBinomialLogitLinkModel::considerCounts()
{
  int j, k;

  initializeCountStorage();
  
  for ( j = 0; j < number_of_groups; j++ )
  {
    max_trials->Val( j ) = trials[ j ]->max();

    n_equal[j] = new CVector( ((int) max_trials->Val( j )) );
    n_larger_equal[j] = new CVector( ((int) max_trials->Val( j )) );
    
    for ( k = ((int) max_trials->Val( j )); k >= 1; k-- )
    {
      n_equal[ j ]->Val( k - 1 ) = trials[ j ]->howManyEqual( ((double) k) );
      if ( k < ((int) max_trials->Val( j )) )
      {
        n_larger_equal[ j ]->Val( k - 1 ) = n_equal[ j ]->Val( k - 1 ) + n_larger_equal[ j ]->Val( k );
      }
      else
      {
        n_larger_equal[ j ]->Val( k - 1 ) = n_equal[ j ]->Val( k - 1 );
      }
    }//end for k
  }//end for j

  if( xi_type == 1 ){  // common xi parameter for all groups
  	// Pool the counts into a single max_trials value and single n_larger_equal vector.
  	// First find the overall max trials
  	double mxTrial = max_trials->max();
  	// Create a vector that contains the total numbers of each trials value
  	CVector *nLargerEqual = new CVector( (int) mxTrial );
  	nLargerEqual->setToZero();
  	for( j = 0; j < number_of_groups; j++ ){
  	  for( k = 0; k < n_larger_equal[ j ]->Len(); k++ ){
  		  nLargerEqual->Val( k ) = nLargerEqual->Val( k ) + n_larger_equal[ j ]->Val( k );
  	  }
  	}
  	max_trials->Val( 0 ) = mxTrial;
        if (n_larger_equal[0] != NULL) {
          delete n_larger_equal[0] ;
        }
  	n_larger_equal[ 0 ] = nLargerEqual;
  } 


}//end



void BayesianBinomialLogitLinkModel::metropolisHastingsUpdateModel( DistributionParameter * par_vals )
{

  if ( par_vals->length() == 1 )
  {
    metropolisHastingsUpdateModel( ((int) par_vals->getScalar()) );
  }
  else
  {
    CVector tmpvec = par_vals->getVector().subVector( 1, par_vals->length() - 1 );
    metropolisHastingsUpdateModel( ((int) par_vals->getVector().Val( 0 )), tmpvec );
  }
}//end



void BayesianBinomialLogitLinkModel::metropolisHastingsUpdateModel( int j, CVector & mu )
{
  int i, k, total;
  double xi, sum1, sum2, wx, w1, w2;

  if( xi_type == 1 )  // common overdispersion parameter case
    xi = exp( log_xi[ 0 ]->lastDraw().getScalar() );
  else if( xi_type == 0 )   // group-specific overdispersion parameter case
    xi = exp( log_xi[ j ]->lastDraw().getScalar() );
  else          // no overdispersion case
  	xi = 1.0;  // This is just a placeholder; the calculations using xi never get used

  #ifdef DEBUGGLM
    printf("mhUpdateModel: enter j = %d, xi = %f, mu = \n", j, xi  );  
    mu.Print();  fflush(stdout);
  #endif

  for ( i = 0; i < mu_linear[ j ]->Len(); i++ )
  {
    wx = exp( mu.Val( i ) );
    mu_linear[ j ]->Val( i ) = wx / ( 1 + wx );
    w1 = xi * mu_linear[ j ]->Val( i );
    w2 = xi * ( 1 - mu_linear[ j ]->Val( i ) );

    //update diagonal covariance for hessians of beta and gamma
    sum1 = 0;
    total = ((int) counts[ j ]->Val( i ));
    if ( total > 0 )
    {
      for ( k = 0; k < total - 1; k++ )
      {
        sum1 += ( 1 / ( ( k + w1 ) * ( k + w1 ) ) );
      }
    }
    b1[ j ]->Val( i ) = sum1;

    sum2 = 0;
    total = ((int) ( trials[ j ]->Val( i ) - counts[ j ]->Val( i ) ) );
    if ( total > 0 )
    {
      for ( k = 0; k < total - 1; k++ )
      {
        sum2 += ( 1 / ( ( k + w2 ) * ( k + w2 ) ) );
      }
    }
    b2[ j ]->Val( i ) = sum2;

    d[ j ]->Val( i )  = xi * xi * ( sum1 + sum2 );
  }

  //printf("mhUpdateModel: mu_linear = and d = \n");
  //mu_linear[j]->Print();
  //d[j]->Print();

  //update (gradient and) hessian of logXi
  if ( ( computed_starting_hessian_logXi->Val( j ) == 0 ) && 
  	( xi_type != 2 ) ) {  // check for overdispersion parameters
  
    computeHessianForLogXi( j, xi );
  }

}//end


void BayesianBinomialLogitLinkModel::metropolisHastingsUpdateModel( int j )
{
  int i, k, total;
  double xi, sum1, sum2, w1, w2;

  if( xi_type == 1 )  // common overdispersion parameter case
    xi = exp( log_xi[ 0 ]->lastDraw().getScalar() );
  else if( xi_type == 0 )   // group-specific overdispersion parameter case
  	xi = exp( log_xi[ j ]->lastDraw().getScalar() );
  else          // no overdispersion case
  	xi = 1.0;  // this is just a placeholder; the calculations using xi never
  	           // get used 

  #ifdef DEBUGGLM
    printf("mhUpdateModel: enter j = %d, xi = %f, mu_linear = \n", j, xi  );  
    mu_linear[j]->Print();  fflush(stdout);
  #endif
  for ( i = 0; i < mu_linear[ j ]->Len(); i++ )
  {
    w1 = xi * mu_linear[ j ]->Val( i );
    w2 = xi * ( 1 - mu_linear[ j ]->Val( i ) );

    //update diagonal covariance for hessians of beta and gamma
    sum1 = 0;
    total = ((int) counts[ j ]->Val( i ));
    if ( total > 0 )
    {
      for ( k = 0; k < total - 1; k++ )
      {
        sum1 += ( 1 / ( ( k + w1 ) * ( k + w1 ) ) );
      }
    }
    b1[ j ]->Val( i ) = sum1;

    sum2 = 0;
    total = ((int) ( trials[ j ]->Val( i ) - counts[ j ]->Val( i ) ) );
    if ( total > 0 )
    {
      for ( k = 0; k < total - 1; k++ )
      {
        sum2 += ( 1 / ( ( k + w2 ) * ( k + w2 ) ) );
      }
    }
    b2[ j ]->Val( i ) = sum2;

    d[ j ]->Val( i )  = xi * xi * ( sum1 + sum2 );
  }

  //printf("mhUpdateModel: d = \n");
  //d[j]->Print();

  //update (gradient and) hessian of logXi
  if ( ( computed_starting_hessian_logXi->Val( j ) == 0 ) && 
  	( xi_type != 2 ) ) {  // check for overdispersion parameters
  
    computeHessianForLogXi( j, xi );
  }

}//end


CVector BayesianBinomialLogitLinkModel::computeWeights( int var_index, DistributionParameter ** predictors, CVector & b_star, CVector & b_0 )
{
  int i, k, total;
  double xi, wx, proposal_mu_linear, proposal_b2, proposal_b1, sum1, sum2;

  CVector proposal_d( trials[ var_index ]->Len() );

  CVector bdiff( b_star - b_0 );
  CVector xb( predictors[ var_index + 2 ]->getMatrix() * bdiff );
   
  if( xi_type == 1 )  // common overdispersion param. case
    xi = exp( log_xi[ 0 ]->lastDraw().getScalar() );
  else if( xi_type == 0 )   // group-specific overdispersion param. case
    xi = exp( log_xi[ var_index ]->lastDraw().getScalar() );
  else     // no overdispersion param. case.  Crease a pretend overdispersion parameter for this calculation.
  	xi = 1;     // Need a less arbitrary choice here 
  	
  for ( i = 0; i < trials[ var_index ]->Len(); i++ )
  {
    wx = exp( xb.Val( i ) );
    proposal_mu_linear = ( mu_linear[ var_index ]->Val( i ) * wx ) / ( 1 - mu_linear[ var_index ]->Val( i ) * ( 1 - wx ) );

    proposal_b2 = xi * ( 1 - proposal_mu_linear );
    proposal_b1 = xi * proposal_mu_linear;

    //update diagonal covariance for hessians of beta and gamma
    sum1 = 0;
    total = ((int) counts[ var_index ]->Val( i ));
    if ( total > 0 )
    {
      for ( k = 0; k < total - 1; k++ )
      {
        sum1 += ( 1 / ( ( k + proposal_b1 ) * ( k + proposal_b1 ) ) );
      }
    }
    
    sum2 = 0;
    total = ((int) ( trials[ var_index ]->Val( i ) - counts[ var_index ]->Val( i ) ) );
    if ( total > 0 )
    {
      for ( k = 0; k < total - 1; k++ )
      {
        sum2 += ( 1 / ( ( k + proposal_b2 ) * ( k + proposal_b2 ) ) );
      }
    }

    proposal_d.Val( i )  = xi * xi * ( sum1 + sum2 );
  }//end for i

  return (proposal_d);

}//end




double BayesianBinomialLogitLinkModel::logXiGradient( int j, double xi )
{
  int i, k, total;
  double sum1, sum2, sum3, sum4, w1, w2, gradient;

  //update gradient of logXi
  sum1 = 0.0;
  for ( k = 0; k < ((int) max_trials->Val(j)); k++ )
  {
    sum1 += ( n_larger_equal[j]->Val(k) / ( xi + k ) );
  }

  sum4 = 0.0;
  for ( i = 0; i < trials[j]->Len(); i++ )
  {
    w1 = xi * mu_linear[ j ]->Val( i );
    w2 = xi * ( 1 - mu_linear[ j ]->Val( i ) );

    sum2 = 0;
    total = ((int) counts[ j ]->Val( i ));
    if ( total > 0 )
    {
      for ( k = 0; k < total - 1; k++ )
      {
        sum2 += ( 1 / ( k + w1 ) );
      }
      sum2 *= mu_linear[ j ]->Val( i );
    }
    
    sum3 = 0;
    total = ((int) ( trials[ j ]->Val( i ) - counts[ j ]->Val( i ) ) );
    if ( total > 0 )
    {
      for ( k = 0; k < total - 1; k++ )
      {
        sum3 += ( 1 / ( k + w2 ) );
      }
      sum3 *= ( 1 - mu_linear[ j ]->Val( i ) );
    }

    sum4 += ( sum2 + sum3 );
  }

  //printf("in logXi gradient: xi = %f, sum1 = %f, sum2 = %f, sum3 = %f,  sum4 = %f\n", xi, sum1, sum2, sum3, sum4 ); fflush(stdout);

  gradient = ( sum4 - sum1 ) * xi;

  return (gradient);

}//end



void BayesianBinomialLogitLinkModel::computeHessianForLogXi(  int j, double xi )
{
  if( xi_type == 2 ){
    printf( "computeHessianForLogXi: Attempt to compute Hessian for nonexistent parameter xi[%d]\n", j );  fflush(stdout);
  	MESSAGE "" ERROR;
  } 
  
  int i, k;
  double sum1, sum2, sum3;

  //printf( "in compute Hessian for logXi: j = %d, xi = %f\n", j, xi ); fflush(stdout);

  sum1 = 0.0;
  for ( k = 0; k < ((int) max_trials->Val(j)); k++ )
  {
    sum1 += ( n_larger_equal[j]->Val(k)  / ( ( xi + k ) * ( xi + k ) ) );
  }
  
  sum2 = 0.0;
  sum3 = 0.0;
  for ( i = 0; i < trials[j]->Len(); i++ )
  {
    sum2 += ( mu_linear[j]->Val(i) * mu_linear[j]->Val(i) *  b1[j]->Val(i) );
    sum3 += ( (1 - mu_linear[j]->Val(i)) * (1 - mu_linear[j]->Val(i)) *  b2[j]->Val(i) );
  }

  //printf("in compute Hessian logXi: sum1 = %f, sum2 = %f, sum3 = %f \n", sum1, sum2, sum3 ); fflush(stdout);
  //printf( "in compute Hessian logXi: b[%d] = \n", j);
  //b[j]->Print();
  double z0;
  if( xi_type == 1 )   // common overdispersion parameter case
  	z0 = xi_z0[ 0 ];
  else                 // group-specific overdispersion parameter case
  	z0 = xi_z0[ j ];
  log_xi_hessian->Val(j) = ( sum1 - sum2 - sum3 ) * xi * xi - 2 * ( ( xi * z0 ) / ( ( xi + z0 ) * ( xi + z0 ) ) );
  log_xi_hessian->Val(j) += logXiGradient( j, xi );

  computed_starting_hessian_logXi->Val( j ) = 1;

}//end



CMatrix BayesianBinomialLogitLinkModel::metropolisHastingsProposalHessian( int n_pars, 
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
            printf( " BayesianBinomialLogitLinkModel::metropolisHastingsProposalHessian: Wrong number of parameters for Gamma [%d].\n", n_pars );
          }
        }//end gamma
      }//end if not null
      else
      {
        printf( " BayesianBinomialLogitLinkModel::metropolisHastingsProposalHessian: Null array.\n" );
      }
    }//end if par_vector
    else
    {
      printf( " BayesianBinomialLogitLinkModel::metropolisHastingsProposalHessian: Null array.\n" );
    }
  }//end if > 0
  else
  {
    printf( " BayesianBinomialLogitLinkModel::metropolisHastingsProposalHessian: No parameters passed.\n" );
  }

  CMatrix null_hessian( 1, 1 );
  null_hessian.setToZero();
  return ( null_hessian );
}//end



CVector BayesianBinomialLogitLinkModel::metropolisHastingsMean( int n_pars, 
	DistributionParameter ** par_vals, DistributionParameter ** par_vector )
{
  //not used

  CVector null_vector( 1 );
  null_vector.setToZero();
  return ( null_vector );

}//end




CMatrix BayesianBinomialLogitLinkModel::metropolisHastingsHessian( int n_pars, DistributionParameter ** par_vals )
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
          printf( " BayesianBinomialLogitLinkModel::metropolisHastingsHessian: Wrong number of parameters for Gamma [%d].\n", n_pars );
        }
      }//end gamma
    }//end if not null
    else
    {
      printf( " BayesianBinomialLogitLinkModel::metropolisHastingsHessian: Null array.\n" );
    }
  }//end if > 0
  else
  {
    printf( " BayesianBinomialLogitLinkModel::metropolisHastingsHessian: No parameters passed.\n" );
  }

  CMatrix null_hessian( 1, 1 );
  null_hessian.setToZero();
  return ( null_hessian );
}//end


// This function never gets used
double BayesianBinomialLogitLinkModel::hessianForlogXi( int i )
{
  if( xi_type == 2 ){
    printf( "hessianForLogXi: Attempt to compute Hessian for nonexistent parameter xi[%d]\n", i );  fflush(stdout);
  	MESSAGE "" ERROR;
  } 
  if( xi_type == 1 )  // common overdispersion param. case
    return ( log_xi_hessian->Val(0) );
  else if( xi_type == 0 )    //  group-specific overdispersion param. case
    return ( log_xi_hessian->Val(i) );

}//end


// If (xi_type == 2), so that there is no overdispersion parameter, this function never gets called.
void BayesianBinomialLogitLinkModel::metropolisHastingsUpdateLogXi( int j )
{
  double xi, denom, var_log_xi;

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
    denom = log_xi_hessian->Val(j);
    if ( denom < 0.0 && ( denom <= - 0.1 ) )   // the upper bound for denom has been changed from -DELTA
    {
      var_log_xi = - 1.0 / denom;
    }
    else 
    {
      double w;

      //too small or negative variance: use gradient estimate of variance
      w =  logXiGradient( j, xi )  - ( ( 2 * xi ) / ( xi + xi_z0[ j ] ) );
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
    logXi_variance->Val( j ) = var_log_xi;
    computed_starting_logXi_variance->Val( j ) = 1;

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
    printf( "BayesianBinomialLogitLinkModel::metropolisHastingsUpdateLogXi: Keeping previous covariance. Current one is %f\n", 
      xi_var.Val(0,0) );
    fflush( stdout );
  }

  #ifdef DEBUGGLM
    printf( "log xi proposal mean = %f, var = %f\n", xi_mean.Val(0), xi_var.Val(0,0) );  fflush(stdout);
  #endif
  
}//end


void BayesianBinomialLogitLinkModel::simulationsToArray( double * simul_output, 
	int start_index, CMatrix * simulated_mu )
{
  int j, s, start_dim, simulations_to_keep;

  simulations_to_keep = simulated_mu->Row();
  start_dim = start_index;
  if( xi_type == 1 ){  // common overdispersion param. case
    for ( s = 0; s < simulations_to_keep; s++ )
    {
      simul_output[ start_dim + s ] = simulated_xi->Val( s, 0 );
    }  	
    start_dim += simulations_to_keep;
  } else if( xi_type == 0){   // group-specific overdispersion param. case
    for ( j = 0; j < number_of_groups; j++ )
    {
      for ( s = 0; s < simulations_to_keep; s++ )
      {
        simul_output[ start_dim + j * simulations_to_keep + s ] = simulated_xi->Val( s, j );
        //printf( "simul xi[%d] = %f\n", s, simulated_xi->Val( s, j ) ); fflush(stdout);
      }
    }
    start_dim += number_of_groups * simulations_to_keep;
  }
  // If there is overdispersion, sample the individual thetas
  if( xi_type != 2 ) {
    generateThetas( simulated_mu );
    //keep thetas
    for ( j = 0; j < number_of_observations; j++ )
    {
      for ( s = 0; s < simulations_to_keep; s++ )
      {
        simul_output[ start_dim + j * simulations_to_keep + s ] = simulated_theta->Val( s, j );
      }
    }
  }

}//end


// If (xi_type == 2) then this function should never be called since there is no 
// overdispersion.
void BayesianBinomialLogitLinkModel::generateThetas( CMatrix * simulated_mu )
{
	if( xi_type == 2 ){
		printf("BayesianBinomialLogitLinkModel::generateThetas: There are no theta parameters to sample\n");
		return;
	}
  int i, j, k, s, simulations_to_keep;
  double xi, wx, wmu, w1, w2;

  simulations_to_keep = simulated_mu->Row();
  k = 0;
  for ( j = 0; j < number_of_groups; j++ )
  {
    //get alpha and beta
    for ( i = 0; i < trials[j]->Len(); i++ )
    {
      for ( s = 0; s < simulations_to_keep; s++ )
      {
      	if( xi_type == 1 )  // common overdispersion param. case
          xi = simulated_xi->Val( s, 0 );
        else               // group-specific overdispersion case
        	xi = simulated_xi->Val( s, j );
        wx = exp( simulated_mu->Val( s, k ) );
        wmu = wx / ( 1 + wx );
        w1 = counts[ j ]->Val( i ) + xi * wmu;
        w2 = trials[ j ]->Val( i ) - counts[ j ]->Val( i ) + xi * ( 1 - wmu );

        simulated_theta->Val( s, k ) = rbeta( w1, w2 );

        //printf( "generate thetas: k = %d, s = %d, xi = %f, wx = %f, wmu = %f, w1 = %f, w2 = %f, theta = %f\n", k, s, xi, wx, wmu, w1, w2, simulated_theta->Val( s, k ) ); fflush(stdout);

      }//end for s
      k++;
    }//end for i
  }//end for j

  //printf("Thetas are\n" );
  //simulated_theta->Print();

}//end



void BayesianBinomialLogitLinkModel::keepSimulation( int simul_number )
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



void BayesianBinomialLogitLinkModel::createOutput( int sim_to_keep )
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
    simulated_theta = new CMatrix( number_of_simulations, number_of_observations );
    
}//end


void BayesianBinomialLogitLinkModel::initializeTemporaryStructures()
{
  //nothing to do here
  //except computing starting (and perhaps fixed) hessian for xi
  
  

}//end


void BayesianBinomialLogitLinkModel::dataAugmentationInitialDraws()
{
  //nothing to do here
}//end


// If (xi_type == 2), so that there is no overdispersion parameter, this function never gets called
void BayesianBinomialLogitLinkModel::drawVariableFromProposal( int var_index )
{
	if( xi_type == 2 ){
		printf("BayesianBinomialLogitLinkModel::drawVariableFromProposal: There is no parameter to draw\n");
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
    printf( " BayesianBinomialLogitLinkModel::drawVariableFromProposal: Wrong argument index [%d].\n", var_index );
  }
}//end


// If (xi_type == 2), so that there is no overdispersion parameter, this function never gets called
void BayesianBinomialLogitLinkModel::updateVariableForProposal( int var_index )
{
	if( xi_type == 2 ){
		printf("BayesianBinomialLogitLinkModel::updateVariableForProposal: There is no parameter to update\n");
		return;
	}
  if ( var_index < number_of_variables )
  {
    metropolisHastingsUpdateLogXi( var_index );
  }
  else
  {
    printf( " BayesianBinomialLogitLinkModel::updateVariableForProposal: Wrong argument index [%d].\n", var_index );
  }
}//end



// If (xi_type == 2), so that there is no overdispersion parameter, this function never 
// gets called.
void BayesianBinomialLogitLinkModel::setCurrentVariableFromProposal( int var_index )
{
  if( xi_type == 2 ){
		printf("BayesianBinomialLogitLinkModel::updateVariableForProposal: There is no parameter to update\n");
		return;
	}
  if ( var_index < number_of_variables )
  {
  	// We have changed the value of xi so we need to update b_denom, b, etc.
  	// If we have updated the common xi parameter, then we need to update
  	// b1, b2, d, etc. for all groups, and otherwise we just need to update them for the
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
    printf( " BayesianBinomialLogitLinkModel::setCurrentVariableFromProposal: Wrong argument index [%d].\n", var_index );
  }
}//end


// If (xi_type == 2), so that there is no overdispersion parameter, this function never 
// gets called.
void BayesianBinomialLogitLinkModel::keepCurrentVariable( int var_index )
{
	if( xi_type == 2 ){
		printf("BayesianBinomialLogitLinkModel::keepCurrentVariable: There is no parameter to keep\n");
		return;
	}
  if ( var_index < number_of_variables )
  {
  	#ifdef DEBUG1
  	  printf("rejected proposed change in xi.  restoring xi[var_index] = %f\n", 
  	    exp( log_xi[ var_index ]->mean().getScalar() ) );  fflush(stdout);
  	#endif
  	
#ifdef FIX1
    DistributionParameter tmpvec = log_xi[ var_index ]->mean();
    log_xi[ var_index ]->setLastDraw( tmpvec );
#else
    log_xi[ var_index ]->setLastDraw( log_xi[ var_index ]->mean() );
#endif
//    log_xi[ var_index ]->setLastDraw( log_xi[ var_index ]->mean() );
  }
  else
  {
    printf( " BayesianBinomialLogitLinkModel::keepCurrentVariable: Wrong argument index [%d].\n", var_index );
  }
}//end


// Calculates the log ratio target density for a proposed change in xi.  The ratio is
// P(new) / P(old).
// If (xi_type == 2), so that there is no overdispersion parameter, this function never 
// gets called.
double BayesianBinomialLogitLinkModel::logRatioTargetDensity( int var_index )
{
	if( xi_type == 2 ){
		MESSAGE "BayesianBinomialLogitLinkModel::logRatioTargetDensity: There is no parameter to calculate\n" ERROR;
	}
	if( (xi_type == 1) && (var_index > 0) )  // common overdispersion parameter case
		MESSAGE "logRatioTargetDensity: var_index out of bounds" ERROR;

  int i, j, k, total;
  double ratio, sum1, sum2, sum3, sum4, xi, xi_star, w1, w2, w1_star, w2_star;

  ratio = 0;
  if ( var_index < number_of_variables )
  {
    xi_star = exp( log_xi[ var_index ]->lastDraw().getScalar() );
    xi = exp( log_xi[ var_index ]->mean().getScalar() );
#ifdef DEBUG1
  printf("xi = %f \n", xi );
  printf("xi_star = %f \n", xi_star );
  printf("xi_z0[var_index] = %f \n", xi_z0[var_index] );
#endif
    sum1 = 0.0;
    for ( k = 0; k < ((int) max_trials->Val( var_index )); k++ )
    {
      sum1 += n_larger_equal[ var_index ]->Val( k ) * log( ( xi_star + k ) / ( xi + k ) );
    }

    sum4 = 0.0;
#ifdef DEBUG1
    printf("n_larger_equal = \n");
    n_larger_equal[ var_index ]->Print();
    printf("mu_linear = \n");
    mu_linear[ var_index ]->Print();
    printf("max_trials = %f \n", max_trials->Val( var_index ) );
#endif
    if( xi_type == 0 ){
	    for ( i = 0; i < trials[ var_index ]->Len(); i++ )
	    {
	      w1 = xi * mu_linear[ var_index ]->Val( i );
	      w2 = xi * ( 1 - mu_linear[ var_index ]->Val( i ) );
	      w1_star = xi_star * mu_linear[ var_index ]->Val( i );
	      w2_star = xi_star * ( 1 - mu_linear[ var_index ]->Val( i ) );
	
	      sum2 = 0;
	      total = ((int) counts[ var_index ]->Val( i ));
	      if ( total > 0 )
	      {
	        for ( k = 0; k < total; k++ )   // Changed by dwoodard from (total-1)
	        {
	          sum2 += log( ( k + w1_star ) / ( k + w1 ) );
	        }
	      }
	    
	      sum3 = 0;
	      total = ((int) ( trials[ var_index ]->Val( i ) - counts[ var_index ]->Val( i ) ) );
	      if ( total > 0 )
	      {
	        for ( k = 0; k < total; k++ )    // Changed by dwoodard from (total-1)
	        {
	          sum3 += log( ( k + w2_star ) / ( k + w2 ) );
	        }
	      }
	
	      sum4 += ( sum2 + sum3 );
	    }
    } else {
    	for ( j = 0; j < number_of_groups; j++ ){
		    for ( i = 0; i < trials[ j ]->Len(); i++ )
		    {
		      w1 = xi * mu_linear[ j ]->Val( i );
		      w2 = xi * ( 1 - mu_linear[ j ]->Val( i ) );
		      w1_star = xi_star * mu_linear[ j ]->Val( i );
		      w2_star = xi_star * ( 1 - mu_linear[ j ]->Val( i ) );
		
		      sum2 = 0;
		      total = ((int) counts[ j ]->Val( i ));
		      if ( total > 0 )
		      {
		        for ( k = 0; k < total; k++ )   // Changed by dwoodard from (total-1)
		        {
		          sum2 += log( ( k + w1_star ) / ( k + w1 ) );
		        }
		      }
		    
		      sum3 = 0;
		      total = ((int) ( trials[ j ]->Val( i ) - counts[ j ]->Val( i ) ) );
		      if ( total > 0 )
		      {
		        for ( k = 0; k < total; k++ )    // Changed by dwoodard from (total-1)
		        {
		          sum3 += log( ( k + w2_star ) / ( k + w2 ) );
		        }
		      }
		
		      sum4 += ( sum2 + sum3 );
		    }
		  }
    }
    
    ratio = sum4 - sum1 - 2 * log( ( xi_star + xi_z0[ var_index ] ) / ( xi + xi_z0[ var_index ] ) );
  }
  else
  {
    printf( " BayesianBinomialLogitLinkModel::logRatioTargetDensity: Wrong argument index [%d].\n", var_index );
  }
#ifdef DEBUG1
  printf("sum1 = %f, sum4 = %f \n", sum1, sum4);
  printf("ratio = %f \n", ratio);
#endif

  return ( ratio );

}//end


// Calculates the log likelihood ratio for a proposed change in beta (the random effects)
// or gamma (the fixed effects), for the case where there is overdispersion.
double BayesianBinomialLogitLinkModel::logRatioTargetDensity( DistributionParameter ** predictors, 
	CVector & b_star, CVector & b_0 )
{
  int i, j, k, total, var_index;
  double sum1, sum2, sum3, ratio, xi, wx, w1, w2, w1_star, w2_star, mu_star;

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
    if( xi_type == 1 )     // common overdispersion param. case
    	xi = exp( log_xi[ 0 ]->lastDraw().getScalar() );
    else if( xi_type == 0)    // group-specific overdispersion case
      xi = exp( log_xi[ var_index ]->lastDraw().getScalar() );

    sum1 = 0.0;
    for ( i = 0; i < trials[ var_index ]->Len(); i++ )
    {
      wx = exp( xb.Val( i ) );
      mu_star = ( mu_linear[ var_index ]->Val( i ) * wx ) / ( 1 - mu_linear[ var_index ]->Val( i ) * ( 1 - wx ) );

      w1 = xi * mu_linear[ var_index ]->Val( i );
      w2 = xi * ( 1 - mu_linear[ var_index ]->Val( i ) );

      w1_star = xi * mu_star;
      w2_star = xi * ( 1 - mu_star );

      sum2 = 0;
      total = ((int) counts[ var_index ]->Val( i ));
      if ( total > 0 )
      {
        for ( k = 0; k < total - 1; k++ )
        {
          sum2 += log( ( k + w1_star ) / ( k + w1 ) );
        }
      }
    
      sum3 = 0;
      total = ((int) ( trials[ var_index ]->Val( i ) - counts[ var_index ]->Val( i ) ) );
      if ( total > 0 )
      {
        for ( k = 0; k < total - 1; k++ )
        {
          sum3 += log( ( k + w2_star ) / ( k + w2 ) );
        }
      }

      sum1 += ( sum2 + sum3 );
    }

    ratio = sum1;

#ifdef DEBUG1
    printf("\nlogRatio: var_index = %d, xi = %f, likelihood ratio = %f\n", var_index, xi, ratio );
    printf("b_star and b_0 = \n" );
    b_star.Print();
    b_0.Print();
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

      if( xi_type == 1 )     // common overdispersion param. case
        xi = exp( log_xi[ 0 ]->lastDraw().getScalar() );
      else if( xi_type == 0)    // group-specific overdispersion case
        xi = exp( log_xi[ j ]->lastDraw().getScalar() );
      sum1 = 0.0;
      for ( i = 0; i < trials[ j ]->Len(); i++ )
      {
        wx = exp( xb.Val( i ) );
        mu_star = ( mu_linear[ j ]->Val( i ) * wx ) / ( 1 - mu_linear[ j ]->Val( i ) * ( 1 - wx ) );

        w1 = xi * mu_linear[ j ]->Val( i );
        w2 = xi * ( 1 - mu_linear[ j ]->Val( i ) );

        w1_star = xi * mu_star;
        w2_star = xi * ( 1 - mu_star );

        sum2 = 0;
        total = ((int) counts[ j ]->Val( i ));
        if ( total > 0 )
        {
          for ( k = 0; k < total - 1; k++ )
          {
            sum2 += log( ( k + w1_star ) / ( k + w1 ) );
          }
        }
    
        sum3 = 0;
        total = ((int) ( trials[ j ]->Val( i ) - counts[ j ]->Val( i ) ) );
        if ( total > 0 )
        {
          for ( k = 0; k < total - 1; k++ )
          {
            sum3 += log( ( k + w2_star ) / ( k + w2 ) );
          }
        }

        sum1 += ( sum2 + sum3 );
      }//end for i

      ratio += sum1;
    }//end for j
  }//end if
  else
  {
    printf( " BayesianBinomialLogitLinkModel::logRatioTargetDensity: Wrong argument type [%d].\n", ((int) predictors[0]->getScalar())  );
  }
  return ( ratio );
}//end



// Calculates the log likelihood ratio for a proposed change in beta (the random effects)
// or gamma (the fixed effects), for the case where there is no overdispersion.
double BayesianBinomialLogitLinkModel::logRatioTargetNoOverdisperse( DistributionParameter ** predictors, 
	CVector & b_star, CVector & b_0 )
{
  int i, j, var_index;
  double sum, ratio, wx, yxb, exb_0, exb_star;

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
    yxb = (*counts[ var_index ]) * xb;
    
    sum = 0.0;
    for ( i = 0; i < trials[ var_index ]->Len(); i++ )
    {
      if( trials[ var_index ]->Val( i ) > 0 ){
        wx = exp( xb.Val( i ) );
        exb_0 = mu_linear[ var_index ]->Val( i ) / ( 1 - mu_linear[ var_index ]->Val( i ) );
        exb_star = exb_0 * exp( xb.Val( i ) );
      	sum += trials[ var_index ]->Val( i ) * log( (1 + exb_star) / (1 + exb_0) );
      }
    }

    ratio = yxb - sum;

#ifdef DEBUG1
    printf("\nlogRatio: var_index = %d, likelihood ratio = %f\n", var_index, ratio );
    printf("b_star and b_0 = \n" );
    b_star.Print();
    b_0.Print();
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
      yxb = (*counts[ j ]) * xb;
      #ifdef DEBUGGAMMA
        printf("counts[ j ] \n");  counts[ j ]->Print();  
        printf("xb = \n");  xb.Print();  printf("yxb = %f\n", yxb );  
        fflush(stdout);
      #endif
      
      sum = 0.0;
      for ( i = 0; i < trials[ j ]->Len(); i++ )
      {
      	if( trials[ j ]->Val( i ) > 0 ){
        	wx = exp( xb.Val( i ) );
        	exb_0 = mu_linear[ j ]->Val( i ) / ( 1 - mu_linear[ j ]->Val( i ) );
          exb_star = exb_0 * exp( xb.Val( i ) );
      		sum += trials[ j ]->Val( i ) * log( (1 + exb_star) / (1 + exb_0) );
        }
      }//end for i

      ratio += yxb - sum;
    }//end for j
  }//end if
  else
  {
    printf( " BayesianBinomialLogitLinkModel::logRatioTargetDensity: Wrong argument type [%d].\n", ((int) predictors[0]->getScalar())  );
  }
  return ( ratio );
}//end



// Calculates the log ratio target density for a proposed change in beta (or gamma?)
double BayesianBinomialLogitLinkModel::logRatioTargetDensity( int n_pars, 
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
    printf( " BayesianBinomialLogitLinkModel::logRatioTargetDensity: No arguments passed.\n" );
  }

  return ( ratio );

}//end

// Calculate the ratio of the proposal density q(xi_star | xi) to 
// q(xi | xi_star).  The proposal is normal on the log-xi space, so 
// the proposal on the xi space is not normal.  We must account for
// the asymmetry, so the log ratio is not zero.
double BayesianBinomialLogitLinkModel::logRatioProposal( int var_index )
{	
	if( xi_type == 2 ){
		MESSAGE "BayesianBinomialLogitLinkModel::logRatioProposal: There is no parameter to calculate\n" ERROR;
	}
	if( (xi_type == 1) && (var_index > 0) )
		MESSAGE "logRatioProposal: var_index out of bounds" ERROR;
		
  double xi_star = exp( log_xi[ var_index ]->lastDraw().getScalar() ),
    xi = exp( log_xi[ var_index ]->mean().getScalar() );

	double ratio = log( xi / xi_star );
	
	return( ratio );
}//end
