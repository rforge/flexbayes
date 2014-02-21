#include "R.h"
#include "Rmath.h"

#include "Const.h"
#include "DistributionExtended.h"


double MATH_PI = M_PI;

/* InvChisqDistribution implementation */

InvChisqDistribution::InvChisqDistribution()
{
  dfreedom_init = 1;
  scale_init = 1;
  dfreedom_updated = 1;
  scale_updated = 1;
  last_item_drawn = 1;
}//end


InvChisqDistribution::InvChisqDistribution( double nu0 )
{
  dfreedom_init = nu0;
  scale_init = 1;
  dfreedom_updated = nu0;
  scale_updated = 1;
  last_item_drawn = 1;
}//end


InvChisqDistribution::~InvChisqDistribution()
{
  //do nothing
}//end


InvChisqDistribution & InvChisqDistribution::operator =(const InvChisqDistribution & distr)
{
  if ( this == &distr )
  {
    return *this;
  }

  dfreedom_init = distr.dfreedom_init;
  scale_init = distr.scale_init;
  dfreedom_updated = distr.dfreedom_updated;
  scale_updated = distr.scale_updated;
  last_item_drawn = distr.last_item_drawn;

  return *this;

}//end



InvChisqDistribution::InvChisqDistribution( const InvChisqDistribution & distr )
{
  dfreedom_init = distr.dfreedom_init;
  scale_init = distr.scale_init;
  dfreedom_updated = distr.dfreedom_updated;
  scale_updated = distr.scale_updated;
  last_item_drawn = distr.last_item_drawn;
}//end copy constructor


Distribution * InvChisqDistribution::clone()
{

  InvChisqDistribution * distr = new InvChisqDistribution( (*this) );
  return distr;
}


void InvChisqDistribution::setDegreesOfFreedom( double nu0 )
{
  dfreedom_init = nu0;
  dfreedom_updated = nu0;
}//end


void InvChisqDistribution::setScale( double sigma02 )
{
  scale_init = sigma02;
  scale_updated = scale_init ;
  last_item_drawn = scale_init;

}//end



void InvChisqDistribution::update( DistributionParameter * par_list, int list_size )
{

  if ( list_size == 2)
  {
    //first parameter is additional_df, second parameter is additional_scale
    update( par_list[0].getScalar(), par_list[1].getScalar() );
  }
  else
  {
    Rprintf( "InvChisqDistribution::update: wrong number of paramaters [%d].\n", list_size );
  }
}//end


void InvChisqDistribution::update( double additional_df, double additional_scale )
{
  if ( dfreedom_init > 0 )
  {
    dfreedom_updated = dfreedom_init + additional_df;
    scale_updated = ( dfreedom_init * scale_init + additional_scale )/ dfreedom_updated;
  }
  else
  {
    dfreedom_updated = additional_df;
    scale_updated = additional_scale / dfreedom_updated;
  }
}//end


void InvChisqDistribution::setLastDraw( DistributionParameter & last_draw )
{

  setLastDraw( last_draw.getScalar() );
 
}//end


void InvChisqDistribution::setLastDraw( double last_draw )
{
  last_item_drawn = last_draw;

}//end


DistributionParameter InvChisqDistribution::draw()
{
  DistributionParameter last_draw ( drawOneItem() );

  return last_draw;
}//end


double InvChisqDistribution::drawOneItem()
{
  double x;

  if ( dfreedom_updated >= NON_INFORMATIVE_SIGMA_DF )
  {
    x = rchisq( dfreedom_updated );
    last_item_drawn = dfreedom_updated * scale_updated / x;
  }
  else
  {
    Rprintf( "InvChisqDistribution::drawOneItem: degrees of freedom [%f] parameter is too small. Using instead [%f] df.\n", dfreedom_updated, NON_INFORMATIVE_SIGMA_DF );
    x = rchisq( NON_INFORMATIVE_SIGMA_DF );
    last_item_drawn = NON_INFORMATIVE_SIGMA_DF * scale_updated / x;
  }

  return last_item_drawn;
}//end


DistributionParameter InvChisqDistribution::lastDraw()
{
  DistributionParameter last_draw( last_item_drawn );

  return last_draw;
}//end
  

DistributionParameter InvChisqDistribution::mode()
{
  double val;

  val = ( dfreedom_updated / ( dfreedom_updated + 2.0 ) ) * scale_updated;
  DistributionParameter mode_val( val );

  return mode_val;
}//end


DistributionParameter InvChisqDistribution::mean()
{
  double val;

  if ( dfreedom_updated > 2 + DELTA )
  {
    val = ( dfreedom_updated / ( dfreedom_updated - 2.0 ) ) * scale_updated;
  }
  else
  {
    Rprintf( "InvChisqDistribution::mean: degrees of freedom [%f < 2.0]. Mean is infinity.\n", dfreedom_updated );
    val = 1.0 / DELTA;
  }

  DistributionParameter mean_val( val );

  return mean_val;
}//end


DistributionParameter InvChisqDistribution::variance()
{
  double val;
  
  if ( dfreedom_updated > 4 + DELTA )
  {
    val = dfreedom_updated / ( dfreedom_updated - 2.0 );
    val = ( 2.0 * val * val / ( dfreedom_updated - 4.0 ) ) * ( scale_updated * scale_updated );
  }
  else
  {
    Rprintf( "InvChisqDistribution::variance: degrees of freedom [%f < 4.0]. Variance is infinity.\n", dfreedom_updated );
    val = 1.0 / DELTA;
  }

  DistributionParameter var( val );

  return var;
}//end


double InvChisqDistribution::logDensity( DistributionParameter & value )
{
  double vD2, val, x;

  x = value.getScalar();
  vD2 = dfreedom_updated / 2.0;
  val = vD2 * log( vD2 ) - lgammafn( vD2 ) + vD2 * log( scale_updated )
      - ( vD2 + 1.0 ) * log( x ) - ( vD2 * scale_updated / x );

  return val;
}//end





/* NormalDistribution implementation */

NormalDistribution::NormalDistribution( int dim_vector )
{
  dim = dim_vector;
  non_informative = false;
  mean_init = new CVector( dim );
  cov_init = new CMatrix( dim, dim );

  very_first_mean = new CVector( dim );
  very_first_cov = new CMatrix( dim, dim );
  very_first_inv_cov = new CMatrix( dim, dim );

  inv_cov_init = new CMatrix( dim, dim );
  inv_cov_beta_init = new CVector( dim );
  mean_updated = new CVector( dim );
  inv_determinant_init = 0;

  cholesky_invCov_updated = new CVector( ( dim * ( dim + 1 ) ) / 2  );

  last_item_drawn = new CVector( dim );
}//end
  

NormalDistribution::~NormalDistribution()
{
  delete very_first_mean;
  delete very_first_cov;
  delete very_first_inv_cov;
  delete mean_init;
  delete cov_init;
  delete inv_cov_init;
  delete inv_cov_beta_init;
  delete mean_updated;
  delete cholesky_invCov_updated;
  delete last_item_drawn;
}//end


void NormalDistribution::clean()
{
  delete very_first_mean;
  delete very_first_cov;
  delete very_first_inv_cov;
  delete mean_init;
  delete cov_init;
  delete inv_cov_init;
  delete inv_cov_beta_init;
  delete mean_updated;
  delete cholesky_invCov_updated;
  delete last_item_drawn;
}//end



NormalDistribution & NormalDistribution::operator =(const NormalDistribution & distr)
{
  if ( this == &distr )
  {
    return *this;
  }

  clean();

  dim = distr.dim;
  non_informative = distr.non_informative;

  mean_init = new CVector( dim );
  (*mean_init) = ( *(distr.mean_init) );

  very_first_mean = new CVector( dim );
  (*very_first_mean) = ( *(distr.very_first_mean) );

  cov_init = new CMatrix( dim, dim );
  (*cov_init) = ( *(distr.cov_init) );

  very_first_cov = new CMatrix( dim, dim );
  (*very_first_cov) = ( *(distr.very_first_cov) );

  very_first_inv_cov = new CMatrix( dim, dim );
  (*very_first_inv_cov) = ( *(distr.very_first_inv_cov) );

  inv_cov_init = new CMatrix( dim, dim );
  (*inv_cov_init) = ( *(distr.inv_cov_init) );

  inv_cov_beta_init = new CVector( dim );
  (*inv_cov_beta_init) = ( *(distr.inv_cov_beta_init) );
  inv_determinant_init = distr.inv_determinant_init;

  mean_updated = new CVector( dim );
  (*mean_updated) = ( *(distr.mean_updated) );

  cholesky_invCov_updated = new CVector( ( dim * ( dim + 1 ) ) / 2  );
  (*cholesky_invCov_updated) = ( *(distr.cholesky_invCov_updated) );

  last_item_drawn = new CVector( dim );
  (*last_item_drawn) = ( *(distr.last_item_drawn) );

  return *this;
}//end



NormalDistribution::NormalDistribution( const NormalDistribution & distr )
{
  dim = distr.dim;
  non_informative = distr.non_informative;

  mean_init = new CVector( dim );
  (*mean_init) = ( *(distr.mean_init) );

  very_first_mean = new CVector( dim );
  (*very_first_mean) = ( *(distr.very_first_mean) );

  cov_init = new CMatrix( dim, dim );
  (*cov_init) = ( *(distr.cov_init) );

  very_first_cov = new CMatrix( dim, dim );
  (*very_first_cov) = ( *(distr.very_first_cov) );

  very_first_inv_cov = new CMatrix( dim, dim );
  (*very_first_inv_cov) = ( *(distr.very_first_inv_cov) );

  inv_cov_init = new CMatrix( dim, dim );
  (*inv_cov_init) = ( *(distr.inv_cov_init) );

  inv_cov_beta_init = new CVector( dim );
  (*inv_cov_beta_init) = ( *(distr.inv_cov_beta_init) );
  inv_determinant_init = distr.inv_determinant_init;

  mean_updated = new CVector( dim );
  (*mean_updated) = ( *(distr.mean_updated) );

  cholesky_invCov_updated = new CVector( ( dim * ( dim + 1 ) ) / 2  );
  (*cholesky_invCov_updated) = ( *(distr.cholesky_invCov_updated) );

  last_item_drawn = new CVector( dim );
  (*last_item_drawn) = ( *(distr.last_item_drawn) );

}//end copy constructor


Distribution * NormalDistribution::clone()
{
  NormalDistribution * distr = new NormalDistribution ( (*this) );
  return distr;
}//end


void NormalDistribution::setMean( CVector * mean )
{
  (*mean_init) = (*mean);
  (*very_first_mean) = (*mean);
  (*mean_updated) = (*mean);
  (*last_item_drawn) = (*mean);
}//end


void NormalDistribution::setCovariance( CMatrix * cov ) throw( rtErr )
{
  non_informative = false;

  (*cov_init) = (*cov);

  (*very_first_cov) = (*cov);

  //update initial inverse covariance
  (*inv_cov_init) = cov_init->inverse();

  (*very_first_inv_cov) = (*inv_cov_init);

  //update initial product cov * beta
  (*inv_cov_beta_init) = (*inv_cov_init) * (*mean_init);


  //ready from random draw
  try
  {
    (*cholesky_invCov_updated) = inv_cov_init->choleskyDecomposition();
    inv_determinant_init = invCovDeterminant();
  }
  catch( rtErr choldcError )
  {
    Rprintf( "NormalDistribution::setCovariance: Covariance matrix does not seem to be symmetric positive definite.\n" );
    char the_error[] = "NormalDistribution::setCovariance: Covariance matrix does not seem to be symmetric positive definite.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

  //printf("normal: cov init is \n"); cov_init->Print();
  //printf(" inv cov init is \n"); inv_cov_init->Print();
  //printf(" variance is \n"); covariance().Print();
  //printf(" inv conv init inverse is \n"); inv_cov_init->inverse().Print();

}//end

void NormalDistribution::setCovariance( CMatrix * cov, CVector * cholesky_invCov )
{
  non_informative = false;

  (*cov_init) = (*cov);

  (*very_first_cov) = (*cov);

  //ready from random draw
  (*cholesky_invCov_updated) = (*cholesky_invCov);
  inv_determinant_init = invCovDeterminant();

  //update initial inverse covariance
  CMatrix inverse_cholLT( cov->Row(), cov->Col() );
  inverse_cholLT.assignInverseOfLowerTriangular( (*cholesky_invCov) );
  (*inv_cov_init) = inverse_cholLT.T() * inverse_cholLT;

  (*very_first_inv_cov) = (*inv_cov_init);

  //update initial product cov * beta
  (*inv_cov_beta_init) = (*inv_cov_init) * (*mean_init);

}//end



void NormalDistribution::scaleInitialMean( double scaleMult )
{
  (*mean_init) = (*very_first_mean) * scaleMult;
  //update initial product cov * beta
  (*inv_cov_beta_init) = (*inv_cov_init) * (*mean_init);

  (*mean_updated) = (*mean_init);
  (*last_item_drawn) = (*mean_init);
}//end 



void NormalDistribution::updateInitialMean( CVector * mean )
{
  (*mean_init) = (*mean);
  //update initial product cov * beta
  (*inv_cov_beta_init) = (*inv_cov_init) * (*mean_init);

  (*mean_updated) = (*mean);
  (*last_item_drawn) = (*mean);
}//end


void NormalDistribution::updateInitialCovariance( CMatrix * cov ) throw( rtErr )
{
  non_informative = false;

  (*cov_init) = (*cov);

  //update initial inverse covariance
  (*inv_cov_init) = cov_init->inverse();

  //update initial product cov * beta
  (*inv_cov_beta_init) = (*inv_cov_init) * (*mean_init);

  //ready from random draw
  try
  {
    (*cholesky_invCov_updated) = inv_cov_init->choleskyDecomposition();
    inv_determinant_init = invCovDeterminant();
  }
  catch( rtErr choldcError )
  {
    Rprintf( "NormalDistribution::updateInitialCovariance: Covariance matrix does not seem to be symmetric positive definite.\n" );

    char the_error[] = "NormalDistribution::updateInitialCovariance: Covariance matrix does not seem to be symmetric positive definite.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

}//end


void NormalDistribution::updateInitialCovariance( CMatrix * cov, CVector * cholesky_invCov )
{
  non_informative = false;

  (*cov_init) = (*cov);

  //ready from random draw
  (*cholesky_invCov_updated) = (*cholesky_invCov);
  inv_determinant_init = invCovDeterminant();

  //update initial inverse covariance
  CMatrix inverse_cholLT( cov->Row(), cov->Col() );
  inverse_cholLT.assignInverseOfLowerTriangular( (*cholesky_invCov) );
  (*inv_cov_init) = inverse_cholLT.T() * inverse_cholLT;

  //( (*inv_cov_init) * (*cov_init) ).Print();


  //update initial product cov * beta
  (*inv_cov_beta_init) = (*inv_cov_init) * (*mean_init);

}//end




void NormalDistribution::updateVeryFirstCovariance( CMatrix & cov )
{
  (*very_first_cov) = cov;
  (*very_first_inv_cov) = cov.inverse();
}//end



void NormalDistribution::setNonInformative()
{
  non_informative = true;
  mean_init->setToZero();
  (*very_first_mean) = (*mean_init);
  (*mean_updated) = (*mean_init);
  (*last_item_drawn) = (*mean_init);

  CMatrix cov( dim, dim );
  cov.setDiagonal( 1.0 );
  //setCovariance( &cov );

  (*cov_init) = cov;

  (*very_first_cov) = cov;

  //update initial inverse covariance
  (*inv_cov_init) = cov_init->inverse();

  (*very_first_inv_cov) = (*inv_cov_init);

  //update initial product cov * beta
  (*inv_cov_beta_init) = (*inv_cov_init) * (*mean_init);


  //ready from random draw
  (*cholesky_invCov_updated) = inv_cov_init->choleskyDecomposition();
  inv_determinant_init = invCovDeterminant();

  
}//end


void NormalDistribution::setLastDraw( DistributionParameter & last_beta )
{
  setLastDraw( last_beta.getVector() );
}//end


void NormalDistribution::setLastDraw( const CVector & last_beta )
{
  (*last_item_drawn) = last_beta;
}//end


void NormalDistribution::scaleLastDraw( double scale )
{
  if ( scale > 0 )
  {
    CVector transformed_vector( dim );

    transformed_vector = ( (*last_item_drawn) * scale ) + ( (*mean_updated) * (1.0 - scale) );    
    (*last_item_drawn) = transformed_vector;
  }
  else
  {
    Rprintf( " NormalDistribution::scaleLastDraw: negative or zero scale [%f] is not allowed.", scale );
  }
}


void NormalDistribution::updateScale( double scale )
{
  if ( scale > 0 )
  {
    cholesky_invCov_updated->multiplyByScalar( 1.0/ scale );
  }
  else
  {
    Rprintf( " NormalDistribution::updateScale: negative or zero scale is not allowed." );
  }
}//end


void NormalDistribution::update( DistributionParameter * par_list, int list_size )
{

  if ( list_size == 2 )
  {
    //posterior update formulas when beta is a normal
#ifdef FIX1
    CVector tmpvec = par_list[0].getVector();
    CMatrix tmpmat = par_list[1].getMatrix();
    update( tmpvec , tmpmat );
#else
    update( par_list[0].getVector(), par_list[1].getMatrix() );
#endif
//    update( par_list[0].getVector(), par_list[1].getMatrix() );
  }
  else if ( list_size == 3 )
  {
    //posterior update formulas when beta is a Student-t (mix of normal and Inv-Chisq)
#ifdef FIX1
    CVector tmpvec = par_list[1].getVector();
    CMatrix tmpmat = par_list[2].getMatrix();
    update( par_list[0].getScalar() , tmpvec , tmpmat );
#else
    update( par_list[0].getScalar(), par_list[1].getVector(), par_list[2].getMatrix() );
#endif
//    update( par_list[0].getScalar(), par_list[1].getVector(), par_list[2].getMatrix() );
  }

}//end


//posterior update formulas when beta is a normal
void NormalDistribution::update( CVector & mean_shift, CMatrix & cov_shift )
{
  //printf("update beta non informative is %d\n", non_informative);

  if ( non_informative )
  {
    //printf( "update beta: cov shift is \n");
    //cov_shift.Print();
    //printf( "update beta: mean shift is \n");
    //mean_shift.Print();    

    (*cholesky_invCov_updated) = cov_shift.choleskyDecomposition();

    (*mean_updated) = cholesky_invCov_updated->asCholeskyDecompositionSolve( mean_shift );
  }
  else
  {
    //printf( "cov shift is \n" ); fflush(stdout);  cov_shift.Print(); fflush(stdout);
    //printf( "inv cov init is \n" ); fflush(stdout);  inv_cov_init->Print(); fflush(stdout);
    //printf( "mean shift is \n" ); fflush(stdout);  mean_shift.Print(); fflush(stdout);
    //printf( "inv cov beta init is \n" ); fflush(stdout);  inv_cov_beta_init->Print(); fflush(stdout);

    (*cholesky_invCov_updated) = ( cov_shift + (*inv_cov_init) ).choleskyDecomposition();
#ifdef DEBUG2
  Rprintf("Updating normal.  inv_cov_init = \n");
  inv_cov_init->Print();
#endif

#ifdef FIX1
    CVector tmpvec = (*inv_cov_beta_init) +  mean_shift;
    (*mean_updated) = cholesky_invCov_updated->asCholeskyDecompositionSolve( tmpvec );
#else
    (*mean_updated) = cholesky_invCov_updated->asCholeskyDecompositionSolve( (*inv_cov_beta_init) +  mean_shift );
#endif
//    (*mean_updated) = cholesky_invCov_updated->asCholeskyDecompositionSolve( (*inv_cov_beta_init) +  mean_shift );
  }

  //true parameters are
  //(*cov_updated) = ( cov_shift + (*inv_cov_init) ).Inverse();
  //(*mean_updated) = cov_updated * ( (*inv_cov_beta_init) +  mean_shift );

  (*last_item_drawn) = (*mean_updated);
}//end

//posterior update formulas when beta is a Student-t (mix of normal and Inv-Chisq)
void NormalDistribution::update( double precision_tau2, CVector & mean_shift, CMatrix & cov_shift )
{
  //prior must be proper
  (*cholesky_invCov_updated) = ( cov_shift + ( precision_tau2 * (*inv_cov_init) ) ).choleskyDecomposition();

#ifdef FIX1
  CVector tmpvec = ( (*inv_cov_beta_init) * precision_tau2 ) +  mean_shift;
  (*mean_updated) = cholesky_invCov_updated->asCholeskyDecompositionSolve( tmpvec );
#else
  (*mean_updated) = cholesky_invCov_updated->asCholeskyDecompositionSolve( ( (*inv_cov_beta_init) * precision_tau2 ) +  mean_shift );
#endif
//  (*mean_updated) = cholesky_invCov_updated->asCholeskyDecompositionSolve( ( (*inv_cov_beta_init) * precision_tau2 ) +  mean_shift );

  //true parameters are
  //(*cov_updated) = ( cov_shift + ( precision_tau2 * (*inv_cov_init) ) ).Inverse();
  //(*mean_updated) = cov_updated * ( ( precision_tau2 * (*inv_cov_beta_init) ) +  mean_shift );

  (*last_item_drawn) = (*mean_updated);
}//end



/* Calculate the posterior mean and covariance matrix for the regression coefficients
in a multiple linear regression model.  This in fact is used to update the random effect vector, conditional
on the random effect variance and other parameters.  precision_tau2 is the random effect
precision vector, while mean_shift = X'Y / sigma^2 and cov_shift = X'X / sigma^2.  Recall that
if beta | tauVec^2 ~ N( beta_0, diag(tauVec^2) ) and Y | beta, sigma^2 ~ N( X'beta, sigma^2 ), then
beta | tauVec^2, Y ~ N( beta_hat, Sigma_beta ) where Sigma_beta = [X'X/sigma^2 + diag(tauVec^2)^-1 ]^-1
and beta_hat = Sigma_beta * [ X'Y/sigma^2 + diag(tauVec^2)^-1 * beta_0 ].   Of course, in the middle
of a linear model beta_0 and tau^2 have the appropriate meanings, namely the mean structure
of the random effects and the random effect variance, respectively.  A good reference
for this calculation is p.356 in Gelman, Carlin, Stern, and Rubin (2004). */
void NormalDistribution::update( CVector & precision_tau2, CVector & mean_shift, CMatrix & cov_shift )
{
  //prior must be proper
  // Here *inv_cov_init is always the identity matrix.
  CMatrix tmpmat = (*inv_cov_init);
  for( int i=0; i < precision_tau2.Len(); i++ ){
  	tmpmat.Val(i,i) = tmpmat.Val(i,i) * precision_tau2.Val(i);
  }
  // cholesky_invCov_updated is the Cholesky decomposition of Sigma_beta^-1
  (*cholesky_invCov_updated) = ( cov_shift + tmpmat ).choleskyDecomposition();

  // Here inv_cov_beta_init is beta_0 in the above description.
  CVector tmpvec = ( tmpmat * (*inv_cov_beta_init) ) +  mean_shift;
  // The following conceptually finds Sigma_beta and multiplies
  // by tmpvec
  (*mean_updated) = cholesky_invCov_updated->asCholeskyDecompositionSolve( tmpvec );

  (*last_item_drawn) = (*mean_updated);
}//end


CMatrix NormalDistribution::covariance()
{
  CMatrix Sigma( dim, dim );
  CMatrix invCholesky( dim, dim );

  //printf(" cholesky_invCov_updated is \n"); fflush(stdout); 
  //cholesky_invCov_updated->Print(); fflush(stdout);

  invCholesky.assignInverseOfLowerTriangular( (*cholesky_invCov_updated)  );
  Sigma = invCholesky.T() * invCholesky;
  
  /*
  int i, ibase;  
  CVector b( dim );

  for ( ibase = 0; ibase < dim; ibase++ )
  {
    //get ibase-th canonical base vector
    for ( i = 0; i < dim; i++ )
    {
      b.Val( i ) = 0;
    }
    b.Val( ibase ) = 1;

    invCholesky.setColumn( ibase, cholesky_invCov_updated->asUpperTriangularSolve( b ) );
  }
  Sigma = invCholesky * invCholesky.T();
  */

  return Sigma;
}//end


double NormalDistribution::invCovDeterminant()
{

  return cholesky_invCov_updated->asCholeskyDecompositionDeterminant();
}//end


DistributionParameter NormalDistribution::draw()
{
#ifdef FIX1
  CVector tmpvec = drawOneItem();
  DistributionParameter val( tmpvec );
#else
  DistributionParameter val( drawOneItem() );
#endif
//  DistributionParameter val( drawOneItem() );

  return val;
}//end


CVector NormalDistribution::drawOneItem()
{
  int i;
  CVector standard_normal( dim );
  
  for ( i = 0; i < dim; i++ )
  {
    standard_normal.Val(i) = norm_rand();
  }

  //draw
  //CVector lTmean( cholesky_invCov_updated->asUpperTriangularMultiplyRight( *mean_updated ) );

  (*last_item_drawn) =  (*mean_updated) + cholesky_invCov_updated->asUpperTriangularSolve( standard_normal );

#ifdef DEBUG2
  Rprintf("normal draw.  Mean = \n");
  mean_updated->Print();
  Rprintf("standard normal draw = \n");
  standard_normal.Print();
  Rprintf("transformed draw = \n");
  last_item_drawn->Print();
#endif

  return (*last_item_drawn);
}//end



CVector NormalDistribution::drawItemFromPrior()
{
  int i;
  CVector standard_normal( dim );
  
  for ( i = 0; i < dim; i++ )
  {
    standard_normal.Val(i) = norm_rand();
  }

  //draw
  //CVector lTmean( cholesky_invCov_updated->asUpperTriangularMultiplyRight( *mean_updated ) );

  (*mean_updated) = (*mean_init);

  try
  {
    (*cholesky_invCov_updated) = inv_cov_init->choleskyDecomposition();
    inv_determinant_init = invCovDeterminant();
  }
  catch( rtErr choldcError )
  {
    Rprintf( "NormalDistribution::drawItemFromPrior: Covariance matrix does not seem to be symmetric positive definite.\n" );

    char the_error[] = "NormalDistribution::DrawItemFromPrior: Covariance matrix does not seem to be symmetric positive definite.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

  (*last_item_drawn) =  (*mean_updated) + cholesky_invCov_updated->asUpperTriangularSolve( standard_normal );

  return (*last_item_drawn);
}//end


DistributionParameter NormalDistribution::lastDraw()
{
#ifdef FIX1
  CVector tmpvec = lastItemDrawn();
  DistributionParameter val( tmpvec );
#else
  DistributionParameter val( lastItemDrawn() );
#endif
//  DistributionParameter val( lastItemDrawn() );

  return val;
}//end


DistributionParameter NormalDistribution::mode()
{
  DistributionParameter val( (*mean_updated) );

  return val;
}//end


DistributionParameter NormalDistribution::mean()
{
  DistributionParameter val( (*mean_updated) );

  return val;
}//end


DistributionParameter NormalDistribution::variance()
{
#ifdef FIX1
  CMatrix tmpmat = covariance();
  DistributionParameter val( tmpmat );
#else
  DistributionParameter val( covariance() );
#endif
//  DistributionParameter val( covariance() );

  return val;
}//end


double NormalDistribution::logDensity( DistributionParameter & value )
{
  double val, mD;
  CVector v( dimension() );
  CVector x( dimension() );

  x = value.getVector();

#ifdef FIX1
  CVector tmpvec = x - (*mean_updated);
  v = cholesky_invCov_updated->asUpperTriangularMultiplyRight( tmpvec );
#else
  v = cholesky_invCov_updated->asUpperTriangularMultiplyRight( x - (*mean_updated) );
#endif
//  v = cholesky_invCov_updated->asUpperTriangularMultiplyRight( x - (*mean_updated) );
  mD = v * v;

  val = - ( (double) dimension() ) * log( 2.0 * MATH_PI ) + log( invCovDeterminant() ) - mD;
  val *= 0.5;

  return val;
}//end






/* StudentTDistribution implementation */


StudentTDistribution::StudentTDistribution()
{
  inv_chisq_dist = NULL;
  normal_dist = NULL;
}//end


StudentTDistribution::StudentTDistribution( int dim_vector, double dfreedom )
{
  normal_dist = new NormalDistribution( dim_vector );
  inv_chisq_dist = new InvChisqDistribution( dfreedom );
  inv_chisq_dist->setScale( 1.0 );
}//end


StudentTDistribution::~StudentTDistribution()
{

  if ( normal_dist != NULL )
  {
    delete normal_dist;
  }

  if ( inv_chisq_dist != NULL )
  {
    delete inv_chisq_dist;
  }
}//end


void StudentTDistribution::clean()
{
  if ( normal_dist != NULL )
  {
    delete normal_dist;
  }

  if ( inv_chisq_dist != NULL )
  {
    delete inv_chisq_dist;
  }

}//end



StudentTDistribution & StudentTDistribution::operator =(const StudentTDistribution & distr )
{
  if ( this == &distr )
  {
    return *this;
  }

  clean();

  normal_dist = new NormalDistribution( ( *(distr.normal_dist) ) );
  inv_chisq_dist = new InvChisqDistribution( ( *(distr.inv_chisq_dist) ) );


  return *this;
}//end



StudentTDistribution::StudentTDistribution( const StudentTDistribution & distr )
{

  normal_dist = new NormalDistribution( ( *(distr.normal_dist) ) );
  inv_chisq_dist = new InvChisqDistribution( ( *(distr.inv_chisq_dist) ) );

}//end copy constructor


Distribution * StudentTDistribution::clone()
{
  StudentTDistribution * distr = new StudentTDistribution( (*this) );
  return distr;
}//end


double StudentTDistribution::logDensity( DistributionParameter & value )
{
  double val, mD;
  CVector v( dimension() );
  CVector x( dimension() );

  x = value.getVector();

#ifdef FIX1
  CVector tmpvec = x - mean().getVector();
  v = normal_dist->choleskyInvCov().asUpperTriangularMultiplyRight( tmpvec );
#else
  v = normal_dist->choleskyInvCov().asUpperTriangularMultiplyRight( x - mean().getVector() );
#endif
//  v = normal_dist->choleskyInvCov().asUpperTriangularMultiplyRight( x - mean().getVector() );
  mD = v * v;
  mD /= ( (double) degreesOfFreedom() );

  val = lgammafn( 0.5 * ( (double) degreesOfFreedom() + dimension() ) )
       - lgammafn( 0.5 * ( (double) degreesOfFreedom() ) );
  val -= ( 0.5 * ( (double) dimension() ) * log( ( (double) degreesOfFreedom() ) / MATH_PI ) );
  val += 0.5 * log( normal_dist->invCovDeterminant() );
  val -= 0.5* ( (double) degreesOfFreedom() + dimension() ) * log( 1 + mD );

  return val;
}//end


CVector StudentTDistribution::drawOneItem()
{
  double inv_chisq_draw = sqrt( inv_chisq_dist->drawOneItem() );
  CVector v( normal_dist->drawOneItem() - normal_dist->mean().getVector() );
  v.multiplyByScalar( inv_chisq_draw );
  return v + normal_dist->mean().getVector();
}//end


DistributionParameter StudentTDistribution::draw()
{
#ifdef FIX1
  CVector tmpvec = drawOneItem();
  DistributionParameter val( tmpvec );
#else
  DistributionParameter val( drawOneItem() );
#endif
//  DistributionParameter val( drawOneItem() );

  return val;
}//end


void StudentTDistribution::update( DistributionParameter * par_list, int list_size )
{
  double beta_dist, precision_tau2;

  //first update inv chisq
  beta_dist = M_D( normal_dist->initialInverseCovariance(), lastItemDrawn() - normal_dist->initialMean() );
  inv_chisq_dist->update( (double) dimension(), beta_dist );
  //draw a tau2
  //and then update normal
  precision_tau2 = 1.0 / inv_chisq_dist->draw().getScalar();
#ifdef FIX1
  CVector tmpvec = par_list[0].getVector();
  CMatrix tmpmat = par_list[1].getMatrix();
  normal_dist->update( precision_tau2 , tmpvec , tmpmat );
#else
  normal_dist->update( precision_tau2, par_list[0].getVector(), par_list[1].getMatrix() );
#endif
//  normal_dist->update( precision_tau2, par_list[0].getVector(), par_list[1].getMatrix() );
  
}//end


void StudentTDistribution::update( CVector & mean_shift, CMatrix & cov_shift )
{
  double beta_dist, precision_tau2;

  //first update inv chisq
  beta_dist = M_D( normal_dist->initialInverseCovariance(), lastItemDrawn() - normal_dist->initialMean() );
  inv_chisq_dist->update( (double) dimension(), beta_dist );
  //draw a tau2
  //and then update normal
  precision_tau2 = 1.0 / inv_chisq_dist->draw().getScalar();
  normal_dist->update( precision_tau2, mean_shift, cov_shift );
  
}//end


CVector StudentTDistribution::lastItemDrawn()
{
  return normal_dist->lastItemDrawn();
}//end


DistributionParameter StudentTDistribution::lastDraw()
{
  return normal_dist->lastDraw();
}//end


void StudentTDistribution::setLastDraw( DistributionParameter & par )
{
  normal_dist->setLastDraw( par );
}//end


void StudentTDistribution::setLastDraw( CVector & par )
{
  normal_dist->setLastDraw( par );
}//end


DistributionParameter StudentTDistribution::mode()
{
  return normal_dist->mode();
}//end


DistributionParameter StudentTDistribution::mean()
{
  return normal_dist->mean();
}//end


DistributionParameter StudentTDistribution::variance()
{
  double df_ratio;
  CMatrix cov( dimension(), dimension() );

  cov = normal_dist->covariance();

  //printf("T: var: df = %f, scale matrix is \n", degreesOfFreedom() );
  //cov.Print();

  if ( degreesOfFreedom() > 2.0 + DELTA )
  {
    df_ratio = degreesOfFreedom() / ( degreesOfFreedom() - 2.0 );
    cov.multiplyByScalar( df_ratio );
  }
  else
  {
    Rprintf( "StudentTDistribution::variance: degrees of freedom [%f < 2.0]. Variance is infinity.\n", degreesOfFreedom() );
    cov.multiplyByScalar( 1.0 / DELTA );
  }

  DistributionParameter var( cov );

  return var;
}//end






/* ConstrainedInvChisqDistribution implementation */

ConstrainedInvChisqDistribution::ConstrainedInvChisqDistribution()
{
  interval = NULL;
  left_bounded_only = false;
  right_bounded_only = false;
  initial_dfreedom = 0;
  initial_scale = 1.0;
  dfreedom = 0;
  scale = 1.0;
  last_item_drawn = scale;

}//end


ConstrainedInvChisqDistribution::ConstrainedInvChisqDistribution( double lower, double upper )
{
  setInterval( lower, upper );
  initial_dfreedom = 0;
  initial_scale = 1.0;
  dfreedom = 0;
  scale = 1.0;
}//end



ConstrainedInvChisqDistribution::ConstrainedInvChisqDistribution( double bound, char * type )
{

  setInterval( bound, type );

  if ( !strcmp( type, "right" ) )
  {
    last_item_drawn = bound / 2.0;
  }
  else
  {
    //unbounded case
    if ( bound > 0 )
    {
      last_item_drawn = bound + 0.5;
    }
    else
    {
      last_item_drawn = 1.0;
    }
  } 
}//end


void ConstrainedInvChisqDistribution::setInterval( double lower, double upper )
{
  interval = new CVector( 2 );
  interval->Val(0) = lower;
  interval->Val(1) = upper;
  left_bounded_only = false;
  right_bounded_only = false;

  last_item_drawn = ( upper - lower ) / 2.0;

}//end


void ConstrainedInvChisqDistribution::setInterval( double bound, char * type )
{
  interval = new CVector( 1 );
  interval->Val(0) = bound;

  if ( !strcmp( type, "left" ) )
  {
    left_bounded_only = true;
    right_bounded_only = false;
  }
  else if ( !strcmp( type, "right" ) )
  {
    right_bounded_only = true;
    left_bounded_only = false;
  }

  if ( !strcmp( type, "right" ) )
  {
    last_item_drawn = bound / 2.0;
  }
  else
  {
    //unbounded case
    if ( bound > 0 )
    {
      last_item_drawn = bound + 0.5;
    }
    else
    {
      last_item_drawn = 1.0;
    }
  } 
}//end



ConstrainedInvChisqDistribution::~ConstrainedInvChisqDistribution()
{
  if ( interval != NULL )
  {
    delete interval;
  }
}//end


void ConstrainedInvChisqDistribution::clean()
{
  if ( interval != NULL )
  {
    delete interval;
    interval = NULL;
  }

}//en


void ConstrainedInvChisqDistribution::setDegreesOfFreedom( double nu0 )
{
  initial_dfreedom = nu0;
  dfreedom = nu0;
}//end


double ConstrainedInvChisqDistribution::upperBound()
{
  double bound;

  if ( !left_bounded_only )
  {
    if ( right_bounded_only )
    {
      bound = interval->Val(0);
    }
    else
    {
      bound = interval->Val(1);
    }  
  }
  else
  {
    bound = -1;
  }

  return bound;
}//end


double ConstrainedInvChisqDistribution::lowerBound()
{
  double bound;

  if ( left_bounded_only || !right_bounded_only )
  {
    bound = interval->Val(0);
  }
  else
  {
    bound = 0;
  }

  return bound;
}//end



void ConstrainedInvChisqDistribution::setValidLastItemDrawn( double sigma02 )
{
  if ( left_bounded_only )
  {
    if ( sigma02 >= interval->Val(0) )
    {
      last_item_drawn = sigma02;
    }
  }
  else if ( right_bounded_only )
  {
    if ( sigma02 <= interval->Val(0) )
    {
      last_item_drawn = sigma02;
    }
  }
  else
  {
    last_item_drawn = sigma02;
  }
}//



void ConstrainedInvChisqDistribution::setScale( double sigma02 )
{
  initial_scale = sigma02;
  scale = initial_scale;
  setValidLastItemDrawn( sigma02 );
}//end


void ConstrainedInvChisqDistribution::update( double additional_df, double additional_scale )
{
  if ( initial_dfreedom > 0 )
  {
    dfreedom = initial_dfreedom + additional_df;
    scale = ( initial_dfreedom * initial_scale + additional_scale )/ dfreedom;
  }
  else
  {
    dfreedom = additional_df;
    scale = additional_scale / dfreedom;
  }
}//end


void ConstrainedInvChisqDistribution::update( DistributionParameter * par_list, int list_size )
{
  if ( list_size == 2)
  {
    //first parameter is additional_df, second parameter is additional_scale
    update( par_list[0].getScalar(), par_list[1].getScalar() );
  }
  else
  {
    Rprintf( "ConstrainedInvChisqDistribution::update: wrong number of paramaters [%d].\n", list_size );
  }
}//end



ConstrainedInvChisqDistribution & ConstrainedInvChisqDistribution::operator=(const ConstrainedInvChisqDistribution & distr )
{
  if ( this == &distr )
  {
    return *this;
  }

  clean();

  if ( distr.interval != NULL )
  {
    interval = new CVector( ( *(distr.interval) ) );
  }

  left_bounded_only = distr.left_bounded_only;
  right_bounded_only = distr.right_bounded_only;
  initial_dfreedom = distr.initial_dfreedom;
  initial_scale = distr.initial_scale;
  dfreedom = distr.dfreedom;
  scale = distr.scale;
  last_item_drawn = distr.last_item_drawn;

  return *this;
}//end



ConstrainedInvChisqDistribution::ConstrainedInvChisqDistribution( const ConstrainedInvChisqDistribution & distr )
{
  clean();

  if ( distr.interval != NULL )
  {
    interval = new CVector( ( *(distr.interval) ) );
  }

  left_bounded_only = distr.left_bounded_only;
  right_bounded_only = distr.right_bounded_only;
  initial_dfreedom = distr.initial_dfreedom;
  initial_scale = distr.initial_scale;
  dfreedom = distr.dfreedom;
  scale = distr.scale;
  last_item_drawn = distr.last_item_drawn;
}//end


Distribution * ConstrainedInvChisqDistribution::clone()
{
  ConstrainedInvChisqDistribution * distr = new ConstrainedInvChisqDistribution( (*this) );
  return distr;
}//end


double ConstrainedInvChisqDistribution::CDFInvChisq( double val )
{
  double cdf_value, x;

  //printf( "CDFInvChisq entered: val = %f\n", val ); fflush(stdout);

  if ( val <= 0 )
  {
    cdf_value = 0.0;
  }
  else
  {
    x = ( dfreedom * scale ) / val;
    cdf_value = 1.0 - pchisq( x, dfreedom, 1, 0 );
  }

  //printf( "CDFInvChisq exiting: cdf_val = %f\n", cdf_value ); fflush(stdout);

  return cdf_value;
}//end



double ConstrainedInvChisqDistribution::inverseCDFInvChisq( double val ) throw( rtErr )
{
  double inv_cdf_value;

  if ( val < 1.0 && val > 0 )
  {
    inv_cdf_value = ( dfreedom * scale ) / qchisq( 1.0 - val, dfreedom, 1, 0 );
  }
  else if ( val == 1.0 )
  {
    inv_cdf_value = VERY_LARGE_NUMBER;
  }
  else if ( val == 0.0 )
  {
    inv_cdf_value = 0.0;
  }
  else
  {
    Rprintf( "ConstrainedInvChisqDistribution::inverseCDFInvChisq: value provided [%f] is not valid.\n", val );
    char the_error[] = "ConstrainedInvChisqDistribution::inverseCDFInvChisq: value provided is not valid.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

  //printf( "CDFinverseCDFInvChisq entered: inv cdf val = %f\n", inv_cdf_value ); fflush(stdout);

  return inv_cdf_value;
}//end


double ConstrainedInvChisqDistribution::drawOneItem()
{
  double U, val, sample, inv_chisq_at_lower_bound, inv_chisq_at_upper_bound, area_at_bounds;

  if ( right_bounded_only )
  {
    inv_chisq_at_upper_bound = CDFInvChisq( upperBound() );
    inv_chisq_at_lower_bound = 0;
    area_at_bounds = inv_chisq_at_upper_bound;
  }
  else if ( left_bounded_only )
  {
    inv_chisq_at_lower_bound = CDFInvChisq( lowerBound() );
    area_at_bounds = 1.0 - inv_chisq_at_lower_bound;
  }
  else
  {
    inv_chisq_at_upper_bound = CDFInvChisq( upperBound() );
    inv_chisq_at_lower_bound = CDFInvChisq( lowerBound() );
    area_at_bounds =  inv_chisq_at_upper_bound - inv_chisq_at_lower_bound;
  }

  U = runif( 0.0, 1.0 );
  val = U * area_at_bounds + inv_chisq_at_lower_bound;

  sample = inverseCDFInvChisq( val );

  last_item_drawn = sample;

  return sample;  
}//end



DistributionParameter ConstrainedInvChisqDistribution::draw()
{
  double sample;

  sample = drawOneItem();
  DistributionParameter val( sample );
  last_item_drawn = sample;

  return val;

}//end


DistributionParameter ConstrainedInvChisqDistribution::lastDraw()
{
  DistributionParameter val( last_item_drawn );

  return val;

}//end


void ConstrainedInvChisqDistribution::setLastDraw( double last_draw )
{
  last_item_drawn = last_draw;
}//end


void ConstrainedInvChisqDistribution::setLastDraw( DistributionParameter & last_draw )
{
  last_item_drawn = last_draw.getScalar();
}//end


DistributionParameter ConstrainedInvChisqDistribution::mode()
{
  InvChisqDistribution inv_chisq( dfreedom );
  inv_chisq.setScale( scale );

  return inv_chisq.mode();
}//end


DistributionParameter ConstrainedInvChisqDistribution::mean()
{
  InvChisqDistribution inv_chisq( dfreedom );
  inv_chisq.setScale( scale );

  return inv_chisq.mean();
}//end


DistributionParameter ConstrainedInvChisqDistribution::variance()
{
  InvChisqDistribution inv_chisq( dfreedom );
  inv_chisq.setScale( scale );

  return inv_chisq.variance();
}//end


double ConstrainedInvChisqDistribution::logDensity( DistributionParameter & value )
{
  InvChisqDistribution inv_chisq( dfreedom );
  inv_chisq.setScale( scale );

  return inv_chisq.logDensity( value );
}//end




/*
  implementation of Wishart Distribution
*/

WishartDistribution::WishartDistribution()
{
  dim = 0;
  initial_dfreedom = 0;
  initial_w = NULL;
  inv_initial_w = NULL;
  cholesky_initial_w = NULL;
  dfreedom = 0;

  inv_w = NULL;
  cholesky_updated_inv_w = NULL;
  scale_determinant = 0;
  last_item_drawn = NULL;
}//end


WishartDistribution::WishartDistribution( double df )
{
  dim = 0;
  initial_dfreedom = df;
  initial_w = NULL;
  inv_initial_w = NULL;
  cholesky_initial_w = NULL;
  dfreedom = df;

  inv_w = NULL;
  cholesky_updated_inv_w = NULL;
  scale_determinant = 0;
  last_item_drawn = NULL;
}//end


WishartDistribution::WishartDistribution( double df, CMatrix * v )
{
  dim = v->Col();
  initial_dfreedom = df;
  initial_w = new CMatrix( (*v) );
  //update initial inverse covariance
  inv_initial_w = new CMatrix( initial_w->inverse() );

  //ready from random draw
  cholesky_initial_w = new CVector( initial_w->choleskyDecomposition() );

  dfreedom = df;

  inv_w = new CMatrix( (*inv_initial_w)  );
  cholesky_updated_inv_w = new CVector( inv_w->choleskyDecomposition() );
  scale_determinant = scaleDeterminant();
  last_item_drawn = new CMatrix( (*v) );
}//end


WishartDistribution::~WishartDistribution()
{
  clean();
}//end


void WishartDistribution::clean()
{
  if ( initial_w != NULL )
  {
    delete initial_w;
    initial_w = NULL;
  }

  if ( inv_initial_w != NULL )
  {
    delete inv_initial_w;
    inv_initial_w = NULL;
  }

  if ( cholesky_initial_w != NULL )
  {
    delete cholesky_initial_w;
    cholesky_initial_w = NULL;
  }

  if ( inv_w != NULL )
  {
    delete inv_w;
    inv_w = NULL;
  }

  if ( cholesky_updated_inv_w != NULL )
  {
    delete cholesky_updated_inv_w;
    cholesky_updated_inv_w = NULL;
  }

  if ( last_item_drawn != NULL )
  {
    delete last_item_drawn;
    last_item_drawn = NULL;
  } 

  dim = 0;
  initial_dfreedom = 0;
  scale_determinant = 0;
  dfreedom = 0;
}//end



void WishartDistribution::cleanVectorsAndMatrices()
{
  if ( initial_w != NULL )
  {
    delete initial_w;
    initial_w = NULL;
  }

  if ( inv_initial_w != NULL )
  {
    delete inv_initial_w;
    inv_initial_w = NULL;
  }

  if ( cholesky_initial_w != NULL )
  {
    delete cholesky_initial_w;
    cholesky_initial_w = NULL;
  }

  if ( inv_w != NULL )
  {
    delete inv_w;
    inv_w = NULL;
  }

  if ( cholesky_updated_inv_w != NULL )
  {
    delete cholesky_updated_inv_w;
    cholesky_updated_inv_w = NULL;
  }

  if ( last_item_drawn != NULL )
  {
    delete last_item_drawn;
    last_item_drawn = NULL;
  } 
}//end



void WishartDistribution::setDegreesOfFreedom( double df )
{
  initial_dfreedom = df;
  dfreedom = df;
}//end


void WishartDistribution::setScaleMatrix( CMatrix * v )
{
  cleanVectorsAndMatrices();

  dim = v->Col();
  initial_w = new CMatrix( (*v) );
  //update initial inverse covariance
  inv_initial_w = new CMatrix( initial_w->inverse() );

  //ready from random draw
  cholesky_initial_w = new CVector( initial_w->choleskyDecomposition() );

  inv_w = new CMatrix( (*inv_initial_w)  );
  cholesky_updated_inv_w = new CVector( inv_w->choleskyDecomposition() );
  scale_determinant = scaleDeterminant();
  last_item_drawn = new CMatrix( (*v) );
}//end  


double WishartDistribution::scaleDeterminant()
{
  return ( 1.0 / cholesky_updated_inv_w->asCholeskyDecompositionDeterminant() );
}//end


WishartDistribution::WishartDistribution( const WishartDistribution & distr )
{
  clean();
  
  dim = distr.dim;
  initial_dfreedom = distr.initial_dfreedom;
  initial_w = new CMatrix( ( *(distr.initial_w) ) );
  //update initial inverse covariance
  inv_initial_w = new CMatrix( ( *(distr.inv_initial_w) ) );

  //ready from random draw
  cholesky_initial_w = new CVector( ( *(distr.cholesky_initial_w) ) );

  dfreedom = distr.dfreedom;
  inv_w = new CMatrix( ( *(distr.inv_w) ) );
  cholesky_updated_inv_w = new CVector( ( *(distr.cholesky_updated_inv_w) ) );
  scale_determinant = distr.scale_determinant;
  last_item_drawn = new CMatrix( ( *(distr.last_item_drawn) ) );
}//end  


WishartDistribution & WishartDistribution::operator =(const WishartDistribution & distr )
{

  if ( this == &distr )
  {
    return *this;
  }

  clean();

  dim = distr.dim;
  initial_dfreedom = distr.initial_dfreedom;
  initial_w = new CMatrix( ( *(distr.initial_w) ) );
  //update initial inverse covariance
  inv_initial_w = new CMatrix( ( *(distr.inv_initial_w) ) );

  //ready from random draw
  cholesky_initial_w = new CVector( ( *(distr.cholesky_initial_w) ) );

  dfreedom = distr.dfreedom;
  inv_w = new CMatrix( ( *(distr.inv_w) ) );
  cholesky_updated_inv_w = new CVector( ( *(distr.cholesky_updated_inv_w) ) );
  scale_determinant = distr.scale_determinant;
  last_item_drawn = new CMatrix( ( *(distr.last_item_drawn) ) );

  return *this;
}//end


Distribution * WishartDistribution::clone()
{

  WishartDistribution * distr = new WishartDistribution( (*this) );
  return distr;
}//end



double WishartDistribution::logDensity( DistributionParameter & value )
{
  int i;
  double density, value_term, exp_term, const_term, scale_term;
  CMatrix val( value.getMatrix() );

  if ( val.Col() == dim && val.Row() == dim && ( val - val.T() ).isZero( 0.1 * DELTA ) )
  {
    value_term = val.Det();
    value_term = 0.5 * ( dfreedom - (dim + 1) ) * log( value_term );

    exp_term = ( (*inv_w) * val ).trace();
    exp_term /= 2.0;

    const_term = 0.5 * dim * dfreedom * log(2.0) + 0.25 * dim * (dim - 1) * log( MATH_PI );
    for ( i = 0; i < dim; i++ )
    {
      const_term += lgammafn( 0.5 * ( dfreedom - i ) );
    }

    scale_term = 0.5 * dfreedom * log( scale_determinant );

    density = - const_term - scale_term + value_term - exp_term;
  }
  else
  {
    density = 0;
  }

  return density;
}//end



CMatrix WishartDistribution::drawOneItem()
{
  int i, j;

  CMatrix trg( dim, dim );
  CMatrix U( dim, dim );

  for ( i = 0; i < dim; i++ )
  {
    trg.Val( i, i ) = sqrt( rchisq( dfreedom - i  ) );
    for ( j = 0; j < i; j++ )
    {
      trg.Val( i, j ) = norm_rand();
    }
  }

  //a wishart(dfreedom, Identity)
  U = trg * trg.T();

  trg.assignInverseOfLowerTriangular( (*cholesky_updated_inv_w) );

  (*last_item_drawn) = trg.T() * ( U * trg );

  return (*last_item_drawn);

}//end



DistributionParameter WishartDistribution::draw()
{

#ifdef FIX1
  CMatrix tmpmat = drawOneItem();
  DistributionParameter val( tmpmat );
#else
  DistributionParameter val( drawOneItem() );
#endif
//  DistributionParameter val( drawOneItem() );

  return val;
}//end



void WishartDistribution::update( double additional_df, CMatrix & scale_shift )
{
  if ( initial_dfreedom > 0 )
  {
    dfreedom = initial_dfreedom + additional_df;
    (*inv_w) = (*inv_initial_w) + scale_shift;
  }
  else
  {
    dfreedom = additional_df;
    (*inv_w) = scale_shift;
  }

  (*cholesky_updated_inv_w) = inv_w->choleskyDecomposition();
}//end



void WishartDistribution::update( DistributionParameter * par_list, int list_size )
{

  if ( list_size == 2)
  {
    //first parameter is additional_df, second parameter is additional_scale
#ifdef FIX1
    CMatrix tmpmat = par_list[1].getMatrix();
    update( par_list[0].getScalar() , tmpmat );
#else
    update( par_list[0].getScalar(), par_list[1].getMatrix() );
#endif
//    update( par_list[0].getScalar(), par_list[1].getMatrix() );
  }
  else
  {
    Rprintf( "WishartDistribution::update: wrong number of paramaters [%d].\n", list_size );
  }

}//end



CMatrix WishartDistribution::scaleMatrix()
{
  CMatrix w( dim, dim );
  CMatrix invCholesky( dim, dim );

  invCholesky.assignInverseOfLowerTriangular( (*cholesky_updated_inv_w)  );
  w = invCholesky.T() * invCholesky;

  return w;
}//end


CMatrix WishartDistribution::lastItemDrawn()
{
  return (*last_item_drawn);
}//end


DistributionParameter WishartDistribution::lastDraw()
{
  DistributionParameter val( (*last_item_drawn) );
  return val;
}//end


void WishartDistribution::setLastDraw( DistributionParameter & par )
{
  (*last_item_drawn) = par.getMatrix();
}//end


DistributionParameter WishartDistribution::mode()
{
  //this is the mean not necessarily the mode
#ifdef FIX1
  CMatrix tmpmat = dfreedom * scaleMatrix();
  DistributionParameter val( tmpmat );
#else
  DistributionParameter val( dfreedom * scaleMatrix() );
#endif
// DistributionParameter val( dfreedom * scaleMatrix() );
 return val; 
}//end



DistributionParameter WishartDistribution::mean()
{
#ifdef FIX1
  CMatrix tmpmat = dfreedom * scaleMatrix();
  DistributionParameter val( tmpmat );
#else
  DistributionParameter val( dfreedom * scaleMatrix() );
#endif
//  DistributionParameter val( dfreedom * scaleMatrix() );
  return val;
}//end



DistributionParameter WishartDistribution::variance()
{
  //this is the scale matrix not necessarily the variance
#ifdef FIX1
  CMatrix tmpmat = scaleMatrix();
  DistributionParameter val( tmpmat );
#else
  DistributionParameter val( scaleMatrix() );
#endif
//  DistributionParameter val( scaleMatrix() );
  return val;
}//end




/*
  implementing InvWishartDistribution as the inverse of a WishartDistribution
*/


InvWishartDistribution::InvWishartDistribution()
{
  dim = 0;
  initial_dfreedom = 0;
  initial_w = NULL;
  inv_initial_w = NULL;
  cholesky_initial_w = NULL;
  dfreedom = 0;

  inv_w = NULL;
  cholesky_updated_inv_w = NULL;
  scale_determinant = 0;
  last_item_drawn = NULL;
  inv_lastDraw = NULL;
}//end


InvWishartDistribution::InvWishartDistribution( double df )
{
  dim = 0;
  initial_dfreedom = df;
  initial_w = NULL;
  inv_initial_w = NULL;
  cholesky_initial_w = NULL;
  dfreedom = df;

  inv_w = NULL;
  cholesky_updated_inv_w = NULL;
  scale_determinant = 0;
  last_item_drawn = NULL;
  inv_lastDraw = NULL;
}//end


InvWishartDistribution::InvWishartDistribution( double df, CMatrix * v )
{
  dim = v->Col();
  initial_dfreedom = df;
  initial_w = new CMatrix( (*v) );
  //update initial inverse covariance
  inv_initial_w = new CMatrix( initial_w->inverse() );

  //ready from random draw
  cholesky_initial_w = new CVector( initial_w->choleskyDecomposition() );

  dfreedom = df;

  inv_w = new CMatrix( (*inv_initial_w)  );
  cholesky_updated_inv_w = new CVector( inv_w->choleskyDecomposition() );
  scale_determinant = scaleDeterminant();
  last_item_drawn = new CMatrix( ((double) (1.0 / dfreedom)) * (*inv_w) );
  inv_lastDraw = new CMatrix( last_item_drawn->inverse() );
  
}//end


InvWishartDistribution::~InvWishartDistribution()
{
  clean();
}//end


void InvWishartDistribution::clean()
{
  if ( initial_w != NULL )
  {
    delete initial_w;
    initial_w = NULL;
  }

  if ( inv_initial_w != NULL )
  {
    delete inv_initial_w;
    inv_initial_w = NULL;
  }

  if ( cholesky_initial_w != NULL )
  {
    delete cholesky_initial_w;
    cholesky_initial_w = NULL;
  }

  if ( inv_w != NULL )
  {
    delete inv_w;
    inv_w = NULL;
  }

  if ( cholesky_updated_inv_w != NULL )
  {
    delete cholesky_updated_inv_w;
    cholesky_updated_inv_w = NULL;
  }

  if ( last_item_drawn != NULL )
  {
    delete last_item_drawn;
    last_item_drawn = NULL;
  } 

  if ( inv_lastDraw != NULL )
  {
    delete inv_lastDraw;
    inv_lastDraw = NULL;
  } 


  dim = 0;
  initial_dfreedom = 0;
  scale_determinant = 0;
  dfreedom = 0;
}//end



void InvWishartDistribution::cleanVectorsAndMatrices()
{
  if ( initial_w != NULL )
  {
    delete initial_w;
    initial_w = NULL;
  }

  if ( inv_initial_w != NULL )
  {
    delete inv_initial_w;
    inv_initial_w = NULL;
  }

  if ( cholesky_initial_w != NULL )
  {
    delete cholesky_initial_w;
    cholesky_initial_w = NULL;
  }

  if ( inv_w != NULL )
  {
    delete inv_w;
    inv_w = NULL;
  }

  if ( cholesky_updated_inv_w != NULL )
  {
    delete cholesky_updated_inv_w;
    cholesky_updated_inv_w = NULL;
  }

  if ( last_item_drawn != NULL )
  {
    delete last_item_drawn;
    last_item_drawn = NULL;
  } 

  if ( inv_lastDraw != NULL )
  {
    delete inv_lastDraw;
    inv_lastDraw = NULL;
  } 
}//end



void InvWishartDistribution::setDegreesOfFreedom( double df )
{
  initial_dfreedom = df;
  dfreedom = df;
}//end


void InvWishartDistribution::setScaleMatrix( CMatrix * v )
{
  cleanVectorsAndMatrices();

  dim = v->Col();
  initial_w = new CMatrix( (*v) );
  //update initial inverse covariance
  inv_initial_w = new CMatrix( initial_w->inverse() );

  //ready from random draw
  cholesky_initial_w = new CVector( initial_w->choleskyDecomposition() );

  inv_w = new CMatrix( (*inv_initial_w)  );
  cholesky_updated_inv_w = new CVector( inv_w->choleskyDecomposition() );
  scale_determinant = scaleDeterminant();
  last_item_drawn = new CMatrix( ((double) (1.0 / dfreedom)) * (*inv_w) );
  inv_lastDraw = new CMatrix( last_item_drawn->inverse() );
}//end  


double InvWishartDistribution::scaleDeterminant()
{
  return ( 1.0 / cholesky_updated_inv_w->asCholeskyDecompositionDeterminant() );
}//end


InvWishartDistribution::InvWishartDistribution( const InvWishartDistribution & distr )
{
  clean();
  
  dim = distr.dim;
  initial_dfreedom = distr.initial_dfreedom;
  initial_w = new CMatrix( ( *(distr.initial_w) ) );
  //update initial inverse covariance
  inv_initial_w = new CMatrix( ( *(distr.inv_initial_w) ) );

  //ready from random draw
  cholesky_initial_w = new CVector( ( *(distr.cholesky_initial_w) ) );

  dfreedom = distr.dfreedom;
  inv_w = new CMatrix( ( *(distr.inv_w) ) );
  cholesky_updated_inv_w = new CVector( ( *(distr.cholesky_updated_inv_w) ) );
  scale_determinant = distr.scale_determinant;
  last_item_drawn = new CMatrix( ( *(distr.last_item_drawn) ) );
  inv_lastDraw = new CMatrix( ( *(distr.inv_lastDraw ) ) );
}//end  


InvWishartDistribution & InvWishartDistribution::operator =(const InvWishartDistribution & distr )
{

  if ( this == &distr )
  {
    return *this;
  }

  clean();

  dim = distr.dim;
  initial_dfreedom = distr.initial_dfreedom;
  initial_w = new CMatrix( ( *(distr.initial_w) ) );
  //update initial inverse covariance
  inv_initial_w = new CMatrix( ( *(distr.inv_initial_w) ) );

  //ready from random draw
  cholesky_initial_w = new CVector( ( *(distr.cholesky_initial_w) ) );

  dfreedom = distr.dfreedom;
  inv_w = new CMatrix( ( *(distr.inv_w) ) );
  cholesky_updated_inv_w = new CVector( ( *(distr.cholesky_updated_inv_w) ) );
  scale_determinant = distr.scale_determinant;
  last_item_drawn = new CMatrix( ( *(distr.last_item_drawn) ) );
  inv_lastDraw = new CMatrix( ( *(distr.inv_lastDraw) ) );

  return *this;
}//end


Distribution * InvWishartDistribution::clone()
{

  InvWishartDistribution * distr = new InvWishartDistribution( (*this) );
  return distr;
}//end



double InvWishartDistribution::logDensity( DistributionParameter & value )
{
  int i;
  double density, value_term, exp_term, const_term, scale_term;
  CMatrix val( value.getMatrix() );

  if ( val.Col() == dim && val.Row() == dim && ( val - val.T() ).isZero( 0.1 * DELTA ) )
  {
    value_term = val.Det();
    value_term = 0.5 * ( dfreedom + (dim + 1) ) * log( value_term );

    exp_term = ( scaleMatrix() * val.inverse() ).trace();
    exp_term /= 2.0;

    const_term = 0.5 * dim * dfreedom * log(2.0) + 0.25 * dim * (dim - 1) * log( MATH_PI );
    for ( i = 0; i < dim; i++ )
    {
      const_term += lgammafn( 0.5 * ( dfreedom - i ) );
    }

    scale_term = 0.5 * dfreedom * log( scale_determinant );

    density = - const_term + scale_term - value_term - exp_term; 
  }
  else
  {
    density = 0;
  }

  return density;
}//end



CMatrix InvWishartDistribution::drawOneItem() throw( rtErr )
{
  int i;

  //printf("InvWishart: draw started. dim is %d.\n", dim ); fflush(stdout);

  int k = (int) dfreedom;

  if ( k < dim )
  {
    Rprintf( "InvWishartDistribution::drawOneItem: degrees of freedom are smaller than dimension.\n" );
    char the_error[] = "InvWishartDistribution::drawOneItem: degrees of freedom are smaller than dimension.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

  
  CMatrix drawW( dim, dim );
  CVector x( dim );
  CMatrix auxX( 1, dim );

  x.setToZero();
  drawW.setToZero();
  NormalDistribution snd( dim );
  snd.setMean( &x );
  CMatrix covX( scaleMatrix() );
  snd.setCovariance( &covX );

  for ( i = 0; i < k; i++ )
  {
    snd.draw();
#ifdef FIX1
    CVector tmpvec = snd.lastItemDrawn();
    auxX.setRow( 0 , tmpvec );
#else
    auxX.setRow( 0, snd.lastItemDrawn() );
#endif
//    auxX.setRow( 0, snd.lastItemDrawn() );
#ifdef FIX1
    CMatrix tmpmat = auxX.xTransposedX();
    drawW.add( tmpmat );
#else
    drawW.add( auxX.xTransposedX() );
#endif
//    drawW.add( auxX.xTransposedX() );
  }

  (*inv_lastDraw) = drawW;
  (*last_item_drawn) = drawW.inverse();

  //printf( "invWishart: draw: scaleMatrix * dfreedom is \n" );
  //covX.Print();

  //printf( "invWishart: draw: inv last draw is \n" ); fflush(stdout);
  //inv_lastDraw->Print(); fflush(stdout);


  //printf( "invWishart: draw: last draw is \n" ); fflush(stdout);
  //last_item_drawn->Print(); fflush(stdout);  
  

  /*
  int j;
  CMatrix trg( dim, dim );
  CMatrix U_inv( dim, dim );
  CVector vTrg( ( ( dim * ( dim + 1 ) ) / 2 ) );

  for ( i = 0; i < dim; i++ )
  {
    vTrg.setLowerTriangularVal( i, i, sqrt( S_rchisq( dfreedom - i  )  ) );
    for ( j = 0; j < i; j++ )
    {
      vTrg.setLowerTriangularVal( i, j, norm_rand( S_evaluator ) );
    }
  }

  printf( "vTrg is \n" ); fflush(stdout);
  vTrg.Print(); fflush(stdout);

  //a wishart(dfreedom, Identity)
  //U = vTrg * vTrg.T();
  //But need U.inv as well as original

  CMatrix cholW( dim, dim );
  CMatrix uTrg( dim, dim );
  CMatrix U( dim, dim );

  cholW.assignLowerTriangular( (*cholesky_updated_inv_w) );
  uTrg.assignLowerTriangular( vTrg );
  U = uTrg * uTrg.T();
  (*inv_lastDraw) = cholW * ( U * cholW.T() );

  printf( "invWishart: draw: inv last draw is \n" ); fflush(stdout);
  inv_lastDraw->Print(); fflush(stdout);

  //get the inverse
  trg.assignInverseOfLowerTriangular( vTrg );
  U_inv = trg.T() * trg;

  trg.assignInverseOfLowerTriangular( (*cholesky_updated_inv_w) );

  //this is a draw from an  InvWishart( dfreedom, Scale )
  (*last_item_drawn) = trg.T() * ( U_inv * trg );

  printf( "invWishart: draw: last draw is \n" ); fflush(stdout);
  last_item_drawn->Print(); fflush(stdout);  

  CMatrix iden( dim, dim );
  iden = (*last_item_drawn) * ( (*inv_lastDraw) );

  printf( "Wishart draw * inv is \n" ); fflush(stdout);
  iden.Print(); fflush(stdout);

  //(*inv_lastDraw) = last_item_drawn->inverse();
  */

  return (*last_item_drawn);

}//end



DistributionParameter InvWishartDistribution::draw()
{

#ifdef FIX1
  CMatrix tmpmat = drawOneItem();
  DistributionParameter val( tmpmat );
#else
  DistributionParameter val( drawOneItem() );
#endif
//  DistributionParameter val( drawOneItem() );

  return val;
}//end



void InvWishartDistribution::update( double additional_df, CMatrix & scale_shift )
{
  if ( initial_dfreedom > 0 )
  {
    dfreedom = initial_dfreedom + additional_df;
    (*inv_w) = (*inv_initial_w) + scale_shift;

    //printf("invWishart: inv_w is\n" ); fflush(stdout);
    //inv_w->Print(); fflush(stdout);
  }
  else
  {
    dfreedom = additional_df;
    (*inv_w) = scale_shift;
  }

  //printf( "invWishart: will get chol decomp\n" ); fflush(stdout);

  (*cholesky_updated_inv_w) = inv_w->choleskyDecomposition();

  //printf( "invWishart: chol decomp is \n" ); fflush(stdout);
  //cholesky_updated_inv_w->Print(); fflush(stdout);

  //printf("invWishart: out of update\n"); fflush(stdout);
}//end



void InvWishartDistribution::update( DistributionParameter * par_list, int list_size )
{

  if ( list_size == 2)
  {
    //first parameter is additional_df, second parameter is additional_scale
#ifdef FIX1
    CMatrix tmpmat = par_list[1].getMatrix();
    update( par_list[0].getScalar() , tmpmat );
#else
    update( par_list[0].getScalar(), par_list[1].getMatrix() );
#endif
//    update( par_list[0].getScalar(), par_list[1].getMatrix() );
  }
  else
  {
    Rprintf( "InvWishartDistribution::update: wrong number of paramaters [%d].\n", list_size );
  }

}//end



CMatrix InvWishartDistribution::scaleMatrix()
{
  CMatrix w( dim, dim );
  CMatrix invCholesky( dim, dim );

  invCholesky.assignInverseOfLowerTriangular( (*cholesky_updated_inv_w)  );
  w = invCholesky.T() * invCholesky;

  //printf("scaleMatrix: w * inv_w is \n");
  //CMatrix aux( (*inv_w) * w );
  //aux.Print();

  return w;
}//end


CMatrix InvWishartDistribution::lastItemDrawn()
{
  return (*last_item_drawn);
}//end

CMatrix InvWishartDistribution::inverseLastItemDrawn()
{
  return (*inv_lastDraw);
}//end


DistributionParameter InvWishartDistribution::lastDraw()
{
  DistributionParameter val( (*last_item_drawn) );
  return val;
}//end


void InvWishartDistribution::setLastDraw( DistributionParameter & par )
{
  (*last_item_drawn) = par.getMatrix();
  (*inv_lastDraw) = last_item_drawn->inverse();
}//end


DistributionParameter InvWishartDistribution::mode()
{
  double denom;

  //this is the mean not necessarily the mode
  denom = dfreedom - dim - 1;
  if ( denom > 0 )
  {
    denom = 1 / denom;
#ifdef FIX1
    CMatrix tmpmat = denom * scaleMatrix();
    DistributionParameter val( tmpmat );
#else
    DistributionParameter val( denom * scaleMatrix() );
#endif
//    DistributionParameter val( denom * scaleMatrix() );
    return val; 
  }
  else
  {
    Rprintf( "InvWishartDistribution::mode: cannot compute it when degrees of freedom is not larger than the dimension + 1.\n" );
#ifdef FIX1
    CMatrix tmpmat = scaleMatrix();
    DistributionParameter val( tmpmat );
#else
    DistributionParameter val( scaleMatrix() );
#endif
//    DistributionParameter val( scaleMatrix() );
    return val;
  }
}//end



DistributionParameter InvWishartDistribution::mean()
{
  double denom;

  //this is the mean not necessarily the mode
  denom = dfreedom - dim - 1;
  if ( denom > 0 )
  {
    denom = 1 / denom;
#ifdef FIX1
    CMatrix tmpmat = denom * scaleMatrix();
    DistributionParameter val( tmpmat );
#else
    DistributionParameter val( denom * scaleMatrix() );
#endif
//    DistributionParameter val( denom * scaleMatrix() );
    return val; 
  }
  else
  {
    Rprintf( "InvWishartDistribution::mean: cannot compute it when degrees of freedom is not larger than the dimension + 1.\n" );
#ifdef FIX1
    CMatrix tmpmat = scaleMatrix();
    DistributionParameter val( tmpmat );
#else
    DistributionParameter val( scaleMatrix() );
#endif
//    DistributionParameter val( scaleMatrix() );
    return val;
  }

}//end



DistributionParameter InvWishartDistribution::variance()
{
  //this is the scale matrix not necessarily the variance
#ifdef FIX1
  CMatrix tmpmat = scaleMatrix();
  DistributionParameter val( tmpmat );
#else
  DistributionParameter val( scaleMatrix() );
#endif
//  DistributionParameter val( scaleMatrix() );
  return val;
}//end






//Implementation of DuMouchelDistribution 

DuMouchelDistribution::DuMouchelDistribution()
{
  tau0 = 1.0;
  tau02 = 1.0;
  last_item_drawn = 1.0;
}//end


DuMouchelDistribution::DuMouchelDistribution( double tau )
{
  tau0 = tau;
  tau02 = tau * tau;
  last_item_drawn = tau02;
}//end


DuMouchelDistribution::~DuMouchelDistribution()
{
  //do nothing
}//end

void DuMouchelDistribution::setParameter( double tau )
{
  tau0 = tau;
  tau02 = tau * tau;
  last_item_drawn = tau02;
}//end

  
DuMouchelDistribution & DuMouchelDistribution::operator=(const DuMouchelDistribution & distr )
{
  if ( this == &distr )
  {
    return *this;
  }

  tau0 = distr.tau0;
  tau02 = distr.tau02;
  last_item_drawn = distr.last_item_drawn;

  return (*this);
}//end


DuMouchelDistribution::DuMouchelDistribution( const DuMouchelDistribution & distr )
{
  tau0 = distr.tau0;
  tau02 = distr.tau02;
  last_item_drawn = distr.last_item_drawn;
}//end


Distribution * DuMouchelDistribution::clone()
{
  DuMouchelDistribution * distr = new DuMouchelDistribution( (*this) );
  return distr;
}//end



void DuMouchelDistribution::update( DistributionParameter * par_list, int list_size )
{
  tau0 = par_list[0].getScalar();
  tau02 = tau0 * tau0; 
}//end



double DuMouchelDistribution::drawOneItem()
{
  double U, odds_ratio;
  
  U = runif( 0.0, 1.0 );

  if ( U < 1.0 )
  {
    odds_ratio = ( U / ( 1 - U ) );
  }
  else
  {
    odds_ratio = VERY_LARGE_NUMBER;
  }

  last_item_drawn = tau02 * odds_ratio * odds_ratio;

  return last_item_drawn;

}//end


DistributionParameter DuMouchelDistribution::draw()
{
  DistributionParameter val( drawOneItem() );
  return val;
}//end


DistributionParameter DuMouchelDistribution::lastDraw()
{
  DistributionParameter val( last_item_drawn );
  return val;
}//end


void DuMouchelDistribution::setLastDraw( double last_draw )
{
  last_item_drawn = last_draw;
}//end


void DuMouchelDistribution::setLastDraw( DistributionParameter & last_draw )
{
  last_item_drawn = last_draw.getScalar();
}//end



DistributionParameter DuMouchelDistribution:: mode()
{
  DistributionParameter val( 0.0 );
  return val;
}//end


DistributionParameter DuMouchelDistribution::mean()
{
  DistributionParameter val( VERY_LARGE_NUMBER );
  return val;
}//end
 

DistributionParameter DuMouchelDistribution::variance()
{
  DistributionParameter val( VERY_LARGE_NUMBER );
  return val;
}//end
 

double DuMouchelDistribution::logDensity( DistributionParameter & value )
{
  double tau, density;

  tau = sqrt( value.getScalar() );
  density = ( 0.5 * tau0 ) / ( tau * ( tau0 + tau ) * ( tau0 + tau ) );

  return density;
}//end



//Implementation of UniformShrinkageDistribution 

UniformShrinkageDistribution::UniformShrinkageDistribution()
{
  tau02 = 1.0;
  last_item_drawn = 1.0;
}//end


UniformShrinkageDistribution::UniformShrinkageDistribution( double tau2 )
{
  tau02 = tau2;
  last_item_drawn = tau2;
}//end


UniformShrinkageDistribution::~UniformShrinkageDistribution()
{
  //do nothing
}//end


void UniformShrinkageDistribution::setParameter( double tau2 )
{
  tau02 = tau2;
  last_item_drawn = tau2;  
}//end


UniformShrinkageDistribution & UniformShrinkageDistribution::operator=(const UniformShrinkageDistribution & distr )
{
  if ( this == &distr )
  {
    return *this;
  }

  tau02 = distr.tau02;
  last_item_drawn = distr.last_item_drawn;

  return (*this);
}//end


UniformShrinkageDistribution::UniformShrinkageDistribution( const UniformShrinkageDistribution & distr )
{
  tau02 = distr.tau02;
  last_item_drawn = distr.last_item_drawn;
}//end


Distribution * UniformShrinkageDistribution::clone()
{
  UniformShrinkageDistribution * distr = new UniformShrinkageDistribution( (*this) );
  return distr;
}//end


void UniformShrinkageDistribution::update( DistributionParameter * par_list, int list_size )
{
  tau02 = par_list[0].getScalar();
}//end


double UniformShrinkageDistribution::drawOneItem()
{
  double U, odds_ratio;
  
  U = runif( 0.0, 1.0 );

  if ( U < 1.0 )
  {
    odds_ratio = ( U / ( 1 - U ) );
  }
  else
  {
    odds_ratio = VERY_LARGE_NUMBER;
  }

  last_item_drawn = tau02 * odds_ratio;

  return last_item_drawn;

}//end


DistributionParameter UniformShrinkageDistribution::draw()
{
  DistributionParameter val( drawOneItem() );
  return val;
}//end


DistributionParameter UniformShrinkageDistribution::lastDraw()
{
  DistributionParameter val( last_item_drawn );
  return val;
}//end


void UniformShrinkageDistribution::setLastDraw( double last_draw )
{
  last_item_drawn = last_draw;
}//end


void UniformShrinkageDistribution::setLastDraw( DistributionParameter & last_draw )
{
  last_item_drawn = last_draw.getScalar();
}//end



DistributionParameter UniformShrinkageDistribution:: mode()
{
  DistributionParameter val( 0.0 );
  return val;
}//end


DistributionParameter UniformShrinkageDistribution::mean()
{
  DistributionParameter val( VERY_LARGE_NUMBER );
  return val;
}//end
 

DistributionParameter UniformShrinkageDistribution::variance()
{
  DistributionParameter val( VERY_LARGE_NUMBER );
  return val;
}//end
 

double UniformShrinkageDistribution::logDensity( DistributionParameter & value )
{
  double tau2, density;
  
  tau2 = value.getScalar();
  density = tau02 / ( ( tau02 + tau2 ) * ( tau02 + tau2 ) );

  return density;
}//end








/* Proper non Informative posteriors for Hierarchical linear model */

ProperNonInfoPosteriorHLM::ProperNonInfoPosteriorHLM()
{
  //assuming uniform shrinkage
  setPriorDistribution( "uniform_shrinkage" );

  tau0 = 0;
  tau02 = 0;
  dfreedom = 0;
  scale = 1;
  last_item_drawn = 1;

  prob = new CVector( NUMBER_EVALS );
  delta_step = 0;

  mode_x = 0;
  mean_x = 0;
  var_x = 0;
}//end


ProperNonInfoPosteriorHLM::ProperNonInfoPosteriorHLM( double tau2, double df ) throw( rtErr )
{
  if ( tau2 <= 0.0 )
  {
    Rprintf( "ProperNonInfoPosteriorHLM: invalid parameter tau02 [%f].\n", tau2 );
    char the_error[] = "ProperNonInfoPosteriorHLM: invalid parameter tau02.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

  tau02 = tau2;
  tau0 = sqrt( tau02 );

  //assuming uniform shrinkage
  setPriorDistribution( "uniform_shrinkage" );

  dfreedom = df;
  scale = 1;
  last_item_drawn = 1;

  prob = new CVector( NUMBER_EVALS );
  delta_step = 0;

  mode_x = 0;
  mean_x = 0;
  var_x = 0;

}//end



void ProperNonInfoPosteriorHLM::initialize( double tau2, double df ) throw( rtErr )
{
  if ( tau2 <= 0.0 )
  {
    Rprintf( "ProperNonInfoPosteriorHLM: invalid parameter tau02 [%f].\n", tau2 );
    char the_error[] = "ProperNonInfoPosteriorHLM: invalid parameter tau02.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

  tau02 = tau2;
  tau0 = sqrt( tau02 );

  //assuming uniform shrinkage
  setPriorDistribution( "uniform_shrinkage" );

  dfreedom = df;
  scale = 1;
  last_item_drawn = 1;

  prob = new CVector( NUMBER_EVALS );
  delta_step = 0;

  mode_x = 0;
  mean_x = 0;
  var_x = 0;

}//end



ProperNonInfoPosteriorHLM::~ProperNonInfoPosteriorHLM()
{
  clean();
}//end


void ProperNonInfoPosteriorHLM::clean()
{
  dfreedom = 0;
  scale = 1;
  last_item_drawn = 1;
  delta_step = 0;

  tau02 = 0;
  tau0 = 0;

  mode_x = 0;
  mean_x = 0;
  var_x = 0;

  if ( prob != NULL )
  {
    delete prob;
    prob = NULL;
  }

  if ( distr_type != NULL )
  {
    delete [] distr_type;
    distr_type = NULL;
  }

}//end



void ProperNonInfoPosteriorHLM::computeCDF()
{
  int i;
  double inv_chisq_99, total, max_dens, val;

  inv_chisq_99 = ( dfreedom * scale ) / qchisq( 1.0 - 0.99, dfreedom, 1, 0 );

  delta_step = inv_chisq_99 / ((double) NUMBER_EVALS - 2 );

  total = 0;
  prob->Val( 0 ) = 0;

  CVector dens( NUMBER_EVALS );
  dens.Val( 0 ) = 0;
  max_dens = 0;
  mode_x = 0;
  mean_x = 0;
  var_x = 0;

  for ( i = 1; i < NUMBER_EVALS; i++ )
  {
    val = i * delta_step;
    dens.Val( i ) = exp( logDensity( val ) );
    total += dens.Val( i );
    prob->Val( i ) = total;

    //now compute mode
    if ( dens.Val( i ) > max_dens )
    {
      mode_x = val;
      max_dens = dens.Val(i);
    }
  }

  //normalize and as a side effect
  //now get the mean and variance
  for ( i = 1; i < NUMBER_EVALS; i++ )
  {
    prob->Val( i ) /= total;

    dens.Val( i ) =  dens.Val( i ) / total;
    mean_x += i * delta_step * dens.Val( i );
  }

  for ( i = 1; i < NUMBER_EVALS; i++ )
  {
    var_x += ( i * delta_step - mean_x ) * ( i * delta_step - mean_x ) * dens.Val( i );
  }

  //printf( "pnipHLM: mean = %f, mode = %f, var = %f\n", mean_x, mode_x, var_x );

}//end


void ProperNonInfoPosteriorHLM::setScale( double scl )
{
  scale = scl;
  last_item_drawn = scale;
  computeCDF();

}//end


void ProperNonInfoPosteriorHLM::setPriorDistribution( char * type ) throw( rtErr )
{
  if ( !strcmp( type, "uniform_shrinkage" ) )
  {
    distr_type = new char[ 18 ];
    sprintf( distr_type, "uniform_shrinkage" );
  }
  else if ( !strcmp( type, "du_mouchel" ) )
  {
    delete [] distr_type;
    distr_type = new char [ 11 ];
    sprintf( distr_type, "du_mouchel" );
  }
  else
  {
    Rprintf( "ProperNonInfoPosteriorHLM::setDistribution: type [%s] is unknown.\n", type );
    char the_error[] = "ProperNonInfoPosteriorHLM::setDistribution: type is unknown.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

  //printf( "*** setting prior/post to %s\n", distr_type );

}//end



ProperNonInfoPosteriorHLM::ProperNonInfoPosteriorHLM( const ProperNonInfoPosteriorHLM & distr )
{
  dfreedom = distr.dfreedom;
  scale = distr.scale;
  last_item_drawn = distr.last_item_drawn;
  delta_step = distr.delta_step;

  tau02 = distr.tau02;
  tau0 = distr.tau0;

  mode_x = distr.mode_x;
  mean_x = distr.mean_x;
  var_x = distr.var_x;

  if ( distr.prob != NULL )
  {
    prob = new CVector( (*distr.prob) );
  }

  if ( distr.distr_type != NULL && strlen( distr.distr_type ) > 0 )
  {
    distr_type = new char [ strlen( distr.distr_type ) + 1 ];
    sprintf( distr_type, "%s", distr.distr_type );
  }

}//end



ProperNonInfoPosteriorHLM & ProperNonInfoPosteriorHLM::operator=( const ProperNonInfoPosteriorHLM & distr )
{
  if ( this == &distr )
  {
    return *this;
  }

  clean();

  dfreedom = distr.dfreedom;
  scale = distr.scale;
  last_item_drawn = distr.last_item_drawn;
  delta_step = distr.delta_step;

  tau02 = distr.tau02;
  tau0 = distr.tau0;

  mode_x = distr.mode_x;
  mean_x = distr.mean_x;
  var_x = distr.var_x;

  if ( distr.prob != NULL )
  {
    prob = new CVector( (*distr.prob) );
  }

  if ( distr.distr_type != NULL && strlen( distr.distr_type ) > 0 )
  {
    distr_type = new char [ strlen( distr.distr_type ) + 1 ];
    sprintf( distr_type, "%s", distr.distr_type );
  }

  return (*this);
}//end


Distribution * ProperNonInfoPosteriorHLM::clone()
{
  ProperNonInfoPosteriorHLM * distr = new ProperNonInfoPosteriorHLM( (*this) );
  return distr;
}//end



double ProperNonInfoPosteriorHLM::inverseCDF( double p ) throw( rtErr )
{
  bool found;
  int i;
  double q, w1, w2;

  q = 0;

  if ( p < 1.0 && p > 0 )
  {
    found = false;
    i = 0;
    while ( !found && i < NUMBER_EVALS )
    {
      if ( prob->Val(i) >= p )
      {
        found = true;
      }
      else
      {
        i++;
      }
    }//end loop

    if ( found )
    {
      if ( prob->Val(i) == p )
      {
        q = i * delta_step;
      }
      else
      {
        if ( i > 0 )
	{
          w1 = ( p - prob->Val( i - 1 ) ) / ( prob->Val( i ) - prob->Val( i - 1) );
          w2 = ( prob->Val( i ) - p ) / ( prob->Val( i ) - prob->Val( i - 1) );
          q = ( w1 * ( i - 1 ) + w2 * i ) * delta_step;
        }
        else
	{
          Rprintf( "ProperNonInfoPosteriorHLM::inverseCDF: something wrong with input value [%f]\n", p );
          char the_error[] = "ProperNonInfoPosteriorHLM::inverseCDF: something wrong with input value.";
          rtErr runtime_error( the_error );
          throw runtime_error;
        }
      }
    }//found
    else
    {
      q = NUMBER_EVALS * delta_step;
    }
  }
  else if ( p == 1.0 )
  {
    q = NUMBER_EVALS * delta_step;
  }
  else if ( p == 0.0 )
  {
    q = 0.0;
  }
  else
  {
    Rprintf( "ProperNonInfoPosteriorHLM::inverseCDF: value provided [%f] is not valid.\n", p );
    char the_error[] = "ProperNonInfoPosteriorHLM::inverseCDF: value provided is not valid.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }

  return q;

}//end



DistributionParameter ProperNonInfoPosteriorHLM::draw()
{
  double U;
  
  U = runif( 0.0,  1.0);
  DistributionParameter x( inverseCDF( U ) );

  last_item_drawn = x.getScalar();

  return x;

}//end  


void ProperNonInfoPosteriorHLM::update( DistributionParameter * par_list, int list_size )
{
  dfreedom = par_list[0].getScalar();
  setScale( par_list[1].getScalar() );

}//end


DistributionParameter ProperNonInfoPosteriorHLM::lastDraw()
{
  DistributionParameter val( last_item_drawn );
  return val;
}//end


void ProperNonInfoPosteriorHLM::setLastDraw( double par )
{
  last_item_drawn = par;
}//end



void ProperNonInfoPosteriorHLM::setLastDraw( DistributionParameter & par )
{
  last_item_drawn = par.getScalar();
}//end


DistributionParameter ProperNonInfoPosteriorHLM::mode()
{
  DistributionParameter val( mode_x );
  return val;
}//end


DistributionParameter ProperNonInfoPosteriorHLM::mean()
{
  DistributionParameter val( mean_x );
  return val;
}//end

  
DistributionParameter ProperNonInfoPosteriorHLM::variance()
{
  DistributionParameter val( var_x );
  return val;
}//end
  

double ProperNonInfoPosteriorHLM::logDensity( DistributionParameter & value )
{
  double val, dens;

  val = value.getScalar();
  dens = logDensity( val );

  return dens;

}//end


double ProperNonInfoPosteriorHLM::logDensity( double val )
{
  //NOTE: this does not return a normalized logDensity

  double dens, denom;

  if ( !strcmp( distr_type, "uniform_shrinkage" ) )
  {
    denom = tau02 + val;
    dens = log( val / ( denom * denom ) );
  }
  else
  {
    denom = tau0 + sqrt( val );
    dens = log( tau0 / ( denom * denom ) );
  }//end

  dens -= ( 0.5 * dfreedom * log( val ) + 0.5 * ( scale / val ) );

  return dens;

}//end
