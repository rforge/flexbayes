#include "R.h"
#include "Rmath.h"

#include "RejectionSamplingDistributionExtended.h"


CInvChisqByRejSampling::CInvChisqByRejSampling()
{
  inv_chisq = NULL;
  interval = NULL;
  left_bounded_only = false;
  right_bounded_only = false;
}//end


CInvChisqByRejSampling::CInvChisqByRejSampling( double nu0 )
{
  inv_chisq = new InvChisqDistribution( nu0 );
}//end


void CInvChisqByRejSampling::setInterval( double lower, double upper )
{
  interval = new CVector( 2 );
  interval->Val(0) = lower;
  interval->Val(1) = upper;
  left_bounded_only = false;
  right_bounded_only = false;
}//end


void CInvChisqByRejSampling::setInterval( double bound, char * type )
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

}//end


double CInvChisqByRejSampling::upperBound()
{
  double val;
  val = -1;

  if ( interval != NULL )
  {
    if ( !left_bounded_only && !right_bounded_only )
    {
      val = interval->Val( 1 );
    }
    else if ( right_bounded_only )
    {
      val = interval->Val( 0 );
    }
  }

  return val;
}//end



double CInvChisqByRejSampling::lowerBound()
{
  double val;
  val = -1;

  if ( interval != NULL )
  {
    if ( !left_bounded_only && !right_bounded_only )
    {
      val = interval->Val( 0 );
    }
    else if ( left_bounded_only )
    {
      val = interval->Val( 0 );
    }
  }

  return val;
}//end



CInvChisqByRejSampling::~CInvChisqByRejSampling()
{
  if ( inv_chisq != NULL )
  {
    delete inv_chisq;
  }

  if ( interval != NULL )
  {
    delete interval;
  }
}//end


void CInvChisqByRejSampling::clean()
{
  if ( inv_chisq != NULL )
  {
    delete inv_chisq;
  }

  if ( interval != NULL )
  {
    delete interval;
  }
}//end


CInvChisqByRejSampling::CInvChisqByRejSampling( const CInvChisqByRejSampling & distr )
{
  inv_chisq = NULL;
  interval = NULL;

  if ( distr.inv_chisq != NULL )
  {
    inv_chisq = new InvChisqDistribution( ( *(distr.inv_chisq) ) );
  }

  if ( distr.interval != NULL )
  {
    interval = new CVector( ( *(distr.interval) ) );
  }

  left_bounded_only = distr.left_bounded_only;
  right_bounded_only = distr.right_bounded_only;
}//end


CInvChisqByRejSampling & CInvChisqByRejSampling::operator=( const CInvChisqByRejSampling & distr )
{
  if ( this == &distr )
  {
    return *this;
  }

  clean();

  inv_chisq = NULL;
  interval = NULL;

  if ( distr.inv_chisq != NULL )
  {
    inv_chisq = new InvChisqDistribution( ( *(distr.inv_chisq) ) );
  }

  if ( distr.interval != NULL )
  {
    interval = new CVector( ( *(distr.interval) ) );
  }

  left_bounded_only = distr.left_bounded_only;
  right_bounded_only = distr.right_bounded_only;

  return *this;
}//end


double CInvChisqByRejSampling::ratioOfDensities( DistributionParameter & val )
{
  bool in_interval;
  double x, ratio;
  
  x = val.getScalar();
  in_interval = false;

  if ( !left_bounded_only && !right_bounded_only )
  {
    in_interval = ( x >= interval->Val(0) && x <= interval->Val( 1 ) );
  }
  else if ( left_bounded_only )
  {
    in_interval = ( x >= interval->Val(0) );
  }
  else if ( right_bounded_only )
  {
    in_interval = ( x <= interval->Val(0) );
  }

  if ( in_interval )
  {
    ratio = 1;
  }
  else
  {
    ratio = 0;
  }

  return ratio;

}//end


DistributionParameter CInvChisqByRejSampling::draw()
{
  return inv_chisq->draw();
}//end






/* 
  implementing Uniform Shrinkage and DuMouchel for Bayesian Hierarchical Linear Model
*/

ProperNonInfoPosteriorForBHLM::ProperNonInfoPosteriorForBHLM()
{
  distr_type = NULL;
  scale_var = 1.0;
  dfreedom_a = 0;;
  dfreedom_b = 0;
  tau02 = 0;
  tau0 = 0;
  const_tau0_ratio_gammas = 0;
  ratio_gammas_b_over_a = 0;
  ratio_probs_b_over_a = 0;
  ratio_a_over_b = 0;
  inv_chisq_a = NULL;
  inv_chisq_b = NULL;
  min_ratio = NULL;
  last_item_drawn = tau02;
}//end


ProperNonInfoPosteriorForBHLM::ProperNonInfoPosteriorForBHLM( double df_a, double df_b, double tau_0_2 )
{
  //assuming uniform shrinkage
  distr_type = new char[ 18 ];
  sprintf( distr_type, "uniform_shrinkage" );

  dfreedom_a = df_a;;
  dfreedom_b = df_b;
  tau02 = tau_0_2;
  tau0 = sqrt( tau02 );

  inv_chisq_a = new ConstrainedInvChisqDistribution( tau02, "right" );
  inv_chisq_a->setDegreesOfFreedom( df_a );
  inv_chisq_b = new ConstrainedInvChisqDistribution( tau02, "left" );
  inv_chisq_a->setDegreesOfFreedom( df_b );

  //compute log of ratios
  ratio_gammas_b_over_a = lgammafn( 0.5 * df_b ) - lgammafn( 0.5 * df_a );
  const_tau0_ratio_gammas = ( (df_b - df_a) / 2 ) * log( 2 * tau02 ) + ratio_gammas_b_over_a;

  //not yet defined
  ratio_probs_b_over_a = 0;
  ratio_a_over_b = 0;
  scale_var = 1.0;

  min_ratio = new CVector( 2 );
  last_item_drawn = tau02;
}//end



void ProperNonInfoPosteriorForBHLM::initialize( double df_a, double df_b, double tau_0_2 )
{
  //assuming uniform shrinkage
  distr_type = new char[ 18 ];
  sprintf( distr_type, "uniform_shrinkage" );

  dfreedom_a = df_a;;
  dfreedom_b = df_b;
  tau02 = tau_0_2;
  tau0 = sqrt( tau02 );

  inv_chisq_a = new ConstrainedInvChisqDistribution( tau02, "right" );
  inv_chisq_a->setDegreesOfFreedom( df_a );
  inv_chisq_b = new ConstrainedInvChisqDistribution( tau02, "left" );
  inv_chisq_a->setDegreesOfFreedom( df_b );

  //compute log of ratios
  ratio_gammas_b_over_a = lgammafn( 0.5 * df_b ) - lgammafn( 0.5 * df_a );
  const_tau0_ratio_gammas = ( (df_b - df_a) / 2 ) * log( 2 * tau02 ) + ratio_gammas_b_over_a;

  //not yet defined
  ratio_probs_b_over_a = 0;
  ratio_a_over_b = 0;
  scale_var = 1.0;

  min_ratio = new CVector( 2 );
  last_item_drawn = tau02;
}//end


void ProperNonInfoPosteriorForBHLM::setPriorDistribution( char * type )
{
  if ( !strcmp( type, "uniform_shrinkage" ) )
  {
    //do nothing
    //default
  }
  else if ( !strcmp( type, "du_mouchel" ) )
  {
    delete [] distr_type;
    distr_type = new char [ 11 ];
    sprintf( distr_type, "du_mouchel" );
  }
  else
  {
    Rprintf( "ProperNonInfoPosteriorForBHLM::setDistribution: type [%s] is unknown.\n", type );
    // Bug fix by Dawn:  the following was:   exit(1);    which causes problems within SPLUS
    MESSAGE "ProperNonInfoPosteriorForBHLM::setDistribution: type [%s] is unknown.\n" ERROR;
  }
}//end



ProperNonInfoPosteriorForBHLM::ProperNonInfoPosteriorForBHLM( const ProperNonInfoPosteriorForBHLM & distr )
{
  scale_var = distr.scale_var;
  dfreedom_a = distr.dfreedom_a;
  dfreedom_b = distr.dfreedom_b;
  tau02 = distr.tau02;
  tau0 = distr.tau0;
  const_tau0_ratio_gammas = distr.const_tau0_ratio_gammas;
  ratio_gammas_b_over_a = distr.ratio_gammas_b_over_a;
  ratio_probs_b_over_a = distr.ratio_probs_b_over_a;
  ratio_a_over_b = distr.ratio_a_over_b;
  last_item_drawn = distr.last_item_drawn;

  if ( distr.distr_type != NULL && strlen( distr.distr_type ) > 0 )
  {
    distr_type = new char [ strlen( distr.distr_type ) + 1 ];
    sprintf( distr_type, "%s", distr.distr_type );
  }
  else
  {
    distr_type = NULL;
  }

  if ( distr.inv_chisq_a != NULL )
  {
    inv_chisq_a = new ConstrainedInvChisqDistribution( ( *(distr.inv_chisq_a) ) );
  }
  else
  {
    inv_chisq_a = NULL;
  }

  if ( distr.inv_chisq_b != NULL )
  {
    inv_chisq_b = new ConstrainedInvChisqDistribution( ( *(distr.inv_chisq_b) ) );
  }
  else
  {
    inv_chisq_b = NULL;
  }

  if ( distr.min_ratio != NULL )
  {
    min_ratio = new CVector( ( *(distr.min_ratio) ) );
  }

}//end



ProperNonInfoPosteriorForBHLM & ProperNonInfoPosteriorForBHLM::operator=( const ProperNonInfoPosteriorForBHLM & distr )
{
  if ( this == &distr )
  {
    return *this;
  }

  clean();

  scale_var = distr.scale_var;
  dfreedom_a = distr.dfreedom_a;
  dfreedom_b = distr.dfreedom_b;
  tau02 = distr.tau02;
  tau0 = distr.tau0;
  const_tau0_ratio_gammas = distr.const_tau0_ratio_gammas;
  ratio_gammas_b_over_a = distr.ratio_gammas_b_over_a;
  ratio_probs_b_over_a = distr.ratio_probs_b_over_a;
  ratio_a_over_b = distr.ratio_a_over_b;
  last_item_drawn = distr.last_item_drawn;

  if ( distr.distr_type != NULL && strlen( distr.distr_type ) > 0 )
  {
    distr_type = new char [ strlen( distr.distr_type ) + 1 ];
    sprintf( distr_type, "%s", distr.distr_type );
  }

  if ( distr.inv_chisq_a != NULL )
  {
    inv_chisq_a = new ConstrainedInvChisqDistribution( ( *(distr.inv_chisq_a) ) );
  }
  else
  {
    inv_chisq_a = NULL;
  }

  if ( distr.inv_chisq_b != NULL )
  {
    inv_chisq_b = new ConstrainedInvChisqDistribution( ( *(distr.inv_chisq_b) ) );
  }
  else
  {
    inv_chisq_b = NULL;
  }

  if ( distr.min_ratio != NULL )
  {
    min_ratio = new CVector( ( *(distr.min_ratio) ) );
  }

  return *this;
}//end


Distribution * ProperNonInfoPosteriorForBHLM::clone()
{
  ProperNonInfoPosteriorForBHLM * distr = new ProperNonInfoPosteriorForBHLM( (*this) );
  return distr;
}//end



void ProperNonInfoPosteriorForBHLM::clean()
{
  if ( inv_chisq_a != NULL )
  {
    delete inv_chisq_a;
  }

  if ( inv_chisq_b != NULL )
  {
    delete inv_chisq_b;
  }

  if ( min_ratio != NULL )
  {
    delete min_ratio;
  }

  if ( distr_type != NULL )
  {
    delete [] distr_type;
  }

  scale_var = 1.0;
  dfreedom_a = 0;;
  dfreedom_b = 0;
  tau02 = 0;
  tau0 = 0;
  const_tau0_ratio_gammas = 0;
  ratio_gammas_b_over_a = 0;
  ratio_probs_b_over_a = 0;
  ratio_a_over_b = 0;
  inv_chisq_a = NULL;
  inv_chisq_b = NULL;
  min_ratio = NULL;
  distr_type = NULL;
  last_item_drawn = tau02;
}//end
  

ProperNonInfoPosteriorForBHLM::~ProperNonInfoPosteriorForBHLM()
{
  if ( inv_chisq_a != NULL )
  {
    delete inv_chisq_a;
  }

  if ( inv_chisq_b != NULL )
  {
    delete inv_chisq_b;
  }

  if ( min_ratio != NULL )
  {
    delete min_ratio;
  }

  if ( distr_type != NULL )
  {
    delete [] distr_type;
  }

}//end


void ProperNonInfoPosteriorForBHLM::setScale( double val )
{
  double p_a, p_b;

  Rprintf( "PNIP: setScale val = %f\n", val );

  if ( val > 0 )
  {
    inv_chisq_a->setScale( val );
    inv_chisq_b->setScale( val );
   
    ratio_a_over_b = exp( const_tau0_ratio_gammas - ( ( dfreedom_b - dfreedom_a ) / 2 ) * log( val ) );

    p_b = pchisq( val * tau02, dfreedom_b, 1, 0 );
    p_a = 1.0 - pchisq( val / tau02, dfreedom_a, 1, 0 );

    ratio_probs_b_over_a = p_b / p_a;
    ratio_a_over_b *= ratio_probs_b_over_a;

    if ( ratio_a_over_b > 1.0 )
    {
      min_ratio->Val( 0 ) = 1.0 / ratio_a_over_b;
      min_ratio->Val( 1 ) = 1.0;
    }
    else
    {
      min_ratio->Val( 0 ) = 1.0;
      min_ratio->Val( 1 ) = ratio_a_over_b;
    }
  }
  else
  {
    Rprintf( "ProperNonInfoPosteriorForBHLM::setScale: Zero or negative scale [%f] is not valid.\n", val );
    // Bug fix by Dawn:  the following was:    exit(1);  which causes problems within SPLUS
    MESSAGE "ProperNonInfoPosteriorForBHLM::setScale: Zero or negative scale [%f] is not valid.\n" ERROR;
  }

}//end


double ProperNonInfoPosteriorForBHLM::ratioOfDensities( DistributionParameter & par )
{
  int index;
  double val, ratio, var_tau2, sqrt_val;

  val = par.getScalar();
  if ( val <= tau02 )
  {
    index = 0;
    var_tau2 = tau02;
  }
  else
  {
    index = 1;
    var_tau2 = val;
  }
   
  if ( !strcmp( distr_type, "uniform_shrinkage" ) )
  {
    ratio = min_ratio->Val( index ) * ( ( var_tau2 * var_tau2 ) / ( ( tau02 + val ) * ( tau02 + val ) ) );
  }
  else if ( !strcmp( distr_type, "du_mouchel" ) )
  {
    sqrt_val = sqrt( val );
    ratio = min_ratio->Val( index ) * ( var_tau2 / ( ( tau0 + sqrt_val ) * ( tau0 + sqrt_val ) ) );
  }

  return ratio;

}//end



DistributionParameter ProperNonInfoPosteriorForBHLM::draw()
{
  double U, val;

  Rprintf( "PNIP: draw \n");

  U = runif( 0.0, 1.0 );
  if ( U <= 0.5 )
  {
    Rprintf( "PNIP: will draw inv_chisq A \n");
    val = inv_chisq_a->drawOneItem();
  }
  else
  {
    Rprintf( "PNIP: will draw inv_chisq B \n");
    val = inv_chisq_b->drawOneItem();
  }    

  Rprintf( "PNIP: item drawn is %f\n", val );

  last_item_drawn = val;

  DistributionParameter par_val( val );
  return par_val;

}//end


double ProperNonInfoPosteriorForBHLM::logDensity( DistributionParameter & value )
{
  double val, val_a, val_b, val_ba;

  val_a = inv_chisq_a->logDensity( value );
  val_b = inv_chisq_b->logDensity( value );

  val_ba = val_b - val_a;
  if ( val_ba > 0 )
  {
    val = val_b  + log( 1.0  + exp( - val_ba ) );
  }
  else
  {
    val = val_a  + log( 1.0 + exp( val_ba ) );
  }

  val = log( 0.5 ) +  val;
 
  return val;
}//end



void ProperNonInfoPosteriorForBHLM::update( DistributionParameter * par_list, int list_size )
{
  inv_chisq_a->update( par_list, list_size );
  inv_chisq_b->update( par_list, list_size );
}//end



void ProperNonInfoPosteriorForBHLM::update( double add_df, double add_scale )
{

  inv_chisq_a->update( add_df, add_scale );
  inv_chisq_b->update( add_df, add_scale );
}//end



DistributionParameter ProperNonInfoPosteriorForBHLM::lastDraw()
{
  DistributionParameter val( last_item_drawn );
  return val;
}//end



void ProperNonInfoPosteriorForBHLM::setLastDraw( double par )
{
  last_item_drawn = par;
}//end



void ProperNonInfoPosteriorForBHLM::setLastDraw( DistributionParameter & par )
{
  last_item_drawn = par.getScalar();
}//end


DistributionParameter ProperNonInfoPosteriorForBHLM::mode()
{
  //this will not return the true mode
  //just the average of the two modes

  DistributionParameter val( 0.5 * ( inv_chisq_a->mode().getScalar() + inv_chisq_b->mode().getScalar() ) );
  return val;
}//end


DistributionParameter ProperNonInfoPosteriorForBHLM::mean()
{
  DistributionParameter val( 0.5 * ( inv_chisq_a->mean().getScalar() + inv_chisq_b->mean().getScalar() ) );
  return val;
}//end



DistributionParameter ProperNonInfoPosteriorForBHLM::variance()
{
  DistributionParameter val( 0.5 * ( inv_chisq_a->variance().getScalar() + inv_chisq_b->variance().getScalar() ) );
  return val;
}//end




