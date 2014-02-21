#include "R.h"
#include "Rmath.h"

#include "DistributionMixture.h"


extern int randomStart( int n );

DistributionMixture::DistributionMixture()
{
  n_mixtures = 0;
  distr_list = NULL;
  initial_props = NULL;
  props = NULL;
  cumu_props = NULL;
  draw_stats = NULL;
  save_stats = false;
  last_component_drawn = 0;
}//end


DistributionMixture::DistributionMixture( int n_mixs )
{
  if ( n_mixs > 0 )
  {
    n_mixtures = n_mixs;
    distr_list = new Distribution * [ n_mixtures ];
    props = new CVector( n_mixtures );
    props->setToZero();
    initial_props = new CVector( n_mixtures );
    initial_props->setToZero();
    cumu_props = new CVector( n_mixtures );
    cumu_props->setToZero();
    
    save_stats = false;
    draw_stats = NULL;

    last_component_drawn = 0;
  }
  else
  {
    Rprintf( "\nDistributionMixture: Number of mixtures must be positive.\n" );
  }
}//end



void DistributionMixture::saveStatistics( bool val )
{
  save_stats = val;
  if ( draw_stats == NULL && save_stats )
  {
    draw_stats = new CVector( n_mixtures );
    draw_stats->setToZero();
  }
}//end


DistributionMixture::~DistributionMixture()
{
  if ( n_mixtures > 0 )
  {
    delete [] distr_list;
    n_mixtures = 0;
    delete props;
    delete initial_props;
    delete cumu_props;
    if ( draw_stats != NULL )
    {
      delete draw_stats;
    }
    save_stats = false;
    last_component_drawn = 0;
  }
}//end


void DistributionMixture::set( int index, Distribution * distr, double prop )
{
  int i;

  distr_list[ index ] = distr;
  props->Val( index ) = prop;
  initial_props->Val( index ) = prop;

  for ( i = index; i < n_mixtures; i++ )
  {
    cumu_props->Val( i ) += prop;
  }

}//end



DistributionParameter * DistributionMixture::stratifiedDrawFromMixture( int number_draws )
{
  int i, j, k, n_draws, draws_component, start_index, total_draws, missing_draws;

  if ( number_draws <= 0 )
  {
    Rprintf( "DistributionMixture::drawFromMixture: wrong number of  draws [%d] specified. Assuming one draw.\n", number_draws );
    n_draws = 1;
  }
  else
  {
    n_draws = number_draws;
  }

  DistributionParameter * values;
  values = new DistributionParameter [ n_draws ];
 
  start_index = randomStart( n_mixtures );

  total_draws = 0;
  for ( i = 0; i < n_mixtures; i++ )
  {
    j = ( i + start_index ) % n_mixtures;
    draws_component = (int) ( ( (double) n_draws) * props->Val( j ) );
    for ( k = 0; k < draws_component; k++ )
    {
      values[ total_draws ] = distr_list[ j ]->draw();
      total_draws++;
    }
  }

  if ( total_draws < n_draws )
  {
    missing_draws = n_draws - total_draws;
    for ( k = 0; k < missing_draws; k++ )
    {
      values[ total_draws ] = distr_list[ start_index ]->draw();
      total_draws++;
    }
    j = start_index;
  }

  last_component_drawn = j;

  return values;
}//end



void DistributionMixture::updateProportion( int i, double val )
{

  if ( i < n_mixtures && i >= 0 )
  {
    if ( val >= 0 )
    {
      props->Val(i) = val;
    }
    else
    {
      Rprintf( "DistributionMixture::updateProportion: negative proportion [%f] ignored.\n", val );
    }
  }
  else
  {
    Rprintf( "DistributionMixture::updateProportion: index [%d] out of range.\n", i );
  }

}//end


void DistributionMixture::normalizeProportions()
{
  int i;
  double total;

  total = 0;
  for ( i = 0; i < n_mixtures; i++ )
  {
    total += props->Val( i );
  }//end for i

  if ( total > 0 )
  {
    for ( i = 0; i < n_mixtures; i++ )
    {
      props->Val( i ) /= total;
      if ( i == 0 )
      {    
        cumu_props->Val( i ) = props->Val( i );
      }
      else
      {
        cumu_props->Val( i ) = cumu_props->Val( i - 1 ) + props->Val( i );
      }
    }//end for i
  }
  else
  {
    Rprintf( " DistributionMixture::normalizeProportions: zero or negative cumulative probability [%f].\n", total );
  }

}//end


void DistributionMixture::updateLogProportion( int i, double val )
{

  if ( i < n_mixtures && i >= 0 )
  {
    props->Val(i) = val;
  }
  else
  {
    Rprintf( "DistributionMixture::updateProportion: index [%d] out of range.\n", i );
  }

}//end


void DistributionMixture::normalizeLogProportions()
{
  int i;
  double max_prop, total;

  if ( n_mixtures > 1 )
  {
    max_prop = props->Val(0);
    for ( i = 1; i < n_mixtures; i++ )
    {
      if ( props->Val(i) > max_prop )
      {
        max_prop = props->Val(i);
      }
    }//end for i

    total = 0;
    for ( i = 0; i < n_mixtures; i++ )
    {
      if ( props->Val(i) - max_prop >= -200 )
      {
        props->Val(i) = exp( props->Val(i) - max_prop );
        total += props->Val(i);
      }
      else
      {
        props->Val(i) = 0;
      }
    }//end for i

    if ( total > 0 )
    {
      for ( i = 0; i < n_mixtures; i++ )
      {
        props->Val( i ) /= total;
        if ( i == 0 )
        {    
          cumu_props->Val( i ) = props->Val( i );
        }
        else
        {
          cumu_props->Val( i ) = cumu_props->Val( i - 1 ) + props->Val( i );
        }
      }//end for i
    }
    else
    {
      Rprintf( " DistributionMixture::normalizeProportions: zero or negative cumulative probability [%f].\n", total );
    }
  }
  else
  {
    props->Val(0) = 1;
    cumu_props->Val(0) = 1;
  }
}//end



void DistributionMixture::viewCurrentProportions()
{
  props->Print();
}//end


void DistributionMixture::viewInitialProportions()
{
  initial_props->Print();
}//end
   

void DistributionMixture::printDrawingStats()
{
  int i;
  double total;

  if ( draw_stats != NULL )
  {
    Rprintf( "\nDistributionMixture: Drawing Statistics: \n" );
    total = draw_stats->Len() * draw_stats->Mean();
    for ( i = 0; i < n_mixtures; i++ )
    {
      Rprintf( "component[ %d ] = %f (%f out of %f).\n", i, draw_stats->Val(i) / total, draw_stats->Val(i), total );
    }
  
    Rprintf( "\n" );
  }

}//end


void DistributionMixture::printDrawingStats( double * stats )
{
  int i;

  if ( draw_stats != NULL )
  {
    for ( i = 0; i < n_mixtures; i++ )
    {
      stats[i] = draw_stats->Val(i);
    }
  }
}//end


void DistributionMixture::setLastDraw( DistributionParameter & par )
{
  int i;

  for ( i = 0; i < n_mixtures; i++ )
  {
    distr_list[i]->setLastDraw( par );
  }

}//end



void DistributionMixture::setLastDrawForComponent( int comp_index, DistributionParameter & par )
{
  distr_list[ comp_index ]->setLastDraw( par );
}//end



DistributionParameter DistributionMixture::lastDraw()
{
  return distr_list[ last_component_drawn ]->lastDraw();

}//end



DistributionParameter DistributionMixture::lastDrawFromComponent( int comp_index )
{
  return distr_list[ comp_index ]->lastDraw();

}//end



DistributionParameter DistributionMixture::draw()
{
  bool found;
  int i, k, index;
  double random_index;
  DistributionParameter value;

  if ( cumu_props->Val( n_mixtures - 1 ) != 1.0 && cumu_props->Val( n_mixtures - 1 ) > 0.0 )
  {
    for ( i = 0; i < n_mixtures; i++ )
    {
      props->Val( i ) /= cumu_props->Val( n_mixtures - 1 );
    }

    cumu_props->Val(0) = props->Val(0);
    for ( i = 1; i < n_mixtures; i++ )
    {
      cumu_props->Val( i ) = cumu_props->Val( i - 1 ) + props->Val(i);
    }
  }
  else if ( cumu_props->Val( n_mixtures - 1 ) <= 0.0 )
  {
    Rprintf( "DistributionMixture::draw: zero or negative proportions.\n" );
    // Bug fix by Dawn: the following was: exit(1);   which causes problems within S-PLUS.
    MESSAGE "DistributionMixture::draw: zero or negative proportions.\n" ERROR;
  }

  random_index = runif( 0.0, 1.0 );

  found = false;
  k = 0;
  while ( k < n_mixtures && !found )
  {
    if ( random_index <= cumu_props->Val( k ) )
    {
      index = k;
      found = true;
    }
    else
    {
      k++;
    }
  }//end loop

  //printf("mix: index = %d, found = %d\n", index, found );

  if ( !found )
  {
    //this condition should never be reached
    Rprintf( "DistributionMixture::draw: something wrong in the choice of component [random index = %f].\n", random_index );
    index = n_mixtures - 1;  
  }

  value = distr_list[ index ]->draw();

  //printf("mix draw is \n");
  //value.Print();

  last_component_drawn = index;

  if ( save_stats )
  {
    draw_stats->Val( index ) += 1;
  }

  return value;
}//end



DistributionParameter DistributionMixture::drawFromComponent( int comp_index )
{
  DistributionParameter value( distr_list[ comp_index ]->draw() );

  last_component_drawn = comp_index;

  if ( save_stats )
  {
    draw_stats->Val( comp_index ) += 1;
  }

  return value;
}//end



DistributionParameter DistributionMixture::mean()
{
  int i, rows, cols;

  rows = distr_list[0]->mean().getMatrix().Row();
  cols = distr_list[0]->mean().getMatrix().Col();

  CMatrix mu( rows, cols );
  CMatrix partial_mu( rows, cols );

  mu.setToZero();
  for ( i = 0; i < n_mixtures; i++ )
  {
    partial_mu = distr_list[i]->mean().getMatrix();
    partial_mu.multiplyByScalar( props->Val(i) );
    mu.add( partial_mu );
  }

  DistributionParameter mu_par( mu );

  return mu_par;
  
}//end


DistributionMixture::DistributionMixture( const DistributionMixture & distr )
{
  int i;

  n_mixtures = distr.n_mixtures;
  last_component_drawn = distr.last_component_drawn;
  props = new CVector( (*(distr.props)) );
  distr_list = new Distribution * [ n_mixtures ];
  for ( i = 0; i < n_mixtures; i++ )
  {
    distr_list[i] = distr.distr_list[i]->clone();
  }

  save_stats = distr.save_stats;
  if ( distr.draw_stats != NULL )
  {
    draw_stats = new CVector( (*(distr.draw_stats)) );
  }
  else
  {
    draw_stats = NULL;
  }

}//end


Distribution * DistributionMixture::clone()
{
  DistributionMixture * mix = new DistributionMixture( (*this) );
  return mix;
}//end
  

DistributionParameter DistributionMixture::mode()
{

  Rprintf( "DistributionMixture::mode(): Sorry this method is not implemented yet. It returns the mean.\n" );
  return mean();

}//end


DistributionParameter DistributionMixture::variance()
{
  int i, rows, cols;

  rows = distr_list[0]->variance().getMatrix().Row();
  cols = distr_list[0]->variance().getMatrix().Col();

  CMatrix cov( rows, cols );
  CMatrix partial_cov( rows, cols );

  cov.setToZero();
  for ( i = 0; i < n_mixtures; i++ )
  {
    partial_cov = distr_list[i]->variance().getMatrix();
    partial_cov.multiplyByScalar( props->Val(i) );
    cov.add( partial_cov );
  }

  DistributionParameter cov_par( cov );

  return cov_par;
  
}//end


double DistributionMixture::logDensity( DistributionParameter & value )
{
  double SMALL_EXP = -100;

  int i, max_i;
  double max_val, density;
  double * val = new double[ n_mixtures ];

  val[0] = distr_list[0]->logDensity( value );
  max_val = val[0];
  max_i = 0;
  for ( i = 1; i < n_mixtures; i++ )
  {
    val[i] = distr_list[i]->logDensity( value );
    if ( val[i] > max_val )
    {
      max_val = val[i];
      max_i = i;
    }
  }

  density = 0;
  for ( i = 1; i < n_mixtures; i++ )
  {
    val[i] -= max_val;
    if ( val[i] < SMALL_EXP )
    {
      val[i] = 0;
    }
    else
    {
      val[i] = exp( val[i] ) * props->Val(i);
    }
    density += val[i];
  }

  density = log( density ) + max_val;
  return density;

}//end 


void DistributionMixture::update( DistributionParameter * par_list, int list_size )
{

  Rprintf( "DistributionMixture::update: Sorry you need to update each component separately, since it is not known what the components are. \n" );

}//end
