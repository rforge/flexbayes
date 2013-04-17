#include "R.h"
#include "Rmath.h"

#include "RejectionSampling.h"
#include "Vector.h"


RejectionSampling::RejectionSampling()
{
  number_of_samples = 0;
  number_drawings = 0;
  save_stats = false;
  acceptance_stats = 0;
}//end


RejectionSampling::RejectionSampling( int number_samples )
{
  number_of_samples = number_samples;
  save_stats = false;
  number_drawings = 0;
  acceptance_stats = 0;
}//end


RejectionSampling::~RejectionSampling()
{
  //do nothing
}//end


RejectionSampling & RejectionSampling::operator=( const RejectionSampling & rej )
{

  if ( this == &rej )
  {
    return *this;
  }

  number_of_samples = rej.number_of_samples;
  save_stats = rej.save_stats;
  number_drawings = rej.number_drawings;
  acceptance_stats = rej. acceptance_stats;

  return *this;
}//end


RejectionSampling::RejectionSampling(  const RejectionSampling & rej )
{
  number_of_samples = rej.number_of_samples;
  save_stats = rej.save_stats;
  number_drawings = rej.number_drawings;
  acceptance_stats = rej. acceptance_stats;
}//end



DistributionParameter ** RejectionSampling::drawSamples( RejectionSamplingDistribution * densities )
{
  int i;
  bool accepted;
  double U;
  DistributionParameter ** samples;

  samples = new DistributionParameter * [ number_of_samples ];
  for ( i = 0; i < number_of_samples; i++ )
  {
    accepted = false;
    while ( !accepted )
    {
      DistributionParameter Y( densities->draw() );
      U = runif( 0.0, 1.0 );
      if ( U <= densities->ratioOfDensities( Y ) )
      {
        accepted = true;
        samples[ i ] = new DistributionParameter( Y );

        if ( save_stats )
	{
          acceptance_stats++;
        }
      }
      
      if ( save_stats )
      {
        number_drawings++;
      }

    }//end while loop

  }//end for loop

  return (samples);

}//end


DistributionParameter ** RejectionSampling::drawSamples( int n_samples, RejectionSamplingDistribution * densities )
{
  number_of_samples = n_samples;
  return drawSamples( densities );

}//end


DistributionParameter RejectionSampling::drawASample( RejectionSamplingDistribution * densities )
{
  bool accepted;
  double U;
  DistributionParameter * sample;

  accepted = false;
  while ( !accepted )
  {
    DistributionParameter Y( densities->draw() );
    U = runif( 0.0, 1.0 );
    if ( U <= densities->ratioOfDensities( Y ) )
    {
      accepted = true;
      sample = new DistributionParameter( Y );

      if ( save_stats )
      {
        acceptance_stats++;
      }
    }
    
    if ( save_stats )
    {
      number_drawings++;
    }

  }//end while loop

  DistributionParameter par_sample( (*sample) );
  delete sample;

  return par_sample;
}//end



DistributionParameter RejectionSampling::AcceptanceRate()
{
  CVector stats( 3 );

  if ( save_stats )
  {
    if ( number_drawings > 0 )
    {
      stats.Val( 0 ) = acceptance_stats / ( (double) number_drawings );
    }
    else
    {
      stats.Val( 0 ) = 0;
    }
     
    stats.Val( 1 ) = acceptance_stats;
    stats.Val( 2 ) = number_drawings;
  }//end if

  DistributionParameter par_val( stats );
  return par_val;
}//end
