

#ifndef BAYES_DISTRIBUTION_H
#define BAYES_DISTRIBUTION_H

#include "DistributionParameter.h"

class Distribution
{
  public:
    virtual ~Distribution() { };
    virtual Distribution * clone() = 0;
    virtual double logDensity( DistributionParameter & value ) = 0;
    virtual DistributionParameter draw() = 0;
    virtual void update( DistributionParameter * par_list, int list_size ) = 0;
    virtual DistributionParameter lastDraw() = 0;
    virtual void setLastDraw( DistributionParameter & par ) = 0;
    virtual DistributionParameter mode() = 0;
    virtual DistributionParameter mean() = 0;
    virtual DistributionParameter variance() = 0;
};//end class


#endif /* BAYES_DISTRIBUTION_H */
