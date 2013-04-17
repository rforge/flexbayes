#ifndef BAYES_REJECTION_SAMPLING_DENSITIES_H
#define BAYES_REJECTION_SAMPLING_DENSITIES_H

#include "DistributionParameter.h"

/*
  abstract class
  to be used with RejectionSampling class
*/
class RejectionSamplingDensities
{
  public:
    virtual ~RejectionSamplingDensities() { };
    virtual double ratioOfDensities( DistributionParameter & par ) = 0;
    virtual DistributionParamater samplingDensity() = 0;
}; //end RejectionSamplingDensities class

#endif /* BAYES_REJECTION_SAMPLING_DENSITIES_H */
